"""Build the primary-voxel response dataset from Geant4 outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import uproot

from .schema import ENCOUNTER_KEY, PRIMARY_KEY


KEV_IN_JOULES = 1.602176634e-16


def _tree_frame(root: uproot.ReadOnlyDirectory, tree: str, columns: list[str]) -> pd.DataFrame:
    arrays = root[tree].arrays(columns, library="np")
    return pd.DataFrame({name: arrays[name] for name in columns})


def _decode_strings(frame: pd.DataFrame) -> pd.DataFrame:
    for column in frame.columns:
        if frame[column].dtype.kind not in {"O", "S", "U"}:
            continue
        frame[column] = frame[column].map(
            lambda value: value.decode() if isinstance(value, bytes) else value
        )
    return frame


def _aggregate_steps(steps: pd.DataFrame) -> pd.DataFrame:
    keys = ["SeedID", "EventID", "TrackID", "VoxelID"]
    steps = steps.sort_values(keys + ["StepNumber"])
    grouped = steps.groupby(keys, sort=False, observed=True)
    result = grouped.agg(
        ExitEnergy_eV=("KEPost_eV", "last"),
        Edep_eV=("Edep_eV", "sum"),
        TrackLength_nm=("StepLength_nm", "sum"),
        StepCount=("StepNumber", "size"),
        MinEnergy_eV=("KEPost_eV", "min"),
    ).reset_index()
    result["MeanLET_keV_um"] = np.divide(
        result["Edep_eV"],
        result["TrackLength_nm"],
        out=np.zeros(len(result), dtype=float),
        where=result["TrackLength_nm"].to_numpy() > 0,
    )
    return result


def _assign_splits(frame: pd.DataFrame, seed: int = 20260923) -> pd.DataFrame:
    primaries = frame[PRIMARY_KEY].drop_duplicates().copy()
    primaries["split"] = ""
    rng = np.random.default_rng(seed)
    for _, indices in primaries.groupby("case_id", sort=True).groups.items():
        indices = np.asarray(list(indices))
        rng.shuffle(indices)
        n = len(indices)
        train_end = int(round(0.70 * n))
        val_end = int(round(0.85 * n))
        primaries.loc[indices[:train_end], "split"] = "train"
        primaries.loc[indices[train_end:val_end], "split"] = "val"
        primaries.loc[indices[val_end:], "split"] = "test"
    return frame.merge(primaries, on=PRIMARY_KEY, how="left", validate="many_to_one")


def build_case(case: dict, input_dir: Path) -> pd.DataFrame:
    root_path = input_dir / case["transport_root"]
    damage_path = input_dir / case["damage_csv"]
    with uproot.open(root_path) as root:
        primary = _decode_strings(_tree_frame(root, "Ntuples/primary", [
            "EventID", "SeedID", "ParticleType", "InitialEnergy_eV",
            "InitialPosX_nm", "InitialPosY_nm", "InitialPosZ_nm",
            "InitialDirX", "InitialDirY", "InitialDirZ", "VoxelSize_nm",
            "Nx", "Ny", "Nz",
        ]))
        phase = _decode_strings(_tree_frame(root, "Ntuples/phase_space", [
            "EventID", "SeedID", "TrackID", "ParentID", "ParticleType",
            "PositionX_nm", "PositionY_nm", "PositionZ_nm",
            "LocalPositionX_nm", "LocalPositionY_nm", "LocalPositionZ_nm",
            "DirectionX", "DirectionY", "DirectionZ", "KineticEnergy_eV",
            "Time_ns", "VoxelID", "VoxelX", "VoxelY", "VoxelZ",
        ]))
        steps = _tree_frame(root, "Ntuples/steps", [
            "EventID", "SeedID", "TrackID", "StepNumber", "KEPost_eV",
            "Edep_eV", "StepLength_nm", "VoxelID",
        ])

    if phase.duplicated(["SeedID", "EventID", "TrackID", "VoxelID"]).any():
        raise ValueError(f"Repeated track/voxel entries in {root_path}")

    frame = phase.merge(
        primary.drop(columns=["ParticleType"]),
        on=["SeedID", "EventID"],
        how="left",
        validate="many_to_one",
    ).merge(
        _aggregate_steps(steps),
        on=["SeedID", "EventID", "TrackID", "VoxelID"],
        how="left",
        validate="one_to_one",
    )

    displacement = frame[["PositionX_nm", "PositionY_nm", "PositionZ_nm"]].to_numpy(copy=True)
    displacement -= frame[["InitialPosX_nm", "InitialPosY_nm", "InitialPosZ_nm"]].to_numpy()
    initial_direction = frame[["InitialDirX", "InitialDirY", "InitialDirZ"]].to_numpy()
    direction_norm = np.linalg.norm(initial_direction, axis=1, keepdims=True)
    initial_direction = initial_direction / np.maximum(direction_norm, 1e-12)
    longitudinal_nm = np.sum(displacement * initial_direction, axis=1)
    transverse = displacement - longitudinal_nm[:, None] * initial_direction

    frame["case_id"] = case["case_id"]
    frame["Particle"] = case.get("particle", frame["ParticleType"])
    frame["InitialEnergy_MeV"] = frame["InitialEnergy_eV"] / 1e6
    frame["Distance_um"] = longitudinal_nm / 1000.0
    frame["TransverseDisplacement_um"] = np.linalg.norm(transverse, axis=1) / 1000.0
    frame["EntryEnergy_MeV"] = frame["KineticEnergy_eV"] / 1e6
    frame["WorldEntryX_nm"] = frame["PositionX_nm"]
    frame["WorldEntryY_nm"] = frame["PositionY_nm"]
    frame["WorldEntryZ_nm"] = frame["PositionZ_nm"]
    frame["ExitEnergy_MeV"] = frame["ExitEnergy_eV"] / 1e6
    frame["Edep_keV"] = frame["Edep_eV"] / 1e3
    frame = frame.rename(columns={
        "DirectionX": "EntryDirX", "DirectionY": "EntryDirY", "DirectionZ": "EntryDirZ",
        "LocalPositionX_nm": "LocalEntryX_nm", "LocalPositionY_nm": "LocalEntryY_nm",
        "LocalPositionZ_nm": "LocalEntryZ_nm",
    })

    damage = pd.read_csv(damage_path).rename(columns={
        "upstream_seedID": "SeedID", "upstream_eventID": "EventID",
        "upstream_voxelID": "VoxelID", "dose_Gy": "Dose_Gy",
        "total_sb": "TotalSB", "total_ssb": "TotalSSB",
        "total_cssb": "TotalCSSB", "total_dsb": "TotalDSB",
        "total_cdsb": "TotalCDSB", "direct_sb": "DirectSB",
        "indirect_sb": "IndirectSB",
    })
    damage_columns = [
        "SeedID", "EventID", "VoxelID", "Dose_Gy", "TotalSB", "TotalSSB",
        "TotalCSSB", "TotalDSB", "TotalCDSB", "DirectSB", "IndirectSB",
    ]
    if damage.duplicated(["SeedID", "EventID", "VoxelID"]).any():
        raise ValueError(f"Damage rows are not unique in {damage_path}")
    frame = frame.merge(
        damage[damage_columns],
        on=["SeedID", "EventID", "VoxelID"],
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    if not frame["_merge"].eq("both").all():
        missing = int(frame["_merge"].ne("both").sum())
        raise ValueError(f"{missing} phase-space rows have no damage row in {case['case_id']}")
    frame = frame.drop(columns="_merge")
    frame["AnyDamage"] = (frame["TotalSB"] > 0).astype(float)
    frame["AnyDSB"] = (frame["TotalDSB"] > 0).astype(float)
    return frame


def add_dna_response(frame: pd.DataFrame, config: dict, input_dir: Path) -> pd.DataFrame:
    dna_frames = []
    for case in config["cases"]:
        with uproot.open(input_dir / case["dna_root"]) as root:
            dna = _tree_frame(root, "ntuple/EventEdep", [
                "Edep_J", "upstream_seedID", "upstream_eventID",
                "upstream_voxelID", "upstream_trackID",
            ])
        dna = dna.rename(columns={
            "upstream_seedID": "SeedID", "upstream_eventID": "EventID",
            "upstream_voxelID": "VoxelID", "upstream_trackID": "TrackID",
        })
        dna["case_id"] = case["case_id"]
        dna["DNAEdep_keV"] = dna["Edep_J"] / KEV_IN_JOULES
        if dna.duplicated(ENCOUNTER_KEY).any():
            raise ValueError(f"DNA event rows are not unique in {case['dna_root']}")
        dna_frames.append(dna[ENCOUNTER_KEY + ["DNAEdep_keV"]])

    dna_response = frame.merge(
        pd.concat(dna_frames, ignore_index=True),
        on=ENCOUNTER_KEY,
        how="left",
        validate="one_to_one",
        indicator=True,
    )
    if not dna_response["_merge"].eq("both").all():
        missing = int(dna_response["_merge"].ne("both").sum())
        raise ValueError(f"{missing} phase-space rows have no DNA EventEdep row")
    return dna_response.drop(columns="_merge")


def build_dataset(input_dir: Path, manifest: Path, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    config = json.loads(manifest.read_text())
    response = pd.concat(
        [build_case(case, input_dir) for case in config["cases"]],
        ignore_index=True,
    )
    response = _assign_splits(response)
    response = response.sort_values(ENCOUNTER_KEY).reset_index(drop=True)
    ps_response = add_dna_response(response, config, input_dir)

    range_frame = response.groupby(PRIMARY_KEY, as_index=False, observed=True).agg(
        Particle=("Particle", "first"),
        InitialEnergy_MeV=("InitialEnergy_MeV", "first"),
        InitialDirX=("InitialDirX", "first"),
        InitialDirY=("InitialDirY", "first"),
        InitialDirZ=("InitialDirZ", "first"),
        ObservedRange_um=("Distance_um", "max"),
        LastExitEnergy_MeV=("ExitEnergy_MeV", "last"),
        VoxelSize_nm=("VoxelSize_nm", "first"),
        Nz=("Nz", "first"),
        split=("split", "first"),
    )
    range_frame["GridLimited"] = (
        range_frame["ObservedRange_um"]
        >= range_frame["VoxelSize_nm"] * range_frame["Nz"] / 1000.0 - 0.31
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    response.to_csv(output_dir / "response.csv.gz", index=False)
    ps_response.to_csv(output_dir / "ps_response.csv.gz", index=False)
    range_frame.to_csv(output_dir / "range.csv.gz", index=False)
    summary = {
        "encounters": len(response),
        "primaries": len(range_frame),
        "cases": response.groupby("case_id").size().to_dict(),
        "splits": response.groupby("split").size().to_dict(),
        "damage_positive_fraction": float(response["AnyDamage"].mean()),
        "dsb_positive_fraction": float(response["AnyDSB"].mean()),
        "grid_limited_primaries": int(range_frame["GridLimited"].sum()),
        "ps_response_rows": len(ps_response),
        "dna_edep_keV": {
            "median": float(ps_response["DNAEdep_keV"].median()),
            "minimum": float(ps_response["DNAEdep_keV"].min()),
            "maximum": float(ps_response["DNAEdep_keV"].max()),
        },
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    return response, range_frame


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    package_root = Path(__file__).resolve().parents[2]
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=package_root / "configs/data_100.json")
    parser.add_argument("--output-dir", type=Path, default=package_root / "data/processed")
    args = parser.parse_args()
    response, ranges = build_dataset(args.input_dir, args.manifest, args.output_dir)
    print(f"Wrote {len(response):,} encounters and {len(ranges):,} primary ranges to {args.output_dir}")


if __name__ == "__main__":
    main()
