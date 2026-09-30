"""Apply the phase-space baseline to every row in a voxel-transport ROOT file."""

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import uproot


def load_phase_space(path: Path) -> pd.DataFrame:
    columns = [
        "EventID", "SeedID", "TrackID", "ParentID", "ParticleType",
        "KineticEnergy_eV", "DirectionX", "DirectionY", "DirectionZ",
        "LocalPositionX_nm", "LocalPositionY_nm", "LocalPositionZ_nm",
        "PositionX_nm", "PositionY_nm", "PositionZ_nm",
        "VoxelID", "VoxelX", "VoxelY", "VoxelZ",
    ]
    with uproot.open(path) as root:
        arrays = root["Ntuples/phase_space"].arrays(columns, library="np")
        primary_columns = [
            "EventID", "SeedID", "InitialEnergy_eV",
            "InitialPosX_nm", "InitialPosY_nm", "InitialPosZ_nm",
            "InitialDirX", "InitialDirY", "InitialDirZ",
        ]
        primary_arrays = root["Ntuples/primary"].arrays(primary_columns, library="np")
    frame = pd.DataFrame({column: arrays[column] for column in columns})
    primary = pd.DataFrame({column: primary_arrays[column] for column in primary_columns})
    frame = frame.merge(primary, on=["SeedID", "EventID"], how="left", validate="many_to_one")
    frame["Particle"] = frame["ParticleType"].map(
        lambda value: value.decode() if isinstance(value, bytes) else str(value)
    ).replace({"lithium+++": "lithium"})
    frame["EntryEnergy_MeV"] = frame["KineticEnergy_eV"] / 1e6
    frame["InitialEnergy_MeV"] = frame["InitialEnergy_eV"] / 1e6
    displacement = frame[["PositionX_nm", "PositionY_nm", "PositionZ_nm"]].to_numpy(copy=True)
    displacement -= frame[["InitialPosX_nm", "InitialPosY_nm", "InitialPosZ_nm"]].to_numpy()
    direction = frame[["InitialDirX", "InitialDirY", "InitialDirZ"]].to_numpy(copy=True)
    direction /= np.maximum(np.linalg.norm(direction, axis=1, keepdims=True), 1e-12)
    frame["Distance_um"] = np.sum(displacement * direction, axis=1) / 1000.0
    return frame.rename(columns={
        "DirectionX": "EntryDirX", "DirectionY": "EntryDirY", "DirectionZ": "EntryDirZ",
        "LocalPositionX_nm": "LocalEntryX_nm", "LocalPositionY_nm": "LocalEntryY_nm",
        "LocalPositionZ_nm": "LocalEntryZ_nm",
        "PositionX_nm": "WorldEntryX_nm", "PositionY_nm": "WorldEntryY_nm",
        "PositionZ_nm": "WorldEntryZ_nm",
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    artifact = joblib.load(args.model_dir / "models.joblib")
    frame = load_phase_space(args.input_root)
    inputs = frame[artifact["input_columns"]]
    frame["PredictedDNAEdep_keV"] = np.expm1(
        artifact["edep"].predict(inputs)
    ).clip(min=0).ravel()
    frame["P_AnyDamage"] = artifact["any_damage"].predict_proba(inputs)[:, 1]
    frame["P_AnyDSB"] = artifact["any_dsb"].predict_proba(inputs)[:, 1]
    predicted_counts = np.expm1(artifact["counts"].predict(inputs)).clip(min=0)
    for index, target in enumerate(artifact["count_targets"]):
        frame[f"Expected{target}"] = predicted_counts[:, index]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output, index=False)
    print(f"Wrote {len(frame):,} predicted voxel responses to {args.output}")


if __name__ == "__main__":
    main()
