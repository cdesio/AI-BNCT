import argparse
import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import uproot


PRIMARY_COLUMNS = [
    "Energy",
    "PosX_um",
    "PosY_um",
    "PosZ_um",
    "MomX",
    "MomY",
    "MomZ",
    "StopPosX_um",
    "StopPosY_um",
    "StopPosZ_um",
    "TraLen_cell_um",
    "TraLen_chro_um",
    "EventID",
    "distance",
    "SeedID",
]

TRACK_COLUMNS = [
    "EventID",
    "TrackID",
    "StepNumber",
    "KineticEnergy_eV",
    "EnergyDeposited_eV",
    "StepLength_um",
    "PrePosX_um",
    "PrePosY_um",
    "PrePosZ_um",
    "PostPosX_um",
    "PostPosY_um",
    "PostPosZ_um",
    "MomX",
    "MomY",
    "MomZ",
    "SeedID",
]

DAMAGE_COLUMNS = [
    "EventID",
    "EnergyDeposited_eV",
    "BaseDamage",
    "StrandDamage",
    "DirectBreaks",
    "IndirectBreaks",
    "Position_x_um",
    "Position_y_um",
    "Position_z_um",
    "TypeClassificationInt",
    "SourceClassificationInt",
    "SeedID",
]

CLASSIFICATION_COLUMNS = [
    "EventID",
    "None",
    "SSB",
    "SSBp",
    "2SSB",
    "DSB",
    "DSBp",
    "DSBpp",
    "SeedID",
]

CHROMOSOME_HIT_COLUMNS = [
    "EventID",
    "e_chromosome_kev",
    "e_dna_kev",
    "SeedID",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build primary- and step-level CSV datasets from MolecularBNCT ROOT files "
            "that contain Ntuples/primary_tracks."
        )
    )
    parser.add_argument(
        "--input-glob",
        default="/Users/yw18581/work/MolecularBNCT-main/bnct_campaign/*.root",
    )
    parser.add_argument(
        "--primary-output",
        default="data/bnct_tracked_primary_summary.csv",
    )
    parser.add_argument(
        "--step-output",
        default="data/bnct_primary_track_steps.csv",
    )
    parser.add_argument("--max-files", type=int)
    return parser.parse_args()


def parse_metadata_from_filename(path):
    stem = Path(path).stem
    metadata = {
        "ParticleType": pd.NA,
        "InitialEnergy_MeV": np.nan,
        "SourceDistance_um": np.nan,
        "Compartment": pd.NA,
        "ExpectedPrimaries": np.nan,
    }

    old_match = re.match(
        r"(?P<particle>alpha|lithium\+{3})(?P<energy>\d+p\d+)_"
        r"(?P<source_distance>\d+)um_(?P<compartment>[A-Za-z]+)_"
        r"(?P<n_primaries>\d+)$",
        stem,
    )
    if old_match:
        particle_raw = old_match.group("particle")
        metadata.update(
            {
                "ParticleType": "lithium" if particle_raw.startswith("lithium") else particle_raw,
                "InitialEnergy_MeV": float(old_match.group("energy").replace("p", ".")),
                "SourceDistance_um": float(old_match.group("source_distance")),
                "Compartment": old_match.group("compartment"),
                "ExpectedPrimaries": int(old_match.group("n_primaries")),
            }
        )
        return metadata

    new_match = re.match(
        r"(?P<particle>[A-Za-z]+)_(?P<energy>\d+p\d+)_(?P<seed>\d+)_t(?P<time>\d+)$",
        stem,
    )
    if new_match:
        metadata.update(
            {
                "ParticleType": new_match.group("particle"),
                "InitialEnergy_MeV": float(new_match.group("energy").replace("p", ".")),
                "ExpectedPrimaries": 1,
            }
        )
    return metadata


def tree_to_frame(root_file, tree_name, columns, required=True):
    key = f"Ntuples/{tree_name}"
    if key not in root_file:
        if required:
            raise KeyError(f"{key} is missing")
        return pd.DataFrame(columns=columns)

    tree = root_file[key]
    available = set(tree.keys())
    use_columns = [column for column in columns if column in available]
    missing = [column for column in columns if column not in available]
    if missing and required:
        raise KeyError(f"{key} is missing required branches: {missing}")
    if not use_columns:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(tree.arrays(use_columns, library="np"))


def add_primary_geometry(primary_df):
    starts = primary_df[["PosX_um", "PosY_um", "PosZ_um"]].to_numpy(float)
    stops = primary_df[["StopPosX_um", "StopPosY_um", "StopPosZ_um"]].to_numpy(float)
    segments = stops - starts
    segment_lengths_sq = np.sum(segments**2, axis=1)
    segment_lengths = np.sqrt(segment_lengths_sq)

    start_radius = np.linalg.norm(starts, axis=1)
    stop_radius = np.linalg.norm(stops, axis=1)

    t_closest = np.zeros(len(primary_df), dtype=float)
    nonzero = segment_lengths_sq > 0
    t_closest[nonzero] = np.clip(
        -np.sum(starts[nonzero] * segments[nonzero], axis=1) / segment_lengths_sq[nonzero],
        0.0,
        1.0,
    )
    closest_points = starts + t_closest[:, None] * segments
    closest_distance = np.linalg.norm(closest_points, axis=1)

    primary_df = primary_df.copy()
    primary_df["StartDistToOrigin_um"] = start_radius
    primary_df["StopDistToOrigin_um"] = stop_radius
    primary_df["DeltaRadius_um"] = stop_radius - start_radius
    primary_df["TrackLength_um"] = segment_lengths
    primary_df["ClosestApproachT"] = t_closest
    primary_df["ClosestApproachX_um"] = closest_points[:, 0]
    primary_df["ClosestApproachY_um"] = closest_points[:, 1]
    primary_df["ClosestApproachZ_um"] = closest_points[:, 2]
    primary_df["ClosestDistanceToOrigin_um"] = closest_distance
    return primary_df


def add_step_geometry(track_df):
    pre = track_df[["PrePosX_um", "PrePosY_um", "PrePosZ_um"]].to_numpy(float)
    post = track_df[["PostPosX_um", "PostPosY_um", "PostPosZ_um"]].to_numpy(float)
    mid = 0.5 * (pre + post)
    track_df = track_df.copy()
    track_df["StepMidX_um"] = mid[:, 0]
    track_df["StepMidY_um"] = mid[:, 1]
    track_df["StepMidZ_um"] = mid[:, 2]
    track_df["StepPreRadius_um"] = np.linalg.norm(pre, axis=1)
    track_df["StepPostRadius_um"] = np.linalg.norm(post, axis=1)
    track_df["StepMidRadius_um"] = np.linalg.norm(mid, axis=1)
    return track_df


def summarize_tracks(track_df):
    if len(track_df) == 0:
        return pd.DataFrame()
    return (
        track_df.groupby(["EventID", "SeedID"], as_index=False)
        .agg(
            NumTrackSteps=("StepNumber", "size"),
            NumTrackIDs=("TrackID", "nunique"),
            TotalTrackStepLength_um=("StepLength_um", "sum"),
            TotalTrackEnergyDeposited_eV=("EnergyDeposited_eV", "sum"),
            MaxTrackStepEnergyDeposited_eV=("EnergyDeposited_eV", "max"),
            MeanTrackStepEnergyDeposited_eV=("EnergyDeposited_eV", "mean"),
            MinStepMidRadius_um=("StepMidRadius_um", "min"),
            MeanStepMidRadius_um=("StepMidRadius_um", "mean"),
            MaxStepMidRadius_um=("StepMidRadius_um", "max"),
            FirstKineticEnergy_eV=("KineticEnergy_eV", "first"),
            LastKineticEnergy_eV=("KineticEnergy_eV", "last"),
        )
    )


def summarize_damage(damage_df):
    if len(damage_df) == 0:
        return pd.DataFrame(columns=["EventID", "SeedID"])

    positions = damage_df[["Position_x_um", "Position_y_um", "Position_z_um"]].to_numpy(float)
    damage_df = damage_df.copy()
    damage_df["DamageRadius_um"] = np.linalg.norm(positions, axis=1)
    return (
        damage_df.groupby(["EventID", "SeedID"], as_index=False)
        .agg(
            DamageRows=("EnergyDeposited_eV", "size"),
            RecordedEnergyDeposited_eV=("EnergyDeposited_eV", "sum"),
            MaxDamageEnergyDeposited_eV=("EnergyDeposited_eV", "max"),
            SumBaseDamage=("BaseDamage", "sum"),
            SumStrandDamage=("StrandDamage", "sum"),
            SumDirectBreaks=("DirectBreaks", "sum"),
            SumIndirectBreaks=("IndirectBreaks", "sum"),
            MinDamageRadius_um=("DamageRadius_um", "min"),
            MeanDamageRadius_um=("DamageRadius_um", "mean"),
            MaxDamageRadius_um=("DamageRadius_um", "max"),
        )
    )


def build_for_file(path):
    metadata = parse_metadata_from_filename(path)
    root_file = uproot.open(path)

    primary_df = tree_to_frame(root_file, "primary_source", PRIMARY_COLUMNS)
    track_df = tree_to_frame(root_file, "primary_tracks", TRACK_COLUMNS)
    damage_df = tree_to_frame(root_file, "damage", DAMAGE_COLUMNS, required=False)
    class_df = tree_to_frame(root_file, "classification", CLASSIFICATION_COLUMNS, required=False)
    hit_df = tree_to_frame(root_file, "chromosome_hits", CHROMOSOME_HIT_COLUMNS, required=False)

    root_name = os.path.basename(path)
    primary_df = primary_df.reset_index(drop=True)
    primary_df.insert(0, "PrimaryRow_ID", np.arange(len(primary_df), dtype=int))
    primary_df.insert(0, "RootFile", root_name)
    track_df.insert(0, "RootFile", root_name)
    for frame in (primary_df, track_df):
        for key, value in metadata.items():
            frame[key] = value

    primary_df = add_primary_geometry(primary_df)
    track_df = add_step_geometry(track_df)

    track_summary = summarize_tracks(track_df)
    damage_summary = summarize_damage(damage_df)
    class_df = class_df.rename(columns={"None": "NoneCount", "2SSB": "TwoSSB"})

    primary_summary = primary_df.merge(track_summary, on=["EventID", "SeedID"], how="left")
    primary_summary = primary_summary.merge(damage_summary, on=["EventID", "SeedID"], how="left")
    if len(class_df):
        primary_summary = primary_summary.merge(class_df, on=["EventID", "SeedID"], how="left")
    if len(hit_df):
        primary_summary = primary_summary.merge(hit_df, on=["EventID", "SeedID"], how="left")

    zero_fill = [
        "NumTrackSteps",
        "NumTrackIDs",
        "TotalTrackStepLength_um",
        "TotalTrackEnergyDeposited_eV",
        "DamageRows",
        "RecordedEnergyDeposited_eV",
        "SumBaseDamage",
        "SumStrandDamage",
        "SumDirectBreaks",
        "SumIndirectBreaks",
        "NoneCount",
        "SSB",
        "SSBp",
        "TwoSSB",
        "DSB",
        "DSBp",
        "DSBpp",
        "e_chromosome_kev",
        "e_dna_kev",
    ]
    for column in zero_fill:
        if column in primary_summary:
            primary_summary[column] = primary_summary[column].fillna(0.0)

    primary_summary["HasDamage"] = (
        primary_summary[["SSB", "SSBp", "TwoSSB", "DSB", "DSBp", "DSBpp"]]
        .fillna(0.0)
        .sum(axis=1)
        > 0
    ).astype(int)
    primary_summary["DamageSeverity"] = np.select(
        [
            primary_summary[["DSB", "DSBp", "DSBpp"]].fillna(0.0).sum(axis=1) > 0,
            primary_summary["HasDamage"] > 0,
        ],
        [2, 1],
        default=0,
    )
    primary_summary["LogTotalTrackEnergyDeposited_eV"] = np.log1p(
        primary_summary["TotalTrackEnergyDeposited_eV"].fillna(0.0)
    )
    primary_summary["LogRecordedEnergyDeposited_eV"] = np.log1p(
        primary_summary["RecordedEnergyDeposited_eV"].fillna(0.0)
    )
    primary_summary["TrackingAvailable"] = 1
    primary_summary["DamageJoinMethod"] = "event_seed"

    return primary_summary, track_df


def main():
    args = parse_args()
    paths = sorted(glob.glob(args.input_glob))
    if args.max_files is not None:
        paths = paths[: args.max_files]
    if not paths:
        raise FileNotFoundError(f"No ROOT files matched {args.input_glob!r}")

    primary_frames = []
    step_frames = []
    for index, path in enumerate(paths, start=1):
        try:
            primary_summary, step_df = build_for_file(path)
        except Exception as exc:
            print(f"skip {path}: {exc}")
            continue
        primary_frames.append(primary_summary)
        step_frames.append(step_df)
        print(
            f"{index}/{len(paths)} {os.path.basename(path)}: "
            f"{len(primary_summary)} primaries, {len(step_df)} track steps, "
            f"{int(primary_summary['HasDamage'].sum())} damaged"
        )

    if not primary_frames:
        raise RuntimeError("No tracked ROOT datasets were built.")

    primary_out = pd.concat(primary_frames, ignore_index=True)
    step_out = pd.concat(step_frames, ignore_index=True)
    os.makedirs(os.path.dirname(args.primary_output), exist_ok=True)
    os.makedirs(os.path.dirname(args.step_output), exist_ok=True)
    primary_out.to_csv(args.primary_output, index=False)
    step_out.to_csv(args.step_output, index=False)

    print(f"Primary rows written: {len(primary_out)}")
    print(f"Step rows written: {len(step_out)}")
    print(f"Primary output: {args.primary_output}")
    print(f"Step output: {args.step_output}")


if __name__ == "__main__":
    main()
