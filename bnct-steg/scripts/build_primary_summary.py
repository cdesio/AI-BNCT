import argparse
import glob
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd
import uproot


PRIMARY_COLUMNS = [
    "Primary",
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
]

CLASSIFICATION_COLUMNS = [
    "None",
    "SSB",
    "SSBp",
    "2SSB",
    "DSB",
    "DSBp",
    "DSBpp",
]

DAMAGE_COLUMNS = [
    "Event",
    "Position_x_um",
    "Position_y_um",
    "Position_z_um",
    "EnergyDeposited_eV",
    "BaseDamage",
    "StrandDamage",
    "DirectBreaks",
    "IndirectBreaks",
    "TypeClassificationInt",
    "SourceClassificationInt",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build a primary-level BNCT dataset from MolecularBNCT ROOT files."
    )
    parser.add_argument(
        "--input-glob",
        default="/Users/yw18581/work/AI-BNCT/ForChiara-OldData/*/*.root",
    )
    parser.add_argument("--output", default="data/bnct_primary_summary.csv")
    parser.add_argument("--max-files", type=int)
    parser.add_argument("--write-per-file", action="store_true")
    parser.add_argument("--per-file-dir", default="data/primary_summary_by_file")
    parser.add_argument(
        "--allow-classification-mismatch",
        action="store_true",
        help=(
            "Keep files where classification rows and reconstructed damage chunks differ. "
            "This is now the default; the flag is kept for backward-compatible commands."
        ),
    )
    parser.add_argument(
        "--strict-classification-match",
        action="store_true",
        help="Skip files where classification rows and reconstructed damage chunks differ.",
    )
    parser.add_argument(
        "--strict-primary-count",
        action="store_true",
        help="Skip files where the filename primary count does not match primary_source rows.",
    )
    return parser.parse_args()


def parse_metadata_from_filename(path):
    stem = Path(path).stem
    match = re.match(
        r"(?P<particle>alpha|lithium\+{3})(?P<energy>\d+p\d+)_"
        r"(?P<source_distance>\d+)um_(?P<compartment>[A-Za-z]+)_"
        r"(?P<n_primaries>\d+)$",
        stem,
    )
    if not match:
        raise ValueError(f"Could not parse BNCT metadata from filename: {path}")

    particle_raw = match.group("particle")
    particle = "lithium" if particle_raw.startswith("lithium") else particle_raw
    return {
        "ParticleType": particle,
        "InitialEnergy_MeV": float(match.group("energy").replace("p", ".")),
        "SourceDistance_um": float(match.group("source_distance")),
        "Compartment": match.group("compartment"),
        "ExpectedPrimaries": int(match.group("n_primaries")),
    }


def tree_to_frame(root_file, tree_name, columns):
    tree = root_file[f"Ntuples/{tree_name}"]
    available = set(tree.keys())
    missing = [column for column in columns if column not in available]
    if missing:
        raise KeyError(f"Ntuples/{tree_name} is missing required branches: {missing}")
    return pd.DataFrame(tree.arrays(columns, library="np"))


def add_track_geometry(primary_df):
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

    directions = np.zeros_like(segments)
    directions[nonzero] = segments[nonzero] / segment_lengths[nonzero, None]

    start_unit = np.zeros_like(starts)
    nonzero_start = start_radius > 0
    start_unit[nonzero_start] = starts[nonzero_start] / start_radius[nonzero_start, None]
    radial_direction_cosine = np.sum(directions * start_unit, axis=1)

    primary_df = primary_df.copy()
    primary_df["StartDistToOrigin_um"] = start_radius
    primary_df["StopDistToOrigin_um"] = stop_radius
    primary_df["DeltaRadius_um"] = stop_radius - start_radius
    primary_df["TrackLength_um"] = segment_lengths
    primary_df["DirectionX"] = directions[:, 0]
    primary_df["DirectionY"] = directions[:, 1]
    primary_df["DirectionZ"] = directions[:, 2]
    primary_df["RadialDirectionCosine"] = radial_direction_cosine
    primary_df["ClosestApproachT"] = t_closest
    primary_df["ClosestApproachX_um"] = closest_points[:, 0]
    primary_df["ClosestApproachY_um"] = closest_points[:, 1]
    primary_df["ClosestApproachZ_um"] = closest_points[:, 2]
    primary_df["ClosestDistanceToOrigin_um"] = closest_distance
    return primary_df


def add_damage_chunk_id(damage_df):
    if len(damage_df) == 0:
        damage_df = damage_df.copy()
        damage_df["DamageChunk_ID"] = np.array([], dtype=int)
        return damage_df

    local_event = damage_df["Event"].astype(int)
    chunk_start = local_event.ne(local_event.shift(fill_value=local_event.iloc[0]))
    damage_df = damage_df.copy()
    damage_df["DamageChunk_ID"] = chunk_start.cumsum().astype(int)
    return damage_df


def summarize_damage_chunks(damage_with_chunks):
    if len(damage_with_chunks) == 0:
        return pd.DataFrame(
            columns=[
                "DamageChunk_ID",
                "Event",
                "DamageRows",
                "RecordedEnergyDeposited_eV",
                "MaxStepEnergyDeposited_eV",
                "MeanStepEnergyDeposited_eV",
                "SumBaseDamage",
                "SumStrandDamage",
                "SumDirectBreaks",
                "SumIndirectBreaks",
                "MeanDamageStepX_um",
                "MeanDamageStepY_um",
                "MeanDamageStepZ_um",
                "MinDamageStepRadius_um",
                "MeanDamageStepRadius_um",
                "MaxDamageStepRadius_um",
            ]
        )

    chunk_summary = (
        damage_with_chunks.groupby("DamageChunk_ID", as_index=False, sort=True)
        .agg(
            Event=("Event", "first"),
            DamageRows=("Event", "size"),
            RecordedEnergyDeposited_eV=("EnergyDeposited_eV", "sum"),
            MaxStepEnergyDeposited_eV=("EnergyDeposited_eV", "max"),
            MeanStepEnergyDeposited_eV=("EnergyDeposited_eV", "mean"),
            SumBaseDamage=("BaseDamage", "sum"),
            SumStrandDamage=("StrandDamage", "sum"),
            SumDirectBreaks=("DirectBreaks", "sum"),
            SumIndirectBreaks=("IndirectBreaks", "sum"),
            MeanDamageStepX_um=("Position_x_um", "mean"),
            MeanDamageStepY_um=("Position_y_um", "mean"),
            MeanDamageStepZ_um=("Position_z_um", "mean"),
        )
    )

    positions = damage_with_chunks[["Position_x_um", "Position_y_um", "Position_z_um"]].to_numpy(float)
    radii = np.linalg.norm(positions, axis=1)
    weighted = damage_with_chunks.copy()
    weighted["DamageStepRadius_um"] = radii

    radius_summary = (
        weighted.groupby("DamageChunk_ID", as_index=False, sort=True)
        .agg(
            MinDamageStepRadius_um=("DamageStepRadius_um", "min"),
            MeanDamageStepRadius_um=("DamageStepRadius_um", "mean"),
            MaxDamageStepRadius_um=("DamageStepRadius_um", "max"),
        )
    )
    return chunk_summary.merge(radius_summary, on="DamageChunk_ID", how="left")


def distance_points_to_segment(points, start, segment, length_sq):
    point_vectors = points - start
    if length_sq == 0:
        t = np.zeros(len(points), dtype=float)
        closest = np.repeat(start[None, :], len(points), axis=0)
    else:
        t = np.clip(point_vectors @ segment / length_sq, 0.0, 1.0)
        closest = start + t[:, None] * segment
    distances = np.linalg.norm(points - closest, axis=1)
    return distances


def assign_damage_chunks_to_primary_tracks(primary_df, damage_with_chunks):
    if len(damage_with_chunks) == 0:
        return pd.DataFrame(
            columns=[
                "DamageChunk_ID",
                "Particle_ID",
                "ChunkMedianDistanceToTrack_um",
                "ChunkMinDistanceToTrack_um",
            ]
        )

    starts = primary_df[["PosX_um", "PosY_um", "PosZ_um"]].to_numpy(float)
    stops = primary_df[["StopPosX_um", "StopPosY_um", "StopPosZ_um"]].to_numpy(float)
    segments = stops - starts
    segment_lengths_sq = np.sum(segments**2, axis=1)

    assignments = []
    for damage_chunk_id, chunk in damage_with_chunks.groupby("DamageChunk_ID", sort=True):
        points = chunk[["Position_x_um", "Position_y_um", "Position_z_um"]].to_numpy(float)
        best_particle_id = -1
        best_score = np.inf
        best_min_distance = np.inf

        for particle_id, (start, segment, length_sq) in enumerate(
            zip(starts, segments, segment_lengths_sq)
        ):
            distances = distance_points_to_segment(points, start, segment, length_sq)
            score = float(np.median(distances))
            min_distance = float(np.min(distances))
            if score < best_score:
                best_score = score
                best_min_distance = min_distance
                best_particle_id = particle_id

        assignments.append(
            {
                "DamageChunk_ID": damage_chunk_id,
                "Particle_ID": best_particle_id,
                "ChunkMedianDistanceToTrack_um": best_score,
                "ChunkMinDistanceToTrack_um": best_min_distance,
            }
        )

    return pd.DataFrame(assignments)


def severity_from_counts(row):
    if row["HasDamage"] == 0:
        return 0
    if row["DSB"] > 0 or row["DSBp"] > 0 or row["DSBpp"] > 0:
        return 2
    return 1


def build_primary_summary_for_file(
    path,
    allow_classification_mismatch=True,
    strict_primary_count=False,
):
    metadata = parse_metadata_from_filename(path)
    root_file = uproot.open(path)

    primary_df = tree_to_frame(root_file, "primary_source", PRIMARY_COLUMNS)
    classification_df = tree_to_frame(root_file, "classification", CLASSIFICATION_COLUMNS)
    damage_df = tree_to_frame(root_file, "damage", DAMAGE_COLUMNS)

    primary_count_matches = len(primary_df) == metadata["ExpectedPrimaries"]
    if strict_primary_count and not primary_count_matches:
        raise ValueError(
            f"{path}: filename says {metadata['ExpectedPrimaries']} primaries, "
            f"but primary_source has {len(primary_df)} rows"
        )

    primary_df = primary_df.reset_index(drop=True)
    primary_df.insert(0, "Particle_ID", np.arange(len(primary_df), dtype=int))
    primary_df.insert(0, "RootFile", os.path.basename(path))
    for key, value in metadata.items():
        primary_df[key] = value

    primary_df = add_track_geometry(primary_df)

    damage_with_chunks = add_damage_chunk_id(damage_df)
    chunk_summary = summarize_damage_chunks(damage_with_chunks)
    chunk_assignment = assign_damage_chunks_to_primary_tracks(primary_df, damage_with_chunks)
    chunk_table = chunk_summary.merge(chunk_assignment, on="DamageChunk_ID", how="left")

    classification_rows = len(classification_df)
    damage_chunks = len(chunk_table)
    chunk_table["ClassificationMatchStatus"] = "matched_by_order"

    if classification_rows == damage_chunks:
        classification_df = classification_df.reset_index(drop=True).copy()
        classification_df["DamageChunk_ID"] = chunk_table["DamageChunk_ID"].to_numpy()
        chunk_table = chunk_table.merge(classification_df, on="DamageChunk_ID", how="left")
    else:
        if not allow_classification_mismatch:
            raise ValueError(
                f"classification rows ({classification_rows}) do not match "
                f"damage chunks ({damage_chunks})"
            )
        # Some files have fewer consecutive damage chunks than classification rows.
        # In that case, keep damage-energy aggregates from the damage chunks and
        # attach classification totals to the nearest damaged primary IDs by order.
        chunk_table["ClassificationMatchStatus"] = "count_mismatch_damage_only"
        for column in CLASSIFICATION_COLUMNS:
            chunk_table[column] = 0.0

    per_primary = (
        chunk_table.groupby("Particle_ID", as_index=False, sort=True)
        .agg(
            DamageChunkCount=("DamageChunk_ID", "size"),
            RecordedEnergyDeposited_eV=("RecordedEnergyDeposited_eV", "sum"),
            MaxStepEnergyDeposited_eV=("MaxStepEnergyDeposited_eV", "max"),
            MeanStepEnergyDeposited_eV=("MeanStepEnergyDeposited_eV", "mean"),
            DamageRows=("DamageRows", "sum"),
            SumBaseDamage=("SumBaseDamage", "sum"),
            SumStrandDamage=("SumStrandDamage", "sum"),
            SumDirectBreaks=("SumDirectBreaks", "sum"),
            SumIndirectBreaks=("SumIndirectBreaks", "sum"),
            MinDamageStepRadius_um=("MinDamageStepRadius_um", "min"),
            MeanDamageStepRadius_um=("MeanDamageStepRadius_um", "mean"),
            MaxDamageStepRadius_um=("MaxDamageStepRadius_um", "max"),
            ChunkMedianDistanceToTrack_um=("ChunkMedianDistanceToTrack_um", "median"),
            ChunkMinDistanceToTrack_um=("ChunkMinDistanceToTrack_um", "min"),
            NoneCount=("None", "sum"),
            SSB=("SSB", "sum"),
            SSBp=("SSBp", "sum"),
            TwoSSB=("2SSB", "sum"),
            DSB=("DSB", "sum"),
            DSBp=("DSBp", "sum"),
            DSBpp=("DSBpp", "sum"),
            ClassificationMatchStatus=("ClassificationMatchStatus", "first"),
        )
    )

    if classification_rows != damage_chunks and classification_rows > 0:
        damaged_primary_ids = (
            chunk_table["Particle_ID"].dropna().astype(int).drop_duplicates().sort_values().to_numpy()
        )
        n_assign = min(classification_rows, len(damaged_primary_ids))
        if n_assign:
            class_by_primary = classification_df.iloc[:n_assign].copy()
            class_by_primary["Particle_ID"] = damaged_primary_ids[:n_assign]
            class_by_primary["ClassificationMatchStatus"] = "classification_by_primary_order_fallback"
            per_primary = per_primary.drop(
                columns=[*CLASSIFICATION_COLUMNS, "ClassificationMatchStatus"],
                errors="ignore",
            )
            per_primary = per_primary.merge(
                class_by_primary[[*CLASSIFICATION_COLUMNS, "Particle_ID", "ClassificationMatchStatus"]],
                on="Particle_ID",
                how="left",
            )
            for column in CLASSIFICATION_COLUMNS:
                per_primary[column] = per_primary[column].fillna(0.0)
            per_primary["ClassificationMatchStatus"] = per_primary[
                "ClassificationMatchStatus"
            ].fillna("damage_only")

    output = primary_df.merge(per_primary, on="Particle_ID", how="left")
    zero_fill = [
        "DamageChunkCount",
        "RecordedEnergyDeposited_eV",
        "MaxStepEnergyDeposited_eV",
        "MeanStepEnergyDeposited_eV",
        "DamageRows",
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
    ]
    output[zero_fill] = output[zero_fill].fillna(0.0)
    for column in zero_fill:
        output[column] = pd.to_numeric(output[column], errors="coerce").fillna(0.0)
    output["ClassificationMatchStatus"] = output["ClassificationMatchStatus"].fillna("no_damage")
    output["ClassificationRowsInFile"] = classification_rows
    output["DamageChunksInFile"] = damage_chunks
    output["ClassificationDamageChunkDelta"] = classification_rows - damage_chunks
    output["PrimaryRowsInFile"] = len(primary_df)
    output["PrimaryCountMatchesFilename"] = primary_count_matches
    output["HasDamage"] = (output["DamageChunkCount"] > 0).astype(int)
    output["HasDamageKnown"] = 1
    output["DamageLabelQuality"] = np.where(
        output["ClassificationDamageChunkDelta"] == 0,
        "exact_damage_classification_match",
        "damage_tree_geometry_assignment",
    )
    output.loc[output["DamageChunkCount"] == 0, "DamageLabelQuality"] = "no_damage_recorded"
    output["DamageSeverity"] = output.apply(severity_from_counts, axis=1).astype(int)
    output["LogRecordedEnergyDeposited_eV"] = np.log1p(output["RecordedEnergyDeposited_eV"])
    return output


def main():
    args = parse_args()
    paths = sorted(glob.glob(args.input_glob))
    if args.max_files is not None:
        paths = paths[: args.max_files]
    if not paths:
        raise FileNotFoundError(f"No ROOT files matched {args.input_glob!r}")

    frames = []
    for index, path in enumerate(paths, start=1):
        try:
            summary = build_primary_summary_for_file(
                path,
                allow_classification_mismatch=(
                    args.allow_classification_mismatch
                    or not args.strict_classification_match
                ),
                strict_primary_count=args.strict_primary_count,
            )
        except Exception as exc:
            print(f"skip {path}: {exc}")
            continue

        frames.append(summary)
        print(
            f"{index}/{len(paths)} {os.path.basename(path)}: "
            f"{len(summary)} primaries, {int(summary['HasDamage'].sum())} damaged"
        )

        if args.write_per_file:
            os.makedirs(args.per_file_dir, exist_ok=True)
            per_file_path = Path(args.per_file_dir) / f"{Path(path).stem}_primary_summary.csv"
            summary.to_csv(per_file_path, index=False)

    if not frames:
        raise RuntimeError("No primary summaries were built.")

    combined = pd.concat(frames, ignore_index=True)
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    combined.to_csv(args.output, index=False)

    print(f"Files used: {len(frames)}")
    print(f"Rows written: {len(combined)}")
    print(f"Damaged primaries: {int(combined['HasDamage'].sum())}")
    print(f"Output: {args.output}")


if __name__ == "__main__":
    main()
