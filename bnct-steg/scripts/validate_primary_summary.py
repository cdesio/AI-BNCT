import argparse
import json
import os

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


REQUIRED_COLUMNS = [
    "RootFile",
    "Particle_ID",
    "ParticleType",
    "InitialEnergy_MeV",
    "SourceDistance_um",
    "Compartment",
    "ExpectedPrimaries",
    "PosX_um",
    "PosY_um",
    "PosZ_um",
    "StopPosX_um",
    "StopPosY_um",
    "StopPosZ_um",
    "TraLen_cell_um",
    "TraLen_chro_um",
    "StartDistToOrigin_um",
    "StopDistToOrigin_um",
    "TrackLength_um",
    "ClosestApproachT",
    "ClosestDistanceToOrigin_um",
    "DamageChunkCount",
    "RecordedEnergyDeposited_eV",
    "DamageRows",
    "SSB",
    "SSBp",
    "TwoSSB",
    "DSB",
    "DSBp",
    "DSBpp",
    "ClassificationMatchStatus",
    "ClassificationRowsInFile",
    "DamageChunksInFile",
    "ClassificationDamageChunkDelta",
    "HasDamage",
    "DamageSeverity",
]

NONNEGATIVE_COLUMNS = [
    "ExpectedPrimaries",
    "TraLen_cell_um",
    "TraLen_chro_um",
    "StartDistToOrigin_um",
    "StopDistToOrigin_um",
    "TrackLength_um",
    "ClosestDistanceToOrigin_um",
    "DamageChunkCount",
    "RecordedEnergyDeposited_eV",
    "DamageRows",
    "SSB",
    "SSBp",
    "TwoSSB",
    "DSB",
    "DSBp",
    "DSBpp",
]

FINITE_REQUIRED_COLUMNS = [
    "Particle_ID",
    "ExpectedPrimaries",
    "PosX_um",
    "PosY_um",
    "PosZ_um",
    "StopPosX_um",
    "StopPosY_um",
    "StopPosZ_um",
    "TraLen_cell_um",
    "TraLen_chro_um",
    "StartDistToOrigin_um",
    "StopDistToOrigin_um",
    "TrackLength_um",
    "ClosestApproachT",
    "ClosestDistanceToOrigin_um",
    "DamageChunkCount",
    "RecordedEnergyDeposited_eV",
    "DamageRows",
    "SSB",
    "SSBp",
    "TwoSSB",
    "DSB",
    "DSBp",
    "DSBpp",
    "HasDamage",
    "DamageSeverity",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Validate a BNCT primary-level summary CSV.")
    parser.add_argument("path", nargs="?", default="data/bnct_primary_summary.csv")
    parser.add_argument("--output-dir", default="data/primary_summary_validation")
    parser.add_argument("--max-tracklength-diff-um", type=float, default=1e-5)
    parser.add_argument("--max-radius-diff-um", type=float, default=1e-5)
    parser.add_argument("--max-closest-distance-diff-um", type=float, default=1e-5)
    parser.add_argument("--warn-assignment-distance-um", type=float, default=0.5)
    return parser.parse_args()


def add_issue(issues, severity, name, count, detail=""):
    if count:
        issues.append(
            {
                "severity": severity,
                "check": name,
                "count": int(count),
                "detail": detail,
            }
        )


def validate_required_columns(df, issues):
    missing = [column for column in REQUIRED_COLUMNS if column not in df.columns]
    add_issue(issues, "error", "missing_required_columns", len(missing), ", ".join(missing))


def validate_numeric_values(df, issues):
    for column in FINITE_REQUIRED_COLUMNS:
        if column in df.columns:
            add_issue(
                issues,
                "error",
                f"nonfinite_{column}",
                (~np.isfinite(df[column].to_numpy(float))).sum(),
            )

    for column in NONNEGATIVE_COLUMNS:
        if column in df.columns:
            add_issue(
                issues,
                "error",
                f"negative_{column}",
                (df[column] < 0).sum(),
            )


def validate_file_counts(df, issues):
    if "RootFile" not in df.columns or "ExpectedPrimaries" not in df.columns:
        return
    per_file = (
        df.groupby("RootFile", dropna=False)
        .agg(Rows=("Particle_ID", "size"), ExpectedPrimaries=("ExpectedPrimaries", "first"))
        .reset_index()
    )
    mismatched = per_file[per_file["Rows"] != per_file["ExpectedPrimaries"]]
    add_issue(
        issues,
        "error",
        "rows_do_not_match_expected_primaries",
        len(mismatched),
        "; ".join(
            f"{row.RootFile}: rows={row.Rows}, expected={row.ExpectedPrimaries}"
            for row in mismatched.itertuples(index=False)
        ),
    )


def validate_geometry(df, issues, args):
    starts = df[["PosX_um", "PosY_um", "PosZ_um"]].to_numpy(float)
    stops = df[["StopPosX_um", "StopPosY_um", "StopPosZ_um"]].to_numpy(float)
    segments = stops - starts
    segment_lengths = np.linalg.norm(segments, axis=1)
    start_radius = np.linalg.norm(starts, axis=1)
    stop_radius = np.linalg.norm(stops, axis=1)

    track_diff = np.abs(segment_lengths - df["TrackLength_um"].to_numpy(float))
    add_issue(
        issues,
        "error",
        "track_length_mismatch",
        (track_diff > args.max_tracklength_diff_um).sum(),
        f"max diff={track_diff.max():.3g} um",
    )

    start_diff = np.abs(start_radius - df["StartDistToOrigin_um"].to_numpy(float))
    stop_diff = np.abs(stop_radius - df["StopDistToOrigin_um"].to_numpy(float))
    add_issue(
        issues,
        "error",
        "start_radius_mismatch",
        (start_diff > args.max_radius_diff_um).sum(),
        f"max diff={start_diff.max():.3g} um",
    )
    add_issue(
        issues,
        "error",
        "stop_radius_mismatch",
        (stop_diff > args.max_radius_diff_um).sum(),
        f"max diff={stop_diff.max():.3g} um",
    )

    segment_lengths_sq = np.sum(segments**2, axis=1)
    t_closest = np.zeros(len(df), dtype=float)
    nonzero = segment_lengths_sq > 0
    t_closest[nonzero] = np.clip(
        -np.sum(starts[nonzero] * segments[nonzero], axis=1) / segment_lengths_sq[nonzero],
        0.0,
        1.0,
    )
    closest = starts + t_closest[:, None] * segments
    closest_distance = np.linalg.norm(closest, axis=1)
    closest_diff = np.abs(closest_distance - df["ClosestDistanceToOrigin_um"].to_numpy(float))
    add_issue(
        issues,
        "error",
        "closest_distance_mismatch",
        (closest_diff > args.max_closest_distance_diff_um).sum(),
        f"max diff={closest_diff.max():.3g} um",
    )

    add_issue(
        issues,
        "error",
        "closest_approach_t_outside_unit_interval",
        ((df["ClosestApproachT"] < 0) | (df["ClosestApproachT"] > 1)).sum(),
    )


def validate_labels(df, issues):
    has_damage_expected = (df["DamageChunkCount"] > 0).astype(int)
    add_issue(
        issues,
        "error",
        "has_damage_inconsistent_with_damage_chunks",
        (df["HasDamage"].astype(int) != has_damage_expected).sum(),
    )

    severe = ((df["DSB"] > 0) | (df["DSBp"] > 0) | (df["DSBpp"] > 0)).astype(int)
    mild = ((df["HasDamage"] > 0) & (severe == 0)).astype(int)
    severity_expected = severe * 2 + mild
    add_issue(
        issues,
        "error",
        "damage_severity_inconsistent_with_counts",
        (df["DamageSeverity"].astype(int) != severity_expected).sum(),
    )


def validate_assignment_qa(df, issues, args):
    damaged = df[df["HasDamage"] > 0].copy()
    if damaged.empty:
        add_issue(issues, "error", "no_damaged_primaries", 1)
        return

    if "ChunkMedianDistanceToTrack_um" in damaged.columns:
        large = damaged["ChunkMedianDistanceToTrack_um"] > args.warn_assignment_distance_um
        add_issue(
            issues,
            "warning",
            "large_chunk_median_distance_to_track",
            large.sum(),
            f"threshold={args.warn_assignment_distance_um} um",
        )

    if "ClassificationDamageChunkDelta" in df.columns:
        add_issue(
            issues,
            "warning",
            "classification_damage_chunk_delta_nonzero",
            (df["ClassificationDamageChunkDelta"] != 0).sum(),
        )


def write_plots(df, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(df["ClosestDistanceToOrigin_um"], bins=50, alpha=0.8)
    ax.set_xlabel("Closest distance to origin (um)")
    ax.set_ylabel("Primaries")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "closest_distance_hist.pdf"))
    plt.close(fig)

    bins = np.linspace(0, np.nanpercentile(df["ClosestDistanceToOrigin_um"], 99), 10)
    radial = df.copy()
    radial["RadiusBin"] = pd.cut(radial["ClosestDistanceToOrigin_um"], bins=bins, include_lowest=True)
    radial_summary = (
        radial.groupby("RadiusBin", observed=True)
        .agg(
            Primaries=("HasDamage", "size"),
            DamageFraction=("HasDamage", "mean"),
            RadiusMid_um=("ClosestDistanceToOrigin_um", "mean"),
        )
        .reset_index()
    )
    radial_summary.to_csv(os.path.join(output_dir, "observed_radial_damage.csv"), index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(radial_summary["RadiusMid_um"], radial_summary["DamageFraction"], marker="o")
    ax.set_xlabel("Closest distance to origin (um)")
    ax.set_ylabel("Observed damage fraction")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, "observed_radial_damage.pdf"))
    plt.close(fig)


def main():
    args = parse_args()
    df = pd.read_csv(args.path)
    issues = []

    validate_required_columns(df, issues)
    if any(issue["severity"] == "error" and issue["check"] == "missing_required_columns" for issue in issues):
        result = {"rows": int(len(df)), "issues": issues, "passed": False}
    else:
        validate_numeric_values(df, issues)
        validate_file_counts(df, issues)
        validate_geometry(df, issues, args)
        validate_labels(df, issues)
        validate_assignment_qa(df, issues, args)
        write_plots(df, args.output_dir)
        result = {
            "rows": int(len(df)),
            "files": int(df["RootFile"].nunique()),
            "damaged_primaries": int(df["HasDamage"].sum()),
            "damage_fraction": float(df["HasDamage"].mean()),
            "issues": issues,
            "passed": not any(issue["severity"] == "error" for issue in issues),
        }

    os.makedirs(args.output_dir, exist_ok=True)
    with open(os.path.join(args.output_dir, "validation_report.json"), "w") as f:
        json.dump(result, f, indent=4)

    print(json.dumps(result, indent=4))
    print(f"Output: {args.output_dir}")

    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
