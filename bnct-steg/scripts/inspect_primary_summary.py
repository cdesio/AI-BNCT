import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description="Inspect a BNCT primary-level summary CSV.")
    parser.add_argument("path", nargs="?", default="data/bnct_primary_summary.csv")
    return parser.parse_args()


def main():
    args = parse_args()
    df = pd.read_csv(args.path)

    print(f"Rows: {len(df)}")
    print(f"Files: {df['RootFile'].nunique()}")
    print(f"Damaged primaries: {int(df['HasDamage'].sum())} ({df['HasDamage'].mean():.4f})")

    group_cols = ["ParticleType", "InitialEnergy_MeV", "SourceDistance_um", "Compartment"]
    grouped = (
        df.groupby(group_cols, dropna=False)
        .agg(
            Primaries=("HasDamage", "size"),
            DamageFraction=("HasDamage", "mean"),
            MeanClosestDistance_um=("ClosestDistanceToOrigin_um", "mean"),
            MeanRecordedEdep_eV=("RecordedEnergyDeposited_eV", "mean"),
        )
        .reset_index()
        .sort_values(group_cols)
    )
    print("\nDamage fraction by condition:")
    print(grouped.to_string(index=False))

    qa_cols = [
        "ChunkMedianDistanceToTrack_um",
        "ChunkMinDistanceToTrack_um",
        "RecordedEnergyDeposited_eV",
        "ClosestDistanceToOrigin_um",
        "TrackLength_um",
    ]
    print("\nNumeric QA summary:")
    print(df[qa_cols].describe().to_string())


if __name__ == "__main__":
    main()
