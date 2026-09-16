import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from config import make_config
from dataset import load_configured_data, resolve_event_splits


def main():
    config = make_config()
    df = load_configured_data(config, start=0, stop=config.TOTAL_EVENTS)
    total_events, train_split, val_split = resolve_event_splits(config)

    print(f"DATA_MODE: {config.DATA_MODE}")
    print(f"Rows: {len(df)}")
    print(
        "Split sizes: "
        f"train={train_split}, "
        f"validation={val_split - train_split}, "
        f"test={total_events - val_split}"
    )
    for feature in config.BINARY_FEATURE_NAMES:
        print(f"{feature}: positive={int(df[feature].sum())}, fraction={df[feature].mean():.4f}")
    print(df.describe().to_string())


if __name__ == "__main__":
    main()
