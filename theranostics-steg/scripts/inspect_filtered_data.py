import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from config import Config
from dataset import load_theranostics_photons


def main():
    config = Config()
    df = load_theranostics_photons(config, start=0, stop=config.TOTAL_EVENTS)
    print(f"Rows: {len(df)}")
    print(df.describe().to_string())


if __name__ == "__main__":
    main()
