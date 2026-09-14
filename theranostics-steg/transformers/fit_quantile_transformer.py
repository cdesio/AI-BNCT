import os
import pickle
import sys

import numpy as np
from sklearn.preprocessing import QuantileTransformer
from tqdm import tqdm

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
from config import Config  # noqa: E402
from dataset import load_theranostics_photons  # noqa: E402


def main():
    config = Config()
    df = load_theranostics_photons(config, start=0, stop=config.TOTAL_EVENTS)
    total_events = len(df)
    train_split = int(0.7 * total_events)
    val_split = int(0.85 * total_events)
    print(f"Loaded {total_events} filtered photon rows")
    print(
        "Split sizes: "
        f"train={train_split}, "
        f"validation={val_split - train_split}, "
        f"test={total_events - val_split}"
    )
    if len(df) == 0:
        raise RuntimeError("No rows found after theranostics photon filters.")

    quantile_transformers = {}
    for feature in tqdm(config.FEATURE_NAMES, desc="Fitting transformers"):
        data = df[feature].to_numpy().reshape(-1, 1)
        if not np.isfinite(data).all():
            raise RuntimeError(f"Non-finite values found in {feature}")

        transformer = QuantileTransformer(
            n_quantiles=min(10000, len(data)),
            output_distribution="normal",
            random_state=42,
            subsample=None,
        )
        transformer.fit(data)
        quantile_transformers[feature] = transformer

    os.makedirs(os.path.dirname(config.TRANSFORMER_LOC), exist_ok=True)
    with open(config.TRANSFORMER_LOC, "wb") as f:
        pickle.dump(quantile_transformers, f)

    print(f"Saved {config.TRANSFORMER_LOC}")


if __name__ == "__main__":
    main()
