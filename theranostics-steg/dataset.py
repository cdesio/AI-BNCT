import glob
import pickle

import numpy as np
import pandas as pd
import torch
from torch.utils.data import IterableDataset


FILTER_COLUMNS = ["particleName", "boundary", "volumeName"]


def load_theranostics_photons(config, start=0, stop=None):
    paths = sorted(glob.glob(config.DATA_GLOB))
    if config.MAX_FILES is not None:
        paths = paths[: config.MAX_FILES]
    if not paths:
        raise FileNotFoundError(f"No CSV files matched DATA_GLOB={config.DATA_GLOB!r}")

    usecols = list(dict.fromkeys([*FILTER_COLUMNS, *config.FEATURE_NAMES]))
    frames = []
    rows_needed = None if stop is None else stop

    for path in paths:
        df = pd.read_csv(path, usecols=usecols)
        mask = (
            (df["particleName"] == "gamma")
            & (df["boundary"] == "enter")
            & df["volumeName"].str.startswith("physPhotonDetector", na=False)

        )
        df = df.loc[mask, config.FEATURE_NAMES]
        frames.append(df)

        if rows_needed is not None and sum(len(frame) for frame in frames) >= rows_needed:
            break

    data = pd.concat(frames, ignore_index=True)
    data = data.replace([np.inf, -np.inf], np.nan).dropna()
    if stop is not None:
        data = data.iloc[:stop]
    if start:
        data = data.iloc[start:]
    return data.reset_index(drop=True)


def resolve_event_splits(config):
    data = load_theranostics_photons(config, stop=config.TOTAL_EVENTS)
    total_events = len(data)
    if total_events == 0:
        raise ValueError("No rows remain after applying the theranostics photon filters.")

    train_split = int(0.7 * total_events)
    val_split = int(0.85 * total_events)
    if train_split == 0 or val_split <= train_split or total_events <= val_split:
        raise ValueError(
            "Not enough filtered rows to create train, validation, and test splits. "
            f"Found {total_events} rows."
        )

    return total_events, train_split, val_split


class TheranosticsPhotonDataset(IterableDataset):
    def __init__(self, config, start=0, stop=None, input_scaler_flag=False):
        self.config = config
        self.features = config.FEATURE_NAMES
        self.step = config.STEP_SIZE
        self.data = load_theranostics_photons(config, start=start, stop=stop)

        self.input_scaler_flag = input_scaler_flag
        if self.input_scaler_flag:
            with open(config.TRANSFORMER_LOC, "rb") as f:
                self.quantile_transformers = pickle.load(f)

    def __iter__(self):
        values = self.data[self.features].to_numpy(dtype=np.float32)

        if self.input_scaler_flag:
            values = values.copy()
            for i, feature in enumerate(self.features):
                transformer = self.quantile_transformers[feature]
                values[:, i] = transformer.transform(values[:, i].reshape(-1, 1)).ravel()

        for start in range(0, len(values), self.step):
            yield torch.tensor(values[start : start + self.step], dtype=torch.float32)
