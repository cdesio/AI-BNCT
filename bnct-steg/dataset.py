import glob
import pickle

import numpy as np
import pandas as pd
import torch
from torch.utils.data import IterableDataset

from config import BNCT_DAMAGE_INDICATOR_FEATURES, BNCT_RAW_DAMAGE_COLUMNS


def add_damage_indicators(df):
    damage_positive = df[BNCT_RAW_DAMAGE_COLUMNS] > 0
    df = df.copy()
    df["Damage"] = damage_positive.any(axis=1).astype(float)
    for damage_column in BNCT_RAW_DAMAGE_COLUMNS:
        df[f"{damage_column}_Positive"] = damage_positive[damage_column].astype(float)
    return df


def load_bnct_csv_damage(config, start=0, stop=None):
    paths = sorted(glob.glob(config.DATA_GLOB))
    if config.MAX_FILES is not None:
        paths = paths[: config.MAX_FILES]
    if not paths:
        raise FileNotFoundError(f"No CSV files matched DATA_GLOB={config.DATA_GLOB!r}")

    derived_features = set(BNCT_DAMAGE_INDICATOR_FEATURES)
    read_columns = [
        *[feature for feature in config.FEATURE_NAMES if feature not in derived_features],
        *BNCT_RAW_DAMAGE_COLUMNS,
    ]
    read_columns = list(dict.fromkeys(read_columns))

    frames = []
    rows_needed = None if stop is None else stop

    for path in paths:
        df = pd.read_csv(path, usecols=read_columns)
        df = add_damage_indicators(df)
        df = df[config.FEATURE_NAMES]
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


def load_primary_summary(config, start=0, stop=None):
    df = pd.read_csv(config.PRIMARY_SUMMARY_PATH)
    missing = [feature for feature in config.FEATURE_NAMES if feature not in df.columns]
    if missing:
        raise KeyError(
            f"{config.PRIMARY_SUMMARY_PATH} is missing configured features: {missing}"
        )

    df = df[config.FEATURE_NAMES].copy()
    for binary_feature in getattr(config, "BINARY_FEATURE_NAMES", []):
        if binary_feature in df.columns:
            df[binary_feature] = (df[binary_feature].astype(float) > 0.5).astype(float)

    data = df.replace([np.inf, -np.inf], np.nan).dropna()
    if stop is not None:
        data = data.iloc[:stop]
    if start:
        data = data.iloc[start:]
    return data.reset_index(drop=True)


def load_configured_data(config, start=0, stop=None):
    data_mode = getattr(config, "DATA_MODE", "csv_damage")
    if data_mode == "csv_damage":
        return load_bnct_csv_damage(config, start=start, stop=stop)
    if data_mode == "primary_summary":
        return load_primary_summary(config, start=start, stop=stop)
    raise ValueError(f"Unsupported DATA_MODE={data_mode!r}")


def resolve_event_splits(config):
    data = load_configured_data(config, stop=config.TOTAL_EVENTS)
    total_events = len(data)
    if total_events == 0:
        raise ValueError("No BNCT rows remain after finite-value filtering.")

    train_split = int(0.7 * total_events)
    val_split = int(0.85 * total_events)
    if train_split == 0 or val_split <= train_split or total_events <= val_split:
        raise ValueError(
            "Not enough BNCT CSV rows to create train, validation, and test splits. "
            f"Found {total_events} rows."
        )

    return total_events, train_split, val_split


class BNCTCsvDamageDataset(IterableDataset):
    def __init__(self, config, start=0, stop=None, input_scaler_flag=False):
        self.config = config
        self.features = config.FEATURE_NAMES
        self.step = config.STEP_SIZE
        self.data = load_configured_data(config, start=start, stop=stop)

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
