"""Train grouped direct-distance range, transport, and damage baselines."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import ExtraTreesClassifier, ExtraTreesRegressor, GradientBoostingRegressor
from sklearn.metrics import (
    average_precision_score,
    mean_absolute_error,
    mean_squared_error,
    roc_auc_score,
)
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .schema import (
    CONDITION_COLUMNS,
    DAMAGE_COUNT_TARGETS,
    DAMAGE_INPUT_COLUMNS,
    TRANSPORT_TARGETS,
)


RANGE_FEATURES = [
    "Particle", "InitialEnergy_MeV", "InitialDirX", "InitialDirY", "InitialDirZ"
]

LOG_TRANSPORT_TARGETS = {
    "EntryEnergy_MeV", "ExitEnergy_MeV", "Edep_keV", "TrackLength_nm",
    "StepCount", "MeanLET_keV_um", "TransverseDisplacement_um",
}


def transform_transport(values: pd.DataFrame | np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float).copy()
    for index, name in enumerate(TRANSPORT_TARGETS):
        if name in LOG_TRANSPORT_TARGETS:
            array[:, index] = np.log1p(np.clip(array[:, index], 0, None))
    return array


def inverse_transport(values: np.ndarray) -> np.ndarray:
    array = np.asarray(values, dtype=float).copy()
    for index, name in enumerate(TRANSPORT_TARGETS):
        if name in LOG_TRANSPORT_TARGETS:
            array[:, index] = np.expm1(array[:, index]).clip(min=0)
    return array


def make_preprocessor(columns: list[str]) -> ColumnTransformer:
    categorical = [column for column in columns if column == "Particle"]
    numeric = [column for column in columns if column not in categorical]
    return ColumnTransformer([
        ("particle", OneHotEncoder(handle_unknown="ignore", sparse_output=False), categorical),
        ("numeric", StandardScaler(), numeric),
    ])


def make_regressor(features: list[str], seed: int) -> Pipeline:
    return Pipeline([
        ("features", make_preprocessor(features)),
        ("model", ExtraTreesRegressor(
            n_estimators=300, min_samples_leaf=3, max_features=0.9,
            n_jobs=-1, random_state=seed,
        )),
    ])


def make_classifier(features: list[str], seed: int) -> Pipeline:
    return Pipeline([
        ("features", make_preprocessor(features)),
        ("model", ExtraTreesClassifier(
            n_estimators=300, min_samples_leaf=4, class_weight="balanced",
            max_features=0.9, n_jobs=-1, random_state=seed,
        )),
    ])


def _binary_metrics(y_true: np.ndarray, probability: np.ndarray) -> dict[str, float]:
    result = {"average_precision": float(average_precision_score(y_true, probability))}
    if len(np.unique(y_true)) == 2:
        result["roc_auc"] = float(roc_auc_score(y_true, probability))
    return result


def _regression_metrics(y_true: np.ndarray, prediction: np.ndarray, names: list[str]) -> dict:
    return {
        name: {
            "mae": float(mean_absolute_error(y_true[:, index], prediction[:, index])),
            "rmse": float(mean_squared_error(y_true[:, index], prediction[:, index]) ** 0.5),
        }
        for index, name in enumerate(names)
    }


def train_models(data_dir: Path, output_dir: Path, seed: int = 20260923) -> dict:
    response = pd.read_csv(data_dir / "response.csv.gz")
    ranges = pd.read_csv(data_dir / "range.csv.gz")
    train = response[response["split"] == "train"].copy()
    test = response[response["split"] == "test"].copy()
    range_train = ranges[ranges["split"] == "train"].copy()
    range_test = ranges[ranges["split"] == "test"].copy()

    range_models = {}
    range_metrics = {}
    for quantile in (0.05, 0.5, 0.95):
        model = Pipeline([
            ("features", make_preprocessor(RANGE_FEATURES)),
            ("model", GradientBoostingRegressor(
                loss="quantile", alpha=quantile, n_estimators=150,
                min_samples_leaf=4, random_state=seed,
            )),
        ])
        model.fit(range_train[RANGE_FEATURES], range_train["ObservedRange_um"])
        range_models[f"q{int(quantile * 100):02d}"] = model
        prediction = model.predict(range_test[RANGE_FEATURES])
        range_metrics[f"q{int(quantile * 100):02d}_mae"] = float(
            mean_absolute_error(range_test["ObservedRange_um"], prediction)
        )

    transport = make_regressor(CONDITION_COLUMNS, seed)
    transport.fit(train[CONDITION_COLUMNS], transform_transport(train[TRANSPORT_TARGETS]))
    transport_prediction = inverse_transport(transport.predict(test[CONDITION_COLUMNS]))

    damage_any = make_classifier(DAMAGE_INPUT_COLUMNS, seed)
    damage_dsb = make_classifier(DAMAGE_INPUT_COLUMNS, seed + 1)
    damage_counts = make_regressor(DAMAGE_INPUT_COLUMNS, seed + 2)
    direct_damage_any = make_classifier(CONDITION_COLUMNS, seed + 3)
    direct_damage_dsb = make_classifier(CONDITION_COLUMNS, seed + 4)
    direct_damage_counts = make_regressor(CONDITION_COLUMNS, seed + 5)
    damage_any.fit(train[DAMAGE_INPUT_COLUMNS], train["AnyDamage"])
    damage_dsb.fit(train[DAMAGE_INPUT_COLUMNS], train["AnyDSB"])
    damage_counts.fit(train[DAMAGE_INPUT_COLUMNS], np.log1p(train[DAMAGE_COUNT_TARGETS]))
    direct_damage_any.fit(train[CONDITION_COLUMNS], train["AnyDamage"])
    direct_damage_dsb.fit(train[CONDITION_COLUMNS], train["AnyDSB"])
    direct_damage_counts.fit(
        train[CONDITION_COLUMNS], np.log1p(train[DAMAGE_COUNT_TARGETS])
    )

    any_probability = damage_any.predict_proba(test[DAMAGE_INPUT_COLUMNS])[:, 1]
    dsb_probability = damage_dsb.predict_proba(test[DAMAGE_INPUT_COLUMNS])[:, 1]
    count_prediction = np.expm1(damage_counts.predict(test[DAMAGE_INPUT_COLUMNS])).clip(min=0)
    direct_any_probability = direct_damage_any.predict_proba(test[CONDITION_COLUMNS])[:, 1]
    direct_dsb_probability = direct_damage_dsb.predict_proba(test[CONDITION_COLUMNS])[:, 1]
    direct_count_prediction = np.expm1(
        direct_damage_counts.predict(test[CONDITION_COLUMNS])
    ).clip(min=0)

    chained = test[CONDITION_COLUMNS].copy()
    chained_transport = transport_prediction
    for index, target in enumerate(TRANSPORT_TARGETS):
        chained[target] = chained_transport[:, index]
    chained_any_probability = damage_any.predict_proba(chained[DAMAGE_INPUT_COLUMNS])[:, 1]
    chained_dsb_probability = damage_dsb.predict_proba(chained[DAMAGE_INPUT_COLUMNS])[:, 1]

    metrics = {
        "rows": {"train": len(train), "test": len(test)},
        "primaries": {
            "train": int(train[["case_id", "SeedID", "EventID"]].drop_duplicates().shape[0]),
            "test": int(test[["case_id", "SeedID", "EventID"]].drop_duplicates().shape[0]),
        },
        "range": range_metrics,
        "transport": _regression_metrics(
            test[TRANSPORT_TARGETS].to_numpy(), transport_prediction, TRANSPORT_TARGETS
        ),
        "damage_oracle_transport": {
            "any_damage": _binary_metrics(test["AnyDamage"].to_numpy(), any_probability),
            "any_dsb": _binary_metrics(test["AnyDSB"].to_numpy(), dsb_probability),
            "counts": _regression_metrics(
                test[DAMAGE_COUNT_TARGETS].to_numpy(), count_prediction, DAMAGE_COUNT_TARGETS
            ),
        },
        "damage_direct_distance": {
            "any_damage": _binary_metrics(
                test["AnyDamage"].to_numpy(), direct_any_probability
            ),
            "any_dsb": _binary_metrics(test["AnyDSB"].to_numpy(), direct_dsb_probability),
            "counts": _regression_metrics(
                test[DAMAGE_COUNT_TARGETS].to_numpy(),
                direct_count_prediction,
                DAMAGE_COUNT_TARGETS,
            ),
        },
        "damage_direct_distance_chained": {
            "any_damage": _binary_metrics(
                test["AnyDamage"].to_numpy(), chained_any_probability
            ),
            "any_dsb": _binary_metrics(test["AnyDSB"].to_numpy(), chained_dsb_probability),
        },
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    artifact = {
        "range_models": range_models,
        "transport": transport,
        "damage_any": damage_any,
        "damage_dsb": damage_dsb,
        "damage_counts": damage_counts,
        "direct_damage_any": direct_damage_any,
        "direct_damage_dsb": direct_damage_dsb,
        "direct_damage_counts": direct_damage_counts,
        "condition_columns": CONDITION_COLUMNS,
        "transport_targets": TRANSPORT_TARGETS,
        "damage_input_columns": DAMAGE_INPUT_COLUMNS,
        "damage_count_targets": DAMAGE_COUNT_TARGETS,
        "seed": seed,
    }
    joblib.dump(artifact, output_dir / "models.joblib")
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=20260923)
    args = parser.parse_args()
    metrics = train_models(args.data_dir, args.output_dir, args.seed)
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
