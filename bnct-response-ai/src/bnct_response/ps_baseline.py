"""Train the phase-space-entry to DNA-voxel-response baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, mean_absolute_error, mean_squared_error, roc_auc_score

from .baseline import make_classifier, make_regressor
from .schema import PS_HISTORY_INPUT_COLUMNS, PS_INPUT_COLUMNS, PS_RESPONSE_COUNT_TARGETS


def _binary_metrics(truth: pd.Series, probability: np.ndarray) -> dict[str, float]:
    return {
        "roc_auc": float(roc_auc_score(truth, probability)),
        "average_precision": float(average_precision_score(truth, probability)),
    }


def _regression_metrics(truth: np.ndarray, prediction: np.ndarray, names: list[str]) -> dict:
    return {
        name: {
            "mae": float(mean_absolute_error(truth[:, index], prediction[:, index])),
            "rmse": float(mean_squared_error(truth[:, index], prediction[:, index]) ** 0.5),
        }
        for index, name in enumerate(names)
    }


def train_ps_baseline(
    data_dir: Path, output_dir: Path, input_mode: str = "local", seed: int = 20260930
) -> dict:
    frame = pd.read_csv(data_dir / "ps_response.csv.gz")
    train = frame[frame["split"] == "train"].copy()
    test = frame[frame["split"] == "test"].copy()
    input_columns = PS_INPUT_COLUMNS if input_mode == "local" else PS_HISTORY_INPUT_COLUMNS

    edep = make_regressor(input_columns, seed)
    any_damage = make_classifier(input_columns, seed + 1)
    any_dsb = make_classifier(input_columns, seed + 2)
    counts = make_regressor(input_columns, seed + 3)

    edep.fit(train[input_columns], np.log1p(train["DNAEdep_keV"]))
    any_damage.fit(train[input_columns], train["AnyDamage"])
    any_dsb.fit(train[input_columns], train["AnyDSB"])
    counts.fit(train[input_columns], np.log1p(train[PS_RESPONSE_COUNT_TARGETS]))

    edep_prediction = np.expm1(edep.predict(test[input_columns])).clip(min=0).reshape(-1, 1)
    damage_probability = any_damage.predict_proba(test[input_columns])[:, 1]
    dsb_probability = any_dsb.predict_proba(test[input_columns])[:, 1]
    count_prediction = np.expm1(counts.predict(test[input_columns])).clip(min=0)

    metrics = {
        "rows": {"train": len(train), "test": len(test)},
        "input_mode": input_mode,
        "input_columns": input_columns,
        "primaries": {
            "train": int(train[["case_id", "SeedID", "EventID"]].drop_duplicates().shape[0]),
            "test": int(test[["case_id", "SeedID", "EventID"]].drop_duplicates().shape[0]),
        },
        "dna_edep": _regression_metrics(
            test[["DNAEdep_keV"]].to_numpy(), edep_prediction, ["DNAEdep_keV"]
        ),
        "any_damage": _binary_metrics(test["AnyDamage"], damage_probability),
        "any_dsb": _binary_metrics(test["AnyDSB"], dsb_probability),
        "damage_counts": _regression_metrics(
            test[PS_RESPONSE_COUNT_TARGETS].to_numpy(),
            count_prediction,
            PS_RESPONSE_COUNT_TARGETS,
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump({
        "edep": edep,
        "any_damage": any_damage,
        "any_dsb": any_dsb,
        "counts": counts,
        "input_columns": input_columns,
        "input_mode": input_mode,
        "count_targets": PS_RESPONSE_COUNT_TARGETS,
        "seed": seed,
    }, output_dir / "models.joblib")
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input-mode", choices=["local", "history"], default="local")
    parser.add_argument("--seed", type=int, default=20260930)
    args = parser.parse_args()
    print(json.dumps(
        train_ps_baseline(args.data_dir, args.output_dir, args.input_mode, args.seed), indent=2
    ))


if __name__ == "__main__":
    main()
