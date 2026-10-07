"""Train a two-step phase-space-entry to DNA-damage baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    f1_score,
    mean_absolute_error,
    mean_squared_error,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.linear_model import LogisticRegression

from .baseline import make_classifier, make_regressor
from .schema import PS_HISTORY_INPUT_COLUMNS, PS_INPUT_COLUMNS, PS_RESPONSE_COUNT_TARGETS


def _regression_metrics(truth: np.ndarray, prediction: np.ndarray) -> dict:
    return {
        name: {
            "mae": float(mean_absolute_error(truth[:, index], prediction[:, index])),
            "rmse": float(mean_squared_error(truth[:, index], prediction[:, index]) ** 0.5),
        }
        for index, name in enumerate(PS_RESPONSE_COUNT_TARGETS)
    }


def _f1_threshold(truth: np.ndarray, probability: np.ndarray) -> float:
    thresholds = np.linspace(0.05, 0.95, 181)
    scores = np.asarray([f1_score(truth, probability >= threshold) for threshold in thresholds])
    return float(thresholds[int(np.argmax(scores))])


def _predict_counts(model, frame: pd.DataFrame, input_columns: list[str]) -> np.ndarray:
    return np.expm1(model.predict(frame[input_columns])).clip(min=0)


def _probability_logit(probability: np.ndarray) -> np.ndarray:
    clipped = np.clip(probability, 1e-6, 1 - 1e-6)
    return np.log(clipped / (1 - clipped)).reshape(-1, 1)


def calibrated_probability(calibrator: LogisticRegression, probability: np.ndarray) -> np.ndarray:
    return calibrator.predict_proba(_probability_logit(probability))[:, 1]


def train_ps_two_step(
    data_dir: Path, output_dir: Path, input_mode: str = "local", seed: int = 20261007
) -> dict:
    frame = pd.read_csv(data_dir / "ps_response.csv.gz")
    train = frame[frame["split"] == "train"].copy()
    validation = frame[frame["split"] == "val"].copy()
    test = frame[frame["split"] == "test"].copy()
    damaged_train = train[train["AnyDamage"] == 1].copy()
    input_columns = PS_INPUT_COLUMNS if input_mode == "local" else PS_HISTORY_INPUT_COLUMNS

    edep = make_regressor(input_columns, seed)
    damage_gate = make_classifier(input_columns, seed + 1)
    positive_counts = make_regressor(input_columns, seed + 2)

    edep.fit(train[input_columns], np.log1p(train["DNAEdep_keV"]))
    damage_gate.fit(train[input_columns], train["AnyDamage"])
    positive_counts.fit(
        damaged_train[input_columns], np.log1p(damaged_train[PS_RESPONSE_COUNT_TARGETS])
    )

    validation_raw_probability = damage_gate.predict_proba(validation[input_columns])[:, 1]
    calibrator = LogisticRegression(random_state=seed + 3)
    calibrator.fit(
        _probability_logit(validation_raw_probability), validation["AnyDamage"].to_numpy()
    )
    validation_probability = calibrated_probability(calibrator, validation_raw_probability)
    threshold = _f1_threshold(validation["AnyDamage"].to_numpy(), validation_probability)

    edep_prediction = np.expm1(edep.predict(test[input_columns])).clip(min=0)
    raw_damage_probability = damage_gate.predict_proba(test[input_columns])[:, 1]
    damage_probability = calibrated_probability(calibrator, raw_damage_probability)
    damage_decision = damage_probability >= threshold
    conditional_counts = _predict_counts(positive_counts, test, input_columns)
    gated_counts = conditional_counts.copy()
    gated_counts[~damage_decision] = 0

    truth_damage = test["AnyDamage"].to_numpy(dtype=int)
    truth_counts = test[PS_RESPONSE_COUNT_TARGETS].to_numpy()
    positive_test = truth_damage == 1
    metrics = {
        "model": "two_step_extra_trees",
        "input_mode": input_mode,
        "input_columns": input_columns,
        "rows": {
            "train": len(train),
            "damaged_train": len(damaged_train),
            "validation": len(validation),
            "test": len(test),
        },
        "primaries": {
            split: int(part[["case_id", "SeedID", "EventID"]].drop_duplicates().shape[0])
            for split, part in [("train", train), ("validation", validation), ("test", test)]
        },
        "dna_edep": {
            "mae": float(mean_absolute_error(test["DNAEdep_keV"], edep_prediction)),
            "rmse": float(mean_squared_error(test["DNAEdep_keV"], edep_prediction) ** 0.5),
        },
        "damage_gate": {
            "threshold": threshold,
            "roc_auc": float(roc_auc_score(truth_damage, damage_probability)),
            "average_precision": float(average_precision_score(truth_damage, damage_probability)),
            "brier_score": float(brier_score_loss(truth_damage, damage_probability)),
            "f1": float(f1_score(truth_damage, damage_decision)),
            "precision": float(precision_score(truth_damage, damage_decision)),
            "recall": float(recall_score(truth_damage, damage_decision)),
            "true_positive_fraction": float(truth_damage.mean()),
            "predicted_positive_fraction": float(damage_decision.mean()),
        },
        "damage_counts": _regression_metrics(truth_counts, gated_counts),
        "conditional_counts_on_damaged_test": _regression_metrics(
            truth_counts[positive_test], conditional_counts[positive_test]
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    joblib.dump({
        "edep": edep,
        "damage_gate": damage_gate,
        "damage_calibrator": calibrator,
        "positive_counts": positive_counts,
        "damage_threshold": threshold,
        "input_columns": input_columns,
        "input_mode": input_mode,
        "count_targets": PS_RESPONSE_COUNT_TARGETS,
        "seed": seed,
    }, output_dir / "models.joblib")
    (output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    identifiers = [
        "case_id", "SeedID", "EventID", "TrackID", "VoxelID", "Particle",
        "EntryEnergy_MeV", "Distance_um", "DNAEdep_keV", "AnyDamage",
    ] + PS_RESPONSE_COUNT_TARGETS
    predictions = test[identifiers].reset_index(drop=True)
    predictions["PredictedDNAEdep_keV"] = edep_prediction
    predictions["RawP_AnyDamage"] = raw_damage_probability
    predictions["P_AnyDamage"] = damage_probability
    predictions["PredictedAnyDamage"] = damage_decision.astype(int)
    for index, target in enumerate(PS_RESPONSE_COUNT_TARGETS):
        predictions[f"Conditional{target}"] = conditional_counts[:, index]
        predictions[f"Predicted{target}"] = gated_counts[:, index]
    predictions.to_csv(output_dir / "test_predictions.csv.gz", index=False)
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--input-mode", choices=["local", "history"], default="local")
    parser.add_argument("--seed", type=int, default=20261007)
    args = parser.parse_args()
    print(json.dumps(
        train_ps_two_step(args.data_dir, args.output_dir, args.input_mode, args.seed), indent=2
    ))


if __name__ == "__main__":
    main()
