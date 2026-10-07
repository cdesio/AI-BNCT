"""Evaluate calibrated-gate plus positive-damage StEG on held-out primaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from .ps_two_step_generate import generate_two_step
from .schema import PS_GENERATIVE_COUNT_TARGETS


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--gate-model-dir", type=Path, required=True)
    parser.add_argument("--steg-model-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=20261007)
    args = parser.parse_args()

    test = pd.read_csv(args.data_dir / "ps_response.csv.gz")
    test = test[test["split"] == "test"].reset_index(drop=True)
    generated = generate_two_step(
        test, args.gate_model_dir, args.steg_model_dir, 1,
        args.device, args.batch_size, args.seed,
    )
    truth_damage = test["AnyDamage"].to_numpy(dtype=int)
    probability = generated["GeneratedP_AnyDamage"].to_numpy()
    metrics = {
        "rows": len(test),
        "damage_gate": {
            "roc_auc": float(roc_auc_score(truth_damage, probability)),
            "average_precision": float(average_precision_score(truth_damage, probability)),
            "brier_score": float(brier_score_loss(truth_damage, probability)),
            "true_positive_fraction": float(truth_damage.mean()),
            "sampled_positive_fraction": float(generated["GeneratedAnyDamage"].mean()),
        },
        "outputs": {},
    }
    for target in PS_GENERATIVE_COUNT_TARGETS:
        truth = test[target].to_numpy()
        sample = generated[f"Generated{target}"].to_numpy()
        metrics["outputs"][target] = {
            "true_mean": float(truth.mean()),
            "generated_mean": float(sample.mean()),
            "true_zero_fraction": float((truth == 0).mean()),
            "generated_zero_fraction": float((sample == 0).mean()),
            "sample_mae": float(np.mean(np.abs(truth - sample))),
            "raw_generated_mean": float(generated[f"RawGenerated{target}"].mean()),
            "raw_generated_zero_fraction": float(
                (generated[f"RawGenerated{target}"] == 0).mean()
            ),
        }

    keep = [
        "case_id", "SeedID", "EventID", "TrackID", "VoxelID", "Particle",
        "EntryEnergy_MeV", "Distance_um", "DNAEdep_keV", "AnyDamage", "AnyDSB",
    ] + PS_GENERATIVE_COUNT_TARGETS
    result = pd.concat(
        [test[keep].reset_index(drop=True), generated.filter(regex="^(Generated|RawGenerated)")],
        axis=1,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_dir / "test_generated.csv.gz", index=False)
    (args.output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
