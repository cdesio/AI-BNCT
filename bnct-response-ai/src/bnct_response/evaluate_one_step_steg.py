"""Evaluate a one-step phase-space StEG on held-out primaries."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch

from .schema import PS_GENERATIVE_COUNT_TARGETS
from .steg import ConditionalStEG, generate, inverse_steg_output


def generate_one_step(
    frame: pd.DataFrame,
    model_dir: Path,
    device_name: str,
    batch_size: int,
    seed: int,
) -> pd.DataFrame:
    """Generate one response sample for every supplied phase-space row."""
    torch.manual_seed(seed)
    np.random.seed(seed)
    checkpoint = torch.load(model_dir / "model.pt", map_location="cpu", weights_only=False)
    metadata = checkpoint["metadata"]
    if metadata["stage"] != "ps_response":
        raise ValueError("The supplied checkpoint must use the local ps_response stage")
    transformers = joblib.load(model_dir / "transformers.joblib")
    condition = transformers["condition"].transform(
        frame[metadata["condition_columns"]]
    ).astype(np.float32)

    device = torch.device(device_name)
    model = ConditionalStEG(
        metadata["output_dim"], metadata["condition_dim"],
        width=metadata.get("width", 256), layers=metadata.get("layers", 5),
    ).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()

    generated_batches = []
    for start in range(0, len(frame), batch_size):
        batch = torch.from_numpy(condition[start:start + batch_size]).to(device)
        generated_batches.append(generate(
            model, batch, metadata["output_dim"], metadata["diffusion_steps"], device,
        ).cpu().numpy())
    scaled = np.concatenate(generated_batches, axis=0)
    values = inverse_steg_output(scaled, transformers["output"], metadata)
    outputs = pd.DataFrame(values, columns=metadata["output_columns"])
    outputs["DNAEdep_keV"] = outputs["DNAEdep_keV"].clip(lower=0)
    for target in PS_GENERATIVE_COUNT_TARGETS:
        outputs[target] = outputs[target].clip(lower=0).round().astype(int)

    raw_counts = outputs[PS_GENERATIVE_COUNT_TARGETS].copy()
    outputs["TotalCDSB"] = np.minimum(outputs["TotalCDSB"], outputs["TotalDSB"])
    minimum_breaks = outputs["TotalSSB"] + outputs["TotalCSSB"] + 2 * outputs["TotalDSB"]
    outputs["TotalSB"] = np.maximum(outputs["TotalSB"], minimum_breaks)
    outputs["AnyDamage"] = (outputs["TotalSB"] > 0).astype(int)
    outputs["AnyDSB"] = (outputs["TotalDSB"] > 0).astype(int)
    return pd.concat(
        [raw_counts.add_prefix("RawGenerated"), outputs.add_prefix("Generated")],
        axis=1,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--seed", type=int, default=20261008)
    args = parser.parse_args()

    test = pd.read_csv(args.data_dir / "ps_response.csv.gz")
    test = test[test["split"] == "test"].reset_index(drop=True)
    generated = generate_one_step(
        test, args.model_dir, args.device, args.batch_size, args.seed
    )

    metrics = {
        "rows": len(test),
        "damage_occurrence": {
            "true_fraction": float(test["AnyDamage"].mean()),
            "generated_fraction": float(generated["GeneratedAnyDamage"].mean()),
        },
        "dna_edep": {
            "true_mean": float(test["DNAEdep_keV"].mean()),
            "generated_mean": float(generated["GeneratedDNAEdep_keV"].mean()),
            "sample_mae": float(np.mean(np.abs(
                test["DNAEdep_keV"].to_numpy()
                - generated["GeneratedDNAEdep_keV"].to_numpy()
            ))),
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
        }

    keep = [
        "case_id", "SeedID", "EventID", "TrackID", "VoxelID", "Particle",
        "EntryEnergy_MeV", "Distance_um", "DNAEdep_keV", "AnyDamage", "AnyDSB",
    ] + PS_GENERATIVE_COUNT_TARGETS
    result = pd.concat([test[keep], generated], axis=1)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output_dir / "test_generated.csv.gz", index=False)
    (args.output_dir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
