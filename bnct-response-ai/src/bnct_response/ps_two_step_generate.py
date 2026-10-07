"""Generate voxel responses with a calibrated gate and positive-damage StEG."""

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch

from .ps_predict import load_phase_space
from .ps_two_step import calibrated_probability
from .schema import PS_GENERATIVE_COUNT_TARGETS
from .steg import ConditionalStEG, generate


def generate_two_step(
    frame: pd.DataFrame,
    gate_dir: Path,
    steg_dir: Path,
    samples_per_row: int,
    device_name: str,
    batch_size: int,
    seed: int,
) -> pd.DataFrame:
    if samples_per_row < 1:
        raise ValueError("samples_per_row must be positive")
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    gate = joblib.load(gate_dir / "models.joblib")
    checkpoint = torch.load(steg_dir / "model.pt", map_location="cpu", weights_only=False)
    metadata = checkpoint["metadata"]
    if metadata["stage"] != "ps_positive_damage":
        raise ValueError("The StEG checkpoint must use stage ps_positive_damage")
    transformers = joblib.load(steg_dir / "transformers.joblib")

    repeated = frame.loc[frame.index.repeat(samples_per_row)].reset_index(drop=True)
    repeated["ResponseSample"] = np.tile(np.arange(samples_per_row), len(frame))
    raw_probability = gate["damage_gate"].predict_proba(
        repeated[gate["input_columns"]]
    )[:, 1]
    probability = calibrated_probability(gate["damage_calibrator"], raw_probability)
    positive = rng.random(len(repeated)) < probability

    edep = np.expm1(gate["edep"].predict(repeated[gate["input_columns"]])).clip(min=0)
    counts = np.zeros((len(repeated), len(PS_GENERATIVE_COUNT_TARGETS)), dtype=int)
    positive_indices = np.flatnonzero(positive)
    if len(positive_indices):
        device = torch.device(device_name)
        model = ConditionalStEG(metadata["output_dim"], metadata["condition_dim"]).to(device)
        model.load_state_dict(checkpoint["state_dict"])
        model.eval()
        for start in range(0, len(positive_indices), batch_size):
            indices = positive_indices[start:start + batch_size]
            condition = transformers["condition"].transform(
                repeated.iloc[indices][metadata["condition_columns"]]
            ).astype(np.float32)
            scaled = generate(
                model, torch.from_numpy(condition).to(device), metadata["output_dim"],
                metadata["diffusion_steps"], device,
            ).cpu().numpy()
            values = transformers["output"].inverse_transform(scaled)
            counts[indices] = np.clip(np.rint(values), 0, None).astype(int)

    outputs = pd.DataFrame(counts, columns=PS_GENERATIVE_COUNT_TARGETS)
    outputs["TotalCDSB"] = np.minimum(outputs["TotalCDSB"], outputs["TotalDSB"])
    minimum_breaks = outputs["TotalSSB"] + outputs["TotalCSSB"] + 2 * outputs["TotalDSB"]
    outputs["TotalSB"] = np.maximum(outputs["TotalSB"], minimum_breaks)
    outputs.loc[positive, "TotalSB"] = np.maximum(outputs.loc[positive, "TotalSB"], 1)
    outputs["AnyDamage"] = positive.astype(int)
    outputs["AnyDSB"] = (outputs["TotalDSB"] > 0).astype(int)
    outputs["DNAEdep_keV"] = edep
    outputs["P_AnyDamage"] = probability
    outputs["RawP_AnyDamage"] = raw_probability
    return pd.concat([repeated.reset_index(drop=True), outputs.add_prefix("Generated")], axis=1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("--gate-model-dir", type=Path, required=True)
    parser.add_argument("--steg-model-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples-per-row", type=int, default=1)
    parser.add_argument("--batch-size", type=int, default=4096)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--seed", type=int, default=20261007)
    args = parser.parse_args()
    result = generate_two_step(
        load_phase_space(args.input_root), args.gate_model_dir, args.steg_model_dir,
        args.samples_per_row, args.device, args.batch_size, args.seed,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"Wrote {len(result):,} generated voxel responses to {args.output}")


if __name__ == "__main__":
    main()
