"""Generate joint DNA-voxel responses from phase-space rows using conditional StEG."""

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch

from .ps_predict import load_phase_space
from .schema import PS_GENERATIVE_COUNT_TARGETS
from .steg import ConditionalStEG, generate


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_root", type=Path)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples-per-row", type=int, default=1)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--seed", type=int, default=20260930)
    args = parser.parse_args()
    if args.samples_per_row < 1:
        parser.error("--samples-per-row must be positive")

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    checkpoint = torch.load(args.model_dir / "model.pt", map_location="cpu", weights_only=False)
    metadata = checkpoint["metadata"]
    if not metadata["stage"].startswith("ps_response"):
        raise ValueError("The supplied StEG checkpoint is not a ps_response model")
    transformers = joblib.load(args.model_dir / "transformers.joblib")
    frame = load_phase_space(args.input_root)
    repeated = frame.loc[frame.index.repeat(args.samples_per_row)].reset_index(drop=True)
    repeated["ResponseSample"] = np.tile(np.arange(args.samples_per_row), len(frame))
    condition = transformers["condition"].transform(
        repeated[metadata["condition_columns"]]
    ).astype(np.float32)

    device = torch.device(args.device)
    model = ConditionalStEG(metadata["output_dim"], metadata["condition_dim"]).to(device)
    model.load_state_dict(checkpoint["state_dict"])
    model.eval()
    generated_scaled = generate(
        model,
        torch.from_numpy(condition).to(device),
        metadata["output_dim"],
        metadata["diffusion_steps"],
        device,
    ).cpu().numpy()
    generated = transformers["output"].inverse_transform(generated_scaled)
    outputs = pd.DataFrame(generated, columns=metadata["output_columns"])
    outputs["DNAEdep_keV"] = outputs["DNAEdep_keV"].clip(lower=0)
    for target in PS_GENERATIVE_COUNT_TARGETS:
        outputs[target] = outputs[target].clip(lower=0).round().astype(int)
    outputs["TotalCDSB"] = np.minimum(outputs["TotalCDSB"], outputs["TotalDSB"])
    minimum_breaks = outputs["TotalSSB"] + outputs["TotalCSSB"] + 2 * outputs["TotalDSB"]
    outputs["TotalSB"] = np.maximum(outputs["TotalSB"], minimum_breaks)
    outputs["AnyDamage"] = (outputs["TotalSB"] > 0).astype(int)
    outputs["AnyDSB"] = (outputs["TotalDSB"] > 0).astype(int)

    provenance = [
        "EventID", "SeedID", "TrackID", "ParentID", "VoxelID",
        "VoxelX", "VoxelY", "VoxelZ", "ResponseSample",
    ]
    result = pd.concat([repeated[provenance], outputs], axis=1)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"Wrote {len(result):,} generated voxel responses to {args.output}")


if __name__ == "__main__":
    main()
