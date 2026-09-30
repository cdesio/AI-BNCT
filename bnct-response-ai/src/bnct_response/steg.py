"""Conditional tabular diffusion (StEG) for transport or damage responses."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
from sklearn.preprocessing import QuantileTransformer
from torch import nn
from torch.utils.data import DataLoader, TensorDataset

from .baseline import make_preprocessor
from .schema import (
    CONDITION_COLUMNS,
    DAMAGE_INPUT_COLUMNS,
    DAMAGE_TARGETS,
    PS_INPUT_COLUMNS,
    PS_HISTORY_INPUT_COLUMNS,
    PS_RESPONSE_TARGETS,
    TRANSPORT_TARGETS,
)


class ResidualBlock(nn.Module):
    def __init__(self, width: int, context_dim: int):
        super().__init__()
        self.norm = nn.LayerNorm(width)
        self.context = nn.Linear(context_dim, width * 2)
        self.net = nn.Sequential(nn.Linear(width, width), nn.SiLU(), nn.Linear(width, width))

    def forward(self, x: torch.Tensor, context: torch.Tensor) -> torch.Tensor:
        scale, shift = self.context(context).chunk(2, dim=-1)
        hidden = self.norm(x) * (1 + scale) + shift
        return x + self.net(hidden)


class ConditionalStEG(nn.Module):
    def __init__(self, output_dim: int, condition_dim: int, width: int = 256, layers: int = 5):
        super().__init__()
        time_dim = 32
        self.time_dim = time_dim
        self.input = nn.Linear(output_dim, width)
        self.condition = nn.Sequential(
            nn.Linear(condition_dim + time_dim, width), nn.SiLU(), nn.Linear(width, width)
        )
        self.blocks = nn.ModuleList([ResidualBlock(width, width) for _ in range(layers)])
        self.output = nn.Linear(width, output_dim)

    def time_embedding(self, t: torch.Tensor) -> torch.Tensor:
        half = self.time_dim // 2
        frequencies = torch.exp(
            -math.log(10_000) * torch.arange(half, device=t.device) / max(half - 1, 1)
        )
        angles = t[:, None] * frequencies[None, :]
        return torch.cat([angles.sin(), angles.cos()], dim=-1)

    def forward(self, noisy: torch.Tensor, t: torch.Tensor, condition: torch.Tensor) -> torch.Tensor:
        context = self.condition(torch.cat([condition, self.time_embedding(t)], dim=-1))
        hidden = self.input(noisy)
        for block in self.blocks:
            hidden = block(hidden, context)
        return self.output(hidden)


def _schedule(steps: int, device: torch.device) -> tuple[torch.Tensor, ...]:
    beta = torch.linspace(1e-4, 0.02, steps, device=device)
    alpha = 1.0 - beta
    alpha_bar = torch.cumprod(alpha, dim=0)
    return beta, alpha, alpha_bar


@torch.no_grad()
def generate(
    model: ConditionalStEG,
    condition: torch.Tensor,
    output_dim: int,
    diffusion_steps: int,
    device: torch.device,
) -> torch.Tensor:
    beta, alpha, alpha_bar = _schedule(diffusion_steps, device)
    value = torch.randn(len(condition), output_dim, device=device)
    for step in reversed(range(diffusion_steps)):
        t = torch.full((len(condition),), step / max(diffusion_steps - 1, 1), device=device)
        noise_prediction = model(value, t, condition)
        coefficient = beta[step] / torch.sqrt(1 - alpha_bar[step])
        mean = (value - coefficient * noise_prediction) / torch.sqrt(alpha[step])
        if step:
            value = mean + torch.sqrt(beta[step]) * torch.randn_like(value)
        else:
            value = mean
    return value


def train_steg(
    data_dir: Path,
    output_dir: Path,
    stage: str,
    epochs: int,
    batch_size: int,
    diffusion_steps: int,
    device_name: str,
    seed: int,
) -> dict:
    torch.manual_seed(seed)
    np.random.seed(seed)
    data_file = "ps_response.csv.gz" if stage.startswith("ps_response") else "response.csv.gz"
    frame = pd.read_csv(data_dir / data_file)
    train = frame[frame["split"] == "train"].copy()
    validation = frame[frame["split"] == "val"].copy()
    if stage == "transport":
        condition_columns = CONDITION_COLUMNS
        output_columns = TRANSPORT_TARGETS
    elif stage == "damage":
        condition_columns = DAMAGE_INPUT_COLUMNS
        output_columns = DAMAGE_TARGETS
    elif stage == "ps_response":
        condition_columns = PS_INPUT_COLUMNS
        output_columns = PS_RESPONSE_TARGETS
    elif stage == "ps_response_history":
        condition_columns = PS_HISTORY_INPUT_COLUMNS
        output_columns = PS_RESPONSE_TARGETS
    else:
        raise ValueError(
            "stage must be 'transport', 'damage', 'ps_response' or 'ps_response_history'"
        )

    condition_transformer = make_preprocessor(condition_columns)
    train_condition = condition_transformer.fit_transform(train[condition_columns]).astype(np.float32)
    validation_condition = condition_transformer.transform(validation[condition_columns]).astype(np.float32)

    output_transformer = QuantileTransformer(
        n_quantiles=min(1000, len(train)), output_distribution="normal",
        subsample=None, random_state=seed,
    )
    train_output = output_transformer.fit_transform(train[output_columns]).astype(np.float32)
    validation_output = output_transformer.transform(validation[output_columns]).astype(np.float32)

    device = torch.device(device_name)
    model = ConditionalStEG(
        output_dim=len(output_columns), condition_dim=train_condition.shape[1]
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=2e-4, weight_decay=1e-5)
    _, _, alpha_bar = _schedule(diffusion_steps, device)
    loader = DataLoader(
        TensorDataset(torch.from_numpy(train_condition), torch.from_numpy(train_output)),
        batch_size=batch_size, shuffle=True, drop_last=False,
    )
    validation_condition_tensor = torch.from_numpy(validation_condition).to(device)
    validation_output_tensor = torch.from_numpy(validation_output).to(device)

    history = []
    best_validation = float("inf")
    best_state = None
    for epoch in range(1, epochs + 1):
        model.train()
        losses = []
        for condition_batch, output_batch in loader:
            condition_batch = condition_batch.to(device)
            output_batch = output_batch.to(device)
            step = torch.randint(0, diffusion_steps, (len(output_batch),), device=device)
            noise = torch.randn_like(output_batch)
            cumulative = alpha_bar[step, None]
            noisy = torch.sqrt(cumulative) * output_batch + torch.sqrt(1 - cumulative) * noise
            normalized_step = step.float() / max(diffusion_steps - 1, 1)
            prediction = model(noisy, normalized_step, condition_batch)
            loss = nn.functional.mse_loss(prediction, noise)
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            losses.append(float(loss.detach().cpu()))

        model.eval()
        with torch.no_grad():
            step = torch.randint(0, diffusion_steps, (len(validation_output_tensor),), device=device)
            noise = torch.randn_like(validation_output_tensor)
            cumulative = alpha_bar[step, None]
            noisy = (
                torch.sqrt(cumulative) * validation_output_tensor
                + torch.sqrt(1 - cumulative) * noise
            )
            prediction = model(
                noisy, step.float() / max(diffusion_steps - 1, 1), validation_condition_tensor
            )
            validation_loss = float(nn.functional.mse_loss(prediction, noise).cpu())
        record = {
            "epoch": epoch,
            "train_loss": float(np.mean(losses)),
            "validation_loss": validation_loss,
        }
        history.append(record)
        print(json.dumps(record))
        if validation_loss < best_validation:
            best_validation = validation_loss
            best_state = {name: value.detach().cpu() for name, value in model.state_dict().items()}

    output_dir.mkdir(parents=True, exist_ok=True)
    metadata = {
        "stage": stage,
        "condition_columns": condition_columns,
        "output_columns": output_columns,
        "condition_dim": int(train_condition.shape[1]),
        "output_dim": len(output_columns),
        "diffusion_steps": diffusion_steps,
        "epochs": epochs,
        "seed": seed,
        "best_validation_loss": best_validation,
    }
    torch.save({"state_dict": best_state, "metadata": metadata}, output_dir / "model.pt")
    joblib.dump(
        {"condition": condition_transformer, "output": output_transformer},
        output_dir / "transformers.joblib",
    )
    (output_dir / "history.json").write_text(json.dumps(history, indent=2) + "\n")
    (output_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    model.load_state_dict(best_state)
    preview_count = min(512, len(validation))
    preview_indices = np.random.default_rng(seed).choice(
        len(validation), size=preview_count, replace=False
    )
    generated_scaled = generate(
        model, validation_condition_tensor[preview_indices], len(output_columns),
        diffusion_steps, device,
    ).cpu().numpy()
    generated = output_transformer.inverse_transform(generated_scaled)
    if stage.startswith("ps_response"):
        generated[:, output_columns.index("DNAEdep_keV")] = np.clip(
            generated[:, output_columns.index("DNAEdep_keV")], 0, None
        )
        for column in output_columns:
            if column != "DNAEdep_keV":
                index = output_columns.index(column)
                generated[:, index] = np.clip(np.rint(generated[:, index]), 0, None)
        dsb_index = output_columns.index("TotalDSB")
        cdsb_index = output_columns.index("TotalCDSB")
        generated[:, cdsb_index] = np.minimum(generated[:, cdsb_index], generated[:, dsb_index])
        sb_index = output_columns.index("TotalSB")
        ssb_index = output_columns.index("TotalSSB")
        cssb_index = output_columns.index("TotalCSSB")
        minimum_breaks = (
            generated[:, ssb_index] + generated[:, cssb_index] + 2 * generated[:, dsb_index]
        )
        generated[:, sb_index] = np.maximum(generated[:, sb_index], minimum_breaks)
    descriptor_columns = ["case_id", "SeedID", "EventID", "Distance_um"]
    if stage.startswith("ps_response"):
        descriptor_columns = [
            "case_id", "SeedID", "EventID", "Particle", "EntryEnergy_MeV", "Distance_um"
        ]
    preview = validation.iloc[preview_indices][descriptor_columns + output_columns].reset_index(
        drop=True
    )
    for index, column in enumerate(output_columns):
        values = generated[:, index]
        if column in DAMAGE_TARGETS:
            values = np.clip(values, 0, None)
        preview[f"generated_{column}"] = values
    preview.to_csv(output_dir / "validation_preview.csv.gz", index=False)
    return metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--stage",
        choices=["transport", "damage", "ps_response", "ps_response_history"],
        required=True,
    )
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--diffusion-steps", type=int, default=100)
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--seed", type=int, default=20260923)
    args = parser.parse_args()
    metadata = train_steg(
        args.data_dir, args.output_dir, args.stage, args.epochs,
        args.batch_size, args.diffusion_steps, args.device, args.seed,
    )
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
