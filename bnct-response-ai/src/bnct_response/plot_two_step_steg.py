"""Plot held-out diagnostics for the two-step positive-damage StEG."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TARGETS = ["TotalSB", "TotalSSB", "TotalCSSB", "TotalDSB", "TotalCDSB"]
LABELS = ["SB", "SSB", "CSSB", "DSB", "CDSB"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluation-dir", type=Path, required=True)
    parser.add_argument("--steg-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.evaluation_dir / "test_generated.csv.gz")
    history = pd.DataFrame(json.loads((args.steg_dir / "history.json").read_text()))

    fig, axes = plt.subplots(2, 3, figsize=(13, 7.5), constrained_layout=True)
    axes[0, 0].plot(history["epoch"], history["train_loss"], alpha=0.55, label="Training")
    axes[0, 0].plot(history["epoch"], history["validation_loss"], linewidth=2,
                    label="Validation")
    axes[0, 0].set_title("Positive-damage StEG convergence")
    axes[0, 0].set_xlabel("Epoch")
    axes[0, 0].set_ylabel("Noise-prediction loss")
    axes[0, 0].legend(frameon=False)
    axes[0, 0].grid(alpha=0.2)

    for ax, target, label in zip(axes.flat[1:], TARGETS, LABELS):
        truth = frame[target].to_numpy()
        generated = frame[f"Generated{target}"].to_numpy()
        upper = int(max(np.quantile(truth, 0.995), np.quantile(generated, 0.995), 1))
        bins = np.arange(-0.5, upper + 1.5)
        ax.hist(truth, bins=bins, density=True, histtype="step", linewidth=2,
                color="#202124", label="Simulation")
        ax.hist(generated, bins=bins, density=True, histtype="step", linewidth=2,
                color="#2878B5", label="Generated")
        ax.set_title(label)
        ax.set_xlabel("Count per voxel")
        ax.set_ylabel("Probability density")
        ax.grid(alpha=0.15)
    axes[0, 1].legend(frameon=False)
    fig.savefig(args.output_dir / "two_step_steg_distributions.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)

    true_positive = frame[frame["AnyDamage"] == 1]
    generated_positive = frame[frame["GeneratedAnyDamage"] == 1]
    true_zero = [(true_positive[target] == 0).mean() for target in TARGETS]
    generated_zero = [
        (generated_positive[f"Generated{target}"] == 0).mean() for target in TARGETS
    ]
    true_mean = [true_positive[target].mean() for target in TARGETS]
    generated_mean = [
        generated_positive[f"Generated{target}"].mean() for target in TARGETS
    ]

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5), constrained_layout=True)
    x = np.arange(len(TARGETS))
    width = 0.36
    axes[0].bar(x - width / 2, true_zero, width, color="#202124", label="Simulation")
    axes[0].bar(x + width / 2, generated_zero, width, color="#2878B5", label="Generated")
    axes[0].set_xticks(x, LABELS)
    axes[0].set_ylim(0, 1)
    axes[0].set_ylabel("Zero fraction within damaged voxels")
    axes[0].set_title("Conditional sparsity")
    axes[0].legend(frameon=False)
    axes[1].bar(x - width / 2, true_mean, width, color="#202124", label="Simulation")
    axes[1].bar(x + width / 2, generated_mean, width, color="#2878B5", label="Generated")
    axes[1].set_xticks(x, LABELS)
    axes[1].set_ylabel("Mean count within damaged voxels")
    axes[1].set_title("Conditional count means")
    fig.savefig(args.output_dir / "two_step_steg_positive_diagnostics.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote figures to {args.output_dir}")


if __name__ == "__main__":
    main()
