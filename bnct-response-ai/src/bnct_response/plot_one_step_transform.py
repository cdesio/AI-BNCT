"""Compare long one-step StEG runs with quantile and log-standard outputs."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TARGETS = ["TotalSB", "TotalSSB", "TotalCSSB", "TotalDSB", "TotalCDSB"]
LABELS = ["SB", "SSB", "CSSB", "DSB", "CDSB"]
RUNS = [
    ("Quantile", "ps_one_step_long_quantile", "#2878B5"),
    ("log1p + standardise", "ps_one_step_long_log", "#D95F02"),
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-dir", type=Path, default=Path("runs"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("runs/one_step_transform_figures"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    records = []
    for label, folder, colour in RUNS:
        metrics = json.loads((args.runs_dir / f"{folder}_eval" / "metrics.json").read_text())
        history = pd.DataFrame(json.loads((args.runs_dir / folder / "history.json").read_text()))
        records.append((label, metrics, history, colour))

    truth = records[0][1]
    truth_means = [truth["outputs"][target]["true_mean"] for target in TARGETS]
    truth_zeros = [truth["outputs"][target]["true_zero_fraction"] for target in TARGETS]
    x = np.arange(len(TARGETS))
    fig, axes = plt.subplots(2, 2, figsize=(13, 8.5), constrained_layout=True)

    axes[0, 0].plot(x, truth_means, color="#202124", marker="o", linewidth=2,
                    label="Simulation")
    axes[0, 1].plot(x, truth_zeros, color="#202124", marker="o", linewidth=2,
                    label="Simulation")
    width = 0.32
    for index, (label, metrics, _, colour) in enumerate(records):
        offset = (index - 0.5) * width
        means = [metrics["outputs"][target]["generated_mean"] for target in TARGETS]
        zeros = [metrics["outputs"][target]["generated_zero_fraction"] for target in TARGETS]
        axes[0, 0].bar(x + offset, means, width, color=colour, label=label)
        axes[0, 1].bar(x + offset, zeros, width, color=colour, label=label)
    for ax in axes[0]:
        ax.set_xticks(x, LABELS)
        ax.grid(axis="y", alpha=0.2)
        ax.legend(frameon=False)
    axes[0, 0].set_title("Generated count means")
    axes[0, 0].set_ylabel("Mean count per held-out voxel")
    axes[0, 1].set_title("Generated zero fractions")
    axes[0, 1].set_ylabel("Zero fraction")

    error_names = ["Mean-count", "Zero-fraction", "Occurrence"]
    summary = []
    for _, metrics, _, _ in records:
        summary.append([
            np.mean([abs(metrics["outputs"][target]["generated_mean"]
                             - metrics["outputs"][target]["true_mean"])
                     for target in TARGETS]),
            np.mean([abs(metrics["outputs"][target]["generated_zero_fraction"]
                             - metrics["outputs"][target]["true_zero_fraction"])
                     for target in TARGETS]),
            abs(metrics["damage_occurrence"]["generated_fraction"]
                - metrics["damage_occurrence"]["true_fraction"]),
        ])
    error_x = np.arange(len(error_names))
    for index, ((label, _, _, colour), values) in enumerate(zip(records, summary)):
        axes[1, 0].bar(error_x + (index - 0.5) * width, values, width,
                       color=colour, label=label)
    axes[1, 0].set_xticks(error_x, error_names)
    axes[1, 0].set_title("Physical distribution errors (lower is better)")
    axes[1, 0].legend(frameon=False)
    axes[1, 0].grid(axis="y", alpha=0.2)

    for label, _, history, colour in records:
        axes[1, 1].plot(history["epoch"], history["validation_loss"],
                        color=colour, label=label)
    axes[1, 1].set_title("Validation diffusion loss")
    axes[1, 1].set_xlabel("Epoch")
    axes[1, 1].set_ylabel("Noise-prediction loss")
    axes[1, 1].legend(frameon=False)
    axes[1, 1].grid(alpha=0.2)

    fig.suptitle("Long one-step local StEG: output-transform comparison", fontsize=15)
    fig.savefig(args.output_dir / "one_step_transform_comparison.png", dpi=180,
                bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote figures to {args.output_dir}")


if __name__ == "__main__":
    main()
