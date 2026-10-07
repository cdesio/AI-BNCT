"""Compare positive-damage StEG parameter-sweep results."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


RUNS = {
    "A: q50 small": "ps_positive_sweep_A_q50_small",
    "B: q100 small": "ps_positive_sweep_B_q100_small",
    "C: q100 medium": "ps_positive_sweep_C_q100_medium",
    "D: log100 small": "ps_positive_sweep_D_log100_small",
}
TARGETS = ["TotalSB", "TotalSSB", "TotalCSSB", "TotalDSB", "TotalCDSB"]
LABELS = ["SB", "SSB", "CSSB", "DSB", "CDSB"]
COLOURS = ["#6C757D", "#A7ADB2", "#2878B5", "#D95F02"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-dir", type=Path, default=Path("runs"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("runs/positive_steg_sweep_figures"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    metrics = {
        label: json.loads((args.runs_dir / f"{folder}_eval" / "metrics.json").read_text())
        for label, folder in RUNS.items()
    }
    histories = {
        label: pd.DataFrame(json.loads((args.runs_dir / folder / "history.json").read_text()))
        for label, folder in RUNS.items()
    }
    truth_means = [next(iter(metrics.values()))["outputs"][t]["true_mean"] for t in TARGETS]
    truth_zeros = [
        next(iter(metrics.values()))["outputs"][t]["true_zero_fraction"] for t in TARGETS
    ]

    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    x = np.arange(len(TARGETS))
    width = 0.18
    axes[0, 0].plot(x, truth_means, color="#202124", marker="o", linewidth=2,
                    label="Simulation")
    axes[0, 1].plot(x, truth_zeros, color="#202124", marker="o", linewidth=2,
                    label="Simulation")
    for index, (label, result) in enumerate(metrics.items()):
        offset = (index - 1.5) * width
        means = [result["outputs"][t]["generated_mean"] for t in TARGETS]
        zeros = [result["outputs"][t]["generated_zero_fraction"] for t in TARGETS]
        axes[0, 0].bar(x + offset, means, width, color=COLOURS[index], label=label)
        axes[0, 1].bar(x + offset, zeros, width, color=COLOURS[index], label=label)
    for ax in axes[0]:
        ax.set_xticks(x, LABELS)
        ax.grid(axis="y", alpha=0.2)
    axes[0, 0].set_title("Generated count means")
    axes[0, 0].set_ylabel("Mean count per test voxel")
    axes[0, 0].legend(frameon=False, fontsize=8)
    axes[0, 1].set_title("Generated zero fractions")
    axes[0, 1].set_ylabel("Zero fraction")
    axes[0, 1].set_ylim(0.7, 1.0)

    mean_errors = []
    zero_errors = []
    correlation_errors = []
    for label, folder in RUNS.items():
        result = metrics[label]
        mean_errors.append(np.mean([
            abs(result["outputs"][t]["generated_mean"] - result["outputs"][t]["true_mean"])
            for t in TARGETS
        ]))
        zero_errors.append(np.mean([
            abs(result["outputs"][t]["generated_zero_fraction"]
                - result["outputs"][t]["true_zero_fraction"])
            for t in TARGETS
        ]))
        frame = pd.read_csv(args.runs_dir / f"{folder}_eval" / "test_generated.csv.gz")
        truth_corr = frame[frame["AnyDamage"] == 1][TARGETS].corr().to_numpy()
        generated_corr = frame[frame["GeneratedAnyDamage"] == 1][
            [f"Generated{t}" for t in TARGETS]
        ].corr().to_numpy()
        correlation_errors.append(np.mean(np.abs(truth_corr - generated_corr)))
    labels = list(RUNS)
    summary_x = np.arange(len(labels))
    axes[1, 0].bar(summary_x - width, mean_errors, width, label="Mean-count error")
    axes[1, 0].bar(summary_x, zero_errors, width, label="Zero-fraction error")
    axes[1, 0].bar(summary_x + width, correlation_errors, width,
                   label="Correlation error")
    axes[1, 0].set_xticks(summary_x, [label.split(":")[0] for label in labels])
    axes[1, 0].set_title("Physical distribution errors (lower is better)")
    axes[1, 0].legend(frameon=False, fontsize=8)
    axes[1, 0].grid(axis="y", alpha=0.2)

    for index, (label, history) in enumerate(histories.items()):
        axes[1, 1].plot(history["epoch"], history["validation_loss"],
                        color=COLOURS[index], label=label)
    axes[1, 1].set_title("Validation diffusion loss")
    axes[1, 1].set_xlabel("Epoch")
    axes[1, 1].set_ylabel("Noise-prediction loss")
    axes[1, 1].legend(frameon=False, fontsize=8)
    axes[1, 1].grid(alpha=0.2)

    fig.suptitle("Positive-damage StEG parameter sweep", fontsize=15)
    fig.savefig(args.output_dir / "positive_steg_sweep_comparison.png", dpi=180,
                bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote figures to {args.output_dir}")


if __name__ == "__main__":
    main()
