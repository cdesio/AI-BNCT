"""Plot one-step versus two-step Extra Trees comparisons."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


COLOURS = ["#5B6770", "#A7ADB2", "#2878B5", "#D95F02"]
COUNT_TARGETS = ["TotalSB", "TotalSSB", "TotalCSSB", "TotalDSB", "TotalCDSB"]


def _read(path: Path) -> dict:
    return json.loads((path / "metrics.json").read_text())


def plot_metric_comparison(
    baseline: Path,
    baseline_history: Path,
    two_step: Path,
    two_step_history: Path,
    output: Path,
) -> None:
    records = [
        ("One-step local", _read(baseline), False),
        ("One-step history", _read(baseline_history), False),
        ("Two-step local", _read(two_step), True),
        ("Two-step history", _read(two_step_history), True),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True)

    x = np.arange(len(COUNT_TARGETS))
    width = 0.19
    for index, (label, metrics, _) in enumerate(records):
        values = [metrics["damage_counts"][target]["mae"] for target in COUNT_TARGETS]
        axes[0].bar(x + (index - 1.5) * width, values, width, label=label,
                    color=COLOURS[index])
    axes[0].set_xticks(x, [name.replace("Total", "") for name in COUNT_TARGETS])
    axes[0].set_ylabel("Mean absolute error")
    axes[0].set_title("Damage-count accuracy")
    axes[0].legend(frameon=False, fontsize=8)

    edep = []
    for _, metrics, is_two_step in records:
        value = metrics["dna_edep"]["mae"] if is_two_step else metrics["dna_edep"]["DNAEdep_keV"]["mae"]
        edep.append(value)
    bars = axes[1].bar(np.arange(4), edep, color=COLOURS)
    axes[1].bar_label(bars, fmt="%.2f", padding=2)
    axes[1].set_xticks(np.arange(4), [item[0].replace(" ", "\n") for item in records])
    axes[1].set_ylabel("MAE (keV)")
    axes[1].set_title("DNA energy deposition")

    width = 0.35
    x = np.arange(4)
    roc = [
        metrics["damage_gate" if is_two_step else "any_damage"]["roc_auc"]
        for _, metrics, is_two_step in records
    ]
    ap = [
        metrics["damage_gate" if is_two_step else "any_damage"]["average_precision"]
        for _, metrics, is_two_step in records
    ]
    axes[2].bar(x - width / 2, roc, width, label="ROC-AUC", color="#2878B5")
    axes[2].bar(x + width / 2, ap, width, label="Average precision", color="#E68613")
    axes[2].set_xticks(x, [item[0].replace(" ", "\n") for item in records])
    axes[2].set_ylim(0.85, 1.0)
    axes[2].set_title("Damage occurrence")
    axes[2].legend(frameon=False)

    fig.suptitle("One-step and two-step Extra Trees on the 1k dataset", fontsize=15)
    fig.savefig(output / "two_step_model_comparison.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_gate_profiles(two_step: Path, two_step_history: Path, output: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), constrained_layout=True, sharey=True)
    for ax, (label, path, colour) in zip(
        axes,
        [("Local", two_step, "#2878B5"), ("History-aware", two_step_history, "#D95F02")],
    ):
        frame = pd.read_csv(path / "test_predictions.csv.gz")
        frame["distance_bin"] = pd.cut(frame["Distance_um"], bins=18, duplicates="drop")
        profile = frame.groupby("distance_bin", observed=True).agg(
            distance=("Distance_um", "mean"),
            simulated=("AnyDamage", "mean"),
            probability=("P_AnyDamage", "mean"),
            gated=("PredictedAnyDamage", "mean"),
        )
        ax.plot(profile["distance"], profile["simulated"], color="#202124", marker="o",
                label="Simulation")
        ax.plot(profile["distance"], profile["probability"], color=colour, linestyle="--",
                marker="o", label="Mean probability")
        ax.plot(profile["distance"], profile["gated"], color=colour, linestyle=":",
                marker="s", label="Positive gate fraction")
        ax.set_title(f"{label} two-step model")
        ax.set_xlabel("Distance along primary trajectory (µm)")
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("Damage fraction")
    axes[0].legend(frameon=False)
    fig.savefig(output / "two_step_damage_vs_distance.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", type=Path, default=Path("runs"))
    parser.add_argument("--output", type=Path, default=Path("runs/ps_two_step_figures_1k"))
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    baseline = args.runs / "ps_baseline_1k"
    baseline_history = args.runs / "ps_baseline_history_1k"
    two_step = args.runs / "ps_two_step_1k"
    two_step_history = args.runs / "ps_two_step_history_1k"
    plot_metric_comparison(baseline, baseline_history, two_step, two_step_history, args.output)
    plot_gate_profiles(two_step, two_step_history, args.output)
    print(f"Wrote figures to {args.output}")


if __name__ == "__main__":
    main()
