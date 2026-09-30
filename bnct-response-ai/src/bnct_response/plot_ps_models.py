"""Plot preliminary comparisons of the four phase-space response models."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


COLOURS = {"truth": "#202124", "local": "#2878B5", "history": "#D95F02"}


def _baseline_predictions(frame: pd.DataFrame, model_dir: Path) -> pd.DataFrame:
    bundle = joblib.load(model_dir / "models.joblib")
    test = frame[frame["split"] == "test"].copy()
    inputs = test[bundle["input_columns"]]
    test["prediction"] = np.expm1(bundle["edep"].predict(inputs)).clip(min=0)
    return test


def _response_profile(ax: plt.Axes, frame: pd.DataFrame, prediction: str, title: str) -> None:
    markers = {"alpha": "o", "lithium": "s"}
    for particle, raw in frame.groupby("Particle"):
        raw = raw.copy()
        raw["energy_bin"] = pd.cut(raw["EntryEnergy_MeV"], bins=12, duplicates="drop")
        part = raw.groupby("energy_bin", observed=True).agg(
            energy=("EntryEnergy_MeV", "mean"),
            truth=("DNAEdep_keV", "mean"),
            predicted=(prediction, "mean"),
        )
        marker = markers.get(str(particle).lower(), "o")
        ax.plot(part["energy"], part["truth"], color=COLOURS["truth"], marker=marker,
                linestyle="-")
        ax.plot(part["energy"], part["predicted"], color="#2878B5", marker=marker,
                linestyle="--")
    ax.set_title(title)
    ax.set_xlabel("Voxel entry energy (MeV)")
    ax.set_ylabel("Mean DNA energy deposition (keV)")
    ax.grid(alpha=0.2)


def _steg_label(label: str, run_dir: Path) -> str:
    metadata = json.loads((run_dir / "metadata.json").read_text())
    return f"StEG: {label} ({metadata['epochs']} epochs)"


def plot_overview(
    data: pd.DataFrame, runs: Path, local_steg: Path, history_steg: Path, output: Path
) -> None:
    variants = [
        ("Extra Trees: local", _baseline_predictions(data, runs / "ps_baseline"), "prediction"),
        ("Extra Trees: history-aware", _baseline_predictions(data, runs / "ps_baseline_history"), "prediction"),
        (_steg_label("local", local_steg), pd.read_csv(local_steg / "validation_preview.csv.gz"), "generated_DNAEdep_keV"),
        (_steg_label("history-aware", history_steg), pd.read_csv(history_steg / "validation_preview.csv.gz"), "generated_DNAEdep_keV"),
    ]
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, (title, frame, prediction) in zip(axes.flat, variants):
        _response_profile(ax, frame, prediction, title)
    handles = [
        Line2D([], [], color=COLOURS["truth"], linestyle="-", label="Simulation"),
        Line2D([], [], color="#2878B5", linestyle="--", label="Model"),
        Line2D([], [], color="#666666", marker="o", linestyle="None", label="Alpha"),
        Line2D([], [], color="#666666", marker="s", linestyle="None", label="Lithium"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False)
    fig.suptitle("Voxel response versus phase-space entry energy", fontsize=15, y=0.98)
    fig.subplots_adjust(top=0.90, bottom=0.12, hspace=0.38, wspace=0.22)
    fig.savefig(output / "four_model_response.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_baseline_metrics(runs: Path, output: Path) -> None:
    records = []
    for label, folder in [("Local", "ps_baseline"), ("History-aware", "ps_baseline_history")]:
        metrics = json.loads((runs / folder / "metrics.json").read_text())
        records.append((label, metrics))
    names = ["Edep MAE\n(keV)", "TotalSB MAE", "Damage ROC-AUC", "DSB ROC-AUC"]
    values = {
        label: [
            metrics["dna_edep"]["DNAEdep_keV"]["mae"],
            metrics["damage_counts"]["TotalSB"]["mae"],
            metrics["any_damage"]["roc_auc"],
            metrics["any_dsb"]["roc_auc"],
        ]
        for label, metrics in records
    }
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), constrained_layout=True)
    x = np.arange(2)
    width = 0.34
    for offset, metric_index in [(-width / 2, 0), (width / 2, 1)]:
        bars = axes[0].bar(x + offset, [v[metric_index] for v in values.values()], width,
                           label=names[metric_index])
        axes[0].bar_label(bars, fmt="%.2f", padding=2, fontsize=9)
    axes[0].set_xticks(x, values.keys())
    axes[0].set_ylabel("Mean absolute error (lower is better)")
    axes[0].legend(frameon=False)
    for offset, metric_index in [(-width / 2, 2), (width / 2, 3)]:
        bars = axes[1].bar(x + offset, [v[metric_index] for v in values.values()], width,
                           label=names[metric_index])
        axes[1].bar_label(bars, fmt="%.3f", padding=2, fontsize=9)
    axes[1].set_xticks(x, values.keys())
    axes[1].set_ylim(0.85, 1.0)
    axes[1].set_ylabel("ROC-AUC (higher is better)")
    axes[1].legend(frameon=False)
    fig.suptitle("Extra Trees baseline: held-out primary test set", fontsize=14)
    fig.savefig(output / "baseline_metrics.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_steg_distributions(local_steg: Path, history_steg: Path, output: Path) -> None:
    variants = [("Local", local_steg), ("History-aware", history_steg)]
    targets = [("DNAEdep_keV", "DNA energy deposition (keV)"), ("TotalSB", "Total strand breaks"),
               ("TotalDSB", "Total double-strand breaks")]
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
    for row, (label, run_dir) in enumerate(variants):
        frame = pd.read_csv(run_dir / "validation_preview.csv.gz")
        for col, (target, xlabel) in enumerate(targets):
            ax = axes[row, col]
            truth = frame[target].to_numpy()
            generated = frame[f"generated_{target}"].to_numpy()
            upper = max(float(np.quantile(truth, 0.99)), float(np.quantile(generated, 0.99)), 1.0)
            bins = np.linspace(0, upper, 25) if target == "DNAEdep_keV" else np.arange(-0.5, np.ceil(upper) + 1.5)
            ax.hist(truth, bins=bins, density=True, histtype="step", linewidth=2,
                    color=COLOURS["truth"], label="Simulation")
            ax.hist(generated, bins=bins, density=True, histtype="step", linewidth=2,
                    color=COLOURS["local" if row == 0 else "history"], label="Generated")
            ax.set_title(f"{label} StEG" if col == 0 else "")
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Density" if col == 0 else "")
            ax.grid(alpha=0.15)
            if row == 0 and col == 2:
                ax.legend(frameon=False)
    fig.suptitle("StEG generated-versus-simulated distributions", fontsize=14)
    fig.savefig(output / "steg_distributions.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def plot_steg_training(local_steg: Path, history_steg: Path, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for label, run_dir, colour in [
        ("Local", local_steg, COLOURS["local"]),
        ("History-aware", history_steg, COLOURS["history"]),
    ]:
        history = pd.DataFrame(json.loads((run_dir / "history.json").read_text()))
        ax.plot(history["epoch"], history["train_loss"], color=colour, alpha=0.45,
                label=f"{label}: training")
        ax.plot(history["epoch"], history["validation_loss"], color=colour, linewidth=2,
                label=f"{label}: validation")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Diffusion noise-prediction loss")
    ax.set_title("StEG training convergence")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    fig.savefig(output / "steg_training_history.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("data/processed/ps_response.csv.gz"))
    parser.add_argument("--runs", type=Path, default=Path("runs"))
    parser.add_argument("--output", type=Path, default=Path("runs/ps_model_figures"))
    parser.add_argument("--local-steg-dir", type=Path)
    parser.add_argument("--history-steg-dir", type=Path)
    args = parser.parse_args()
    local_steg = args.local_steg_dir or args.runs / "ps_steg_local_preliminary"
    history_steg = args.history_steg_dir or args.runs / "ps_steg_history_preliminary"
    args.output.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.data)
    plot_overview(data, args.runs, local_steg, history_steg, args.output)
    plot_baseline_metrics(args.runs, args.output)
    plot_steg_distributions(local_steg, history_steg, args.output)
    plot_steg_training(local_steg, history_steg, args.output)
    print(f"Wrote figures to {args.output}")


if __name__ == "__main__":
    main()
