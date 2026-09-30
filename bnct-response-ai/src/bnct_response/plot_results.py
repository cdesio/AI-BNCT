"""Create meeting-ready held-out baseline result plots."""

from __future__ import annotations

import argparse
from pathlib import Path

import joblib
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)

from .baseline import inverse_transport


CASE_LABELS = {
    "alpha_1p47": r"Alpha, 1.47 MeV",
    "alpha_1p78": r"Alpha, 1.78 MeV",
    "lithium_0p84": r"Lithium, 0.84 MeV",
    "lithium_1p01": r"Lithium, 1.01 MeV",
}
COLORS = {
    "observed": "#202124",
    "direct": "#0072B2",
    "oracle": "#D55E00",
    "transport": "#009E73",
}


def _style() -> None:
    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "legend.fontsize": 8,
        "figure.titlesize": 13,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "grid.alpha": 0.18,
        "grid.linewidth": 0.7,
        "savefig.bbox": "tight",
    })


def _save(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=220, facecolor="white")
    fig.savefig(output_dir / f"{stem}.pdf", facecolor="white")
    plt.close(fig)


def _binned(frame: pd.DataFrame, columns: list[str], bins: int = 8) -> pd.DataFrame:
    edges = np.linspace(frame["Distance_um"].min(), frame["Distance_um"].max(), bins + 1)
    work = frame.copy()
    work["distance_bin"] = pd.cut(
        work["Distance_um"], edges, include_lowest=True, duplicates="drop"
    )
    grouped = work.groupby("distance_bin", observed=True)
    rows = []
    for _, group in grouped:
        row = {"Distance_um": group["Distance_um"].mean(), "n": len(group)}
        for column in columns:
            row[column] = group[column].mean()
            row[f"{column}_se"] = group[column].std(ddof=1) / np.sqrt(len(group))
        rows.append(row)
    return pd.DataFrame(rows)


def _prepare_predictions(test: pd.DataFrame, artifact: dict) -> pd.DataFrame:
    result = test.copy()
    condition = result[artifact["condition_columns"]]
    transport = inverse_transport(artifact["transport"].predict(condition))
    for index, target in enumerate(artifact["transport_targets"]):
        result[f"pred_{target}"] = transport[:, index]

    result["pred_AnyDamage_direct"] = artifact["direct_damage_any"].predict_proba(condition)[:, 1]
    result["pred_AnyDSB_direct"] = artifact["direct_damage_dsb"].predict_proba(condition)[:, 1]
    local = result[artifact["damage_input_columns"]]
    result["pred_AnyDamage_oracle"] = artifact["damage_any"].predict_proba(local)[:, 1]
    result["pred_AnyDSB_oracle"] = artifact["damage_dsb"].predict_proba(local)[:, 1]
    return result


def plot_transport_distance(frame: pd.DataFrame, output_dir: Path) -> None:
    cases = list(CASE_LABELS)
    targets = [
        ("EntryEnergy_MeV", "Entry energy (MeV)"),
        ("Edep_keV", "Deposited energy (keV)"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(13.2, 6.1), sharex="col")
    for column, case in enumerate(cases):
        case_frame = frame[frame["case_id"] == case]
        for row, (target, ylabel) in enumerate(targets):
            summary = _binned(case_frame, [target, f"pred_{target}"])
            ax = axes[row, column]
            ax.errorbar(
                summary["Distance_um"], summary[target], yerr=summary[f"{target}_se"],
                fmt="o", ms=4, capsize=2, color=COLORS["observed"], label="Simulation",
            )
            ax.plot(
                summary["Distance_um"], summary[f"pred_{target}"], marker="s", ms=3,
                color=COLORS["transport"], lw=1.6, label="Extra Trees",
            )
            if row == 0:
                ax.set_title(CASE_LABELS[case])
            if column == 0:
                ax.set_ylabel(ylabel)
            if row == 1:
                ax.set_xlabel("Distance along initial direction (µm)")
    axes[0, 0].legend(frameon=False, loc="best")
    fig.suptitle("Held-out transport response versus distance")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, output_dir, "transport_vs_distance")


def plot_damage_distance(frame: pd.DataFrame, output_dir: Path) -> None:
    outcomes = [
        ("AnyDamage", "Any damage"),
        ("AnyDSB", "Any DSB"),
    ]
    fig, axes = plt.subplots(4, 2, figsize=(10.5, 11.8), sharey=True)
    for row, case in enumerate(CASE_LABELS):
        case_frame = frame[frame["case_id"] == case]
        columns = [
            "AnyDamage", "AnyDSB", "pred_AnyDamage_direct", "pred_AnyDSB_direct",
            "pred_AnyDamage_oracle", "pred_AnyDSB_oracle",
        ]
        summary = _binned(case_frame, columns)
        for column, (observed, label) in enumerate(outcomes):
            ax = axes[row, column]
            ax.errorbar(
                summary["Distance_um"], summary[observed],
                yerr=summary[f"{observed}_se"], fmt="o", ms=4, capsize=2,
                color=COLORS["observed"], label="Observed fraction",
            )
            ax.plot(
                summary["Distance_um"], summary[f"pred_{observed}_direct"],
                color=COLORS["direct"], lw=1.8,
                label="Damage model: initial conditions + distance",
            )
            ax.plot(
                summary["Distance_um"], summary[f"pred_{observed}_oracle"],
                color=COLORS["oracle"], lw=1.4, ls="--", alpha=0.85,
                label="Damage model + simulated local tracking",
            )
            if row == 0:
                ax.set_title(label)
            if column == 0:
                ax.set_ylabel(f"{CASE_LABELS[case]}\nObserved fraction / probability")
            if row == len(CASE_LABELS) - 1:
                ax.set_xlabel("Distance along initial direction (µm)")
            ax.set_ylim(-0.04, 1.04)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, ncol=3, loc="lower center")
    fig.suptitle("Damage probability per encountered voxel versus distance")
    fig.tight_layout(rect=(0, 0.055, 1, 0.97))
    _save(fig, output_dir, "damage_vs_distance")


def plot_classifier_diagnostics(frame: pd.DataFrame, output_dir: Path) -> None:
    outcomes = [
        ("AnyDamage", "Any damage"),
        ("AnyDSB", "Any DSB"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.0))
    for row, (target, label) in enumerate(outcomes):
        truth = frame[target].to_numpy()
        probability = frame[f"pred_{target}_direct"].to_numpy()
        false_positive, true_positive, _ = roc_curve(truth, probability)
        precision, recall, _ = precision_recall_curve(truth, probability)
        observed, predicted = calibration_curve(truth, probability, n_bins=8, strategy="quantile")
        auc = roc_auc_score(truth, probability)
        ap = average_precision_score(truth, probability)

        axes[row, 0].plot(false_positive, true_positive, color=COLORS["direct"], lw=2)
        axes[row, 0].plot([0, 1], [0, 1], color="#777777", ls=":")
        axes[row, 0].set_title(f"{label}: ROC-AUC = {auc:.3f}")
        axes[row, 0].set_xlabel("False-positive rate")
        axes[row, 0].set_ylabel("True-positive rate")

        axes[row, 1].plot(recall, precision, color=COLORS["direct"], lw=2)
        axes[row, 1].axhline(truth.mean(), color="#777777", ls=":")
        axes[row, 1].set_title(f"{label}: AP = {ap:.3f}")
        axes[row, 1].set_xlabel("Recall")
        axes[row, 1].set_ylabel("Precision")

        axes[row, 2].plot(predicted, observed, "o-", color=COLORS["direct"], lw=1.7)
        axes[row, 2].plot([0, 1], [0, 1], color="#777777", ls=":")
        axes[row, 2].set_title(f"{label}: calibration")
        axes[row, 2].set_xlabel("Mean predicted probability")
        axes[row, 2].set_ylabel("Observed fraction")
        for ax in axes[row]:
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1.02)
    fig.suptitle("Direct-distance damage model on held-out primaries")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, output_dir, "damage_classifier_diagnostics")


def plot_transport_parity(frame: pd.DataFrame, output_dir: Path) -> None:
    targets = [
        ("EntryEnergy_MeV", "Entry energy (MeV)"),
        ("ExitEnergy_MeV", "Exit energy (MeV)"),
        ("Edep_keV", "Deposited energy (keV)"),
        ("TrackLength_nm", "Track length (nm)"),
        ("StepCount", "Step count"),
        ("MeanLET_keV_um", "Mean LET (keV/µm)"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(10.8, 7.0))
    for ax, (target, label) in zip(axes.flat, targets):
        observed = frame[target].to_numpy()
        predicted = frame[f"pred_{target}"].to_numpy()
        limit = float(max(np.quantile(observed, 0.995), np.quantile(predicted, 0.995)))
        ax.hexbin(observed, predicted, gridsize=32, mincnt=1, cmap="viridis")
        ax.plot([0, limit], [0, limit], color="#D55E00", ls="--", lw=1.2)
        ax.set_xlim(0, limit)
        ax.set_ylim(0, limit)
        ax.set_xlabel(f"Simulated {label}")
        ax.set_ylabel(f"Predicted {label}")
    fig.suptitle("Transport prediction parity on held-out primaries")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    _save(fig, output_dir, "transport_parity")


def plot_range_support(ranges: pd.DataFrame, artifact: dict, output_dir: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2), sharey=False)
    rng = np.random.default_rng(20260924)
    for ax, particle, title in zip(axes, ["alpha", "lithium"], ["Alpha", "Lithium"]):
        particle_frame = ranges[ranges["Particle"] == particle].copy()
        energies = sorted(particle_frame["InitialEnergy_MeV"].unique())
        for index, energy in enumerate(energies):
            observed = particle_frame[particle_frame["InitialEnergy_MeV"] == energy]
            jitter = rng.normal(0, 0.025, len(observed))
            ax.scatter(
                np.full(len(observed), index) + jitter,
                observed["ObservedRange_um"],
                s=25, color=COLORS["observed"], alpha=0.55,
                label="Held-out primaries" if index == 0 else None,
            )
            model_input = observed.iloc[[0]][
                ["Particle", "InitialEnergy_MeV", "InitialDirX", "InitialDirY", "InitialDirZ"]
            ]
            q05 = artifact["range_models"]["q05"].predict(model_input)[0]
            q50 = artifact["range_models"]["q50"].predict(model_input)[0]
            q95 = artifact["range_models"]["q95"].predict(model_input)[0]
            ax.errorbar(
                index, q50, yerr=[[q50 - q05], [q95 - q50]],
                fmt="D", ms=6, capsize=5, color=COLORS["direct"], lw=1.8,
                label="Predicted median and 5–95% interval" if index == 0 else None,
            )
        ax.set_xticks(range(len(energies)), [f"{energy:.2f}" for energy in energies])
        ax.set_xlabel("Initial energy (MeV)")
        ax.set_ylabel("Maximum observed distance (µm)")
        ax.set_title(title)
        ax.legend(frameon=False, loc="best")
    fig.suptitle("Particle range support on held-out primaries")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    _save(fig, output_dir, "range_support")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--model-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    _style()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.data_dir / "response.csv.gz")
    ranges = pd.read_csv(args.data_dir / "range.csv.gz")
    test = frame[frame["split"] == "test"].copy()
    range_test = ranges[ranges["split"] == "test"].copy()
    artifact = joblib.load(args.model_dir / "models.joblib")
    result = _prepare_predictions(test, artifact)
    plot_range_support(range_test, artifact, args.output_dir)
    plot_transport_distance(result, args.output_dir)
    plot_damage_distance(result, args.output_dir)
    plot_classifier_diagnostics(result, args.output_dir)
    plot_transport_parity(result, args.output_dir)
    print(f"Wrote 5 PNG and 5 PDF figures to {args.output_dir}")


if __name__ == "__main__":
    main()
