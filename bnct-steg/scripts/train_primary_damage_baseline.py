import argparse
import json
import os

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.calibration import calibration_curve
from sklearn.compose import ColumnTransformer
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score,
    brier_score_loss,
    precision_recall_curve,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GroupShuffleSplit
from sklearn.pipeline import Pipeline
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from xgboost import XGBClassifier


CATEGORICAL_FEATURES = [
    "ParticleType",
    "Compartment",
]

NUMERIC_FEATURES = [
    "InitialEnergy_MeV",
    "SourceDistance_um",
    "PosX_um",
    "PosY_um",
    "PosZ_um",
    "MomX",
    "MomY",
    "MomZ",
    "StopPosX_um",
    "StopPosY_um",
    "StopPosZ_um",
    "TraLen_cell_um",
    "TraLen_chro_um",
    "StartDistToOrigin_um",
    "StopDistToOrigin_um",
    "DeltaRadius_um",
    "TrackLength_um",
    "DirectionX",
    "DirectionY",
    "DirectionZ",
    "RadialDirectionCosine",
    "ClosestApproachT",
    "ClosestDistanceToOrigin_um",
]

DAMAGE_TREE_FEATURES = [
    "RecordedEnergyDeposited_eV",
    "LogRecordedEnergyDeposited_eV",
    "DamageRows",
    "SumBaseDamage",
    "SumStrandDamage",
    "SumDirectBreaks",
    "SumIndirectBreaks",
]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Train primary-level BNCT damage prediction baselines."
    )
    parser.add_argument("--input", default="data/bnct_primary_summary.csv")
    parser.add_argument("--output-dir", default="runs/primary_damage_baseline")
    parser.add_argument("--test-size", type=float, default=0.25)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument(
        "--include-damage-tree-features",
        action="store_true",
        help=(
            "Include damage-tree aggregate columns. This is useful for diagnostics, "
            "but it leaks damage information for the current ROOT files."
        ),
    )
    return parser.parse_args()


def make_preprocessor(numeric_features):
    numeric_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
        ]
    )
    categorical_pipe = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="most_frequent")),
            ("onehot", OneHotEncoder(handle_unknown="ignore", sparse_output=False)),
        ]
    )
    return ColumnTransformer(
        transformers=[
            ("numeric", numeric_pipe, numeric_features),
            ("categorical", categorical_pipe, CATEGORICAL_FEATURES),
        ]
    )


def get_feature_names(preprocessor, numeric_features):
    names = list(numeric_features)
    onehot = preprocessor.named_transformers_["categorical"].named_steps["onehot"]
    names.extend(onehot.get_feature_names_out(CATEGORICAL_FEATURES).tolist())
    return names


def build_models(scale_pos_weight, numeric_features, random_state):
    return {
        "logistic_regression": Pipeline(
            steps=[
                ("preprocessor", make_preprocessor(numeric_features)),
                (
                    "classifier",
                    LogisticRegression(
                        class_weight="balanced",
                        max_iter=2000,
                        random_state=random_state,
                    ),
                ),
            ]
        ),
        "random_forest": Pipeline(
            steps=[
                ("preprocessor", make_preprocessor(numeric_features)),
                (
                    "classifier",
                    RandomForestClassifier(
                        n_estimators=400,
                        max_depth=6,
                        min_samples_leaf=10,
                        class_weight="balanced_subsample",
                        random_state=random_state,
                        n_jobs=-1,
                    ),
                ),
            ]
        ),
        "xgboost": Pipeline(
            steps=[
                ("preprocessor", make_preprocessor(numeric_features)),
                (
                    "classifier",
                    XGBClassifier(
                        n_estimators=300,
                        max_depth=3,
                        learning_rate=0.03,
                        subsample=0.9,
                        colsample_bytree=0.9,
                        eval_metric="logloss",
                        scale_pos_weight=scale_pos_weight,
                        random_state=random_state,
                    ),
                ),
            ]
        ),
    }


def evaluate_model(name, model, X_train, y_train, X_test, y_test, output_dir, numeric_features):
    model.fit(X_train, y_train)
    proba = model.predict_proba(X_test)[:, 1]

    metrics = {
        "roc_auc": float(roc_auc_score(y_test, proba)),
        "pr_auc": float(average_precision_score(y_test, proba)),
        "brier_score": float(brier_score_loss(y_test, proba)),
        "test_damage_fraction": float(np.mean(y_test)),
        "mean_predicted_probability": float(np.mean(proba)),
    }

    fpr, tpr, _ = roc_curve(y_test, proba)
    precision, recall, _ = precision_recall_curve(y_test, proba)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].plot(fpr, tpr, label=f"ROC AUC = {metrics['roc_auc']:.3f}")
    axes[0].plot([0, 1], [0, 1], "--", color="gray")
    axes[0].set_xlabel("False positive rate")
    axes[0].set_ylabel("True positive rate")
    axes[0].legend()

    axes[1].plot(recall, precision, label=f"PR AUC = {metrics['pr_auc']:.3f}")
    axes[1].set_xlabel("Recall")
    axes[1].set_ylabel("Precision")
    axes[1].legend()

    frac_pos, mean_pred = calibration_curve(y_test, proba, n_bins=10, strategy="quantile")
    axes[2].plot(mean_pred, frac_pos, marker="o")
    axes[2].plot([0, 1], [0, 1], "--", color="gray")
    axes[2].set_xlabel("Mean predicted probability")
    axes[2].set_ylabel("Observed damage fraction")
    axes[2].set_title(f"Brier = {metrics['brier_score']:.3f}")

    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{name}_diagnostics.pdf"))
    plt.close(fig)

    predictions = X_test[["RootFile", "ParticleType", "InitialEnergy_MeV", "SourceDistance_um", "Compartment", "ClosestDistanceToOrigin_um"]].copy()
    predictions["HasDamage"] = y_test.to_numpy()
    predictions["PredictedDamageProbability"] = proba
    predictions.to_csv(os.path.join(output_dir, f"{name}_test_predictions.csv"), index=False)

    write_feature_importance(name, model, output_dir, numeric_features)
    return metrics


def write_feature_importance(name, model, output_dir, numeric_features):
    preprocessor = model.named_steps["preprocessor"]
    classifier = model.named_steps["classifier"]
    feature_names = get_feature_names(preprocessor, numeric_features)

    if hasattr(classifier, "feature_importances_"):
        importance = classifier.feature_importances_
    elif hasattr(classifier, "coef_"):
        importance = np.abs(classifier.coef_[0])
    else:
        return

    importance_df = (
        pd.DataFrame({"feature": feature_names, "importance": importance})
        .sort_values("importance", ascending=False)
        .reset_index(drop=True)
    )
    importance_df.to_csv(os.path.join(output_dir, f"{name}_feature_importance.csv"), index=False)

    top = importance_df.head(20).iloc[::-1]
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.barh(top["feature"], top["importance"])
    ax.set_xlabel("Importance")
    ax.set_title(f"{name} top features")
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{name}_feature_importance.pdf"))
    plt.close(fig)


def plot_radial_damage(predictions_path, output_dir, model_name):
    df = pd.read_csv(predictions_path)
    bins = np.linspace(0, np.nanpercentile(df["ClosestDistanceToOrigin_um"], 99), 9)
    if len(np.unique(bins)) < 3:
        return
    df["RadiusBin"] = pd.cut(df["ClosestDistanceToOrigin_um"], bins=bins, include_lowest=True)

    radial = (
        df.groupby("RadiusBin", observed=True)
        .agg(
            Primaries=("HasDamage", "size"),
            ObservedDamageFraction=("HasDamage", "mean"),
            PredictedDamageProbability=("PredictedDamageProbability", "mean"),
            RadiusMid_um=("ClosestDistanceToOrigin_um", "mean"),
        )
        .reset_index()
    )
    radial.to_csv(os.path.join(output_dir, f"{model_name}_radial_damage.csv"), index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(radial["RadiusMid_um"], radial["ObservedDamageFraction"], marker="o", label="Observed")
    ax.plot(
        radial["RadiusMid_um"],
        radial["PredictedDamageProbability"],
        marker="o",
        label="Predicted",
    )
    ax.set_xlabel("Closest distance to origin (um)")
    ax.set_ylabel("Damage probability")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(output_dir, f"{model_name}_radial_damage.pdf"))
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    df = pd.read_csv(args.input)
    numeric_features = list(NUMERIC_FEATURES)
    feature_set = "primary_geometry_only"
    if args.include_damage_tree_features:
        numeric_features.extend(DAMAGE_TREE_FEATURES)
        feature_set = "primary_geometry_plus_damage_tree_aggregates"

    required = [*CATEGORICAL_FEATURES, *numeric_features, "HasDamage", "RootFile"]
    missing = [column for column in required if column not in df.columns]
    if missing:
        raise KeyError(f"Input is missing required columns: {missing}")

    df = df[required].replace([np.inf, -np.inf], np.nan)
    y = df["HasDamage"].astype(int)
    X = df.drop(columns=["HasDamage"])

    splitter = GroupShuffleSplit(
        n_splits=1,
        test_size=args.test_size,
        random_state=args.random_state,
    )
    train_idx, test_idx = next(splitter.split(X, y, groups=X["RootFile"]))
    X_train, X_test = X.iloc[train_idx].copy(), X.iloc[test_idx].copy()
    y_train, y_test = y.iloc[train_idx].copy(), y.iloc[test_idx].copy()

    positives = max(int(y_train.sum()), 1)
    negatives = max(int((1 - y_train).sum()), 1)
    models = build_models(
        scale_pos_weight=negatives / positives,
        numeric_features=numeric_features,
        random_state=args.random_state,
    )

    metrics = {
        "feature_set": feature_set,
        "included_numeric_features": numeric_features,
        "rows": int(len(df)),
        "train_rows": int(len(X_train)),
        "test_rows": int(len(X_test)),
        "train_damage_fraction": float(y_train.mean()),
        "test_damage_fraction": float(y_test.mean()),
        "models": {},
    }

    for name, model in models.items():
        model_metrics = evaluate_model(
            name,
            model,
            X_train,
            y_train,
            X_test,
            y_test,
            args.output_dir,
            numeric_features,
        )
        metrics["models"][name] = model_metrics
        plot_radial_damage(
            os.path.join(args.output_dir, f"{name}_test_predictions.csv"),
            args.output_dir,
            name,
        )

    with open(os.path.join(args.output_dir, "metrics.json"), "w") as f:
        json.dump(metrics, f, indent=4)

    print(json.dumps(metrics, indent=4))
    print(f"Output: {args.output_dir}")


if __name__ == "__main__":
    main()
