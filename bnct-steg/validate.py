import torch
import numpy as np
import os
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import train_test_split
from xgboost import XGBClassifier
import pickle
import matplotlib.pyplot as plt
import pandas as pd

from generate import generate_momenta
from plotting import plot_distributions, plot_2d_momenta_pairs, plot_corr_matrix
from diffusion import Diffusion



def run_validation(model, config, loader, epoch, max_batches=20, save_folder=None, plot_distributions_flag=False):
    """
    Run validation on a dataset loader for unconditional 4-momenta generation.
    Compares generated vs real distributions in real space.
    """
    if config.INPUT_SCALER_FLAG:
        with open(config.TRANSFORMER_LOC, 'rb') as f:
                    quantile_transformers = pickle.load(f)
    
    diffusion = Diffusion(config)
                    
    model.eval()

    real_all = []
    fake_all = []

    with torch.no_grad():
        for i, inputs in enumerate(loader):
            if i >= max_batches:
                break
            
            print(f'Batch {i+1}/{max_batches}', end='\r')

            inputs = inputs.to(config.DEVICE)

            # Generate fake momenta
            fake_unscaled, _ = generate_momenta(model, config, diffusion, quantile_transformers=quantile_transformers, batch_size=len(inputs))

            if config.INPUT_SCALER_FLAG:
                real_unscaled = inputs.cpu().numpy()
                for j, target in enumerate(quantile_transformers.keys()):
                    real_unscaled[:,j] = quantile_transformers[target].inverse_transform(real_unscaled[:,j].reshape(-1,1)).flatten()
                real_unscaled = torch.from_numpy(real_unscaled)
            
            else:
                real_unscaled = inputs.cpu()

            real_all.append(real_unscaled)
            fake_all.append(fake_unscaled.cpu())

    # Concatenate all batches
    real_all = torch.cat(real_all).numpy()
    fake_all = torch.cat(fake_all).numpy()

    os.makedirs(f"{save_folder}/data/", exist_ok=True)
    with open(f"{save_folder}/data/real.pkl", 'wb') as f:
        pickle.dump(real_all, f)

    with open(f"{save_folder}/data/fake.pkl", 'wb') as f:
        pickle.dump(fake_all, f)

    # Plot distributions
    if plot_distributions_flag:
        plot_distributions(real_all, fake_all, epoch, suffix="_validation", save_folder=save_folder)
        if getattr(config, "PLOT_2D_FEATURE_PAIRS", False):
            plot_2d_momenta_pairs(real_all, fake_all, epoch, save_folder=save_folder, scaled=False)
        plot_corr_matrix(real_all, fake_all, epoch, save_folder=save_folder)
        # plot_predicted_vs_true(real_all, fake_all, epoch, suffix="_validation", save_folder=save_folder)
        if config.INPUT_SCALER_FLAG:
            real_scaled_all = real_all.copy()
            fake_scaled_all = fake_all.copy()
            for j, target in enumerate(quantile_transformers.keys()):
                    real_scaled_all[:,j] = quantile_transformers[target].transform(real_scaled_all[:,j].reshape(-1,1)).flatten()
                    fake_scaled_all[:,j] = quantile_transformers[target].transform(fake_scaled_all[:,j].reshape(-1,1)).flatten()

            if getattr(config, "PLOT_2D_FEATURE_PAIRS", False):
                plot_2d_momenta_pairs(real_scaled_all, fake_scaled_all, epoch, suffix="_scaled", save_folder=save_folder, scaled=True)

        #     plot_distributions(real_scaled_all, fake_scaled_all, epoch, suffix="_scaled", save_folder=save_folder)

    # Run BDT on real-space 4-momenta
    auc = run_bdt(real_all, fake_all, epoch, config, save_folder=save_folder, plot_distributions_flag=plot_distributions_flag)
    return auc


def run_bdt(real, fake, epoch, config, save_folder=None, plot_distributions_flag=False):
    """
    Train a BDT to distinguish real vs generated 4-momenta.
    """
    X = np.concatenate([real, fake])
    y = np.concatenate([np.ones(len(real)), np.zeros(len(fake))])

    X_df = pd.DataFrame(X, columns=config.FEATURE_NAMES)

    X_train, X_test, y_train, y_test = train_test_split(
        X_df,
        y,
        test_size=0.5,
        random_state=42,
        stratify=y,
    )

    clf = XGBClassifier(max_depth=4, n_estimators=200, learning_rate=0.05)
    clf.fit(X_train, y_train)

    pred = clf.predict_proba(X_test)[:, 1]

    auc = roc_auc_score(y_test, pred)
    fpr, tpr, _ = roc_curve(y_test, pred)

    print(f"BDT AUC: {auc:.3f}")

    if plot_distributions_flag:
        pred_real = pred[y_test == 1]
        pred_fake = pred[y_test == 0]

        # fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 8), gridspec_kw={"height_ratios": [2, 1]})

        # bins = np.linspace(0, 1, 60)
        # ax1.hist(pred_real, bins=bins, density=True, alpha=0.5, label="Real", histtype="stepfilled")
        # ax1.hist(pred_fake, bins=bins, density=True, alpha=0.5, label="Generated", histtype="stepfilled")
        # ax1.set_xlabel("BDT classifier output")
        # ax1.set_ylabel("Density")
        # ax1.set_title(f"Classifier response (epoch {epoch})")
        # ax1.legend()

        # ax2.plot(fpr, tpr, label=f"AUC = {auc:.3f}")
        # ax2.plot([0, 1], [0, 1], "--", color="gray")
        # ax2.set_xlabel("False Positive Rate")
        # ax2.set_ylabel("True Positive Rate")
        # ax2.set_title("ROC Curve")
        # ax2.legend()

        # plt.tight_layout()
        # save_path = f"{save_folder}/bdt_validation_epoch_{epoch}.pdf" if save_folder else f"bdt_validation_epoch_{epoch}.pdf"
        # plt.savefig(save_path)
        # plt.close()

        # Saves plots separately so ROC curve shape is square
        fig, ax1 = plt.subplots(figsize=(6, 6))

        bins = np.linspace(0, 1, 60)
        # bins = 60
        ax1.hist(pred_real, bins=bins, density=True, alpha=0.5,
                label="Real", histtype="stepfilled")
        ax1.hist(pred_fake, bins=bins, density=True, alpha=0.5,
                label="Generated", histtype="stepfilled")

        ax1.set_xlabel("BDT classifier output")
        ax1.set_ylabel("Density")
        ax1.set_title(f"Classifier response (epoch {epoch})")

        ax2 = ax1.twinx()
        ax2.plot(fpr, tpr, color='tab:red', label=f"AUC = {auc:.3f}")
        ax2.plot([0, 1], [0, 1], "--", color="gray")
        ax2.set_ylabel("True Positive Rate")

        ax1.set_ylim(0, None)
        ax2.set_ylim(0, 1)

        lines_1, labels_1 = ax1.get_legend_handles_labels()
        lines_2, labels_2 = ax2.get_legend_handles_labels()
        ax1.legend(lines_1 + lines_2, labels_1 + labels_2, loc="upper left")

        plt.xlim(0, 1)

        fig.tight_layout()
        fig.savefig(f"{save_folder}/bdt_response_epoch_{epoch}.pdf")
        plt.close(fig)

        try:
            import shap
        except ImportError:
            print("SHAP is not installed; skipping SHAP plots.")
        else:
            os.makedirs(f"{save_folder}/shap/", exist_ok=True)

            explainer = shap.TreeExplainer(clf)
            shap_values = explainer.shap_values(X_test)

            plt.figure(figsize=(10, 8))
            shap.summary_plot(shap_values, X_test, show=False)
            plt.tight_layout()
            plt.savefig(f"{save_folder}/shap/shap_summary_epoch_{epoch}.pdf")
            plt.close()

            plt.figure(figsize=(8, 6))
            shap.summary_plot(shap_values, X_test, plot_type="bar", show=False)
            plt.tight_layout()
            plt.savefig(f"{save_folder}/shap/shap_bar_epoch_{epoch}.pdf")
            plt.close()

            importance = np.abs(shap_values).mean(axis=0)
            top_feature = X_test.columns[np.argmax(importance)]

            plt.figure(figsize=(8, 6))
            shap.dependence_plot(
                top_feature,
                shap_values,
                X_test,
                interaction_index="auto",
                show=False,
            )
            plt.tight_layout()
            plt.savefig(f"{save_folder}/shap/shap_dependence_{top_feature}_epoch_{epoch}.pdf")
            plt.close()

            importance_df = (
                pd.DataFrame({
                    "feature": X_test.columns,
                    "importance": np.abs(shap_values).mean(axis=0),
                })
                .sort_values("importance", ascending=False)
            )
            importance_df.to_csv(f"{save_folder}/shap/shap_importance_epoch_{epoch}.csv", index=False)

    return auc
