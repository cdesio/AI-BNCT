import os
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import torch

from config import Config
from colours import cmap_cyan_orange

config = Config()

# 4-momenta feature names
momenta_names = config.FEATURE_NAMES

def get_plot_range(values, nsigma=3):
    """
    Get plotting range excluding extreme outliers.
    """

    mean = np.mean(values)
    std = np.std(values)

    lower = max(np.min(values), mean - nsigma * std)
    upper = min(np.max(values), mean + nsigma * std)

    return lower, upper

# def get_plot_range(values, lower_q=0.5, upper_q=99.5):

#     return np.percentile(
#         values,
#         [lower_q, upper_q]
#     )

def plot_distributions(real, fake, epoch, suffix="", save_folder=None):
    """
    Plot histogram distributions of real vs fake 4-momenta.
    """

    pdf_path = (
        f"{save_folder}/distributions_epoch_{epoch}{suffix}.pdf"
        if save_folder else
        f"distributions_epoch_{epoch}{suffix}.pdf"
    )

    with PdfPages(pdf_path) as pdf:

        for i, name in enumerate(momenta_names):

            real_feat = real[:, i]
            fake_feat = fake[:, i]

            # Convert if torch
            real_feat = real_feat.cpu().numpy() if torch.is_tensor(real_feat) else real_feat
            fake_feat = fake_feat.cpu().numpy() if torch.is_tensor(fake_feat) else fake_feat

            plt.figure(figsize=(6,4))

            lower, upper = get_plot_range(real_feat)

            bins = np.histogram_bin_edges(real_feat, bins=100, range=[lower, upper])

            plt.hist(real_feat, bins=bins, density=True, range=[lower, upper], alpha=0.5, label="Real")
            plt.hist(fake_feat, bins=bins, density=True, range=[lower, upper], alpha=0.5, label="Generated")

            plt.xlabel(name)
            plt.ylabel("Density")
            plt.title(f"{name} {suffix}")
            plt.legend()
            plt.tight_layout()

            pdf.savefig()
            plt.close()

    print(f"Saved distributions to {pdf_path}")


def plot_loss(loss, epoch, save_folder=None):
    plt.figure(figsize=(6,4))
    plt.plot(np.linspace(0, epoch, len(loss)), loss)
    plt.xlabel("Epoch")
    plt.ylabel("Loss")
    plt.title("Loss over time")
    plt.tight_layout()
    save_path = f"{save_folder}/loss_vs_epoch_{epoch}.pdf" if save_folder else f"loss_vs_epoch_{epoch}.pdf"
    plt.savefig(save_path)
    plt.close()

def plot_predicted_vs_true(real, fake, epoch, suffix="", save_folder=None, bins=100): # This I'm pretty sure is useless as we just see the 1d histograms swept out by 90 degrees but will leave in just in case
    """
    For each 4-momenta feature, plot a 2D histogram where:
        x-axis = generated/fake values
        y-axis = real values
    Shows where the model over/under-predicts.
    """
    for i, name in enumerate(momenta_names):
        real_feat = real[:, i]
        fake_feat = fake[:, i]

        real_np = real_feat.cpu().numpy() if isinstance(real_feat, torch.Tensor) else real_feat
        fake_np = fake_feat.cpu().numpy() if isinstance(fake_feat, torch.Tensor) else fake_feat

        # Compute 2D histogram
        H, xedges, yedges = np.histogram2d(fake_np, real_np, bins=bins)

        # Plot heatmap
        plt.figure(figsize=(6,6))
        plt.imshow(H.T, origin='lower', aspect='auto',
                   extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                   cmap='viridis')
        plt.colorbar(label='Counts')
        # Ideal diagonal
        plt.plot([fake_np.min(), fake_np.max()],
                 [fake_np.min(), fake_np.max()], 'r--', linewidth=1)
        plt.xlabel(f"Generated {name}")
        plt.ylabel(f"Real {name}")
        plt.title(f"Generated vs Real {name} {suffix}")
        plt.tight_layout()

        save_path = f"{save_folder}/{name}_pred_vs_true_epoch_{epoch}{suffix}.pdf" if save_folder else f"{name}_pred_vs_true_epoch_{epoch}{suffix}.pdf"
        plt.savefig(save_path)
        plt.close()


def plot_2d_momenta_pairs(real, fake, epoch, suffix="", save_folder=None, scaled=False, bins=100):
    """
    Same as before, but saves ALL plots into a single multi-page PDF.
    """

    particle_names = ["H1", "H2", "L1", "L2"]

    pdf_path = (
        f"{save_folder}/momenta_2d_histograms_epoch_{epoch}{suffix}.pdf"
        if save_folder else
        f"momenta_2d_histograms_epoch_{epoch}{suffix}.pdf"
    )

    with PdfPages(pdf_path) as pdf:

        for p_idx, pname in enumerate(particle_names):
            base = p_idx * 4

            for i in range(4):
                for j in range(i+1, 4):

                    idx_i = base + i
                    idx_j = base + j

                    real_a = real[:, idx_i]
                    real_b = real[:, idx_j]
                    fake_a = fake[:, idx_i]
                    fake_b = fake[:, idx_j]

                    # Convert to numpy
                    real_a = real_a.cpu().numpy() if torch.is_tensor(real_a) else real_a
                    real_b = real_b.cpu().numpy() if torch.is_tensor(real_b) else real_b
                    fake_a = fake_a.cpu().numpy() if torch.is_tensor(fake_a) else fake_a
                    fake_b = fake_b.cpu().numpy() if torch.is_tensor(fake_b) else fake_b

                    if not scaled:

                        # Combine real and generated so ranges are identical
                        all_a = np.concatenate([real_a, fake_a])
                        all_b = np.concatenate([real_b, fake_b])

                        x_min, x_max = get_plot_range(all_a)
                        y_min, y_max = get_plot_range(all_b)

                    else:
                        x_min, x_max = None, None
                        y_min, y_max = None, None


                    hist_range = None
                    if not scaled:
                        hist_range = [
                            [x_min, x_max],
                            [y_min, y_max]
                        ]

                    # Shared bins
                    H_real, xedges, yedges = np.histogram2d(real_a, real_b, bins=bins, range=hist_range)
                    H_fake, _, _ = np.histogram2d(fake_a, fake_b, bins=[xedges, yedges])

                    # Plot
                    fig, ax = plt.subplots(2, 1, figsize=(6, 10))

                    im0 = ax[0].imshow(H_real.T, origin='lower', aspect='auto',
                                       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                                       cmap='viridis')
                    # _, _, _, im0 = ax[0].hist2d(real_a, real_b, bins=bins, range=hist_range, cmap="viridis")
                    fig.colorbar(im0, ax=ax[0])
                    ax[0].set_title(f"{pname} Real: {momenta_names[idx_i]} vs {momenta_names[idx_j]}")
                    ax[0].set_xlabel(momenta_names[idx_i])
                    ax[0].set_ylabel(momenta_names[idx_j])

                    im1 = ax[1].imshow(H_fake.T, origin='lower', aspect='auto',
                                       extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                                       cmap='viridis')
                    # _, _, _, im1 = ax[1].hist2d(fake_a, fake_b, bins=bins, range=hist_range, cmap="viridis")
                    fig.colorbar(im1, ax=ax[1])
                    ax[1].set_title(f"{pname} Gen: {momenta_names[idx_i]} vs {momenta_names[idx_j]}")
                    ax[1].set_xlabel(momenta_names[idx_i])
                    ax[1].set_ylabel(momenta_names[idx_j])

                    if not scaled:
                        ax[0].set_xlim(x_min, x_max)
                        ax[0].set_ylim(y_min, y_max)

                        ax[1].set_xlim(x_min, x_max)
                        ax[1].set_ylim(y_min, y_max)

                    plt.tight_layout()

                    pdf.savefig(fig)
                    plt.close()

    print(f"Saved 2d histograms to {pdf_path}")

def plot_logvar_vs_timestep(config, t, logvar, epoch, suffix="", save_folder=None):
    mean_logvar = np.zeros(config.NOISE_STEPS)

    for i in range(config.NOISE_STEPS):
        mask = (t == i)
        if np.any(mask):
            mean_logvar[i] = logvar[mask].mean()

    plt.plot(range(config.NOISE_STEPS), mean_logvar)
    plt.xlabel("timestep")
    plt.ylabel("mean logvar")
    plt.title(f"Logvar vs Timestep (epoch {epoch})")
    if save_folder is not None:
        plt.savefig(f"{save_folder}/logvar_vs_t_epoch_{epoch}{suffix}.png")
        plt.close()
    else:
        plt.savefig(f"logvar_vs_t_epoch_{epoch}{suffix}.png")
        plt.close()

def plot_BDT_AUC_vs_epoch(BDT_train, BDT_val, epoch, checkpoint_interval, save_folder=None):
    plt.figure(figsize=(6,4))
    plt.plot(np.linspace(checkpoint_interval, (epoch), len(BDT_train)), BDT_train, label="Train")
    plt.plot(np.linspace(checkpoint_interval, (epoch), len(BDT_val)), BDT_val, label='Val')
    plt.axhline(np.min(BDT_val), ls="--", color="r", label=f"Min val AUC: {np.min(BDT_val)}")
    plt.xlabel("Epoch")
    plt.ylabel("BDT AUC")
    plt.title("BDT AUC over time")
    plt.legend()
    plt.tight_layout()
    save_path = f"{save_folder}/BDT_AUC_vs_epoch_{epoch}.pdf" if save_folder else f"loss_vs_epoch_{epoch}.pdf"
    plt.savefig(save_path)
    plt.close()

def plot_corr_matrix(real, fake, epoch, suffix="", save_folder=None):

    if torch.is_tensor(real):
        real = real.cpu().numpy()
    if torch.is_tensor(fake):
        fake = fake.cpu().numpy()

    corr_real = np.corrcoef(real, rowvar=False)
    corr_fake = np.corrcoef(fake, rowvar=False)

    diff = corr_fake - corr_real

    pdf_path = (
        f"{save_folder}/correlation_matrices_epoch_{epoch}{suffix}.pdf"
        if save_folder else
        f"correlation_matrices_epoch_{epoch}{suffix}.pdf"
    )

    with PdfPages(pdf_path) as pdf:
        # Real
        fig, ax = plt.subplots(figsize=(8,6))
        im = ax.imshow(corr_real, vmin=-1, vmax=1, cmap=cmap_cyan_orange)
        plt.colorbar(im, ax=ax)

        ax.set_xticks(range(len(momenta_names)))
        ax.set_yticks(range(len(momenta_names)))
        ax.set_xticklabels(momenta_names, rotation=90)
        ax.set_yticklabels(momenta_names)

        ax.set_title("Correlation Matrix - Real")

        plt.tight_layout()

        pdf.savefig()
        plt.close()

        # Fake
        fig, ax = plt.subplots(figsize=(8,6))
        im = ax.imshow(corr_fake, vmin=-1, vmax=1, cmap=cmap_cyan_orange)
        plt.colorbar(im, ax=ax)

        ax.set_xticks(range(len(momenta_names)))
        ax.set_yticks(range(len(momenta_names)))
        ax.set_xticklabels(momenta_names, rotation=90)
        ax.set_yticklabels(momenta_names)

        ax.set_title("Correlation Matrix - Fake")

        plt.tight_layout()

        pdf.savefig()
        plt.close()

        # Difference
        fig, ax = plt.subplots(figsize=(8,6))
        im = ax.imshow(diff, vmin=-1, vmax=1, cmap=cmap_cyan_orange)
        plt.colorbar(im, ax=ax)

        ax.set_xticks(range(len(momenta_names)))
        ax.set_yticks(range(len(momenta_names)))
        ax.set_xticklabels(momenta_names, rotation=90)
        ax.set_yticklabels(momenta_names)

        ax.set_title("Correlation Difference (Fake - Real)")

        plt.tight_layout()

        pdf.savefig()
        plt.close()
