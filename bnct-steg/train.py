import os
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))
import json
import time
from dataclasses import asdict

import numpy as np
import torch
from torch.utils.data import DataLoader
from tqdm import tqdm

from config import make_config
from dataset import BNCTCsvDamageDataset, resolve_event_splits
from diffusion import Diffusion
from model import MomentaDiffusionModel
from plotting import plot_BDT_AUC_vs_epoch, plot_loss
from run_paths import create_train_run_dir, write_run_metadata
from validate import run_validation


def main():
    time0 = time.time()
    config = make_config()
    total_events, train_split, val_split = resolve_event_splits(config)

    save_folder = create_train_run_dir(config)
    save_folder_loss = f"{save_folder}/loss_over_time"
    save_folder_bdt_auc = f"{save_folder}/BDT_AUC_vs_time"
    os.makedirs(save_folder, exist_ok=True)
    os.makedirs(save_folder_loss, exist_ok=True)
    os.makedirs(save_folder_bdt_auc, exist_ok=True)

    with open(f"{save_folder}/config.json", "w") as f:
        json.dump(asdict(config), f, indent=4)
    write_run_metadata(
        f"{save_folder}/run_metadata.json",
        {
            "run_dir": save_folder,
            "kind": "train",
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "total_events": total_events,
            "train_split": train_split,
            "val_split": val_split,
        },
    )

    train_dataset = BNCTCsvDamageDataset(
        config,
        start=0,
        stop=train_split,
        input_scaler_flag=config.INPUT_SCALER_FLAG,
    )
    val_dataset = BNCTCsvDamageDataset(
        config,
        start=train_split,
        stop=val_split,
        input_scaler_flag=config.INPUT_SCALER_FLAG,
    )

    train_loader = DataLoader(train_dataset, batch_size=config.BATCH_SIZE)
    val_loader = DataLoader(val_dataset, batch_size=config.BATCH_SIZE)

    num_train_iters = int(np.ceil(len(train_dataset.data) / train_dataset.step))
    num_val_iters = int(np.ceil(len(val_dataset.data) / val_dataset.step))

    model = MomentaDiffusionModel(config).to(config.DEVICE)
    optimizer = torch.optim.Adam(model.parameters(), lr=config.LR)
    diffusion = Diffusion(config)

    checkpoint_dir = f"{save_folder}/checkpoints"
    os.makedirs(checkpoint_dir, exist_ok=True)

    losses = []
    bdt_output_train = []
    bdt_output_val = []

    for epoch in range(config.EPOCHS):
        print(f"\n===== Starting epoch {epoch + 1}/{config.EPOCHS} =====")
        model.train()
        running_loss = 0.0

        pbar = tqdm(train_loader, desc=f"Epoch {epoch + 1}", total=num_train_iters, leave=False)
        for i, x0 in enumerate(pbar):
            x0 = x0.to(config.DEVICE)
            loss = diffusion.training_step(model, x0)

            if torch.isnan(loss):
                print(f"NaN detected in loss at batch {i}; skipping batch")
                continue

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            optimizer.step()

            running_loss += loss.item()
            pbar.set_postfix({"avg_loss": f"{running_loss / (i + 1):.4f}"})

        epoch_loss = running_loss / (i + 1)
        losses.append(epoch_loss)
        print(f"Epoch {epoch + 1} finished, average training loss: {epoch_loss:.4f}")

        if (epoch + 1) % config.CHECKPOINT_INTERVAL == 0:
            print(f"Running validation on train_loader for epoch {epoch + 1}...")
            bdt_auc_train = run_validation(
                model,
                config,
                train_loader,
                epoch + 1,
                max_batches=num_train_iters,
                save_folder=save_folder,
                plot_distributions_flag=False,
            )

            print(f"Running validation for epoch {epoch + 1}...")
            bdt_auc = run_validation(
                model,
                config,
                val_loader,
                epoch + 1,
                max_batches=num_val_iters,
                save_folder=save_folder,
                plot_distributions_flag=True,
            )

            checkpoint_path = os.path.join(checkpoint_dir, f"bnct_model_epoch{epoch + 1}.pt")
            torch.save(model.state_dict(), checkpoint_path)
            print(f"Checkpoint saved: {checkpoint_path}")

            bdt_output_train.append(bdt_auc_train)
            bdt_output_val.append(bdt_auc)
            plot_loss(losses, epoch + 1, save_folder=save_folder_loss)
            plot_BDT_AUC_vs_epoch(
                bdt_output_train,
                bdt_output_val,
                epoch + 1,
                config.CHECKPOINT_INTERVAL,
                save_folder=save_folder_bdt_auc,
            )

    if bdt_output_train and bdt_output_val:
        plot_BDT_AUC_vs_epoch(
            bdt_output_train,
            bdt_output_val,
            config.EPOCHS,
            config.CHECKPOINT_INTERVAL,
            save_folder=save_folder_bdt_auc,
        )

    final_model_path = f"{save_folder}/bnct_model_final.pt"
    torch.save(model.state_dict(), final_model_path)
    print(f"Training complete. Final model saved: {final_model_path}")

    elapsed = time.time() - time0
    print(f"Executed in {int(elapsed // 60)}m {elapsed % 60:.2f}s")


if __name__ == "__main__":
    main()
