import os
os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.getcwd(), ".matplotlib"))

import argparse
import time
from dataclasses import asdict

import numpy as np
import torch
from torch.utils.data import DataLoader

from config import Config
from dataset import BNCTCsvDamageDataset, resolve_event_splits
from model import MomentaDiffusionModel
from run_paths import create_eval_run_dir, resolve_train_run_dir, write_run_metadata
from validate import run_validation


def parse_args():
    parser = argparse.ArgumentParser(description="Evaluate a trained BNCT StEG run.")
    parser.add_argument(
        "--run-dir",
        help="Training run directory. Defaults to runs/latest_train_run.txt.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    config = Config()
    train_run_dir = resolve_train_run_dir(config, args.run_dir)
    save_folder = create_eval_run_dir(train_run_dir)
    total_events, train_split, val_split = resolve_event_splits(config)

    test_dataset = BNCTCsvDamageDataset(
        config,
        start=val_split,
        stop=total_events,
        input_scaler_flag=config.INPUT_SCALER_FLAG,
    )
    num_test_iters = int(np.ceil(len(test_dataset.data) / test_dataset.step))
    test_loader = DataLoader(test_dataset, batch_size=config.BATCH_SIZE)

    model = MomentaDiffusionModel(config).to(config.DEVICE)
    model_path = f"{train_run_dir}/bnct_model_final.pt"
    model.load_state_dict(torch.load(model_path, map_location=config.DEVICE))

    write_run_metadata(
        f"{save_folder}/run_metadata.json",
        {
            "run_dir": save_folder,
            "kind": "evaluate",
            "train_run_dir": train_run_dir,
            "model_path": model_path,
            "created_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "total_events": total_events,
            "train_split": train_split,
            "val_split": val_split,
            "config": asdict(config),
        },
    )

    run_validation(
        model,
        config,
        test_loader,
        epoch="final_test",
        max_batches=num_test_iters,
        plot_distributions_flag=True,
        save_folder=save_folder,
    )


if __name__ == "__main__":
    main()
