import json
import os
from datetime import datetime


def create_train_run_dir(config):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = os.path.join(config.RUNS_ROOT, f"{timestamp}_{config.RUN_TAG}")
    os.makedirs(run_dir, exist_ok=False)

    latest_path = os.path.join(config.RUNS_ROOT, "latest_train_run.txt")
    with open(latest_path, "w") as f:
        f.write(run_dir)

    return run_dir


def resolve_train_run_dir(config, requested_run=None):
    if requested_run:
        return requested_run

    latest_path = os.path.join(config.RUNS_ROOT, "latest_train_run.txt")
    if not os.path.exists(latest_path):
        raise FileNotFoundError(
            f"No latest run marker found at {latest_path}. Run train.py first or pass --run-dir."
        )

    with open(latest_path) as f:
        run_dir = f.read().strip()
    if not run_dir:
        raise RuntimeError(f"Latest run marker at {latest_path} is empty.")
    return run_dir


def create_eval_run_dir(train_run_dir):
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    eval_root = os.path.join(train_run_dir, "evaluations")
    eval_dir = os.path.join(eval_root, timestamp)
    os.makedirs(eval_dir, exist_ok=False)
    return eval_dir


def write_run_metadata(path, payload):
    with open(path, "w") as f:
        json.dump(payload, f, indent=4)
