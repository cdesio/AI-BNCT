#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "${PROJECT_DIR}"
export PYTHONPATH="${PROJECT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"

DATA_DIR=${DATA_DIR:-data/processed_1k}
RUNS_DIR=${RUNS_DIR:-runs}
DEVICE=${DEVICE:-cuda}
EPOCHS=${EPOCHS:-800}
MIN_EPOCHS=${MIN_EPOCHS:-400}
PATIENCE=${PATIENCE:-150}

run_variant() {
  local name=$1
  local transform=$2
  local model_dir="${RUNS_DIR}/${name}"

  python -m bnct_response.steg \
    --data-dir "${DATA_DIR}" \
    --output-dir "${model_dir}" \
    --stage ps_response \
    --epochs "${EPOCHS}" \
    --min-epochs "${MIN_EPOCHS}" \
    --early-stopping-patience "${PATIENCE}" \
    --diffusion-steps 100 \
    --batch-size 512 \
    --width 256 \
    --layers 5 \
    --learning-rate 2e-4 \
    --weight-decay 1e-5 \
    --output-transform "${transform}" \
    --seed 20261008 \
    --device "${DEVICE}"

  python -m bnct_response.evaluate_one_step_steg \
    --data-dir "${DATA_DIR}" \
    --model-dir "${model_dir}" \
    --output-dir "${model_dir}_eval" \
    --device "${DEVICE}" \
    --seed 20261008
}

run_variant ps_one_step_long_quantile quantile
run_variant ps_one_step_long_log log_standard

python -m bnct_response.plot_one_step_transform \
  --runs-dir "${RUNS_DIR}" \
  --output-dir "${RUNS_DIR}/one_step_transform_figures"
