#!/usr/bin/env bash
set -euo pipefail

DATA_DIR=${DATA_DIR:-data/processed_1k}
RUNS_DIR=${RUNS_DIR:-runs}
GATE_DIR=${GATE_DIR:-runs/ps_two_step_1k}
DEVICE=${DEVICE:-cuda}

run_variant() {
  local name=$1
  local steps=$2
  local width=$3
  local layers=$4
  local transform=$5
  local model_dir="${RUNS_DIR}/${name}"

  python -m bnct_response.steg \
    --data-dir "${DATA_DIR}" \
    --output-dir "${model_dir}" \
    --stage ps_positive_damage \
    --epochs 300 \
    --min-epochs 80 \
    --early-stopping-patience 30 \
    --diffusion-steps "${steps}" \
    --batch-size 512 \
    --width "${width}" \
    --layers "${layers}" \
    --learning-rate 1e-4 \
    --output-transform "${transform}" \
    --device "${DEVICE}"

  python -m bnct_response.evaluate_two_step_steg \
    --data-dir "${DATA_DIR}" \
    --gate-model-dir "${GATE_DIR}" \
    --steg-model-dir "${model_dir}" \
    --output-dir "${model_dir}_eval" \
    --device "${DEVICE}"
}

run_variant ps_positive_sweep_A_q50_small 50 128 3 quantile
run_variant ps_positive_sweep_B_q100_small 100 128 3 quantile
run_variant ps_positive_sweep_C_q100_medium 100 256 3 quantile
run_variant ps_positive_sweep_D_log100_small 100 128 3 log_standard
