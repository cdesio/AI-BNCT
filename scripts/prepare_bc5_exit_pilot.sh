#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
mode="${1:-prepare}"
case "$mode" in
  prepare|submit|prepare-test|submit-test) ;;
  *) echo "Usage: bash scripts/prepare_bc5_exit_pilot.sh [prepare|submit|prepare-test|submit-test]" >&2; exit 2 ;;
esac
cases=(alpha_1p47 alpha_1p78 lithium_0p84 lithium_1p01)
particles=(alpha alpha lithium lithium)
energies=(1.47 1.78 0.84 1.01)
seed=1234
events=100
if [[ "$mode" == *-test ]]; then
  cases=(alpha_1p47)
  particles=(alpha)
  energies=(1.47)
  events=10
fi
if [[ "$mode" == prepare* ]]; then
  for label in "${cases[@]}"; do
    test ! -e "$ROOT/jobs/exits${events}_${label}_seed${seed}"
  done
  for i in "${!cases[@]}"; do
    python scripts/make_bc5_jobs.py \
      --particle "${particles[$i]}" --energy-mev "${energies[$i]}" \
      --name "exits${events}_${cases[$i]}" --events "$events" --seed "$seed" \
      --save-exits --checkpoint-events 50 --dna-cpus 4 --dna-mem 32GB \
      --upstream-exe "$ROOT/bnct-voxel-ps/build-g4-lithium/bnctVoxelPS" \
      --dna-exe "$ROOT/bnct-dna-simulation/build-g4-lithium/rbe"
  done
else
  for label in "${cases[@]}"; do
    test -f "$ROOT/jobs/exits${events}_${label}_seed${seed}/submit_all.sh"
  done
  for label in "${cases[@]}"; do
    bash "$ROOT/jobs/exits${events}_${label}_seed${seed}/submit_all.sh"
  done
fi
