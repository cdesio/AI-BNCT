# BNCT Response AI

Direct-distance surrogate models for the response of one 300 nm DNA voxel to
one BNCT secondary particle. The scientific query is:

```text
particle + initial energy + direction + distance
    -> local entry state + voxel transport
    -> local DNA damage
```

The model does not generate a full cell or predict unhit voxels. One row is one
primary entering one voxel. Distance is calculated continuously along the
primary's initial direction rather than being represented by a voxel index.
Local entry position and scattered entry direction are transport outputs, not
query inputs. The simple damage baseline marginalizes over those microscopic
variables and uses energy-loss, path-length, step, and LET summaries. Conditional
StEG can later preserve and generate the complete correlated local state.

## Phase-space voxel response

The current primary task replaces Geant4-DNA replay while retaining upstream
Geant4 transport. One phase-space row supplies the local particle state at a
voxel entry:

```text
particle + entry energy + local entry point + entry direction
    -> DNA energy deposition + voxel damage
```

Build the joined phase-space response data and train both approaches:

```bash
PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.data \
  --input-dir ../data_100 --output-dir data/processed

PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.ps_baseline \
  --data-dir data/processed --output-dir runs/ps_baseline

PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.steg \
  --data-dir data/processed --output-dir runs/ps_steg \
  --stage ps_response --epochs 100
```

Apply the baseline directly to a new upstream phase-space ROOT file:

```bash
PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.ps_predict \
  ../data_100/alpha_100_seed1234_ps.root \
  --model-dir runs/ps_baseline --output runs/ps_baseline/alpha_predictions.csv
```

Generate one or more stochastic StEG responses per phase-space row:

```bash
PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.ps_generate \
  ../data_100/alpha_100_seed1234_ps.root \
  --model-dir runs/ps_steg --samples-per-row 10 \
  --output runs/ps_steg/alpha_generated.csv
```

The 100-primary pilot has 10,686 response rows but only 400 independent
primaries. Splits are therefore grouped by upstream primary. The present DNA
ROOT output stores energy deposition and damage hits, but not a complete set of
local DNA-track path summaries; those can be added as targets when future
simulation output records them.

Two conditioning variants are supported:

```text
local:   particle + local entry energy, position and direction
history: local inputs + world entry position, projected distance,
         initial primary energy and direction
```

Train the history-aware variants with:

```bash
PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.ps_baseline \
  --data-dir data/processed --output-dir runs/ps_baseline_history \
  --input-mode history

PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.steg \
  --data-dir data/processed --output-dir runs/ps_steg_history \
  --stage ps_response_history --epochs 100
```

## Data contract

`response.csv.gz` joins the upstream `primary`, `phase_space`, and `steps`
trees to the clustered damage CSV using:

```text
case_id + SeedID + EventID + TrackID + VoxelID
```

The damage CSV is unique at `(SeedID, EventID, VoxelID)` in the current
single-track datasets. Splits are assigned by `(case_id, SeedID, EventID)`, so
voxels from one primary never occur in multiple splits.

`range.csv.gz` contains one row per primary and its maximum observed
longitudinal distance. This is a support model, not a claim of exact stopping
range when a future particle exits the simulated grid.

## Pilot workflow

Run from this directory with the existing `ml` environment:

```bash
PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.data \
  --input-dir ../data_100 --output-dir data/processed

PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.baseline \
  --data-dir data/processed --output-dir runs/baseline

PYTHONPATH=src MPLCONFIGDIR=.matplotlib \
  /opt/anaconda3/envs/ml/bin/python -m bnct_response.plot_results \
  --data-dir data/processed --model-dir runs/baseline \
  --output-dir runs/baseline/figures

PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.steg \
  --data-dir data/processed --output-dir runs/steg_transport --stage transport

PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.steg \
  --data-dir data/processed --output-dir runs/steg_damage --stage damage
```

Sample the chained baseline at a requested condition:

```bash
PYTHONPATH=src /opt/anaconda3/envs/ml/bin/python -m bnct_response.sample \
  --model-dir runs/baseline --particle alpha \
  --initial-energy-mev 1.47 --distance-um 4.0 --samples 20
```

The present pilot has only two initial energies per species and one direction.
Interpolation across a broad energy or angular domain requires a deliberately
sampled simulation campaign.

## Two-step damage model

The two-step Extra Trees model separates sparse damage occurrence from the
positive damage response:

```text
phase-space entry -> damage gate
                  -> no damage: zero counts
                  -> damage: conditional positive-count regressor
```

Build the larger dataset and train local and history-aware variants with:

```bash
PYTHONPATH=src python -m bnct_response.data \
  --input-dir ../data_1k --manifest configs/data_1k.json \
  --output-dir data/processed_1k

PYTHONPATH=src python -m bnct_response.ps_two_step \
  --data-dir data/processed_1k --output-dir runs/ps_two_step_1k \
  --input-mode local

PYTHONPATH=src python -m bnct_response.ps_two_step \
  --data-dir data/processed_1k --output-dir runs/ps_two_step_history_1k \
  --input-mode history

PYTHONPATH=src python -m bnct_response.plot_two_step
```

The gate threshold is selected on validation primaries by maximum F1 and then
held fixed for the test set. Because the Extra Trees classifier uses balanced
class weights, its raw probabilities are calibrated with a logistic model fit
only on validation primaries. `P_AnyDamage` is the calibrated probability and
`RawP_AnyDamage` retains the original classifier score.

### Two-step StEG

Train the local StEG second step only on damaged voxels, preferably on a GPU:

```bash
python -m bnct_response.steg \
  --data-dir data/processed_1k \
  --output-dir runs/ps_steg_positive_damage_1k \
  --stage ps_positive_damage \
  --epochs 200 --diffusion-steps 100 --batch-size 1024 --device cuda
```

Evaluate one stochastic response per held-out test voxel:

```bash
python -m bnct_response.evaluate_two_step_steg \
  --data-dir data/processed_1k \
  --gate-model-dir runs/ps_two_step_1k \
  --steg-model-dir runs/ps_steg_positive_damage_1k \
  --output-dir runs/ps_two_step_steg_eval_1k \
  --device cuda
```

For a new phase-space ROOT file, generate repeated stochastic responses with:

```bash
python -m bnct_response.ps_two_step_generate INPUT_PS.root \
  --gate-model-dir runs/ps_two_step_1k \
  --steg-model-dir runs/ps_steg_positive_damage_1k \
  --output generated.csv --samples-per-row 10 --device cuda
```

Step 1 samples damage occurrence from the calibrated probability. Step 2 is
called only for positive samples and jointly generates the five damage counts.
Energy deposition remains the Extra Trees point prediction in this version.

### Positive-damage StEG parameter sweep

The StEG CLI exposes `--width`, `--layers`, `--learning-rate`,
`--weight-decay`, `--output-transform`, `--min-epochs`, and
`--early-stopping-patience`. Checkpoints store these settings and remain
compatible with older checkpoints that used the default 256-wide, five-layer
network.

Run the four sequential GPU variants with:

```bash
bash scripts/run_positive_steg_sweep.sh
```

The sweep compares 50 versus 100 diffusion steps, two network widths, and
quantile-normal versus `log1p + standardise` outputs. Every model is evaluated
automatically. Evaluation reports both raw generated counts and counts after
the physical consistency constraints, making constraint-induced TotalSB
inflation visible.

The baseline reports three damage views: direct-distance prediction from the
query alone, an oracle diagnostic using true local transport, and a chained
diagnostic using point-predicted transport. The sampler uses the direct-distance
model; conditional StEG is intended to improve the chain by retaining transport
fluctuations.
