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

The baseline reports three damage views: direct-distance prediction from the
query alone, an oracle diagnostic using true local transport, and a chained
diagnostic using point-predicted transport. The sampler uses the direct-distance
model; conditional StEG is intended to improve the chain by retaining transport
fluctuations.
