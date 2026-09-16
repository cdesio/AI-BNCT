# BNCT StEG

This folder is a BNCT-focused copy of the LHCb StEG diffusion model. It keeps the
generic tabular diffusion model and reads the old BNCT CSV files in:

```text
/Users/yw18581/work/AI-BNCT/ForChiara-OldData/CSVfiles/*.csv
```

## Environment

Use the existing conda environment:

```bash
conda activate ml
```

or run commands without activating:

```bash
/opt/anaconda3/bin/conda run -n ml python scripts/inspect_csv_data.py
/opt/anaconda3/bin/conda run -n ml python transformers/fit_quantile_transformer.py
/opt/anaconda3/bin/conda run -n ml python train.py
/opt/anaconda3/bin/conda run -n ml python evaluate.py
```

## Data Contract

Each CSV row is treated as one tabular BNCT damage/track record. The initial
smoke-test features are:

```text
PosX_um
PosY_um
PosZ_um
StopPosX_um
StopPosY_um
StopPosZ_um
Dist_Start_to_Origin
Dist_Stop_to_Start
Damage
```

`Particle_ID` is not used because it is an index within each file, not a physical
quantity to generate.

The raw sparse count columns `SSB`, `DSB`, `DSBp`, and `DSBpp` are read from the
CSV files but are not generated directly in this smoke test. The dataset derives
a single binary indicator instead:

```text
Damage = 1 if any raw damage count is nonzero
```

This avoids asking the unconditional diffusion model to represent mostly-zero
integer count distributions as continuous generated values.

During generation, `Damage` is thresholded back to `0` or `1`.

The CSV filenames contain useful condition information, such as particle type,
energy, source distance, and compartment. Those are not included in this first
unconditional smoke test. A later conditional model should handle them explicitly
rather than asking the unconditional generator to produce categorical labels.

## Workflow

Build the primary-level ROOT-derived dataset:

```bash
python scripts/build_primary_summary.py --output data/bnct_primary_summary.csv
python scripts/inspect_primary_summary.py data/bnct_primary_summary.csv
python scripts/validate_primary_summary.py data/bnct_primary_summary.csv
```

By default, this keeps all primary rows available in `primary_source`, including
files where classification rows do not exactly match reconstructed damage chunks.
Those rows are annotated with `ClassificationDamageChunkDelta`,
`PrimaryCountMatchesFilename`, and `DamageLabelQuality` so they can be filtered
later if needed. To reproduce a stricter build:

```bash
python scripts/build_primary_summary.py \
  --strict-classification-match \
  --strict-primary-count
```

For newer MolecularBNCT ROOT files produced with all-primary tracking enabled,
build both a primary-level summary and a per-step table:

```bash
python scripts/build_tracked_primary_dataset.py \
  --input-glob "/Users/yw18581/work/MolecularBNCT-main/bnct_campaign/*.root" \
  --primary-output data/bnct_tracked_primary_summary.csv \
  --step-output data/bnct_primary_track_steps.csv
```

The tracked ROOT format is preferred for future datasets because
`primary_source`, `primary_tracks`, `damage`, and `classification` can be joined
directly with `EventID` and `SeedID`, avoiding the old damage-chunk matching
workaround.

Train primary-level damage prediction baselines:

```bash
python scripts/train_primary_damage_baseline.py \
  --input data/bnct_primary_summary.csv \
  --output-dir runs/primary_damage_baseline_geometry_only
```

By default, the baseline uses only conditions and primary-source geometry. Damage-tree
aggregate features are excluded because they leak damage information for the
current files. To run an intentionally diagnostic leaky model:

```bash
python scripts/train_primary_damage_baseline.py --include-damage-tree-features
```

Run the StEG-like diffusion smoke test on the primary-level geometry features:

```bash
BNCT_CONFIG=primary_geometry python transformers/fit_quantile_transformer.py

BNCT_CONFIG=primary_geometry \
BNCT_DEVICE=cpu \
BNCT_EPOCHS=1 \
BNCT_CHECKPOINT_INTERVAL=1 \
BNCT_TOTAL_EVENTS=5000 \
BNCT_RUN_TAG=bnct_primary_geometry_smoke \
python train.py

BNCT_CONFIG=primary_geometry \
BNCT_TOTAL_EVENTS=5000 \
python evaluate.py --run-dir runs/<train-run>
```

For a longer GPU run, set `BNCT_DEVICE` to the appropriate PyTorch device in
your environment and increase `BNCT_EPOCHS`, `BNCT_TOTAL_EVENTS`, and optionally
`BNCT_NOISE_STEPS`.

Inspect the available rows and feature ranges:

```bash
python scripts/inspect_csv_data.py
```

Fit per-feature quantile transformers:

```bash
python transformers/fit_quantile_transformer.py
```

Run training:

```bash
python train.py
```

Evaluate the latest training run:

```bash
python evaluate.py
```

Evaluate a specific training run:

```bash
python evaluate.py --run-dir runs/<train-run>
```

## Outputs

Each training call creates a fresh folder:

```text
runs/YYYYMMDD_HHMMSS_bnct_csv_damage_smoke/
```

Each evaluation call writes under the selected training run:

```text
runs/<train-run>/evaluations/YYYYMMDD_HHMMSS/
```

The latest training run is tracked in:

```text
runs/latest_train_run.txt
```

The smoke config uses `DEVICE = "cpu"` and `NOISE_STEPS = 100` so it can run
locally. For a serious GPU run, change `DEVICE` and consider increasing
`NOISE_STEPS`, `EPOCHS`, and possibly model width/depth in `config.py`.
