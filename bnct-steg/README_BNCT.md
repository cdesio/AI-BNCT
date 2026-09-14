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

By default, files are skipped if the number of reconstructed consecutive damage
chunks does not match the number of classification rows. To keep those files
with an explicit fallback assignment status:

```bash
python scripts/build_primary_summary.py --allow-classification-mismatch
```

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
