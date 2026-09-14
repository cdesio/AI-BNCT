# Theranostics StEG Smoke Test

This is a quick StEG adaptation for theranostics photon detector CSV outputs.
It reads the existing CSV files by absolute path and does not copy or move the
source data.

## Features

Rows are filtered to detector-entering gamma photons:

```text
particleName == "gamma"
boundary == "enter"
volumeName starts with "physPhotonDetector"
```

The model uses only continuous columns:

```python
[
    "x_mm",
    "y_mm",
    "z_mm",
    "px",
    "py",
    "pz",
    "kineticEnergy_MeV",
    "time_s",
]
```

## Run

From this folder:

```bash
/opt/anaconda3/bin/conda run -n ml python scripts/inspect_filtered_data.py
/opt/anaconda3/bin/conda run -n ml python transformers/fit_quantile_transformer.py
/opt/anaconda3/bin/conda run -n ml python train.py
```

Defaults are intentionally small: first 5 CSV files, 50k filtered rows, CPU,
10 epochs, and 100 diffusion steps. Change these in `config.py` after the smoke
test passes.

Each `train.py` call creates a fresh timestamped folder under `runs/`, for
example:

```text
runs/20260902_143012_theranostics_photon_smoke/
```

`evaluate.py` evaluates the latest training run by default and writes a fresh
subfolder under that run:

```text
runs/<train-run>/evaluations/<timestamp>/
```

To evaluate a specific run:

```bash
/opt/anaconda3/bin/conda run -n ml python evaluate.py --run-dir runs/<train-run>
```
