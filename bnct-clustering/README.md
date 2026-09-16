# BNCT DNA clustering

This stage converts `bnct-dna-simulation` ROOT output into strand-break counts.
It uses the AlphaGlue continuous-DNA C++ clustering algorithm, which is
effectively the same break-counting algorithm as TAT after removing unused
members. The direct and indirect hit selection follows the TAT workflow:
sum direct energy per sugar and upstream voxel event, then sample a strand break;
sample OH/deoxyribose reactions for indirect breaks. A sugar is counted only
once per source and upstream voxel event.

Each phase-space row entering the DNA simulation is one **DNA replay event**.
The output has one row per `(upstream_seedID, upstream_eventID,
upstream_voxelID)`. Breaks from all DNA replay events in that voxel are
clustered together. Contributing DNA event, particle, primary, track, and
parent IDs are retained as semicolon-separated lists. `ix`, `iy`, `iz` are
decoded from the upstream voxel ID. The per-voxel CSV includes zero-damage
voxel events; `*_by_z.csv` sums counts by Z layer and includes both voxel
event and replay counts. These are per-voxel cluster counts, not a
whole-nucleus damage classification.
`*_cdsb` counts DSB clusters with more than two strand breaks. Hybrid and
mixed source labels come from the AlphaGlue C++ cluster source classification.

## Build and run

The `clustering` Conda environment has Python, pybind11, uproot and scipy:

```bash
cd /Users/yw18581/work/AI-BNCT
/opt/anaconda3/bin/conda run -n clustering cmake -S bnct-clustering -B bnct-clustering/build \
  -DPython3_EXECUTABLE=/opt/anaconda3/envs/clustering/bin/python
cmake --build bnct-clustering/build -j8
/opt/anaconda3/bin/conda run -n clustering python bnct-clustering/run.py \
  bnct-dna-simulation/build/alpha_dna.root --output bnct-clustering/alpha_damage.csv
/opt/anaconda3/bin/conda run -n clustering python bnct-clustering/test_run.py
```

The bundled sugar file is the default. It must match the DNA geometry used
by `rbe`. The grid defaults to 4 x 4 x 40 at 300 nm, matching the upstream
simulation; set `--nx`, `--ny`, `--nz` and `--voxel-size-nm` for another grid.
Use `--seed` to repeat the stochastic damage sampling.
Older DNA ROOT files with `part1_*` branches are readable. In those files,
`upstream_primaryID` is inferred from `part1_particleSource`, since the
separate primary ID was not saved in that schema.

The default `--damage-preset alphaglue` uses 5 to 37.5 eV for direct breaks
and 0.405 for OH/deoxyribose indirect breaks. The optional
`--damage-preset molecular-bnct` uses the MolecularBNCT default 17.5 eV
direct threshold and probability 1.0 for OH strand damage. These presets
change break formation only; MolecularBNCT's `DSB+`/`DSB++` classification is
a distinct, per-damage-record scheme and is not claimed by this output.

The C++ files in `include/` and `src/` were copied from
`AlphaGlue_Collab/Clustering` on 2026-09-16. `pyClustering.cc` started from
the same version and adds explicit complex-DSB counts. Geant4 license headers
are preserved where present.
