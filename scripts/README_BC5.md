# BC5 pipeline jobs

Run the generator from the AI-BNCT checkout on BC5. It creates one run
directory containing `upstream.mac`, `dna.mac`, three independent SLURM jobs,
and `submit_all.sh`. No account or partition is specified.

## One-time builds on BC5

Use the patched lithium Geant4 module for both C++ simulations:

```bash
module use /projects/b56v/software/modulefiles
module load geant4/11.3.0-lithium
source /projects/b56v/software/geant4-v11.3.0-lithium-install/bin/geant4.sh
cmake -S bnct-voxel-ps -B bnct-voxel-ps/build
cmake --build bnct-voxel-ps/build -j4
cmake -S bnct-dna-simulation -B bnct-dna-simulation/build
cmake --build bnct-dna-simulation/build -j4
```

Build the Python clustering extension in the same Python environment used by
the clustering job. That environment needs `numpy`, `scipy`, `uproot`,
`pybind11`, and `cmake`:

```bash
conda activate clustering
cmake -S bnct-clustering -B bnct-clustering/build -DPython3_EXECUTABLE="$(which python)"
cmake --build bnct-clustering/build -j4
```

## Generate and submit

```bash
python scripts/make_bc5_jobs.py --particle alpha --name alpha_1p47 \
  --events 100000 --seed 6069075
./jobs/alpha_1p47_seed6069075/submit_all.sh
```

The alpha template uses 1.47 MeV with a 4 x 4 x 40 grid of 300 nm voxels.
Use `--particle lithium` for the 0.84 MeV Li-7 template. `--events` changes
the upstream `/run/beamOn` count. `--seed` is passed to both simulations and
the stochastic clustering step. The generator's `--project-root` defaults to
the checkout containing the generator; use it if generating scripts for a
different checkout path. `--output-dir` can put the jobs and results elsewhere.

The launcher submits upstream, then DNA with `afterok`, then clustering with
`afterok`. To resume after an upstream result already exists:

```bash
./jobs/alpha_1p47_seed6069075/submit_all.sh dna
```

To rerun only clustering:

```bash
./jobs/alpha_1p47_seed6069075/submit_all.sh clustering
```

Each stage is also independently resubmittable with `sbatch` on its `.sbatch`
file. Successful upstream output is a ROOT file and a binary phase-space file;
the DNA stage reads the binary directly. Clustering writes a per-voxel-event
CSV and a `*_by_z.csv` summary. The job directory contains all results and
`logs/` holds SLURM stdout and stderr. Existing generated directories are not
overwritten by the generator.

Resource defaults are one CPU and 16 GB for upstream, four CPUs and 100 GB
for DNA, and one CPU and 32 GB for clustering. Change these with the
generator's `--*-time`, `--*-mem`, and `--dna-cpus` options before submission.
All requested walltimes are checked against BC5's 24-hour limit. The DNA job
now defaults to 24 hours; choose `--events` from your timing tests so the DNA
stage finishes within that limit with headroom for run-to-run variation.
