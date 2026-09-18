# BC5 pipeline jobs

Run the generator from the AI-BNCT checkout on BC5. It creates one run
directory containing `upstream.mac`, `dna.mac`, four independent SLURM jobs,
and `submit_all.sh`. No account or partition is specified.

## One-time builds on BC5

Use the patched lithium Geant4 module for both C++ simulations:

```bash
module use /projects/b56v/software/modulefiles
module load geant4/11.3.0-lithium
source /projects/b56v/software/geant4-v11.3.0-lithium-install/bin/geant4.sh
G4_CMAKE_FILE=$(find /projects/b56v/software/geant4-v11.3.0-lithium-install \
  -name Geant4Config.cmake -print -quit)
test -n "$G4_CMAKE_FILE"
G4_CMAKE_DIR=${G4_CMAKE_FILE%/*}
cmake -S bnct-voxel-ps -B bnct-voxel-ps/build-g4-lithium \
  -DGeant4_DIR="$G4_CMAKE_DIR" -DWITH_GEANT4_UIVIS=OFF
cmake --build bnct-voxel-ps/build-g4-lithium -j4
cmake -S bnct-dna-simulation -B bnct-dna-simulation/build-g4-lithium \
  -DGeant4_DIR="$G4_CMAKE_DIR" -DWITH_GEANT4_VIS=OFF
cmake --build bnct-dna-simulation/build-g4-lithium -j4
grep '^Geant4_DIR:' bnct-voxel-ps/build-g4-lithium/CMakeCache.txt \
  bnct-dna-simulation/build-g4-lithium/CMakeCache.txt
```

Use fresh build directories after changing Geant4 versions: CMake caches the
previous `Geant4_DIR`. Check both `CMakeCache.txt` files and confirm they point
to the patched installation, not `/projects/b56v/software/conda-envs/geant4`.
The latter requires a newer C++ runtime and can fail to link with
`GLIBCXX_3.4.30` or `GLIBCXX_3.4.32` errors. The batch jobs do not need
visualization; `WITH_GEANT4_VIS=OFF` avoids OpenGL/Qt dependencies.

The generator defaults to executables in `build/`. When using the fresh build
directories above, pass `--upstream-exe` and `--dna-exe` with their absolute
paths to `make_bc5_jobs.py`.

Build the Python clustering extension in the same Python environment used by
the clustering job. That environment needs `numpy`, `scipy`, `uproot`,
`pybind11`, and `cmake`:

```bash
module use /projects/b56v/software/modulefiles
module load clustering/conda
cmake -S bnct-clustering -B bnct-clustering/build -DPython3_EXECUTABLE="$(which python)"
cmake --build bnct-clustering/build -j4
```

## Generate and submit

```bash
python scripts/make_bc5_jobs.py --particle alpha --name alpha_1p47 \
  --events 100000 --seed 6069075 \
  --upstream-exe "$PWD/bnct-voxel-ps/build-g4-lithium/bnctVoxelPS" \
  --dna-exe "$PWD/bnct-dna-simulation/build-g4-lithium/rbe"
./jobs/alpha_1p47_seed6069075/submit_all.sh
```

The alpha template uses 1.47 MeV with a 4 x 4 x 40 grid of 300 nm voxels.
Use `--particle lithium` for the 0.84 MeV Li-7 template. `--events` changes
the upstream `/run/beamOn` count. `--seed` is passed to both simulations and
the stochastic clustering step. The generator's `--project-root` defaults to
the checkout containing the generator; use it if generating scripts for a
different checkout path. `--output-dir` can put the jobs and results elsewhere.
The jobs default to the patched `geant4/11.3.0-lithium` module. If you build
the alpha-only pipeline with Geant4 11.1.3, generate its jobs with matching
`--geant4-module` and `--geant4-sh` values. Pass `--geant4-sh ''` if that module
already sets up the Geant4 runtime. Build and run each executable against the
same Geant4 installation. Lithium runs should use a build with the patched
lithium definitions and models.

The launcher submits upstream, then DNA with `afterok`. DNA replays the binary
phase space in batches of 10,000 records by default. Each finished batch has a
ROOT file and a `.done` marker. When the job approaches its time limit, it
submits another DNA job, which skips completed batches. After all batches
finish, a merge job runs `hadd` and submits clustering. Set batch size with
`--checkpoint-events`; choose a size that normally finishes well within one
DNA job. The merge job loads `root/conda` and requires `hadd`.
Each replay event gets a reproducible seed derived from `--seed` and its
original phase-space record number. With the same executable, input, and
settings, changing the checkpoint size does not change an event's seed.
This seeding scheme changes individual outcomes compared with older DNA runs
that seeded only once per process.

To resume after an upstream result already exists:

```bash
./jobs/alpha_1p47_seed6069075/submit_all.sh dna
```

To rerun only the merge or clustering:

```bash
./jobs/alpha_1p47_seed6069075/submit_all.sh merge
```

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
Jobs request email on `FAIL`, `END`, and `TIME_LIMIT`; delivery depends on
BC5's Slurm mail configuration. Use `--mail-user` to set the recipient.
All requested walltimes are checked against BC5's 24-hour limit. The DNA job
defaults to 24 hours, with 23 hours available for DNA batches. A batch that
cannot finish within one job fails with a request to reduce the checkpoint
size. Completed batches are retained. Use
`--dna-budget-seconds` if the DNA walltime is shorter than 24 hours.
