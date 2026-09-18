# BNCT DNA Voxel Replay

This app is the downstream Geant4-DNA stage for `bnct-voxel-ps`.

It is based on the newer AlphaGlue `DNA-simulation` app, with MolecularBNCT's
patched lithium DNA physics added for BNCT Li-7 replay.  The geometry is the
AlphaGlue continuous 300 nm DNA voxel geometry: a water `chromatinSegment`
filled from sugar and histone binary files.

## What Is Replayed

Every row in `Ntuples/phase_space` from `bnctVoxelPS` should become one primary
event in this DNA simulation.  The particle is launched inside the 300 nm DNA
voxel at the local voxel-entry position and with the recorded direction and
kinetic energy.

Particle IDs in the intermediate binary PS file:

- `3`: alpha
- `53`: BNCT Li-7, replayed as Geant4-DNA `lithium+++`
- `1`: electron
- `2`: gamma
- `11`: positron

## Build

```bash
cd /Users/yw18581/work/AI-BNCT
cmake -S bnct-dna-simulation -B bnct-dna-simulation/build \
  -DCMAKE_PREFIX_PATH=/opt/anaconda3/envs/geant4.1 \
  -DGeant4_DIR=/opt/anaconda3/envs/geant4.1/lib/cmake/Geant4 \
  -DCLHEP_DIR=/opt/anaconda3/envs/geant4.1/lib/CLHEP-2.4.6.2
cmake --build bnct-dna-simulation/build -j8
```

## Use Upstream Binary Phase Space

`bnctVoxelPS` writes the DNA replay binary directly by default:

```bash
cd /Users/yw18581/work/AI-BNCT/bnct-voxel-ps/build
./bnctVoxelPS -mac alpha_300nm_4x4x40.mac -out alpha_ps -seed 1234
```

This produces `alpha_ps.bin`, which can be passed directly to `rbe -in`.

## Optional ROOT Conversion

The converter is kept for older ROOT-only outputs and for debugging. It requires
`uproot` in the Python environment used for conversion.

```bash
python bnct-dna-simulation/convert_bnct_root_to_decay_ps.py \
  bnct-voxel-ps/build/alpha_1000_secondaries.root \
  bnct-dna-simulation/build/alpha_1000_secondaries_ps.bin
```

For a quick smoke test, convert only a few rows:

```bash
python bnct-dna-simulation/convert_bnct_root_to_decay_ps.py \
  bnct-voxel-ps/build/alpha_1000_secondaries.root \
  bnct-dna-simulation/build/smoke_ps.bin \
  --max-rows 10
```

If boundary starts prove fragile, add a small inward displacement along the
recorded direction:

```bash
--nudge-nm 0.01
```

## Run DNA Replay

```bash
cd /Users/yw18581/work/AI-BNCT/bnct-dna-simulation/build
./rbe \
  -mac rbe_PS.in \
  -in /Users/yw18581/work/AI-BNCT/bnct-voxel-ps/build/alpha_ps.bin \
  -out smoke_dna.root \
  -seed 1234
```

By default, `rbe` uses the bundled 300 nm geometry files in `geometryFiles/`.
Pass `-sugar` or `-histone` only when replaying with a different DNA geometry.

For a record range, pass `-start-record N -record-count M`. The replay event
numbers in all output ntuples remain the original zero-based record numbers,
so independently completed ROOT files can be merged with `hadd`. For BNCT
phase-space files, pass `-record-doubles 18` to specify the record width
explicitly. The BC5 job generator runs these ranges in resumable batches and
merges them before clustering.
For phase-space replay, each event's random seed is derived from `-seed` and
its original record number. The seed therefore stays the same across DNA
checkpoints and in an uninterrupted run using this version of the executable.
Earlier results used one seed per process and will not match event by event.

The output ROOT file follows the AlphaGlue layout: `EventEdep`, `Direct`,
`Indirect`, and `Info` ntuples under the `ntuple` directory.

For strand-break clustering and per-Z summaries, see
[`bnct-clustering`](../bnct-clustering/README.md).
