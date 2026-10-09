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

The output ROOT file follows the AlphaGlue layout: `EventEdep`, `Direct`,
`Indirect`, and `Info` ntuples under the `ntuple` directory.

## Optional voxel transitions

Add `--save-exits` to a phase-space replay to write `ntuple/VoxelExit`.
One row is written per event, joined by `EventNum` and the upstream identifiers.
`Outcome=1` means first outward cube crossing; `Outcome=2` means zero kinetic
energy inside the cube; `Outcome=0` means no recognised outcome and must be
excluded from transition training. Missing outcomes have placeholder zero fields.
Positions are local nm, energy is MeV, and directions are dimensionless.
`Particle` records the exiting/stopping charge state. Alpha and lithium charge
exchange replacement tracks are followed as the same physical ion. This assumes
the current single-ion DNA physics; it is not a general nuclear branching tracker.
Recording does not kill tracks or consume random numbers.

Use `-chemOFF --save-exits` for physical-only replay. To recover transitions
matched to old damage, use the original Geant4 version, physics, phase-space
binary, geometry, seed and threading/run configuration. A short truncated replay
can change event seed scheduling relative to a long run; a same-seed smoke test
alone does not establish historical recovery. Compare physical outputs with:

```bash
python compare_replay_physics.py original_dna.root replay_dna.root
```

The comparison reports exact event-deposition and input-record differences.
Chemistry-on/off agreement on a new controlled pair is separate from agreement
with a historical run. Record transitions with chemistry enabled for new joint
exit-and-damage training datasets.

For strand-break clustering and per-Z summaries, see
[`bnct-clustering`](../bnct-clustering/README.md).
