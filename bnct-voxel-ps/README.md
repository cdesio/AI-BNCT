# BNCT Voxel Phase-Space Transport

Minimal Geant4-DNA transport app for BNCT-style alpha and lithium secondaries in
a parameterised grid of water voxels.

The first target is phase-space/tracking generation only. DNA damage simulation
is intentionally left to a downstream app such as MolecularBNCT, AlphaGlue, or
the TAT phase-space workflow.

## Defaults

- Voxel size: 300 nm cube
- Grid: 4 x 4 x 40
- Physical extent: 1.2 um x 1.2 um x 12 um
- Source position: (0, 0, 0)
- Direction: +Z
- Material: water
- Physics: ion-focused Geant4 transport by default; transport-only mode is available for fast geometry visualization

## Build

```bash
cd /Users/yw18581/work/AI-BNCT/bnct-voxel-ps
mkdir -p build
cd build
cmake ..
cmake --build . -j8
```

## Run

```bash
./bnctVoxelPS -mac alpha_300nm_4x4x40.mac -out alpha_300nm_4x4x40 -seed 1234
```

For visualization, run with `-gui` and execute `vis.mac` after initialization,
or include `/control/execute vis.mac` after `/run/initialize` in a macro.
Use `-physics transport` when you only want to inspect geometry/trajectories
without building EM energy-loss tables.

Quick geometry visualization macros are provided:

```bash
./bnctVoxelPS -physics transport -mac alpha_vis.mac -gui -out alpha_vis -seed 1234
./bnctVoxelPS -physics transport -mac lithium_vis.mac -gui -out lithium_vis -seed 1234
```

The ROOT file contains:

- `primary`: one row per event
- `steps`: one row per step inside a voxel
- `phase_space`: one row when a track enters a voxel

The app also writes an AlphaGlue/TAT-compatible binary phase-space file next to
the ROOT file by default:

```bash
./bnctVoxelPS -mac alpha_300nm_4x4x40.mac -out alpha_ps -seed 1234
```

This writes:

- `alpha_ps.root`
- `alpha_ps.bin`

Use `-bin custom_name.bin` to choose a different binary path.  The binary file
contains one 18-double record per replayable phase-space row and can be passed
directly to `bnct-dna-simulation/build/rbe -in`.

## Main Macro Controls

```text
/det/setVoxelSize 300 nm
/det/setGrid 4 4 40
/det/setMaxStep 30 nm
/det/checkOverlaps false
/tracking/killSecondaries false
/tracking/recordSecondaries true
/primary/particle alpha
/primary/energy 1.47 MeV
/primary/position 0 0 0 nm
/primary/direction 0 0 1
```

Secondary tracks are transported by default. A phase-space row is written once
per `(eventID, trackID, voxelID)` when the track enters a voxel, including
crossings between adjacent parameterised voxels.
