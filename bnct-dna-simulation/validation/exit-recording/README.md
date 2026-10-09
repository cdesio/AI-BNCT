# Exit-recording smoke tests

Branch: `codex/dna-voxel-exits`.
Local runtime: Geant4 11.1.3; four worker threads; seed 1234.
Input: first four rows converted from `data_100/alpha_100_seed1234_ps.root`.
These small tests do not establish reproducibility of a historical full run.

## Results

- Build succeeded with `cmake --build bnct-dna-simulation/build -j4`.
- Chemistry on and off, with `--save-exits`: four exit rows each,
  all `Outcome=1`, at the 300 nm cube boundary within geometry tolerance.
- Chemistry-off recording enabled versus disabled: `EventEdep` and
  `PS_data` matched exactly. `VoxelExit` was absent when recording was disabled.
- Chemistry-on versus chemistry-off: input records matched, but deposited
  energies and exit energies/positions differed. Same seed alone was insufficient
  to recover matched transitions in this controlled pair.
- Separate 1 keV alpha started at the cube centre: `Outcome=2`, zero final
  kinetic energy, stopping position inside the cube.
- A lithium stopping fixture failed at particle creation in the local runtime;
  lithium exit/stopping coverage remains unverified here.
- Historical `data_1k` output reports Geant4 11.3, unlike this local runtime.

## Consequence

Do not join chemistry-off replay transitions to existing damage rows on the
assumption that seeds guarantee identical physical tracks. Test recovery using
the original runtime/configuration before accepting such a join. For new matched
training data, run chemistry enabled with `--save-exits`, capturing exits and
damage together.

The ROOT files here are diagnostic smoke-test outputs, not training data.
Use `compare_replay_physics.py` to repeat the physical-output comparison.
