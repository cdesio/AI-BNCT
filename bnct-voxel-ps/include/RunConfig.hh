#ifndef BNCT_RUN_CONFIG_HH
#define BNCT_RUN_CONFIG_HH

#include "globals.hh"

struct RunConfig {
  G4String outputStem = "bnct_voxel_ps";
  G4String binaryOutput;
  G4long seed = 1234;
};

#endif
