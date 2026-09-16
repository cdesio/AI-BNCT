#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "RunConfig.hh"
#include "TrackingOptions.hh"

#include "G4RunManagerFactory.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "Randomize.hh"

#ifdef BNCT_WITH_GEANT4_UIVIS
#include "G4VisExecutive.hh"
#endif

#include <cstdlib>
#include <iostream>

namespace {

void PrintUsage(const char* program)
{
  std::cout << "Usage: " << program
            << " [-mac macro.mac] [-out output_stem] [-seed integer] [-gui]"
            << " [-physics ion|transport] [-bin phase_space.bin]\n";
}

}  // namespace

int main(int argc, char** argv)
{
  G4String macro;
  G4bool gui = false;
  G4String physicsMode = "ion";
  RunConfig config;

  for (G4int i = 1; i < argc; ++i) {
    const G4String arg = argv[i];
    if ((arg == "-h") || (arg == "--help")) {
      PrintUsage(argv[0]);
      return 0;
    } else if (arg == "-mac" && i + 1 < argc) {
      macro = argv[++i];
    } else if (arg == "-out" && i + 1 < argc) {
      config.outputStem = argv[++i];
    } else if (arg == "-bin" && i + 1 < argc) {
      config.binaryOutput = argv[++i];
    } else if (arg == "-seed" && i + 1 < argc) {
      config.seed = std::atol(argv[++i]);
    } else if (arg == "-gui") {
      gui = true;
    } else if (arg == "-physics" && i + 1 < argc) {
      physicsMode = argv[++i];
    } else {
      std::cerr << "Unknown or incomplete argument: " << arg << "\n";
      PrintUsage(argv[0]);
      return 1;
    }
  }

  if (physicsMode != "ion" && physicsMode != "transport") {
    std::cerr << "Unknown physics mode: " << physicsMode << "\n";
    PrintUsage(argv[0]);
    return 1;
  }

  if (config.binaryOutput.empty()) {
    config.binaryOutput = config.outputStem + ".bin";
  }

  CLHEP::HepRandom::setTheSeed(config.seed);

  auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::Serial);
  auto* detector = new DetectorConstruction();
  auto* trackingOptions = new TrackingOptions();
  runManager->SetUserInitialization(detector);
  runManager->SetUserInitialization(new PhysicsList(physicsMode == "transport"));
  runManager->SetUserInitialization(new ActionInitialization(config, detector, trackingOptions));

#ifdef BNCT_WITH_GEANT4_UIVIS
  G4VisExecutive* visManager = nullptr;
  if (gui) {
    visManager = new G4VisExecutive();
    visManager->Initialize();
  }
#endif

  G4UIExecutive* uiExec = nullptr;
  if (gui) {
    uiExec = new G4UIExecutive(argc, argv);
  }

  auto* ui = G4UImanager::GetUIpointer();
  if (!macro.empty()) {
    ui->ApplyCommand("/control/execute " + macro);
  }

  if (uiExec) {
    uiExec->SessionStart();
    delete uiExec;
  } else if (macro.empty()) {
    PrintUsage(argv[0]);
  }

#ifdef BNCT_WITH_GEANT4_UIVIS
  delete visManager;
#endif
  delete runManager;
  delete trackingOptions;
  return 0;
}
