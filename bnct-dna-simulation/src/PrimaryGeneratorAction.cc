//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class

#include "PrimaryGeneratorAction.hh"
#include "Randomize.hh"
#include "G4GeneralParticleSource.hh"
#include "CommandLineParser.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4AnalysisManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4DNAGenericIonsManager.hh"

#include <map>
#include <set>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>



using namespace G4DNAPARSER;
using CLHEP::nanometer;

namespace
{
  G4Mutex messangerInit = G4MUTEX_INITIALIZER;

  // Track bad particle IDs globally (used for logging the warnings at end of run)
  std::set<int> g_badParticleIDs;
  G4Mutex g_badIDsMutex = G4MUTEX_INITIALIZER;
  G4int g_nBadPrimaryWarnings = 0; 
  G4Mutex g_warningCountMutex = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(G4String PS_data)
    : G4VUserPrimaryGeneratorAction(), fpParticleGun(nullptr), fPS_data(PS_data)
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0); // does this do anything?
  if ((IsPhaseSpaceInputActive(parser))||(parser->GetCommandIfActive("-photonPS")))
  {
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);
    fParticleGun->SetParticleEnergy(0);
    fParticleGun->SetParticlePosition(G4ThreeVector(0., 0., 0.));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));
  }
  else
  {
    fpParticleGun = new G4GeneralParticleSource();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  // Log all bad particle IDs when the run is complete
  if (!g_badParticleIDs.empty()) {
    G4cout << "\n========================================" << G4endl;
    G4cout << "SUMMARY: Bad Particle IDs Encountered" << G4endl;
    G4cout << "========================================" << G4endl;
    G4cout << "Total unique unknown particle IDs: " << g_badParticleIDs.size() << G4endl;
    G4cout << "Unknown IDs: ";
    for (auto id : g_badParticleIDs) {
      G4cout << id << " ";
    }
    G4cout << "\n========================================\n" << G4endl;
  }else {

    G4cout << "Great no Bad Particle IDs " << G4endl;

  }
  
  // handle cleanup of particle guns
  if (fParticleGun) delete fParticleGun;
  if (fpParticleGun) delete fpParticleGun;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0); // again, does this do anything?
  if (IsPhaseSpaceInputActive(parser))
  {
    G4long eventNum = anEvent->GetEventID();

    std::ifstream ps_file (fPS_data, std::ifstream::binary);

    if (fDecayPSRecordDoubles == 0)
    {
      ps_file.seekg(0, ps_file.end);
      const std::streamoff fileSize = ps_file.tellg();
      ps_file.seekg(0, ps_file.beg);

      // Support legacy 16-double records and new 18-double records.
      fDecayPSRecordDoubles = 16;
      if (fileSize % static_cast<std::streamoff>(18 * 8) == 0)
      {
        fDecayPSRecordDoubles = 18;
      }
      else if (fileSize % static_cast<std::streamoff>(16 * 8) != 0)
      {
        G4cout << "Warning: PSfile size not divisible by 16 or 18 doubles; assuming 16." << G4endl;
      }
    }

    const size_t recordDoubles = fDecayPSRecordDoubles;
    ps_file.seekg(eventNum * static_cast<std::streamoff>(recordDoubles * 8), ps_file.beg);

    std::vector<double> line(recordDoubles);

    ps_file.read(reinterpret_cast<char *>(line.data()), static_cast<std::streamsize>(recordDoubles * sizeof(double)));

    ps_file.close();
    G4double positionX = line[0];
    G4double positionY = line[1];
    G4double positionZ = line[2];
    G4double momentumX = line[3];
    G4double momentumY = line[4];
    G4double momentumZ = line[5];
    G4double particleEnergy = line[6];
    part1_EventNum = line[7];
    primaryParticle = line[8];
    part1_CopyNum = line[9];
    part1_Time = line[10];
    part1_particleSource = line[11];
    upstreamLocalPosition = G4ThreeVector(positionX, positionY, positionZ);
    upstreamWorldPosition = upstreamLocalPosition;
    upstreamSeedID = -1;
    upstreamTrackID = -1;
    upstreamParentID = -1;
    if (recordDoubles >= 18)
    {
      upstreamWorldPosition = G4ThreeVector(line[12], line[13], line[14]);
      upstreamSeedID = static_cast<G4int>(line[15]);
      upstreamTrackID = static_cast<G4int>(line[16]);
      upstreamParentID = static_cast<G4int>(line[17]);
    }

    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    // G4ParticleDefinition *particle;
    G4ParticleDefinition *particle = nullptr;  // <-- مهم جداً

  
    // Particle configuration map
    // Map for basic particles (e-, gamma, alpha, e+)
    static const std::map<int, std::string> basicParticleMap = {
      {1, "e-"},
      {2, "gamma"},
      {3, "alpha"},
      {11, "e+"}
    };

    // Map for ions: {particleID, {name, Z, A}}
    static const std::map<int, std::tuple<std::string, int, int>> ionMap = {
      {4, {"Rn220", 86, 220}}, 
      {5, {"Po216", 84, 216}}, 
      {6, {"Pb212", 82, 212}},
      {7, {"Bi212", 83, 212}}, 
      {8, {"Tl208", 81, 208}}, // this is duplicated 
      {9, {"Po212", 84, 212}},
      {10, {"Pb208", 82, 208}}, 
      {16, {"Ra226", 88, 226}}, 
      {17, {"Rn222", 86, 222}},
      {18, {"Po218", 84, 218}}, 
      {19, {"Pb214", 82, 214}}, 
      {20, {"Bi214", 83, 214}},
      {21, {"Po214", 84, 214}}, 
      {22, {"Pb210", 82, 210}}, 
      {23, {"Bi210", 83, 210}},
      {24, {"Po210", 84, 210}}, 
      {25, {"Pb206", 82, 206}}, 
      {26, {"Tl210", 81, 210}},
      {27, {"Tl208", 81, 208}},  // this is duplicated 
      {28, {"Hg206", 80, 206}}, 
      {29, {"At218", 85, 218}},
      {30, {"Rn218", 86, 218}}, 
      {31, {"Co60", 27, 60}}, 
      {32, {"Ra223", 88, 223}},
      {33, {"Rn219", 86, 219}}, 
      {34, {"Po215", 84, 215}}, 
      {35, {"Pb211", 82, 211}},
      {36, {"Bi211", 83, 211}}, 
      {37, {"Tl207", 81, 207}}
    };

    // BNCT upstream phase-space extension.  Geant4-DNA lithium models are
    // registered on generic lithium charge states, so use lithium+++ as the
    // replay primary rather than a normal G4IonTable Li7 ion.
    static const std::map<int, std::string> dnaGenericIonMap = {
      {53, "lithium+++"}
    };

    int particleID = static_cast<int>(primaryParticle);

    // First, check if it's a basic particle
    auto basicIt = basicParticleMap.find(particleID);
    if (basicIt != basicParticleMap.end()) {
      primaryName = basicIt->second;
      particle = particleTable->FindParticle(primaryName);
    }
    // Otherwise, check if it's a Geant4-DNA generic ion
    else {
      auto genericIt = dnaGenericIonMap.find(particleID);
      if (genericIt != dnaGenericIonMap.end()) {
        primaryName = "Li7";
        particle = G4DNAGenericIonsManager::Instance()->GetIon(genericIt->second);
      }
      // Otherwise, check if it's a normal ion
      else {
      auto ionIt = ionMap.find(particleID);
      if (ionIt != ionMap.end()) {
        // if using C++17 or later, can bind directly
        auto [name, Z, A] = ionIt->second;

        // This is the same as line 204, however it's compatible with C++ version 11/14, not sure which is used
        // std::string name = std::get<0>(ionIt->second);
        // int Z = std::get<1>(ionIt->second);
        // int A = std::get<2>(ionIt->second);

        primaryName = name;
        particle = G4IonTable::GetIonTable()->GetIon(Z, A, 0);
      }
      // if the particleID is not found in either map, particle remains nullptr so here we handle errors
      else {
        // Should have thread-safe tracking of bad particle IDs
        G4AutoLock lock(&g_badIDsMutex);
        g_badParticleIDs.insert(particleID);
        
        G4AutoLock warnLock(&g_warningCountMutex);
        if (g_nBadPrimaryWarnings < 50) {
          G4Exception("PrimaryGeneratorAction::GeneratePrimaries",
                  "BadPrimaryID",
                  JustWarning,
                  ("Unknown primaryParticle id: " + std::to_string(particleID) + 
                   " at event " + std::to_string(anEvent->GetEventID())).c_str());
      
                   
                   g_nBadPrimaryWarnings++;
                   
                   if (g_nBadPrimaryWarnings == 50) {
            G4cout << "~~~ Any further warnings will be printed as a summary at end of run. ~~~" << G4endl;
          }
        }
        return;  // Skip this event
      }
      }
    }

    // if particle is still null (shouldn't happen if tables are initialised) handle fatal error safely
    if (!particle) {
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries",
              "ParticleCreationFailed",
              FatalException,
              ("Failed to create particle: " + primaryName + 
               " (ID: " + std::to_string(particleID) + ")").c_str());
      return;
    }
    
    // G4cout << primaryName << " position = " << G4ThreeVector(positionX, positionY, positionZ) << " momentum = " << G4ThreeVector(momentumX, momentumY, momentumZ) << " energy = " << particleEnergy << "Ra event number = " << part1_EventNum << " copy number = " << part1_CopyNum << " source particle = " << part1_particleSource <<G4endl;
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticlePosition(G4ThreeVector(positionX, positionY, positionZ));
    fParticleGun->SetParticleEnergy(particleEnergy);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momentumX, momentumY, momentumZ));

    fParticleGun->GeneratePrimaryVertex(anEvent);

    if (parser->GetCommandIfActive("-out") != 0)
    {
      G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();
      analysisManager->FillNtupleIColumn(4, 0, static_cast<G4int>(eventNum));
      analysisManager->FillNtupleIColumn(4, 1, part1_EventNum);
      analysisManager->FillNtupleIColumn(4, 2, part1_CopyNum);
      analysisManager->FillNtupleIColumn(4, 3, part1_particleSource);
      analysisManager->FillNtupleIColumn(4, 4, upstreamSeedID);
      analysisManager->FillNtupleIColumn(4, 5, upstreamTrackID);
      analysisManager->FillNtupleIColumn(4, 6, upstreamParentID);
      analysisManager->FillNtupleDColumn(4, 7, part1_Time / CLHEP::ns);
      analysisManager->FillNtupleDColumn(4, 8, upstreamLocalPosition.x() / nanometer);
      analysisManager->FillNtupleDColumn(4, 9, upstreamLocalPosition.y() / nanometer);
      analysisManager->FillNtupleDColumn(4, 10, upstreamLocalPosition.z() / nanometer);
      analysisManager->FillNtupleDColumn(4, 11, upstreamWorldPosition.x() / nanometer);
      analysisManager->FillNtupleDColumn(4, 12, upstreamWorldPosition.y() / nanometer);
      analysisManager->FillNtupleDColumn(4, 13, upstreamWorldPosition.z() / nanometer);
      analysisManager->AddNtupleRow(4);
    }

  }
  else if (parser->GetCommandIfActive("-photonPS"))
  {
    G4long eventNum = anEvent->GetEventID();


    std::ifstream ps_file (fPS_data, std::ifstream::binary);

    ps_file.seekg(eventNum*8*8, ps_file.beg); // find position of event data. 8 doubles which are each 8 bytes

    double line[8];

    ps_file.read((char *)&line, sizeof line);

    ps_file.close();

    G4double positionX = line[0];
    G4double positionY = line[1];
    G4double positionZ = line[2];
    G4double momentumX = line[3];
    G4double momentumY = line[4];
    G4double momentumZ = line[5];
    G4double particleEnergy = line[6];
    part1_EventNum = line[7];

    // G4cout <<  " position = " << G4ThreeVector(positionX, positionY, positionZ) << " momentum = " << G4ThreeVector(momentumX, momentumY, momentumZ) << " energy = " << particleEnergy << " Co60 event number = " << part1_EventNum <<G4endl;


    G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition *particle = particleTable->FindParticle("e-");
    fParticleGun->SetParticleDefinition(particle);
    fParticleGun->SetParticlePosition(G4ThreeVector(positionX, positionY, positionZ));
    fParticleGun->SetParticleEnergy(particleEnergy);
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(momentumX, momentumY, momentumZ));

    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
  else
  {
    fpParticleGun->GeneratePrimaryVertex(anEvent);
    primaryName = fpParticleGun->GetCurrentSource()->GetParticleDefinition()->GetParticleName();
  }
}
