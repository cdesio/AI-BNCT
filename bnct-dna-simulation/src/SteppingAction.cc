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
#include "SteppingAction.hh"
#include "G4AnalysisManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4ITTrackHolder.hh"
#include "G4Track.hh"
#include <map>
#include "globals.hh"
#include "G4Molecule.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4DNAElastic.hh"
#include "G4DNAElectronSolvation.hh"
#include "DetectorConstruction.hh"
#include "CommandLineParser.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4KDTree.hh"
#include "G4KDTreeResult.hh"
using MapOfDelayedLists =
    std::map<double, std::map<int, G4TrackList *>>;
using namespace G4DNAPARSER;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::SteppingAction(DetectorConstruction* fpDet)
    : G4UserSteppingAction(), fpEventAction(0)
, fpDetector(fpDet)
{
  fpEventAction = (EventAction *)G4EventManager::GetEventManager()->GetUserEventAction();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

SteppingAction::~SteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void SteppingAction::UserSteppingAction(const G4Step *step)
{
  G4double dE = 0;
  G4int flagVolume = 0.;

  const G4String &volumeName = step->GetPreStepPoint()->GetPhysicalVolume()->GetName();

  if (volumeName == "sugar0")
  {
    flagVolume = 1;
  }
  else if (volumeName == "sugar1")
  {
    flagVolume = 2;
  }
  else if (volumeName == "histone")
  {
    flagVolume = 3;
  }
  else if (volumeName == "chromatinSegment")
  {
    flagVolume = 4;
  }
  else if (volumeName == "TrackingVol")
  {
    flagVolume = 5;
  }
  if (flagVolume != 0)
  {

    G4String particleName = step->GetTrack()->GetParticleDefinition()->GetParticleName();
    const PrimaryGeneratorAction *generatorAction = static_cast<const PrimaryGeneratorAction *>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    G4String primaryName = generatorAction->primaryName;

    if ((((particleName == "alpha") || (particleName == "alpha+") || (particleName == "helium")) && (primaryName == "alpha")) || (((particleName == "lithium+++") || (particleName == "lithium++") || (particleName == "lithium+") || (particleName == "lithium")) && (primaryName == "Li7")) || ((primaryName == "e-") && (step->GetTrack()->GetTrackID() == 1)) || (((particleName == "proton") || (particleName == "hydrogen")) && (primaryName == "proton"))||((primaryName == "gamma") && (step->GetTrack()->GetTrackID() == 1))) // add gamma?

    {
      // set flag if left box, so that primary scattering back can be removed as these are counted in part 1
      if ((volumeName == "TrackingVol")&&(step->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="chromatinSegment")&&(fpEventAction->GetLeftBox()==true))
      {
        step->GetTrack()->SetTrackStatus(fStopAndKill);
        return;
      }
      if ((false == fpEventAction->GetStartTrackFound())&&(step->GetPreStepPoint()->GetPhysicalVolume()->GetName()=="chromatinSegment"))
      {
        fpEventAction->SetStartTrackKE(step->GetPreStepPoint()->GetKineticEnergy());
        fpEventAction->SetStartTrackFound();
        fpEventAction->SetStartTrackPos(step->GetPreStepPoint()->GetPosition());
      }
      if ((volumeName == "chromatinSegment")&&(step->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="TrackingVol")&&(fpEventAction->GetStartTrackFound())) //needs to be after enters box because a single step can enter and leave
      {
        fpEventAction->SetLeftBox();
        fpEventAction->SetEndTrackKE(step->GetPostStepPoint()->GetKineticEnergy());
        fpEventAction->SetEndTrackPos(step->GetPostStepPoint()->GetPosition());
      }
      if ((step->GetPostStepPoint()->GetPhysicalVolume()->GetName()=="chromatinSegment")&&(step->GetPostStepPoint()->GetKineticEnergy() == 0))
      {
        fpEventAction->SetStoppedInBox();
        // G4cout << "primary stopped in box" << G4endl;

        fpEventAction->SetEndTrackKE(step->GetPostStepPoint()->GetKineticEnergy());
        fpEventAction->SetEndTrackPos(step->GetPostStepPoint()->GetPosition());
        fpEventAction->AddPathLength(step->GetStepLength());
      }
      if (volumeName=="chromatinSegment")
      {
        fpEventAction->AddPathLength(step->GetStepLength());
      }
    }
  }

  dE = step->GetTotalEnergyDeposit();
  if ((volumeName != "world") && (volumeName != "surroundingMaterial") && (volumeName != "TrackingVol") && (dE != 0))
  {

    fpEventAction->AddEdep(dE);
  }
  if ((flagVolume == 0)||(flagVolume == 5))
  {
    return;
  }

  // remove water molecules created in the sugar volumes or histones
  if ((flagVolume == 1) || (flagVolume == 2) || (flagVolume == 3))
  {
    RemoveTracks(step);
  }

  if ((flagVolume != 3) && (dE > 0))
  {
    CommandLineParser *parser = CommandLineParser::GetParser();
    // Command *command(0);
    if (parser->GetCommandIfActive("-out") == 0)
      return;

    G4AnalysisManager *analysisManager = G4AnalysisManager::Instance();

    G4int eventID = G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetEventID();

    G4ThreeVector prePoint = step->GetPreStepPoint()->GetPosition();
    G4ThreeVector postPoint = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector point = prePoint + G4UniformRand() * (postPoint - prePoint);

    // only save if within 1 nm of sugar (this will allow Rdirect to be varied (typically 0.35nm) but reduce file size)
    auto fPositions0 = ((DetectorConstruction *)fpDetector)->fPositions0;
    auto fPositionsBase0 = ((DetectorConstruction *)fpDetector)->fPositionsBase0;
    auto fPositions1 = ((DetectorConstruction *)fpDetector)->fPositions1;
    auto fPositionsBase1 = ((DetectorConstruction *)fpDetector)->fPositionsBase1;
    auto result0 = fPositions0->NearestInRange(point, 1*nm); 
    auto result1 = fPositions1->NearestInRange(point, 1*nm); 

    if ((result0->GetSize() == 0)&&(result1->GetSize() == 0)) return;

    const PrimaryGeneratorAction *generatorAction = static_cast<const PrimaryGeneratorAction *>(
        G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());

    G4int part1_CopyNum = generatorAction->part1_CopyNum;
    G4double part1_Time = generatorAction->part1_Time;
    G4int part1_particleSource = generatorAction->part1_particleSource;
    G4int part1_EventNum = generatorAction->part1_EventNum;
    G4int upstreamSeedID = generatorAction->upstreamSeedID;
    G4int upstreamTrackID = generatorAction->upstreamTrackID;
    G4int upstreamParentID = generatorAction->upstreamParentID;

    analysisManager->FillNtupleIColumn(1, 0, eventID);
    analysisManager->FillNtupleDColumn(1, 1, dE / eV);
    analysisManager->FillNtupleDColumn(1, 2, point.x() / nanometer);
    analysisManager->FillNtupleDColumn(1, 3, point.y() / nanometer);
    analysisManager->FillNtupleDColumn(1, 4, point.z() / nanometer);
    analysisManager->FillNtupleIColumn(1, 5, particleID[step->GetTrack()->GetParticleDefinition()->GetParticleName()]);
    if (IsPhaseSpaceInputActive(parser))
    {
      analysisManager->FillNtupleIColumn(1, 6, part1_CopyNum);
      analysisManager->FillNtupleDColumn(1, 7, part1_Time + (step->GetTrack()->GetGlobalTime()));
      analysisManager->FillNtupleIColumn(1, 8, part1_particleSource); // save primary name so that DNA damage can be calculated per incoming particle
      analysisManager->FillNtupleIColumn(1, 9, part1_EventNum);
      analysisManager->FillNtupleIColumn(1, 10, upstreamSeedID);
      analysisManager->FillNtupleIColumn(1, 11, upstreamTrackID);
      analysisManager->FillNtupleIColumn(1, 12, upstreamParentID);
    }
    analysisManager->AddNtupleRow(1);
  }
}

void SteppingAction::RemoveTracks(const G4Step *step)
{
  MapOfDelayedLists delayList = G4ITTrackHolder::Instance()->GetDelayedLists();
  for (auto &delayedmap_it : delayList)
  {
    for (auto &trackList : delayedmap_it.second)
    {
      if (nullptr == trackList.second)
      {
        continue;
      }
      G4TrackList::iterator itt = trackList.second->begin();
      G4TrackList::iterator endd = trackList.second->end();
      for (; itt != endd; ++itt)
      {
        G4Track *track = *itt;

        if ((track->GetParentID() !=
             step->GetTrack()->GetTrackID()) ||
            (track->GetPosition() !=
             step->GetPostStepPoint()->GetPosition()))
        {
          continue;
        }
        track->SetTrackStatus(fKillTrackAndSecondaries);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
