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

// SI: This constructor uses speedup options of DNA models

#include "G4EmDNAPhysics_option2_Jinyan.hh"
#include "G4EmDNABuilder.hh"
#include "G4SystemOfUnits.hh"

// ions
#include "G4Alpha.hh"
#include "G4DNAGenericIonsManager.hh"

// utilities
#include "G4EmParameters.hh"
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmBuilder.hh"

// TsDNA-lithium
#include "TsDNAKRIonIonisationScaledModel.hh"
#include "TsDNAKRIonMillerGreenExcitationModel.hh"
#include "TsDNAKRIonBornExcitationModel.hh"
#include "TsDNAChargeIncrease.hh"
#include "TsDNAChargeDecrease.hh"

#include "G4DNAGenericIonsManager.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option2_Jinyan);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics_option2_Jinyan::G4EmDNAPhysics_option2_Jinyan(G4int ver, const G4String& nam)
  : G4EmDNAPhysics_option2(ver, nam)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDNAFast(true);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics_option2_Jinyan::ConstructProcess()
{

   G4EmDNAPhysics_option2::ConstructProcess();

   // lithium
   auto ph = G4PhysicsListHelper::GetPhysicsListHelper();
   auto ionsMgr = G4DNAGenericIonsManager::Instance();
 
   const std::vector<G4String> liNames = {
     "lithium+++", "lithium++", "lithium+", "lithium"
   };
   for (auto& pname : liNames) {
     auto particle = ionsMgr->GetIon(pname);
     if (!particle) continue;
 
     //  Ionisation
     G4DNAIonisation* liIoni =
       new G4DNAIonisation("Li_G4DNAIonisation");
     auto ionModel = new TsDNAKRIonIonisationScaledModel();
     ionModel->SelectStationary(false);
     liIoni->AddEmModel(1, ionModel);
     ph->RegisterProcess(liIoni, particle);
 
     //  Excitation
     G4DNAExcitation* liExc =
       new G4DNAExcitation("Li_G4DNAExcitation");
     // low energy
     auto ex1 = new TsDNAKRIonMillerGreenExcitationModel();
     ex1->SetLowEnergyLimit(70*eV);
     ex1->SetHighEnergyLimit(3.5*MeV);
     ex1->SelectStationary(false);
     liExc->AddEmModel(1, ex1);
     // high energy
     auto ex2 = new TsDNAKRIonBornExcitationModel();
     ex2->SetLowEnergyLimit(3.5*MeV);
     ex2->SetHighEnergyLimit(700*MeV);
     ex2->SelectStationary(false);
     liExc->AddEmModel(2, ex2);
     ph->RegisterProcess(liExc, particle);
 
     //  Charge transfer 
     ph->RegisterProcess(
       new TsDNAChargeIncrease("Li_TsDNAChargeIncrease"),
       particle);
     ph->RegisterProcess(
       new TsDNAChargeDecrease("Li_TsDNAChargeDecrease"),
       particle);
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
