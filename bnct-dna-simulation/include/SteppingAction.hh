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
#pragma once
#include "G4UserSteppingAction.hh"
#include "G4String.hh"
#include <map>
class EventAction;
class G4ParticleDefinition;
class G4VPhysicalVolume;
class DetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(DetectorConstruction* fpDet);
    ~SteppingAction() override;

    void UserSteppingAction(const G4Step* step) override;
    // void Initialize();
private:
  EventAction* fpEventAction;
  DetectorConstruction* fpDetector;

    void RemoveTracks(const G4Step* step);

    std::map<G4String, int> particleID = {
        { "e-", 1 },
        { "gamma", 2 }, 
        { "alpha", 3 },
        { "lithium+++", 53 },
        { "lithium++", 54 },
        { "lithium+", 55 },
        { "lithium", 56 },
        { "Li7", 53 },
        { "Rn220", 4 },
        { "Po216", 5 }, 
        { "Pb212", 6 }, 
        { "Bi212", 7 }, 
        { "Tl208", 8 }, 
        { "Po212", 9 }, 
        { "Pb208", 10 }, 
        { "alpha+", 11 }, 
        { "helium", 12 }, 

        { "At211", 13 }, 
        { "Po211", 14 }, 
        { "Bi207", 15 }, 
        { "Ra226", 16 },
        { "Rn222", 17 },
        { "Po218", 18 },
        { "Pb214", 19 },
        { "Bi214", 20 },
        { "Po214", 21 },
        { "Pb210", 22 },
        { "Bi210", 23 },
        { "Po210", 24 },
        { "Pb206", 25 },
        { "Tl210", 26 },
        { "Tl206", 27 },
        { "Hg206", 28 },
        { "At218", 29 },
        { "Rn218", 30 },

        { "Co60", 31 },



        { "Ra223", 32 },
        { "Rn219", 33 },
        { "Po215", 34 },
        { "Pb211", 35 },
        { "Bi211", 36 },
        { "Tl207", 37 },





        { "Ra225", 38 },
        { "Rn221", 39 },
        { "Po217", 40 },
        { "Pb213", 41 },
        { "Ac225", 42 },
        { "Fr221", 43 },
        { "Ra221", 44 },
        { "Rn217", 45 },
        { "Po213", 46 },
        { "Pb209", 47 },
        { "Bi209", 48 },
        { "Tl205", 49 },
        { "At217", 50 },
        { "Bi213", 51 },
        { "Tl209", 52 },
      





};

};
