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

#ifndef G4EmDNAPhysics_option2_Jinyan_h
#define G4EmDNAPhysics_option2_Jinyan_h 1

#include "G4EmDNAPhysics_option2.hh"
#include "globals.hh"

// 前置声明：在 .cc 中会用到的类
class G4DNAIonisation;
class G4DNAExcitation;
class TsDNAKRIonIonisationScaledModel;
class TsDNAKRIonMillerGreenExcitationModel;
class TsDNAKRIonBornExcitationModel;
class TsDNAChargeIncrease;
class TsDNAChargeDecrease;

// 自定义 Physics 列表，继承自 option2，只添加 lithium 过程
class G4EmDNAPhysics_option2_Jinyan : public G4EmDNAPhysics_option2 {
public:
  explicit G4EmDNAPhysics_option2_Jinyan(
      G4int ver = 1,
      const G4String& name = "G4EmDNAPhysics_option2_Jinyan");
  ~G4EmDNAPhysics_option2_Jinyan() override = default;

  // 调用父类构造标准过程后，插入 lithium 专属过程
  void ConstructProcess() override;
};

#endif // G4EmDNAPhysics_option2_Jinyan_h

