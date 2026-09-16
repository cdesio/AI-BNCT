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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "globals.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "G4LogicalVolume.hh"
#include "G4UnionSolid.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Ellipsoid.hh"
#include "G4PVParameterised.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4NistManager.hh"
#include "G4PhysicalVolumesSearchScene.hh"
#include <G4SystemOfUnits.hh>
#include "G4UserLimits.hh"
#include "G4UnitsTable.hh"
#include "G4LogicalVolumeStore.hh"
#include <fstream>
#include "CommandLineParser.hh"

#include "Randomize.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4RunManager.hh"

#include "G4ProductionCuts.hh"

using namespace G4DNAPARSER;

#define countof(x) (sizeof(x) / sizeof(x[0]))

using namespace std;
using CLHEP::angstrom;
using CLHEP::degree;
using CLHEP::micrometer;
using CLHEP::mm;
using CLHEP::nanometer;

static G4VisAttributes visInvBlue(false, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visInvWhite(false, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visInvPink(false, G4Colour(1.0, 0.0, 1.0));
static G4VisAttributes visInvCyan(false, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visInvRed(false, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visInvGreen(false, G4Colour(0.0, 1.0, 0.0));
static G4VisAttributes visBlue(true, G4Colour(0.0, 0.0, 1.0));
static G4VisAttributes visWhite(true, G4Colour(1.0, 1.0, 1.0));
static G4VisAttributes visBlack(true, G4Colour(0.0, 0.0, 0.0));
static G4VisAttributes visCyan(true, G4Colour(0.0, 1.0, 1.0));
static G4VisAttributes visRed(true, G4Colour(1.0, 0.0, 0.0));
static G4VisAttributes visGreen(true, G4Colour(0.0, 1.0, 0.0));
static G4VisAttributes visGrey(true, G4Colour(0.839216, 0.839216, 0.839216));

namespace
{
G4String AlternateGeometryPath(const G4String& fileName)
{
  const G4String buildPrefix = "../geometryFiles/";
  const G4String sourcePrefix = "geometryFiles/";

  if (fileName.compare(0, buildPrefix.size(), buildPrefix) == 0)
  {
    return sourcePrefix + fileName.substr(buildPrefix.size());
  }
  if (fileName.compare(0, sourcePrefix.size(), sourcePrefix) == 0)
  {
    return buildPrefix + fileName.substr(sourcePrefix.size());
  }
  return "";
}

G4bool OpenGeometryFile(const G4String& fileName, ifstream& file, G4String& openedPath)
{
  openedPath = fileName;
  file.open(openedPath);
  if (file.is_open()) return true;

  const G4String alternatePath = AlternateGeometryPath(fileName);
  if (alternatePath == "") return false;

  file.clear();
  openedPath = alternatePath;
  file.open(openedPath);
  return file.is_open();
}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction() : G4VUserDetectorConstruction()
{
  fBoxSize = 300 * nm;

  fDetectorMessenger = new DetectorMessenger(this);
  placeSugars = false; // only needed for visualisation
  placeHistones = true;

  fPositions0 = new G4KDTree();
  fPositions1 = new G4KDTree();
  fPositionsBase0 = new G4KDTree();
  fPositionsBase1 = new G4KDTree();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{

  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  /***************************************************************************/
  //                               World
  /***************************************************************************/

  G4NistManager *man = G4NistManager::Instance();
  G4Material *waterMaterial = man->FindOrBuildMaterial("G4_WATER");
  G4Material *vacuum = man->FindOrBuildMaterial("G4_Galactic");
  G4Material *air = man->FindOrBuildMaterial("G4_AIR");               // not used but added so that ICRU90 DEDX function works
  G4Material *steel = man->FindOrBuildMaterial("G4_STAINLESS-STEEL"); // not used but added so that ICRU90 DEDX function works

  G4Box *solidWorld{nullptr};
  G4Orb *solidSurroundingMaterial{nullptr};

  solidWorld = new G4Box("world", 1 * mm, 1 * mm, 1 * mm);
  solidSurroundingMaterial = new G4Orb("surroundingMaterial", 1 * um); 
  G4Material *surroundingMaterial;

  CommandLineParser *parser = CommandLineParser::GetParser();
  Command *command(0);
  if ((IsPhaseSpaceInputActive(parser))||(parser->GetCommandIfActive("-photonPS")))
  {
    surroundingMaterial = waterMaterial;
    G4cout << "surrounding = water" << G4endl;
  }
  else
  {
    surroundingMaterial = waterMaterial;
    G4cout << "surrounding = water" << G4endl;

  }

  G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld,
                                                    vacuum,
                                                    "world");

  G4PVPlacement *physiWorld = new G4PVPlacement(0,
                                                G4ThreeVector(),
                                                "world",
                                                logicWorld,
                                                0,
                                                false,
                                                0);

  G4LogicalVolume *logicSurrounding = new G4LogicalVolume(solidSurroundingMaterial,
                                                          surroundingMaterial,
                                                          "surroundingMaterial");

  G4PVPlacement *physiSurrounding = new G4PVPlacement(0,
                                                      G4ThreeVector(),
                                                      logicSurrounding,
                                                      "surroundingMaterial",
                                                      logicWorld,
                                                      0,
                                                      false,
                                                      0);
  logicWorld->SetVisAttributes(&visInvWhite);
  logicSurrounding->SetVisAttributes(&visInvWhite);

  G4double chromatin_X = fBoxSize;
  G4double chromatin_Y = fBoxSize;
  G4double chromatin_Z = fBoxSize;

  G4Box *solidTrackingVol = new G4Box("TrackingVol", 20 * nm + chromatin_X / 2, 20 * nm + chromatin_Y / 2, 20 * nm + chromatin_Z / 2); // volume to track radicals in dkill larger than chromatin
  G4LogicalVolume *logicTrackingVol = new G4LogicalVolume(solidTrackingVol,
                                                          waterMaterial,
                                                          "TrackingVol");
  G4PVPlacement *physiTrackingVol = new G4PVPlacement(0,
                                                      G4ThreeVector(),
                                                      logicTrackingVol,
                                                      "TrackingVol",
                                                      logicSurrounding,
                                                      0,
                                                      false,
                                                      false);
  logicTrackingVol->SetVisAttributes(&visInvWhite);

  G4Box *solidChromatinSegment = new G4Box("chromatinSegment", chromatin_X / 2, chromatin_Y / 2, chromatin_Z / 2);
  G4LogicalVolume *logicChromatinSegment = new G4LogicalVolume(solidChromatinSegment,
                                                               waterMaterial,
                                                               "chromatinSegment");

  chromatinVolume = chromatin_X / m * chromatin_Y / m * chromatin_Z / m; // in m3

  G4PVPlacement *physiChromatinSegment = new G4PVPlacement(0,
                                                           G4ThreeVector(),
                                                           logicChromatinSegment,
                                                           "chromatinSegment",
                                                           logicTrackingVol,
                                                           0,
                                                           0,
                                                           false);
  logicChromatinSegment->SetVisAttributes(&visBlack);

  G4double sugarRadius{2.9 * angstrom};

  G4Orb *solidSugar = new G4Orb("sugar", sugarRadius);
  G4LogicalVolume *logicSugar0 = new G4LogicalVolume(solidSugar,
                                                     waterMaterial,
                                                     "sugar0");

  G4LogicalVolume *logicSugar1 = new G4LogicalVolume(solidSugar,
                                                     waterMaterial,
                                                     "sugar1");

  G4Ellipsoid *solidHistone = new G4Ellipsoid("histone", 2.5 * nanometer, 2.5 * nanometer, 2.5 * nanometer);

  G4LogicalVolume *logicHistone = new G4LogicalVolume(solidHistone,
                                                      waterMaterial,
                                                      "histone");

  const G4String sugarFileName = GetCommandOptionOrDefault(CommandLineParser::GetParser(), "-sugar");

  ifstream f;
  G4String openedSugarFileName;

  if (OpenGeometryFile(sugarFileName, f, openedSugarFileName))
  {
    G4cout << "Sugar file opened OK: " << openedSugarFileName << G4endl;
  }
  else
  {
    throw std::runtime_error("*******SUGAR GEOMETRY FILE NOT FOUND OR OPENED CORRECTY********");
  }

  G4int copynum{0};
  float line[12];

  while (f.read((char *)&line, sizeof line))
  {
    float x1 = line[0], y1 = line[1], z1 = line[2], x2 = line[3], y2 = line[4], z2 = line[5], base_x1 = line[6], base_y1 = line[7], base_z1 = line[8], base_x2 = line[9], base_y2 = line[10], base_z2 = line[11];
    fPositions0->Insert(G4ThreeVector(x1*nm,
                                      y1*nm,
                                      z1*nm));
    fPositions1->Insert(G4ThreeVector(x2*nm,
                                      y2*nm,
                                      z2*nm));
    fPositionsBase0->Insert(G4ThreeVector(base_x1*nm,
                                          base_y1*nm,
                                          base_z1*nm));
    fPositionsBase1->Insert(G4ThreeVector(base_x2*nm,
                                          base_y2*nm,
                                          base_z2*nm));
    if (placeSugars)
    {

      new G4PVPlacement(0,
                        G4ThreeVector(x1*nm,
                                      y1*nm,
                                      z1*nm),
                        logicSugar0,
                        "sugar0",
                        logicChromatinSegment,
                        false,
                        copynum,
                        false);

      new G4PVPlacement(0,
                        G4ThreeVector(x2*nm,
                                      y2*nm,
                                      z2*nm),
                        logicSugar1,
                        "sugar1",
                        logicChromatinSegment,
                        false,
                        copynum,
                        false);
    }
    ++copynum;
  }

  // Histones
  const G4String histoneFileName = GetCommandOptionOrDefault(CommandLineParser::GetParser(), "-histone");

  ifstream f2;
  G4String openedHistoneFileName;

  if (OpenGeometryFile(histoneFileName, f2, openedHistoneFileName))
  {
    G4cout << "Histone file opened OK: " << openedHistoneFileName << G4endl;
  }
  else
  {
    throw std::runtime_error("*******HISTONE GEOMETRY FILE NOT FOUND OR OPENED CORRECTY********");
  }

  float line2[6];

  while (f2.read((char *)&line2, sizeof line2))
  {
    float x1 = line2[0], y1 = line2[1], z1 = line2[2], rotx = line2[3], roty = line2[4], rotz = line2[5];
    fPositions.push_back(G4ThreeVector(x1*nm,
                                       y1*nm,
                                       z1*nm));
    fRotations.push_back(new G4RotationMatrix(roty, // geant convention is y,z,x
                                              rotz,
                                              rotx));
  }

  if (placeHistones)
  {
    for (size_t i = 0; i < fPositions.size(); ++i)
    {
      new G4PVPlacement(fRotations[i],
                        fPositions[i],
                        logicHistone,
                        "histone",
                        logicChromatinSegment,
                        false,
                        i,
                        false);
    }
  }

  numSugar = copynum;

  logicSugar0->SetVisAttributes(&visWhite);
  logicSugar1->SetVisAttributes(&visRed);
  logicHistone->SetVisAttributes(&visGrey);

  G4Region *aRegion = new G4Region("Target");

  G4ProductionCuts *cuts = new G4ProductionCuts();

  G4double defCut = 1 * nanometer;
  cuts->SetProductionCut(defCut, "gamma");
  cuts->SetProductionCut(defCut, "e-");
  cuts->SetProductionCut(defCut, "e+");
  cuts->SetProductionCut(defCut, "proton");

  aRegion->SetProductionCuts(cuts);

  aRegion->AddRootLogicalVolume(logicTrackingVol); // all volumes placed within are included in the region

  return physiWorld;
}
void DetectorConstruction::SetSize(G4double value)
{
  fBoxSize = value;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
