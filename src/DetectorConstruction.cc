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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4PhysicalConstants.hh"  //twopi
#include "G4LogicalSkinSurface.hh"
#include "G4MaterialPropertiesTable.hh"
#include "G4GenericTrap.hh"
#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),fScoringVolume(nullptr),fScoringVolume2(nullptr)
{
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
     delete vacum;
     delete air;
     delete pb;
     delete cu;
     delete cuLess;
     delete pbMore;
     delete kapton;
     delete kaptonLess;
     delete gasMixture;

     delete Drift_cathode_Mat;
     delete Gem_inner_Mat;
     delete Gem_outter_Mat;
     delete Readout_plate_Mat;
     delete Shell_Mat;
     delete Readout_bar_Mat;
     delete Gem_Mat;
     delete world_Mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  G4cout << "Defining the materials" << G4endl;
  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();
  // Build materials
  vacum   = nistManager->FindOrBuildMaterial("G4_Galactic");
  air   = nistManager->FindOrBuildMaterial("G4_AIR");
  pb   = nistManager->FindOrBuildMaterial("G4_Pb");
  fe   = nistManager->FindOrBuildMaterial("G4_Fe");
  w   = nistManager->FindOrBuildMaterial("G4_W");
  cu   = nistManager->FindOrBuildMaterial("G4_Cu");

  G4cout<<"define high density Pb"<<G4endl;
  // define high density Pb
  G4int nComponent;
  G4double fracmass=1.0;
  pbMore   = new G4Material("pbMore",11.34*1*g/cm3,nComponent=1);
  pbMore -> AddMaterial(pb,fracmass);

  G4cout<<"define low density Cu"<<G4endl;
  // define low density Cu
  G4double ratio_outter = 0.77790212; // 1-(3600*pi)/(29400*sqrt(3))
  G4double density = ratio_outter * 8.94 * g/cm3;
  cuLess   = new G4Material("cuLess",density,nComponent=1);
  cuLess -> AddMaterial(cu,fracmass);

  G4cout<<"define Kapton"<<G4endl;
  // Define Kapton material
  density = 1.42 * g/cm3; // Density of Kapton
  kapton = new G4Material("Kapton", density, 3);
  kapton->AddElement(nistManager->FindOrBuildElement("H"), 0.0273);
  kapton->AddElement(nistManager->FindOrBuildElement("C"), 0.7213);
  kapton->AddElement(nistManager->FindOrBuildElement("O"), 0.2514);

  G4double ratio_inner = 0.81183374; // 1-(3600*pi+2500*pi)/2/(29400*sqrt(3))
  density = ratio_inner * 1.42 * g/cm3; // Density of Kapton
  kaptonLess = new G4Material("KaptonLess", density, 3);
  kaptonLess->AddElement(nistManager->FindOrBuildElement("H"), 0.0273);
  kaptonLess->AddElement(nistManager->FindOrBuildElement("C"), 0.7213);
  kaptonLess->AddElement(nistManager->FindOrBuildElement("O"), 0.2514);

  G4cout<<"define mixgass"<<G4endl;
  // Define 70%Ar and 30%CO2 gas mixture
  density = 1.822 * mg/cm3; // Density of the gas mixture
  G4Material* gasMixture = new G4Material("ArCO2GasMixture", density, 3);
  gasMixture->AddElement(nistManager->FindOrBuildElement("Ar"), 0.7);
  gasMixture->AddElement(nistManager->FindOrBuildElement("C"), 0.3/3);
  gasMixture->AddElement(nistManager->FindOrBuildElement("O"), 0.3/3*2);

  G4cout<<"define detector material"<<G4endl;
  // Define detector material
  Drift_cathode_Mat = cu;
  Gem_inner_Mat = kaptonLess;
  Gem_outter_Mat = cuLess;
  Readout_plate_Mat = cu;
  Shell_Mat = kapton;
  Readout_bar_Mat = cu;
  Gem_Mat = gasMixture;
  world_Mat = vacum;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineConstants()
{
  Gem_outter_x=1*m;
  Gem_outter_y=1*m;
  Gem_outter_z=5*um;

  Gem_inner_x=1*m;
  Gem_inner_y=1*m;
  Gem_inner_z=50*um;

  drift_cathode_x=1*m;
  drift_cathode_y=1*m;
  drift_cathode_z=0.1*mm;

  readout_bar_gap_x=210*um;
  readout_bar_x=150*um;
  readout_bar_y=1*m;
  readout_bar_z=0.1*mm;

  readout_plate_x=1*m;
  readout_plate_y=1*m;
  readout_plate_z=0.1*mm;

  gap1 = 4.8*mm;
  gap2 = 2*mm;
  Gem_x=1*m;
  Gem_y=1*m;
  Gem_z=
        gap1
        +num_Gem_inner*gap2
        +(num_Gem_outter*Gem_outter_z+num_Gem_inner*Gem_inner_z)
        +drift_cathode_z
        +readout_bar_z
        +readout_plate_z;

  Box_x=1*m;
  Box_y=1*m;
  Box_z=1*m;

  num_Gem=4;

  experimentalHall_x=1.1*Box_x;
  experimentalHall_y=1.1*Box_y;
  experimentalHall_z=1.1*(Box_z+Gem_z*num_Gem);
 
  Pbplate_x=100*mm;
  Pbplate_y=100*mm;
  Pbplate_z=0.1*mm;

  Pbbox_x=0.5*m;
  Pbbox_y=0.5*m;
  //Pbbox_z=0.5*m;
  //Pbbox_z=20*mm;
  Pbbox_z=1*mm;

/*
  Pbbox_x=20*mm;
  Pbbox_z=0.5*m;
  Pbbox_y=0.5*m;
*/

/*
  Pbbox_x=0.3*m;
  Pbbox_y=0.3*m;
  Pbbox_z=0.3*m;
*/

  Febox_x=0.3*m;
  Febox_y=0.3*m;
  Febox_z=0.3*m;

  Wbox_x=0.3*m;
  Wbox_y=0.3*m;
  Wbox_z=0.3*m;

  Vacbox_x=0.5*m-40*mm;
  Vacbox_y=0.5*m-40*mm;
  Vacbox_z=0.5*m-40*mm;

  pku_box_x = 30*cm;
  pku_box_y = 40*cm;
  pku_box_z = 3*cm;

  pku_bar_x = 3*cm;
  pku_bar_y1 = 20*cm;
  pku_bar_y2 = 40*cm;
  pku_bar_z = 3*cm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume *DetectorConstruction::DefineVolumes()
{

    /*******************************
   * Gem outter plate       *
   *******************************/
  G4VSolid* Gem_outter_box
    = new G4Box("Gem_outter_box",             // World Volume
                Gem_outter_x/2,        // x size
                Gem_outter_y/2,        // y size
                Gem_outter_z/2);       // z size

  G4LogicalVolume* Gem_outter_Log
    = new G4LogicalVolume(Gem_outter_box,
			  Gem_outter_Mat,
                          "Gem_outter_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  Gem_outter_Log->SetVisAttributes(G4VisAttributes::GetInvisible());

    /*******************************
   * Gem inner plate       *
   *******************************/
  G4VSolid* Gem_inner_box
    = new G4Box("Gem_inner_box",             // World Volume
                Gem_inner_x/2,        // x size
                Gem_inner_y/2,        // y size
                Gem_inner_z/2);       // z size

  G4LogicalVolume* Gem_inner_Log
    = new G4LogicalVolume(Gem_inner_box,
			  Gem_inner_Mat,
                          "Gem_inner_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  Gem_inner_Log->SetVisAttributes(G4VisAttributes::GetInvisible());

    /*******************************
   * Drift cathode       *
   *******************************/
  G4VSolid* drift_cathode_box 
    = new G4Box("dricath_box",             // World Volume
                drift_cathode_x/2,        // x size
                drift_cathode_y/2,        // y size
                drift_cathode_z/2);       // z size
  
  G4LogicalVolume* driftcathodeLog 
    = new G4LogicalVolume(drift_cathode_box,
			  Drift_cathode_Mat,
                          "dricathLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  G4VisAttributes* driftcathodeLogAtt = new G4VisAttributes(G4Colour(233/256.,0/256.,0/256,0.7));
  //driftcathodeLog->SetVisAttributes(driftcathodeLogAtt);
  driftcathodeLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    /*******************************
   * Readout bar       *
   *******************************/
  G4VSolid* readout_bar_box 
    = new G4Box("readoutbar_box",             // World Volume
                readout_bar_x/2,        // x size
                readout_bar_y/2,        // y size
                readout_bar_z/2);       // z size

  G4LogicalVolume* readoutbarLog
    = new G4LogicalVolume(readout_bar_box,
			  Readout_bar_Mat,
                          "readoutbarLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  readoutbarLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    /*******************************
   * Readout plate       *
   *******************************/
  G4VSolid* readout_plate_box 
    = new G4Box("readoutplate_box",             // World Volume
                readout_plate_x/2,        // x size
                readout_plate_y/2,        // y size
                readout_plate_z/2);       // z size
  
  G4LogicalVolume* readoutplateLog 
    = new G4LogicalVolume(readout_plate_box,
			  Readout_plate_Mat,
                          "readoutplateLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  G4VisAttributes* readoutplateLogAtt = new G4VisAttributes(G4Colour(233/256.,0/256.,0/256,0.7));
  //readoutplateLog->SetVisAttributes(readoutplateLogAtt);
  readoutplateLog->SetVisAttributes(G4VisAttributes::GetInvisible());


    /*******************************
   * The Gem       *
   *******************************/

  G4VSolid* Gem_box 
    = new G4Box("Gem_box",             // World Volume
                Gem_x/2,        // x size
                Gem_y/2,        // y size
                Gem_z/2);       // z size
  
  G4LogicalVolume* GemLog 
    = new G4LogicalVolume(Gem_box,
			  Gem_Mat,
                          "GemLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  //GemLog->SetVisAttributes(G4VisAttributes::GetInvisible());

    // put Gem outter plate in Gem
    //

  Zcenter = -0.5*Gem_z + drift_cathode_z + gap1 - gap2 - 0.5*Gem_outter_z;
  for (int i=0; i<num_Gem_outter; i++){
    if(i%2==0) Zcenter = Zcenter + gap2 + Gem_outter_z;
    if(i%2==1) Zcenter = Zcenter + Gem_inner_z + Gem_outter_z;
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    Gem_outter_Log,             //its logical volume
                    "Gem_outter",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    1);          //overlaps checking
  }

    // put Gem inner plate in Gem
    //

  Zcenter = -0.5*Gem_z + drift_cathode_z + gap1 - gap2 - Gem_outter_z - 0.5*Gem_inner_z;
  for (int i=0; i<num_Gem_inner; i++){
    Zcenter = Zcenter + gap2 + 2*Gem_outter_z + Gem_inner_z;
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    Gem_inner_Log,             //its logical volume
                    "Gem_inner",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    1);          //overlaps checking
  }

  // put drift cathode in Gem
  //
  Zcenter = -0.5*(Gem_z-drift_cathode_z);
  G4VPhysicalVolume* driftcathodePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    driftcathodeLog,             //its logical volume
                    "driftcathode",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking

  // put readout bar in Gem
  //
  const G4int num_readout_bar=Gem_x/readout_bar_gap_x;

  Zcenter = 0.5*Gem_z-readout_plate_z-0.5*readout_bar_z;
  Xcenter = -0.5*(Gem_x-readout_bar_x)-readout_bar_gap_x;
  for (int i=0; i<num_readout_bar; i++){

    Xcenter = Xcenter + readout_bar_gap_x;
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Xcenter,0,Zcenter),                    //at position
                    readoutbarLog,             //its logical volume
                    "readoutbar",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    0);          //overlaps checking
  }

  // put readout plate in Gem
  //
  Zcenter = 0.5*(Gem_z-readout_plate_z);
  G4VPhysicalVolume* readoutplatePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    readoutplateLog,             //its logical volume
                    "readoutplate",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking

    /*******************************
   * The Pb plate       *
   *******************************/

  G4VSolid* Pbplate_box 
    = new G4Box("Pbplate_box",             // World Volume
                Pbplate_x/2,        // x size
                Pbplate_y/2,        // y size
                Pbplate_z/2);       // z size
  
  G4LogicalVolume* PbplateLog 
    = new G4LogicalVolume(Pbplate_box,
			  pb,
                          "PbplateLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

    /*******************************
   * The Pb box       *
   *******************************/

  G4VSolid* Pbbox_box 
    = new G4Box("Pbbox_box",             // World Volume
                Pbbox_x/2,        // x size
                Pbbox_y/2,        // y size
                Pbbox_z/2);       // z size
  
  G4LogicalVolume* PbboxLog 
    = new G4LogicalVolume(Pbbox_box,
			  pb,
			  //pbMore,
                          "PbboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

    /*******************************
   * The Fe box       *
   *******************************/

  G4VSolid* Febox_box 
    = new G4Box("Febox_box",             // World Volume
                Febox_x/2,        // x size
                Febox_y/2,        // y size
                Febox_z/2);       // z size
  
  G4LogicalVolume* FeboxLog 
    = new G4LogicalVolume(Febox_box,
			  fe,
                          "FeboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

    /*******************************
   * The W box       *
   *******************************/

  G4VSolid* Wbox_box 
    = new G4Box("Wbox_box",             // World Volume
                Wbox_x/2,        // x size
                Wbox_y/2,        // y size
                Wbox_z/2);       // z size
  
  G4LogicalVolume* WboxLog 
    = new G4LogicalVolume(Wbox_box,
			  w,
                          "WboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits


    /*******************************
   * The Vaccum box       *
   *******************************/

  G4VSolid* Vacbox_box 
    = new G4Box("Vacbox_box",             // World Volume
                Vacbox_x/2,        // x size
                Vacbox_y/2,        // y size
                Vacbox_z/2);       // z size
  
  G4LogicalVolume* VacboxLog 
    = new G4LogicalVolume(Vacbox_box,
			  pb,
                          "VacboxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

    /*******************************
   * The Vaccum box PKU       *
   *******************************/

  G4VSolid* Vacbox_box_PKU 
    = new G4Box("Vacbox_box_PKU",             // World Volume
                pku_box_x/2,        // x size
                pku_box_y/2,        // y size
                pku_box_z/2);       // z size
  
  G4LogicalVolume* VacboxPLog 
    = new G4LogicalVolume(Vacbox_box_PKU,
			  vacum,
                          "VacboxPLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  VacboxPLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  G4LogicalVolume* VacboxKLog 
    = new G4LogicalVolume(Vacbox_box_PKU,
			  vacum,
                          "VacboxKLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  VacboxKLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  G4LogicalVolume* VacboxULog 
    = new G4LogicalVolume(Vacbox_box_PKU,
			  vacum,
                          "VacboxULog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  VacboxULog->SetVisAttributes(G4VisAttributes::GetInvisible());


    /*******************************
   * The PKU Bar1 & 2       *
   *******************************/

  G4VSolid* Bar1_box 
    = new G4Box("Bar1_box",             // World Volume
                pku_bar_x/2,        // x size
                pku_bar_y1/2,        // y size
                pku_bar_z/2);       // z size
  
  G4LogicalVolume* Bar1Log 
    = new G4LogicalVolume(Bar1_box,
			  pb,
                          "Bar1Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  G4VSolid* Bar2_box 
    = new G4Box("Bar2_box",             // World Volume
                pku_bar_x/2,        // x size
                pku_bar_y2/2,        // y size
                pku_bar_z/2);       // z size
  
  G4LogicalVolume* Bar2Log 
    = new G4LogicalVolume(Bar2_box,
			  pb,
                          "Bar2Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  // put Bar in PKUBox P
  //
G4RotationMatrix rotmP1;                    //rotation matrix to place modules
rotmP1.rotateZ(90*deg);
G4Transform3D transformP1(rotmP1, G4ThreeVector(0.,pku_bar_y2/2.-pku_bar_x/2.,0.));
  new G4PVPlacement(transformP1,
                    Bar1Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxPLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

G4RotationMatrix rotmP2;                    //rotation matrix to place modules
rotmP2.rotateZ(90*deg);
G4Transform3D transformP2(rotmP1, G4ThreeVector(0.,+pku_bar_x/2.,0.));
  new G4PVPlacement(transformP2,
                    Bar1Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxPLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(pku_bar_y1/2.+pku_bar_x/2.,pku_bar_y2/2.-pku_bar_y1/2.,0),     //at position
                    Bar1Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxPLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-pku_bar_y1/2.-pku_bar_x/2.,0,0),     //at position
                    Bar2Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxPLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

G4RotationMatrix rotmK1;                    //rotation matrix to place modules
rotmK1.rotateZ(-45*deg);
G4Transform3D transformK1(rotmK1, G4ThreeVector(0.,(pku_bar_y1/2.+pku_bar_x/2.)/std::sqrt(2),0.));
new G4PVPlacement(transformK1,                       //rotation+position
            Bar1Log,             //its logical volume
            "PKUbar",                //its name
            VacboxKLog,                //its mother  volume
            false,                   //no boolean operation
            0);               //copy number

G4RotationMatrix rotmK2;                    //rotation matrix to place modules
rotmK2.rotateZ(45*deg);
G4Transform3D transformK2(rotmK2, G4ThreeVector(0.,-(pku_bar_y1/2.+pku_bar_x/2.)/std::sqrt(2),0.));
new G4PVPlacement(transformK2,                       //rotation+position
            Bar1Log,             //its logical volume
            "PKUbar",                //its name
            VacboxKLog,                //its mother  volume
            false,                   //no boolean operation
            0);               //copy number

  G4double a=(pku_bar_y1/2.-pku_bar_x/2.)/std::sqrt(2) + pku_bar_x/2.*std::sqrt(2);
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-a-pku_bar_x/2.,0,0),     //at position
                    Bar2Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxKLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-pku_bar_y1/2.-pku_bar_x/2.,0,0),     //at position
                    Bar2Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxULog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(pku_bar_y1/2.+pku_bar_x/2.,0,0),     //at position
                    Bar2Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxULog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

G4RotationMatrix rotmU1;                    //rotation matrix to place modules
rotmU1.rotateZ(90*deg);
G4Transform3D transformU1(rotmP1, G4ThreeVector(0.,-pku_bar_y2/2.+pku_bar_x/2.,0.));
  new G4PVPlacement(transformU1,
                    Bar1Log,             //its logical volume
                    "PKUbar",                //its name
                    VacboxULog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

    /*******************************
   * The Box       *
   *******************************/

  G4VSolid* Box1_box 
    = new G4Box("Box1_box",             // World Volume
                (Box_x-200*um)/2,        // x size
                (Box_y-200*um)/2,        // y size
                (Box_z-200*um)/2);       // z size
  
  G4LogicalVolume* Box1Log 
    = new G4LogicalVolume(Box1_box,
                          air,
			  //world_Mat,
                          "Box1Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  Box1Log->SetVisAttributes(G4VisAttributes::GetInvisible());

  G4VSolid* Box_box 
    = new G4Box("Box_box",             // World Volume
                Box_x/2,        // x size
                Box_y/2,        // y size
                Box_z/2);       // z size
  
  G4LogicalVolume* BoxLog 
    = new G4LogicalVolume(Box_box,
                          Shell_Mat,
                          //air,
			  //world_Mat,
                          "BoxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  BoxLog->SetVisAttributes(G4VisAttributes::GetInvisible());

  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,0),     //at position
                    Box1Log,             //its logical volume
                    "Box1",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

/*
  // put Pbplate in Box
  //
  G4VPhysicalVolume* PbplatePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    PbplateLog,             //its logical volume
                    "Pbplate",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Vacbox in Pbbox
  //
  G4VPhysicalVolume* VacboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    VacboxLog,             //its logical volume
                    "Vacbox",                //its name
                    PbboxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Vacbox PKU in Pbbox
  //
  G4VPhysicalVolume* VacPboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-1./3.*m,0,0),                    //at position
                    VacboxPLog,             //its logical volume
                    "VacPbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  G4VPhysicalVolume* VacKboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,0),                    //at position
                    VacboxKLog,             //its logical volume
                    "VacKbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  G4VPhysicalVolume* VacUboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(1./3.*m,0,0),                    //at position
                    VacboxULog,             //its logical volume
                    "VacUbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Pbbox in Box
  //
  G4VPhysicalVolume* PbboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    //G4ThreeVector(),                    //at position
                    G4ThreeVector(-1./3.*m,1./3.*m,0),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Pbbox in Box
  //
  G4VPhysicalVolume* PbboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    //G4ThreeVector(),                    //at position
                    G4ThreeVector(0,0,0),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Pbbox in Box
  //
  G4VPhysicalVolume* PbboxPhys1 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,-200*mm),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  G4VPhysicalVolume* PbboxPhys2 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,200*mm),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    1,                       //copy number
                    0);          //overlaps checking
*/

/*
  // put Febox in Box
  //
  G4VPhysicalVolume* FeboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    //G4ThreeVector(),                    //at position
                    G4ThreeVector(0,-1./3.*m,1./3.*m),                    //at position
                    FeboxLog,             //its logical volume
                    "Febox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  // put Wbox in Box
  //
  G4VPhysicalVolume* WboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    //G4ThreeVector(),                    //at position
                    G4ThreeVector(1./3.*m,0,-1./3.*m),                    //at position
                    WboxLog,             //its logical volume
                    "Wbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

    /*******************************
   * The Experimental Hall       *
   *******************************/

  G4VSolid* experimentalHall_box 
    = new G4Box("expHall_box",             // World Volume
                experimentalHall_x/2,        // x size
                experimentalHall_y/2,        // y size
                experimentalHall_z/2);       // z size
  
  G4LogicalVolume* experimentalHallLog 
    = new G4LogicalVolume(experimentalHall_box,
                          //air,
			  world_Mat,
                          "expHallLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  G4VPhysicalVolume* experimentalHallPhys 
    = new G4PVPlacement(0,
                        G4ThreeVector(),   //at (0,0,0)
                        "expHall",
                        experimentalHallLog,
                        0,
                        false, 
                        0);
    
  experimentalHallLog->SetVisAttributes(G4VisAttributes::GetInvisible());

  // put Box in world
  //
  G4VPhysicalVolume* BoxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    BoxLog,             //its logical volume
                    "Box",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking

  // put Gem in world
  //
  Zcenter = -0.5*Box_z-2.5*Gem_z;

  for (int j=0; j<num_Gem; j++){ //in z axis
        if(j!=2) Zcenter += Gem_z;
        else if(j==2) Zcenter += Gem_z+Box_z;
        G4cout<<j<<" Gem with Zcenter = "<<Zcenter<<G4endl;
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    GemLog,             //its logical volume
                    "Gem",                //its name
                    experimentalHallLog,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }
  G4VisAttributes* GemLogAtt = new G4VisAttributes(G4Colour(233/256.,206/256.,238./256,0.1));
  GemLog->SetVisAttributes(GemLogAtt);

  fScoringVolume = GemLog;
  fScoringVolume2 = readoutplateLog;
  fScoringVolume3 = PbboxLog;

  // visualization attributes ------------------------------------------------

  auto visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  experimentalHallLog->SetVisAttributes(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.9));   // LightGray
  BoxLog->SetVisAttributes(visAttributes);


  visAttributes = new G4VisAttributes(G4Colour(0,1,0));
  PbboxLog->SetVisAttributes(visAttributes);

  visAttributes = new G4VisAttributes(G4Colour(0,0,1));
  VacboxLog->SetVisAttributes(visAttributes);

  return experimentalHallPhys;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4cout << "Construt the DetectorGeometry" <<G4endl;
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // define a material
  DefineMaterials();
    
  // define some constant
  DefineConstants();

  // Define volumes
  return DefineVolumes();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
