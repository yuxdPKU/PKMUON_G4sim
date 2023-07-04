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
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4cout << "Construt the DetectorGeometry" <<G4endl;
  
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  //
  // define a material
  //   
    
  G4cout << "Defining the materials" << G4endl;
  // Get nist material manager
  G4NistManager* nistManager = G4NistManager::Instance();
  // Build materials
  G4Material* vacum   = nistManager->FindOrBuildMaterial("G4_Galactic");
  G4Material* air   = nistManager->FindOrBuildMaterial("G4_AIR");
  G4Material* pb   = nistManager->FindOrBuildMaterial("G4_Pb");
  G4Material* cu   = nistManager->FindOrBuildMaterial("G4_Cu");
  G4Material* Water = nistManager->FindOrBuildMaterial("G4_WATER");
  G4Material* fillM = nistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
  G4Material *detM = nistManager->FindOrBuildMaterial("G4_POLYSTYRENE");
  G4int nComponent;
  G4double fracmass=1.0;
  G4Material* pbLess   = new G4Material("pbLess",1E-3*g/cm3,nComponent=1);
  pbLess -> AddMaterial(pb,fracmass);
  G4double temperature = 325000000.*kelvin;
  G4double pressure = 50.*atmosphere;
  G4Material* pbLessHT   = new G4Material("pbLessHT",1E-3*g/cm3,nComponent=1,kStateGas,temperature,pressure);
  pbLessHT -> AddMaterial(pb,fracmass);

G4cout<<"define Kapton"<<G4endl;
  // Define Kapton material
  G4double density = 1.42 * g/cm3; // Density of Kapton
  G4Material* kapton = new G4Material("Kapton", density, 3);
  kapton->AddElement(nistManager->FindOrBuildElement("H"), 0.0273);
  kapton->AddElement(nistManager->FindOrBuildElement("C"), 0.7213);
  kapton->AddElement(nistManager->FindOrBuildElement("O"), 0.2514);

G4cout<<"define mixgass"<<G4endl;
  // Define 70%Ar and 30%CO2 gas mixture
  density = 1.822 * mg/cm3; // Density of the gas mixture
  G4Material* gasMixture = new G4Material("ArCO2GasMixture", density, 3);
  gasMixture->AddElement(nistManager->FindOrBuildElement("Ar"), 0.7);
  gasMixture->AddElement(nistManager->FindOrBuildElement("C"), 0.3/3);
  gasMixture->AddElement(nistManager->FindOrBuildElement("O"), 0.3/3*2);

  G4Material* Drift_cathode_Mat = cu;
  G4Material* Gem_inner_Mat = kapton;
  G4Material* Gem_outter_Mat = cu;
  G4Material* Readout_platte_Mat = cu;
  G4Material* Readout_bar_Mat = cu;
  G4Material* Gem_Mat = gasMixture;
  G4Material* world_Mat = vacum;

  // define some constant
  G4double Xcenter, Ycenter, Zcenter;

  const G4int num_Gem_hole = 8;
  G4double Gem_hole_outter_Diameter=60*um;
  G4double Gem_hole_outter_Length=5*um;
  G4double Gem_Hexagonal_Length=140*um*10*100;//1000x

  G4double Gem_outter_bar1_x=Gem_Hexagonal_Length*num_Gem_hole;
  G4double Gem_outter_bar1_y=0.5*Gem_Hexagonal_Length*std::sqrt(3);
  G4double Gem_outter_bar1_z=5*um;

  G4double Gem_outter_bar2_x=Gem_outter_bar1_x;
  G4double Gem_outter_bar2_y=Gem_outter_bar1_y;
  G4double Gem_outter_bar2_z=Gem_outter_bar1_z;
  G4double Gem_outter_x=Gem_outter_bar1_x;
  G4double Gem_outter_y=Gem_outter_bar1_y*2*num_Gem_hole;//bar1 and bar2
  G4double Gem_outter_z=Gem_outter_bar1_z;

   G4double Gem_hole_inner_rmin1=0*um, Gem_hole_inner_rmin2=0*um;
   G4double Gem_hole_inner_rmax1=30*um, Gem_hole_inner_rmax2=25*um;
   G4double Gem_hole_inner_hz=25*um;
   G4double Gem_hole_inner_phimin=0.*deg, Gem_hole_inner_phimax=360.*deg;

  G4double Gem_inner_bar1_x=Gem_outter_bar1_x;
  G4double Gem_inner_bar1_y=Gem_outter_bar1_y;
  G4double Gem_inner_bar1_z=Gem_hole_inner_hz*2;
  G4double Gem_inner_bar2_x=Gem_inner_bar1_x;
  G4double Gem_inner_bar2_y=Gem_inner_bar1_y;
  G4double Gem_inner_bar2_z=Gem_inner_bar1_z;

  G4double Gem_inner_x=Gem_inner_bar1_x;
  G4double Gem_inner_y=Gem_outter_bar1_y*2*num_Gem_hole;//bar1 and bar2;
  G4double Gem_inner_z=Gem_inner_bar1_z;

  G4double drift_cathode_x=Gem_outter_bar1_x;
  G4double drift_cathode_y=Gem_outter_y;
  G4double drift_cathode_z=0.1*mm;

  G4double readout_bar_x=150*um;
  G4double readout_bar_gap_x=210*um;
  G4double readout_bar_y=Gem_outter_y;
  G4double readout_bar_z=0.1*mm;
  G4double readout_platte_x=Gem_outter_x;
  G4double readout_platte_y=Gem_outter_y;
  G4double readout_platte_z=0.1*mm;

  G4double gap1 = 4.8*mm;
  G4double gap2 = 2*mm;
  G4double Gem_x=Gem_outter_x;
  G4double Gem_y=Gem_outter_y;
  const G4int num_Gem_outter=4*2;
  const G4int num_Gem_inner=4;
  G4double Gem_z=gap1+num_Gem_inner*gap2+(num_Gem_outter*Gem_outter_z+num_Gem_inner*Gem_inner_z)+drift_cathode_z+readout_bar_z+readout_platte_z;

  G4double Box_x=1*m;
  G4double Box_y=1*m;
  G4double Box_z=1*m;

  G4int num_Gem=4;

  G4double experimentalHall_x=1.1*fmax(Gem_outter_x,Box_x);
  G4double experimentalHall_y=1.1*fmax(Gem_outter_y,Box_y);
  G4double experimentalHall_z=1.1*(Box_z+Gem_z*4);
 
  G4double Pbplatte_x=100*mm;
  G4double Pbplatte_y=100*mm;
  G4double Pbplatte_z=0.1*mm;

/*
  G4double Pbbox_x=0.5*m;
  G4double Pbbox_y=0.5*m;
  //G4double Pbbox_z=0.5*m;
  G4double Pbbox_z=20*mm;
  //G4double Pbbox_z=2*mm;
*/

  G4double Pbbox_x=20*mm;
  G4double Pbbox_z=0.5*m;
  G4double Pbbox_y=0.5*m;

  G4double Vacbox_x=0.5*m-40*mm;
  G4double Vacbox_y=0.5*m-40*mm;
  G4double Vacbox_z=0.5*m-40*mm;

    /*******************************
   * Gem hole outter       *
   *******************************/
   G4VSolid* Gem_hole_outter_tub
            = new G4Tubs("Gem_hole_outter_tub",                      //name
                         0*mm, 0.5*Gem_hole_outter_Diameter,       //r1, r2
                         0.5*Gem_hole_outter_Length,               //half-length
                         0., twopi);                    //theta1, theta2

   G4LogicalVolume* Gem_hole_outter_Log
                        = new G4LogicalVolume(Gem_hole_outter_tub,          //solid
                                   Gem_Mat,            //material
                                   "Gem_hole_outter_Log");            //name


    /*******************************
   * Gem outter bar1      *
   *******************************/
  G4VSolid* Gem_outter_bar1_box
    = new G4Box("Gem_outter_bar1_box",             // World Volume
                Gem_outter_bar1_x/2,        // x size
                Gem_outter_bar1_y/2,        // y size
                Gem_outter_bar1_z/2);       // z size

  G4LogicalVolume* Gem_outter_bar1_Log
    = new G4LogicalVolume(Gem_outter_bar1_box,
			  Gem_outter_Mat,
                          "Gem_outter_bar1_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

   //G4VisAttributes* Gem_hole_outter_LogAtt = new G4VisAttributes(G4Colour(0/256.,0/256.,100/256,0.8));

  // put holes within Gem outter bar1
  //
  Xcenter = -0.5*(Gem_outter_bar1_x+Gem_Hexagonal_Length);

  for (int j=0; j<num_Gem_hole; j++){ //in x axis
        Xcenter = Xcenter + Gem_Hexagonal_Length;
        G4cout<<j<<"hole with Xcenter = "<<Xcenter<<G4endl;
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Xcenter,0,0),                    //at position
                    Gem_hole_outter_Log,             //its logical volume
                    "Gem_hole_outter",                //its name
                    Gem_outter_bar1_Log,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }

    /*******************************
   * Gem outter bar2      *
   *******************************/

  G4VSolid* Gem_outter_bar2_box
    = new G4Box("Gem_outter_bar2_box",             // World Volume
                Gem_outter_bar2_x/2,        // x size
                Gem_outter_bar2_y/2,        // y size
                Gem_outter_bar2_z/2);       // z size

  G4LogicalVolume* Gem_outter_bar2_Log
    = new G4LogicalVolume(Gem_outter_bar2_box,
			  Gem_outter_Mat,
                          "Gem_outter_bar2_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  // put holes within Gem outter bar2
  //
  Xcenter = -0.5*(Gem_outter_bar2_x);

  for (int j=0; j<num_Gem_hole-1; j++){ //in x axis
        Xcenter = Xcenter + Gem_Hexagonal_Length;
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Xcenter,0,0),                    //at position
                    Gem_hole_outter_Log,             //its logical volume
                    "Gem_hole_outter",                //its name
                    Gem_outter_bar2_Log,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }

    /*******************************
   * Gem outter platte       *
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

  // put bar within Gem outter
  //
  Ycenter = -0.5*(Gem_outter_y+Gem_outter_bar1_y)-Gem_outter_bar1_y;

  G4int num_Gem_outter_bar1 = 3;
  for (int j=0; j<num_Gem_outter_bar1; j++){
        Ycenter = Ycenter + 2*Gem_outter_bar1_y;
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,Ycenter,0),                    //at position
                    Gem_outter_bar1_Log,             //its logical volume
                    "Gem_outter",                //its name
                    Gem_outter_Log,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }

  Ycenter = -0.5*(Gem_outter_y-Gem_outter_bar2_y)-Gem_outter_bar2_y;

  G4int num_Gem_outter_bar2 = 3;
  for (int j=0; j<num_Gem_outter_bar2; j++){
        Ycenter = Ycenter + 2*Gem_outter_bar2_y;
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,Ycenter,0),                    //at position
                    Gem_outter_bar2_Log,             //its logical volume
                    "Gem_outter",                //its name
                    Gem_outter_Log,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }



    /*******************************
   * Gem hole inner       *
   *******************************/
   G4VSolid* Gem_hole_inner_cons
            = new G4Cons("Gem_hole_inner_cons",                      //name
                Gem_hole_inner_rmin1,Gem_hole_inner_rmax1,
                Gem_hole_inner_rmin2,Gem_hole_inner_rmax2,
                Gem_hole_inner_hz/2,
                Gem_hole_inner_phimin,Gem_hole_inner_phimax);

   G4LogicalVolume* Gem_hole_inner_Log
                        = new G4LogicalVolume(Gem_hole_inner_cons,          //solid
                                   Gem_Mat,            //material
                                   "Gem_hole_inner_Log");            //name

   G4VSolid* Gem_hole_inner_box
            = new G4Box("Gem_hole_inner_box",                      //name
                Gem_hole_inner_rmax1,        // x size
                Gem_hole_inner_rmax1,        // y size
                Gem_hole_inner_hz);       // z size

   G4LogicalVolume* Gem_hole_inner_doubleLog
                        = new G4LogicalVolume(Gem_hole_inner_box,          //solid
                                   Gem_inner_Mat,            //material
                                   "Gem_hole_inner_doubleLog");            //name

  // put cons within holes
  //
  Zcenter = -0.5*Gem_hole_inner_hz;
  new G4PVPlacement(0,                       //no rotation
              G4ThreeVector(0,0,Zcenter),                    //at position
              Gem_hole_inner_Log,             //its logical volume
              "Gem_hole_inner_Log",                //its name
              Gem_hole_inner_doubleLog,                //its mother  volume
              false,                   //no boolean operation
              0);               //copy number

  Zcenter = 0.5*Gem_hole_inner_hz;
  G4RotationMatrix rotm;                    //rotation matrix to place modules
  rotm.rotateX(180*deg);
      G4Transform3D transform(rotm, G4ThreeVector(0.,0.,Zcenter));
  new G4PVPlacement(transform,                       //rotation+position
              Gem_hole_inner_Log,             //its logical volume
              "Gem_hole_inner_Log",                //its name
              Gem_hole_inner_doubleLog,                //its mother  volume
              false,                   //no boolean operation
              1);               //copy number

    /*******************************
   * Gem inner bar1      *
   *******************************/
  G4VSolid* Gem_inner_bar1_box
    = new G4Box("Gem_inner_bar1_box",             // World Volume
                Gem_inner_bar1_x/2,        // x size
                Gem_inner_bar1_y/2,        // y size
                Gem_inner_bar1_z/2);       // z size

  G4LogicalVolume* Gem_inner_bar1_Log
    = new G4LogicalVolume(Gem_inner_bar1_box,
			  Gem_inner_Mat,
                          "Gem_inner_bar1_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  // put holes within Gem inner bar1
  //
  Xcenter = -0.5*(Gem_inner_bar1_x+Gem_Hexagonal_Length);

  for (int j=0; j<num_Gem_hole; j++){ //in x axis
        Xcenter = Xcenter + Gem_Hexagonal_Length;
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Xcenter,0,0),                    //at position
                    Gem_hole_inner_doubleLog,             //its logical volume
                    "Gem_hole_inner_doubleLog",                //its name
                    Gem_inner_bar1_Log,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }

    /*******************************
   * Gem inner bar2      *
   *******************************/
  G4VSolid* Gem_inner_bar2_box
    = new G4Box("Gem_inner_bar2_box",             // World Volume
                Gem_inner_bar2_x/2,        // x size
                Gem_inner_bar2_y/2,        // y size
                Gem_inner_bar2_z/2);       // z size

  G4LogicalVolume* Gem_inner_bar2_Log
    = new G4LogicalVolume(Gem_inner_bar2_box,
			  Gem_inner_Mat,
                          "Gem_inner_bar2_Log",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  // put holes within Gem inner bar2
  //
  Xcenter = -0.5*(Gem_inner_bar2_x);

  for (int j=0; j<num_Gem_hole-1; j++){ //in x axis
        Xcenter = Xcenter + Gem_Hexagonal_Length;
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Xcenter,0,0),                    //at position
                    Gem_hole_inner_doubleLog,             //its logical volume
                    "Gem_hole_inner_doubleLog",                //its name
                    Gem_inner_bar2_Log,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }

    /*******************************
   * Gem inner platte       *
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

  // put bar within Gem inner
  //
  Ycenter = -0.5*(Gem_inner_y+Gem_inner_bar1_y)-Gem_inner_bar1_y;

  G4int num_Gem_inner_bar1 = 3;
  for (int j=0; j<num_Gem_inner_bar1; j++){
        Ycenter = Ycenter + 2*Gem_inner_bar1_y;
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,Ycenter,0),                    //at position
                    Gem_inner_bar1_Log,             //its logical volume
                    "Gem_inner",                //its name
                    Gem_inner_Log,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }

  Ycenter = -0.5*(Gem_inner_y-Gem_inner_bar2_y)-Gem_inner_bar2_y;

  G4int num_Gem_inner_bar2 = 3;
  for (int j=0; j<num_Gem_inner_bar2; j++){
        Ycenter = Ycenter + 2*Gem_inner_bar2_y;
        new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,Ycenter,0),                    //at position
                    Gem_inner_bar2_Log,             //its logical volume
                    "Gem_inner",                //its name
                    Gem_inner_Log,                //its mother  volume
                    false,                   //no boolean operation
                    j+1);               //copy number
  }

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
  
    /*******************************
   * Readout platte       *
   *******************************/
  G4VSolid* readout_platte_box 
    = new G4Box("readoutplatte_box",             // World Volume
                readout_platte_x/2,        // x size
                readout_platte_y/2,        // y size
                readout_platte_z/2);       // z size
  
  G4LogicalVolume* readoutplatteLog 
    = new G4LogicalVolume(readout_platte_box,
			  Readout_platte_Mat,
                          "readoutplatteLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

  G4VisAttributes* readoutplatteLogAtt = new G4VisAttributes(G4Colour(233/256.,0/256.,0/256,0.7));
  //readoutplatteLog->SetVisAttributes(readoutplatteLogAtt);
  readoutplatteLog->SetVisAttributes(G4VisAttributes::GetInvisible());


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

    // put Gem outter platte in Gem
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

    // put Gem inner platte in Gem
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
  G4LogicalVolume* readoutbarLog[num_readout_bar];
  Zcenter = 0.5*Gem_z-readout_platte_z-0.5*readout_bar_z;
  Xcenter = -0.5*(Gem_x-readout_bar_x)-readout_bar_gap_x;
  for (int i=0; i<num_readout_bar; i++){

  G4String name = "readoutbarLog" + std::to_string(i);

  readoutbarLog[i]
    = new G4LogicalVolume(readout_bar_box,
			  Readout_bar_Mat,
                          name,
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits

    Xcenter = Xcenter + readout_bar_gap_x;
    new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(Xcenter,0,Zcenter),                    //at position
                    readoutbarLog[i],             //its logical volume
                    name,                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    i+1,                       //copy number
                    0);          //overlaps checking
  }

  // put readout platte in Gem
  //
  Zcenter = 0.5*(Gem_z-readout_platte_z);
  G4VPhysicalVolume* readoutplattePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(0,0,Zcenter),                    //at position
                    readoutplatteLog,             //its logical volume
                    "readoutplatte",                //its name
                    GemLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    1);          //overlaps checking

    /*******************************
   * The Pb platte       *
   *******************************/

  G4VSolid* Pbplatte_box 
    = new G4Box("Pbplatte_box",             // World Volume
                Pbplatte_x/2,        // x size
                Pbplatte_y/2,        // y size
                Pbplatte_z/2);       // z size
  
  G4LogicalVolume* PbplatteLog 
    = new G4LogicalVolume(Pbplatte_box,
			  pb,
                          "PbplatteLog",
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
                          "PbboxLog",
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
   * The Box       *
   *******************************/

  G4VSolid* Box_box 
    = new G4Box("Box_box",             // World Volume
                Box_x/2,        // x size
                Box_y/2,        // y size
                Box_z/2);       // z size
  
  G4LogicalVolume* BoxLog 
    = new G4LogicalVolume(Box_box,
                          //air,
			  world_Mat,
                          "BoxLog",
                          0,               //opt: fieldManager
                          0,               //opt: SensitiveDetector
                          0);              //opt: UserLimits
  BoxLog->SetVisAttributes(G4VisAttributes::GetInvisible());

/*
  // put Pbplatte in Box
  //
  G4VPhysicalVolume* PbplattePhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    PbplatteLog,             //its logical volume
                    "Pbplatte",                //its name
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
  // put Pbbox in Box
  //
  G4VPhysicalVolume* PbboxPhys 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking
*/

  // put Pbbox in Box
  //
  G4VPhysicalVolume* PbboxPhys1 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(-200*mm,0,0),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    0);          //overlaps checking

  G4VPhysicalVolume* PbboxPhys2 
    = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(200*mm,0,0),                    //at position
                    PbboxLog,             //its logical volume
                    "Pbbox",                //its name
                    BoxLog,                //its mother  volume
                    false,                   //no boolean operation
                    1,                       //copy number
                    0);          //overlaps checking



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
  fScoringVolume2 = readoutplatteLog;
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
