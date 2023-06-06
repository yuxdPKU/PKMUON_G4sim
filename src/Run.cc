//
// 2020.5.8 by siguang wang (siguang@pku.edu.cn  PKU)
//


#include <TTree.h>
#include "Run.hh"
#include <TTree.h>
#include "TVector3.h"
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TMath.h>
#include <TRandom3.h>
#include "TClonesArray.h"
#include <iostream>
#include "RunMessenger.hh"
#include "G4SystemOfUnits.hh"

static Run*  instance = 0;

Run::Run()
  :  fRunMessenger(0)
{
   totEdep=0;
  _GunEng = 0;
  rootFileName = "test.root";
  fRunMessenger = new RunMessenger(this);
}

Run::~Run() {
 
  delete fRunMessenger;
}

Run*  Run::GetInstance()
{
  if (instance == 0 )
    {
      instance = new Run();
    }
  return instance;
}


void Run::initTree() {
  std::cout << "Init Tree File." << std::endl;
  _file = new TFile(rootFileName.c_str(), "RECREATE");
  _tree = new TTree("T1", "Simple Out Tree");

  _tree->Branch("GunEng", &_GunEng, "GunEng/D");
  _tree->Branch("totEdep", &totEdep, "totEdep/D");
  _tree->Branch("cosTh", &cosTh, "cosTh/D");
  
  _tree->Branch("X0", &X0, "X0/D"); //source Position
  _tree->Branch("Y0", &Y0, "Y0/D"); //Source Position
  _tree->Branch("Z0", &Z0, "Z0/D"); //Source Position
  _tree->Branch("px0", &px0, "px0/D"); //initial px
  _tree->Branch("py0", &py0, "py0/D"); //initial py
  _tree->Branch("pz0", &pz0, "pz0/D"); //initial pz

  _tree->Branch("vX",&vX); //Edep X
  _tree->Branch("vY",&vY); //Edep Y
  _tree->Branch("vZ",&vZ); //Edep Z

  _tree->Branch("vEdep",&vEdep); //Edep Energy
  _tree->Branch("vTrkID",&vTrkID); //TrackID
 
/*
  _tree->Branch("X1", &X1, "X1/D"); //final Position
  _tree->Branch("Y1", &Y1, "Y1/D"); //final Position
  _tree->Branch("Z1", &Z1, "Z1/D"); //final Position
  _tree->Branch("px1", &px1, "px1/D"); //final px
  _tree->Branch("py1", &py1, "py1/D"); //final py
  _tree->Branch("pz1", &pz1, "pz1/D"); //final pz
*/

  _tree->Branch("vReadoutTrkid",&vReadoutTrkid); //Trk id in readout bar
  _tree->Branch("vReadoutEdep",&vReadoutEdep); //Edep Energy in readout bar
  _tree->Branch("vReadoutPosX",&vReadoutPosX); //Edep X in readout bar
  _tree->Branch("vReadoutPosY",&vReadoutPosY); //Edep Y in readout bar
  _tree->Branch("vReadoutPosZ",&vReadoutPosZ); //Edep Z in readout bar

  _tree->Branch("vPx",&vPx); //Pvx
  _tree->Branch("vPy",&vPy); //Pvy
  _tree->Branch("vPz",&vPz); //Pvz

  totEdep = 0;
  _GunEng =0;
  cosTh=-2;
  status=false;
  X0=0;
  Y0=0;
  Z0=0;
  px0=0;
  py0=0;
  pz0=0;
  
/*
  X1=0;
  Y1=0;
  Z1=0;
  px1=0;
  py1=0;
  pz1=0;
*/

  vX.clear();
  vY.clear();
  vZ.clear(); 

  vEdep.clear();
  vTrkID.clear();
 
  vReadoutTrkid.clear();
  vReadoutEdep.clear();
  vReadoutPosX.clear();
  vReadoutPosY.clear();
  vReadoutPosZ.clear();
  vPx.clear();
  vPy.clear();
  vPz.clear();

  std::vector<Double_t>().swap(vEdep);
  std::vector<Double_t>().swap(vX);
  std::vector<Double_t>().swap(vY);
  std::vector<Double_t>().swap(vZ);
  std::vector<Int_t>().swap(vTrkID);
  std::vector<Int_t>().swap(vReadoutTrkid);
  std::vector<Double_t>().swap(vReadoutEdep);
  std::vector<Double_t>().swap(vReadoutPosX);
  std::vector<Double_t>().swap(vReadoutPosY);
  std::vector<Double_t>().swap(vReadoutPosZ);
  std::vector<Double_t>().swap(vPx);
  std::vector<Double_t>().swap(vPy);
  std::vector<Double_t>().swap(vPz);
  
  G4cout << "ROOT Name = " << rootFileName << G4endl;

}

void Run::AddEnergy1(double Eng1) {
  if(Eng1>0)  totEdep += Eng1 ;
}



void Run::SetPos(G4double inX,G4double inY,G4double inZ){
  X0 = inX;
  Y0 = inY;
  Z0 = inZ;
}

void Run::SetPxyz(G4double inX,G4double inY,G4double inZ){
  px0 = inX;
  py0 = inY;
  pz0 = inZ;
}

void Run::SetPosF(G4double fiX,G4double fiY,G4double fiZ){
  X1 = fiX;
  Y1 = fiY;
  Z1 = fiZ;
}

void Run::SetPxyzF(G4double fiX,G4double fiY,G4double fiZ){
  px1 = fiX;
  py1 = fiY;
  pz1 = fiZ;
}

void Run::ClearAll(){
  totEdep = 0;
  _GunEng =0;
  cosTh=-2;
  status=false;

  X0=0;
  Y0=0;
  Z0=0;
  
  px0=0;
  py0=0;
  pz0=0;
  
  X1=0;
  Y1=0;
  Z1=0;
  
  px1=0;
  py1=0;
  pz1=0;

  vX.clear();
  vY.clear();
  vZ.clear(); 
  vEdep.clear();
  vTrkID.clear();

  vReadoutTrkid.clear();
  vReadoutEdep.clear();
  vReadoutPosX.clear();
  vReadoutPosY.clear();
  vReadoutPosZ.clear();
  vPx.clear();
  vPy.clear();
  vPz.clear();

  std::vector<Double_t>().swap(vEdep);
  std::vector<Double_t>().swap(vX);
  std::vector<Double_t>().swap(vY);
  std::vector<Double_t>().swap(vZ);
  std::vector<Int_t>().swap(vTrkID);
  std::vector<Int_t>().swap(vReadoutTrkid);
  std::vector<Double_t>().swap(vReadoutEdep);
  std::vector<Double_t>().swap(vReadoutPosX);
  std::vector<Double_t>().swap(vReadoutPosY);
  std::vector<Double_t>().swap(vReadoutPosZ);
  std::vector<Double_t>().swap(vPx);
  std::vector<Double_t>().swap(vPy);
  std::vector<Double_t>().swap(vPz);
 
}
void  Run::AddX(G4double inX){
  vX.push_back(inX);
}
void  Run::AddY(G4double inY){
  vY.push_back(inY);
}
void  Run::AddZ(G4double inZ){
  vZ.push_back(inZ);
}
void  Run::AddEdep(G4double inE){
  vEdep.push_back(inE);
}

void  Run::AddTrkID(G4int iTrkID){
  vTrkID.push_back(iTrkID);
}

void  Run::AddReadoutEdepX(G4double X){
  vReadoutPosX.push_back(X);
}

void  Run::AddReadoutEdepY(G4double Y){
  vReadoutPosY.push_back(Y);
}

void  Run::AddReadoutEdepZ(G4double Z){
  vReadoutPosZ.push_back(Z);
}

void  Run::AddReadoutEdep(G4double E){
  vReadoutEdep.push_back(E);
}

void  Run::AddPx(G4double X){
  vPx.push_back(X);
}

void  Run::AddPy(G4double Y){
  vPy.push_back(Y);
}

void  Run::AddPz(G4double Z){
  vPz.push_back(Z);
}

void  Run::AddReadoutTrkid(G4int i){
  vReadoutTrkid.push_back(i);
}

void Run::Fill(){
  _tree->Fill();
/*
  G4cout<<"this is Run::Fill function: size of vReadoutPosX is "<<vReadoutPosX.size()<<G4endl;
  for(int i=0; i<vReadoutPosX.size(); i++){
  G4cout<<"trk "<<i<<G4endl;
  G4cout<<"x = "<<vReadoutPosX[i]<<G4endl;
  G4cout<<"y = "<<vReadoutPosY[i]<<G4endl;
  G4cout<<"z = "<<vReadoutPosZ[i]<<G4endl;
  G4cout<<"e = "<<vReadoutEdep[i]<<G4endl;
  }
*/
  ClearAll();
}

void Run::saveTree(){
  _file->cd();
  _tree->Write("T1",TObject::kOverwrite);
  _file->Close();
}