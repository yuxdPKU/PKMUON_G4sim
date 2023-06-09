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
// Previous authors: G. Guerrieri, S. Guatelli and M. G. Pia, INFN Genova, Italy
// Authors (since 2007): S. Guatelli, University of Wollongong, Australia
// 
#include "EventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"
#include "RunAction.hh"
#include "G4RunManager.hh"
#include "Run.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction()

{ 
  
}
 
EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{

 Run::GetInstance()->ClearAll();
 G4PrimaryVertex *pv = evt->GetPrimaryVertex(0);
 G4double kinE =  pv-> GetPrimary()->GetKineticEnergy();
 Run::GetInstance()->RecPrimPartEng(kinE/MeV);
 Run::GetInstance()->SetPos(pv->GetX0()/mm,pv->GetY0()/mm,pv->GetZ0()/mm);
 G4ThreeVector mom =  pv->GetPrimary()->GetMomentum();
 Run::GetInstance()->SetPxyz(mom.x()/MeV,mom.y()/MeV,mom.z()/MeV);
 
}
 
void EventAction::EndOfEventAction(const G4Event* evt)
{  //if(Run::GetInstance()->GetStatus()==true){
    //Run::GetInstance()->Fill();
  //}
    Run::GetInstance()->Fill();
}

