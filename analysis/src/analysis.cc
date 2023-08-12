#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sys/time.h>

#include "../include/PoCA.h"

using namespace std;
extern double Z1;
extern double Z2;
extern double Z3;
extern double Z4;

void analysis(){
TRandom *rand = new TRandom();
double sigma=100*0.001;

//TFile *f = new TFile("../../build/root_file/box.root","");
//TFile *f = new TFile("../../build/root_file/pbfew.root","");
TFile *f = new TFile("../../build/root_file/pku.root","");
TTree *t = (TTree*)f->Get("T1");

std::vector<double> *vReadoutPosX=0; //Edep X in readout bar
std::vector<double> *vReadoutPosY=0; //Edep Y in readout bar
std::vector<double> *vReadoutPosZ=0; //Edep Z in readout bar
std::vector<double> *vReadoutEdep=0; //Edep Energy in readout bar
std::vector<double> *vReadoutE=0; //Total Energy in readout bar
std::vector<int> *vReadoutTrkid=0; //Trk id in readout bar
std::vector<int> *vReadoutTrkparentid=0; //Trk parent id in readout bar
std::vector<double> *vPx=0;
std::vector<double> *vPy=0;
std::vector<double> *vPz=0;

std::vector<double> *vPbPosX=0;
std::vector<double> *vPbPosY=0;
std::vector<double> *vPbPosZ=0;
std::vector<int> *vPbTrkid=0;

t->SetBranchAddress("vReadoutPosX",&vReadoutPosX);
t->SetBranchAddress("vReadoutPosY",&vReadoutPosY);
t->SetBranchAddress("vReadoutPosZ",&vReadoutPosZ);
t->SetBranchAddress("vReadoutEdep",&vReadoutEdep);
t->SetBranchAddress("vReadoutE",&vReadoutE);
t->SetBranchAddress("vReadoutTrkid",&vReadoutTrkid);
t->SetBranchAddress("vReadoutTrkparentid",&vReadoutTrkparentid);
t->SetBranchAddress("vPx",&vPx);
t->SetBranchAddress("vPy",&vPy);
t->SetBranchAddress("vPz",&vPz);
t->SetBranchAddress("vPbPosX",&vPbPosX);
t->SetBranchAddress("vPbPosY",&vPbPosY);
t->SetBranchAddress("vPbPosZ",&vPbPosZ);
t->SetBranchAddress("vPbTrkid",&vPbTrkid);
//t->Print();

//new file and new tree
//TFile * fn = new TFile("../root/mu1GeV_box.root","recreate");
//TFile * fn = new TFile("../root/mu1GeV_pbfew.root","recreate");
TFile * fn = new TFile("../root/mu1GeV_pku.root","recreate");
TTree * tn = new TTree("T1","tree");
tn = t->CloneTree(0);

double muonE[4];
double rec_x[4], rec_y[4], rec_z[4];
double rec_x_smear[4], rec_y_smear[4];
double angle, angle_smear;
double poca_x, poca_y, poca_z, dca;
int poca_status;
double poca_x_smear, poca_y_smear, poca_z_smear, dca_smear;
int poca_status_smear;
double angle_gem12;
tn->Branch("muonE", muonE, "muonE[4]/D");
tn->Branch("rec_x", rec_x, "rec_x[4]/D");
tn->Branch("rec_y", rec_y, "rec_y[4]/D");
tn->Branch("rec_z", rec_z, "rec_z[4]/D");
tn->Branch("rec_x_smear", rec_x_smear, "rec_x_smear[4]/D");
tn->Branch("rec_y_smear", rec_y_smear, "rec_y_smear[4]/D");
tn->Branch("angle", &angle, "angle/D");
tn->Branch("angle_smear", &angle_smear, "angle_smear/D");
tn->Branch("poca_x", &poca_x, "poca_x/D");
tn->Branch("poca_y", &poca_y, "poca_y/D");
tn->Branch("poca_z", &poca_z, "poca_z/D");
tn->Branch("dca", &dca, "dca/D");
tn->Branch("poca_status", &poca_status, "poca_status/I");
tn->Branch("poca_x_smear", &poca_x_smear, "poca_x_smear/D");
tn->Branch("poca_y_smear", &poca_y_smear, "poca_y_smear/D");
tn->Branch("poca_z_smear", &poca_z_smear, "poca_z_smear/D");
tn->Branch("dca_smear", &dca_smear, "dca_smear/D");
tn->Branch("poca_status_smear", &poca_status_smear, "poca_status_smear/I");
tn->Branch("angle_gem12", &angle_gem12, "angle_gem12/D");

double gem_z_sim[4] = {Z1,Z2,Z3,Z4};

struct timeval start;
struct timeval end;
unsigned long timer;

gettimeofday(&start, NULL); // 计时开始

int nevent = t->GetEntries();
for(int ievent=0; ievent<nevent; ievent++){
        if (ievent % (int)(nevent / 10) == 0) cout << "Processing progress: " << ievent / (int)(nevent / 10) << "0%" << endl;
        t->GetEntry(ievent);

        // initialization
        for(int igem=0; igem<4; igem++){
                rec_x[igem]=0;
                rec_y[igem]=0;
                rec_z[igem]=0;
        }
        double Etot[4]={0,0,0,0};

        TVector3 *pmom[4];
        double nhit = vReadoutTrkid->size();
        //cout<<"nhit = "<<nhit<<endl;
        for(int itrk=0; itrk<nhit; itrk++){
                // muon track id
                // trk id = 1 <==> parent id = 0
                //if((*vReadoutTrkid)[itrk]==1) continue;
                if((*vReadoutTrkid)[itrk]!=1) continue;
                // cout<<"vReadoutTrkid = "<<(*vReadoutTrkid)[itrk]<<endl;
                //cout<<"trkid = "<<(*vReadoutTrkid)[itrk]<<" , parentid = "<<(*vReadoutTrkparentid)[itrk]<<endl;
                //if((*vReadoutTrkparentid)[itrk]==0) if((*vReadoutTrkid)[itrk]!=1) cout<<"what!!!"<<endl;

                // gem layer number
                int igem=0;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[0])<1 ) igem=0;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[1])<1 ) igem=1;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[2])<1 ) igem=2;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[3])<1 ) igem=3;
                // cout<<"igem = "<<igem<<endl;

                // energy weight position
                muonE[igem]=(*vReadoutE)[itrk];
                Etot[igem]+=(*vReadoutEdep)[itrk];
                rec_x[igem]+=(*vReadoutEdep)[itrk]*(*vReadoutPosX)[itrk];
                rec_y[igem]+=(*vReadoutEdep)[itrk]*(*vReadoutPosY)[itrk];
                //rec_z[igem]+=(*vReadoutEdep)[itrk]*(*vReadoutPosZ)[itrk];
                rec_z[igem]+=(*vReadoutEdep)[itrk]*gem_z_sim[igem];
                pmom[igem] = new TVector3((*vPx)[itrk],(*vPy)[itrk],(*vPz)[itrk]);
        }
        for(int igem=0; igem<4; igem++){
        rec_x[igem] = rec_x[igem]/Etot[igem];
        rec_y[igem] = rec_y[igem]/Etot[igem];
        rec_z[igem] = rec_z[igem]/Etot[igem];
        rec_x_smear[igem] = rec_x[igem]+rand->Gaus(0,sigma);
        rec_y_smear[igem] = rec_y[igem]+rand->Gaus(0,sigma);
        }

        for(int igem=0; igem<4; igem++){
        //cout<<"gem layer "<<igem<<": x="<<rec_x[igem]<<" , y="<<rec_y[igem]<<" , z="<<rec_z[igem]<<endl;
        //cout<<"gem layer "<<igem<<": muon E ="<<muonE[igem]<<endl;
        }

        // Edep in Pb box
        //cout<<"vPbTrkid->size() = "<<vPbTrkid->size()<<endl;
        for(int itrk=0; itrk<vPbTrkid->size(); itrk++){
        if((*vPbTrkid)[itrk]!=1) continue;
        //cout<<"Pb pos x="<<(*vPbPosX)[itrk]<<" , y="<<(*vPbPosY)[itrk]<<" , z="<<(*vPbPosZ)[itrk]<<endl;
        }

        // angle <incoming vector, outcoming vector>
        TVector3 *Pos1 = new TVector3(rec_x[0],rec_y[0],rec_z[0]);
        TVector3 *Pos2 = new TVector3(rec_x[1],rec_y[1],rec_z[1]);
        TVector3 *Pos3 = new TVector3(rec_x[2],rec_y[2],rec_z[2]);
        TVector3 *Pos4 = new TVector3(rec_x[3],rec_y[3],rec_z[3]);
        TVector3 Veci = *Pos2 - *Pos1;
        TVector3 Veco = *Pos4 - *Pos3;
        angle = Veci.Angle(Veco) * 180 / M_PI;
        //cout<<"angle = "<<angle/3.1415926*180.<<endl;
        //angle = cal_ang(rec_x[0],rec_y[0],rec_z[0],rec_x[1],rec_y[1],rec_z[1],rec_x[2],rec_y[2],rec_z[2],rec_x[3],rec_y[3],rec_z[3]); // old function

        TVector3 *Pos1_smear = new TVector3(rec_x_smear[0],rec_y_smear[0],rec_z[0]);
        TVector3 *Pos2_smear = new TVector3(rec_x_smear[1],rec_y_smear[1],rec_z[1]);
        TVector3 *Pos3_smear = new TVector3(rec_x_smear[2],rec_y_smear[2],rec_z[2]);
        TVector3 *Pos4_smear = new TVector3(rec_x_smear[3],rec_y_smear[3],rec_z[3]);
        TVector3 Veci_smear = *Pos2_smear - *Pos1_smear;
        TVector3 Veco_smear = *Pos4_smear - *Pos3_smear;
        angle_smear = Veci_smear.Angle(Veco_smear) * 180 / M_PI;
        //angle_smear = cal_ang(rec_x_smear[0],rec_y_smear[0],rec_z[0],rec_x_smear[1],rec_y_smear[1],rec_z[1],rec_x_smear[2],rec_y_smear[2],rec_z[2],rec_x_smear[3],rec_y_smear[3],rec_z[3]);
        angle_gem12 = pmom[0]->Angle(*pmom[1]);


        // calculate PoCA, DCA

        if(CheckPoCAStatus({rec_x[0],rec_y[0],rec_z[0]},
                           {rec_x[1],rec_y[1],rec_z[1]},
                           {rec_x[2],rec_y[2],rec_z[2]},
                           {rec_x[3],rec_y[3],rec_z[3]}))
        {
                poca_status = 1;
                V3 poca = GetPoCAPoint({rec_x[0],rec_y[0],rec_z[0]},
                                       {rec_x[1],rec_y[1],rec_z[1]},
                                       {rec_x[2],rec_y[2],rec_z[2]},
                                       {rec_x[3],rec_y[3],rec_z[3]});
                poca_x = poca.x;
                poca_y = poca.y;
                poca_z = poca.z;
                dca = GetDCA({rec_x[0],rec_y[0],rec_z[0]},
                             {rec_x[1],rec_y[1],rec_z[1]},
                             {rec_x[2],rec_y[2],rec_z[2]},
                             {rec_x[3],rec_y[3],rec_z[3]});
        }
        else poca_status = 0;


        // old function
        /*
        struct Line l1 = {{rec_x[0],rec_y[0],rec_z[0]}, {rec_x[1],rec_y[1],rec_z[1]}};
        struct Line l2 = {{rec_x[2],rec_y[2],rec_z[2]}, {rec_x[3],rec_y[3],rec_z[3]}};
        struct Line lper = line_perpendicular(l1,l2);
        poca_x = (lper.p1.x + lper.p2.x) / 2;
        poca_y = (lper.p1.y + lper.p2.y) / 2;
        poca_z = (lper.p1.z + lper.p2.z) / 2;
        dca = sqrt((lper.p1.x-lper.p2.x)*(lper.p1.x-lper.p2.x)+(lper.p1.y-lper.p2.y)*(lper.p1.y-lper.p2.y)+(lper.p1.z-lper.p2.z)*(lper.p1.z-lper.p2.z));
        poca_status = lper.status;
        */

        V3 poca_smear = GetPoCAPoint({rec_x_smear[0],rec_y_smear[0],rec_z[0]},
                                     {rec_x_smear[1],rec_y_smear[1],rec_z[1]},
                                     {rec_x_smear[2],rec_y_smear[2],rec_z[2]},
                                     {rec_x_smear[3],rec_y_smear[3],rec_z[3]});
        poca_x_smear = poca_smear.x;
        poca_y_smear = poca_smear.y;
        poca_z_smear = poca_smear.z;
        dca_smear = GetDCA({rec_x_smear[0],rec_y_smear[0],rec_z[0]},
                           {rec_x_smear[1],rec_y_smear[1],rec_z[1]},
                           {rec_x_smear[2],rec_y_smear[2],rec_z[2]},
                           {rec_x_smear[3],rec_y_smear[3],rec_z[3]});

        if(CheckPoCAStatus({rec_x_smear[0],rec_y_smear[0],rec_z[0]},
                           {rec_x_smear[1],rec_y_smear[1],rec_z[1]},
                           {rec_x_smear[2],rec_y_smear[2],rec_z[2]},
                           {rec_x_smear[3],rec_y_smear[3],rec_z[3]})) poca_status_smear = 1;
        else poca_status_smear = 0;

        tn->Fill();
}
gettimeofday(&end, NULL); // 计时结束
timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
printf("time = %ld us\n", timer);

fn->cd();
fn->Write();
fn->Close();
f->Close();

}

