#pragma cling add_include_path("/opt/homebrew/opt/boost/include")
#pragma cling add_library_path("/opt/homebrew/opt/boost/lib")
#pragma cling load("libboost_filesystem.dylib")
#pragma cling load("libboost_system.dylib")

#include <iostream>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <sys/time.h>

using namespace std;

void analysis(){
TRandom *rand = new TRandom();
double sigma=100*0.001;//100um

TFile *f = new TFile("../../build/root_file/all.root","");
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

double muonE[2];
double rec_x[2], rec_y[2], rec_z[2];
double px[2], py[2], pz[2];
double rec_x_smear[2], rec_y_smear[2];
double Eout, Pxout, Pyout, Pzout, Xout, Yout, Zout;
double angle, angle_smear;
double poca_x, poca_y, poca_z, dca;
int poca_status;
double poca_x_smear, poca_y_smear, poca_z_smear, dca_smear;
int poca_status_smear;
double angle_gem12;

double Z1 = 0.;
double Z2 = 13.34;
double gem_z_sim[2] = {Z1,Z2};

struct timeval start;
struct timeval end;
unsigned long timer;

gettimeofday(&start, NULL); // 计时开始

ofstream outfile("UpperGEMmuon.dat");
outfile.setf(ios::fixed,ios::floatfield); outfile.precision(6);

int number = 0;
int nevent = t->GetEntries();
for(int ievent=0; ievent<nevent; ievent++){
        if (ievent % (int)(nevent / 10) == 0) cout << "Processing progress: " << ievent / (int)(nevent / 10) << "0%" << endl;
        t->GetEntry(ievent);

        // initialization
        for(int igem=0; igem<2; igem++){
                rec_x[igem]=0;
                rec_y[igem]=0;
                rec_z[igem]=0;
                px[igem]=0;
                py[igem]=0;
                pz[igem]=0;
        }
        double Etot[2]={0,0};
        bool gemstatus[2]={false};
        bool outstatus=false;

        TVector3 *pmom[2];
        double nhit = vReadoutTrkid->size();
        if(nhit < 3) continue;
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
                int igem=-1;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[0])<0.1 ) igem=0;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[1])<0.1 ) igem=1;
                if( fabs((*vReadoutPosZ)[itrk]-13.44)<0.01){
                        Eout=(*vReadoutE)[itrk];
                        Pxout=(*vPx)[itrk];
                        Pyout=(*vPy)[itrk];
                        Pzout=(*vPz)[itrk];
                        Xout=(*vReadoutPosX)[itrk];
                        Yout=(*vReadoutPosY)[itrk];
                        Zout=(*vReadoutPosZ)[itrk];
                        outstatus=true;
                        continue;
                }
                if(igem==-1) continue;
                // cout<<"igem = "<<igem<<endl;

                // energy weight position
                muonE[igem]=(*vReadoutE)[itrk];
                Etot[igem]+=(*vReadoutEdep)[itrk];
                rec_x[igem]+=(*vReadoutEdep)[itrk]*(*vReadoutPosX)[itrk];
                rec_y[igem]+=(*vReadoutEdep)[itrk]*(*vReadoutPosY)[itrk];
                //rec_z[igem]+=(*vReadoutEdep)[itrk]*(*vReadoutPosZ)[itrk];
                rec_z[igem]+=(*vReadoutEdep)[itrk]*gem_z_sim[igem];
                px[igem]=(*vPx)[itrk];
                py[igem]=(*vPy)[itrk];
                pz[igem]=(*vPz)[itrk];
                pmom[igem] = new TVector3((*vPx)[itrk],(*vPy)[itrk],(*vPz)[itrk]);
                gemstatus[igem]=true;
        }

        bool hasFalse=false;
        for (bool element : gemstatus) {
            if (!element) {
                hasFalse = true;
                break;
            }
        }
        if (hasFalse && !outstatus) {
            //cout<<"Don't have signals in all GEM detector"<<endl;
            continue;
        }

        for(int igem=0; igem<2; igem++){
        rec_x[igem] = rec_x[igem]/Etot[igem];
        rec_y[igem] = rec_y[igem]/Etot[igem];
        rec_z[igem] = rec_z[igem]/Etot[igem];
        rec_x_smear[igem] = rec_x[igem]+rand->Gaus(0,sigma);
        rec_y_smear[igem] = rec_y[igem]+rand->Gaus(0,sigma);
        }
        if(fabs(rec_x[0])<1/sqrt(2)/2*1000 && fabs(rec_y[0])<1/sqrt(2)/2*1000) number++;

        for(int igem=0; igem<2; igem++){
        //cout<<"gem layer "<<igem<<": x="<<rec_x[igem]<<" , y="<<rec_y[igem]<<" , z="<<rec_z[igem]<<endl;
        //cout<<"gem layer "<<igem<<": muon E ="<<muonE[igem]<<endl;
        }

        //outfile<<ievent<<" "<<rec_x[1]<<" "<<rec_y[1]<<" "<<rec_z[1]<<" "<<px[1]<<" "<<py[1]<<" "<<pz[1]<<" "<<muonE[1]<<endl;
        outfile<<ievent<<" "<<Xout<<" "<<Yout<<" "<<Zout<<" "<<Pxout<<" "<<Pyout<<" "<<Pzout<<" "<<Eout<<endl;

}
cout<<"In X:[-1/2sqrt2,1/2sqrt2] m and Y:[-1/2sqrt2,1/2sqrt2] m, number = "<<number<<endl;

gettimeofday(&end, NULL); // 计时结束
timer = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
printf("time = %ld us\n", timer);

f->Close();

}

