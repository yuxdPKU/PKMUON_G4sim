#include "dca.h"

void analysis(){
TRandom *rand = new TRandom();
double sigma=100*0.001;

TFile *f = new TFile("../build/root_file/all.root","");
TTree *t = (TTree*)f->Get("T1");

std::vector<double> *vReadoutPosX=0; //Edep X in readout bar
std::vector<double> *vReadoutPosY=0; //Edep Y in readout bar
std::vector<double> *vReadoutPosZ=0; //Edep Z in readout bar
std::vector<double> *vReadoutEdep=0; //Edep Energy in readout bar
std::vector<int> *vReadoutTrkid=0; //Trk id in readout bar
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
t->SetBranchAddress("vReadoutTrkid",&vReadoutTrkid);
t->SetBranchAddress("vPx",&vPx);
t->SetBranchAddress("vPy",&vPy);
t->SetBranchAddress("vPz",&vPz);
t->SetBranchAddress("vPbPosX",&vPbPosX);
t->SetBranchAddress("vPbPosY",&vPbPosY);
t->SetBranchAddress("vPbPosZ",&vPbPosZ);
t->SetBranchAddress("vPbTrkid",&vPbTrkid);
//t->Print();

//new file and new tree
TFile * fn = new TFile("mu1GeV_new.root","recreate");
TTree * tn = new TTree("T1","tree");
tn = t->CloneTree(0);

double rec_x[4], rec_y[4], rec_z[4];
double rec_x_smear[4], rec_y_smear[4];
double angle, angle_smear;
double dca_x, dca_y, dca_z, dca;
double dca_x_smear, dca_y_smear, dca_z_smear, dca_smear;
double angle_gem12;
tn->Branch("rec_x", rec_x, "rec_x[4]/D");
tn->Branch("rec_y", rec_y, "rec_y[4]/D");
tn->Branch("rec_z", rec_z, "rec_z[4]/D");
tn->Branch("rec_x_smear", rec_x_smear, "rec_x_smear[4]/D");
tn->Branch("rec_y_smear", rec_y_smear, "rec_y_smear[4]/D");
tn->Branch("angle", &angle, "angle/D");
tn->Branch("angle_smear", &angle_smear, "angle_smear/D");
tn->Branch("dca_x", &dca_x, "dca_x/D");
tn->Branch("dca_y", &dca_y, "dca_y/D");
tn->Branch("dca_z", &dca_z, "dca_z/D");
tn->Branch("dca", &dca, "dca/D");
tn->Branch("dca_x_smear", &dca_x_smear, "dca_x_smear/D");
tn->Branch("dca_y_smear", &dca_y_smear, "dca_y_smear/D");
tn->Branch("dca_z_smear", &dca_z_smear, "dca_z_smear/D");
tn->Branch("dca_smear", &dca_smear, "dca_smear/D");
tn->Branch("angle_gem12", &angle_gem12, "angle_gem12/D");

double gem_z_sim[4] = {-513.34, -500, 513.24, 526.58};
int nevent = t->GetEntries();
for(int ievent=0; ievent<nevent; ievent++){
cout<<endl<<"running event "<<ievent<<endl;
        t->GetEntry(ievent);
        for(int igem=0; igem<4; igem++){
        rec_x[igem]=0;
        rec_y[igem]=0;
        rec_z[igem]=0;
        }
        double Etot[4]={0,0,0,0};
TVector3 *pmom[4];
//cout<<"nhit = "<<vReadoutTrkid->size()<<endl;
        for(int itrk=0; itrk<vReadoutTrkid->size(); itrk++){
if((*vReadoutTrkid)[itrk]!=1) continue;

                int igem=0;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[0])<1 ) igem=0;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[1])<1 ) igem=1;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[2])<1 ) igem=2;
                if( fabs((*vReadoutPosZ)[itrk]-gem_z_sim[3])<1 ) igem=3;
                Etot[igem]+=(*vReadoutEdep)[itrk];
                rec_x[igem]+=(*vReadoutEdep)[itrk]*(*vReadoutPosX)[itrk];
                rec_y[igem]+=(*vReadoutEdep)[itrk]*(*vReadoutPosY)[itrk];
                rec_z[igem]+=(*vReadoutEdep)[itrk]*(*vReadoutPosZ)[itrk];
                //pmom[igem] = TVector3(vPx[itrk],vPy[itrk],vPz[itrk]);
                pmom[igem] = new TVector3((*vPx)[itrk],(*vPy)[itrk],(*vPz)[itrk]);
//        cout<<"igem = "<<igem<<endl;
//        cout<<"vReadoutTrkid = "<<(*vReadoutTrkid)[itrk]<<endl;
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
        }

        // Edep in Pb box
        //cout<<"vPbTrkid->size() = "<<vPbTrkid->size()<<endl;
        for(int itrk=0; itrk<vPbTrkid->size(); itrk++){
        if((*vPbTrkid)[itrk]!=1) continue;
        //cout<<"Pb pos x="<<(*vPbPosX)[itrk]<<" , y="<<(*vPbPosY)[itrk]<<" , z="<<(*vPbPosZ)[itrk]<<endl;
        }

        angle = cal_ang(rec_x[0],rec_y[0],rec_z[0],rec_x[1],rec_y[1],rec_z[1],rec_x[2],rec_y[2],rec_z[2],rec_x[3],rec_y[3],rec_z[3]);
        //cout<<"angle = "<<angle<<endl;
        angle_smear = cal_ang(rec_x_smear[0],rec_y_smear[0],rec_z[0],rec_x_smear[1],rec_y_smear[1],rec_z[1],rec_x_smear[2],rec_y_smear[2],rec_z[2],rec_x_smear[3],rec_y_smear[3],rec_z[3]);
        angle_gem12 = pmom[0]->Angle(*pmom[1]);

    struct Line l1 = {{rec_x[0],rec_y[0],rec_z[0]}, {rec_x[1],rec_y[1],rec_z[1]}};
    struct Line l2 = {{rec_x[2],rec_y[2],rec_z[2]}, {rec_x[3],rec_y[3],rec_z[3]}};
    struct Line lper = line_perpendicular(l1,l2);
    dca_x = (lper.p1.x + lper.p2.x) / 2;
    dca_y = (lper.p1.y + lper.p2.y) / 2;
    dca_z = (lper.p1.z + lper.p2.z) / 2;
    dca = sqrt((lper.p1.x-lper.p2.x)*(lper.p1.x-lper.p2.x)+(lper.p1.y-lper.p2.y)*(lper.p1.y-lper.p2.y)+(lper.p1.z-lper.p2.z)*(lper.p1.z-lper.p2.z));
    //cout<<"dca_x = "<<dca_x<<" , dca_y = "<<dca_y<<" , dca_z = "<<dca_z<<endl;
    //cout<<"dca = "<<dca<<endl;

    struct Line l1_smear = {{rec_x_smear[0],rec_y_smear[0],rec_z[0]}, {rec_x_smear[1],rec_y_smear[1],rec_z[1]}};
    struct Line l2_smear = {{rec_x_smear[2],rec_y_smear[2],rec_z[2]}, {rec_x_smear[3],rec_y_smear[3],rec_z[3]}};
    struct Line lper_smear = line_perpendicular(l1_smear,l2_smear);
    dca_x_smear = (lper_smear.p1.x + lper_smear.p2.x) / 2;
    dca_y_smear = (lper_smear.p1.y + lper_smear.p2.y) / 2;
    dca_z_smear = (lper_smear.p1.z + lper_smear.p2.z) / 2;
    dca_smear = sqrt((lper_smear.p1.x-lper_smear.p2.x)*(lper_smear.p1.x-lper_smear.p2.x)+(lper_smear.p1.y-lper_smear.p2.y)*(lper_smear.p1.y-lper_smear.p2.y)+(lper_smear.p1.z-lper_smear.p2.z)*(lper_smear.p1.z-lper_smear.p2.z));

        tn->Fill();
}

fn->cd();
fn->Write();
fn->Close();
f->Close();

}

