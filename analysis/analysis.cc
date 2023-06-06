using namespace std;

#include <iostream>
#include <cmath>
#include <vector>

struct point {
    double x;
    double y;
    double z;
};

struct Line {
    struct point p1;
    struct point p2;
};

double cal_ang(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4);
struct Line line_perpendicular(struct Line l1, struct Line l2);
struct point cal_normal_vector(struct Line l1, struct Line l2);
struct point normalization(struct point p1);

void analysis(){
TRandom *rand = new TRandom();
double sigma=100*0.001;

TFile *f = new TFile("./build/mu1GeV.root","");
TTree *t = (TTree*)f->Get("T1");

std::vector<double> *vReadoutPosX=0; //Edep X in readout bar
std::vector<double> *vReadoutPosY=0; //Edep Y in readout bar
std::vector<double> *vReadoutPosZ=0; //Edep Z in readout bar
std::vector<double> *vReadoutEdep=0; //Edep Energy in readout bar
std::vector<int> *vReadoutTrkid=0; //Trk id in readout bar
std::vector<double> *vPx=0;
std::vector<double> *vPy=0;
std::vector<double> *vPz=0;

t->SetBranchAddress("vReadoutPosX",&vReadoutPosX);
t->SetBranchAddress("vReadoutPosY",&vReadoutPosY);
t->SetBranchAddress("vReadoutPosZ",&vReadoutPosZ);
t->SetBranchAddress("vReadoutEdep",&vReadoutEdep);
t->SetBranchAddress("vReadoutTrkid",&vReadoutTrkid);
t->SetBranchAddress("vPx",&vPx);
t->SetBranchAddress("vPy",&vPy);
t->SetBranchAddress("vPz",&vPz);
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
        t->GetEntry(ievent);
        for(int igem=0; igem<4; igem++){
        rec_x[igem]=0;
        rec_y[igem]=0;
        rec_z[igem]=0;
        }
        double Etot[4]={0,0,0,0};
TVector3 *pmom[4];
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
        }
        for(int igem=0; igem<4; igem++){
        rec_x[igem] = rec_x[igem]/Etot[igem];
        rec_y[igem] = rec_y[igem]/Etot[igem];
        rec_z[igem] = rec_z[igem]/Etot[igem];
        rec_x_smear[igem] = rec_x[igem]+rand->Gaus(0,sigma);
        rec_y_smear[igem] = rec_y[igem]+rand->Gaus(0,sigma);
        }
        angle = cal_ang(rec_x[0],rec_y[0],rec_z[0],rec_x[1],rec_y[1],rec_z[1],rec_x[2],rec_y[2],rec_z[2],rec_x[3],rec_y[3],rec_z[3]);
        angle_smear = cal_ang(rec_x_smear[0],rec_y_smear[0],rec_z[0],rec_x_smear[1],rec_y_smear[1],rec_z[1],rec_x_smear[2],rec_y_smear[2],rec_z[2],rec_x_smear[3],rec_y_smear[3],rec_z[3]);
        angle_gem12 = pmom[0]->Angle(*pmom[1]);

    struct Line l1 = {{rec_x[0],rec_y[0],rec_z[0]}, {rec_x[1],rec_y[1],rec_z[1]}};
    struct Line l2 = {{rec_x[2],rec_y[2],rec_z[2]}, {rec_x[3],rec_y[3],rec_z[3]}};
    struct Line lper = line_perpendicular(l1,l2);
    dca_x = (lper.p1.x + lper.p2.x) / 2;
    dca_y = (lper.p1.y + lper.p2.y) / 2;
    dca_z = (lper.p1.z + lper.p2.z) / 2;
    dca = sqrt((lper.p1.x-lper.p2.x)*(lper.p1.x-lper.p2.x)+(lper.p1.y-lper.p2.y)*(lper.p1.y-lper.p2.y)+(lper.p1.z-lper.p2.z)*(lper.p1.z-lper.p2.z));

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

double cal_ang(double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3, double z3, double x4, double y4, double z4){
   // 计算两条直线的向量
   double v1x = x2 - x1;
   double v1y = y2 - y1;
   double v1z = z2 - z1;
   double v2x = x4 - x3;
   double v2y = y4 - y3;
   double v2z = z4 - z3;

   // 计算两条直线的点积和向量长度
   double dot_product = v1x * v2x + v1y * v2y + v1z * v2z;
   double length1 = sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
   double length2 = sqrt(v2x * v2x + v2y * v2y + v2z * v2z);

   // 计算夹角
   double angle = acos(dot_product / (length1 * length2)) * 180 / M_PI;
   return angle;
}

struct Line line_perpendicular(struct Line l1, struct Line l2) {
    double x1 = l1.p1.x;
    double y1 = l1.p1.y;
    double z1 = l1.p1.z;
    double x2 = l1.p2.x;
    double y2 = l1.p2.y;
    double z2 = l1.p2.z;
    double x3 = l2.p1.x;
    double y3 = l2.p1.y;
    double z3 = l2.p1.z;
    double x4 = l2.p2.x;
    double y4 = l2.p2.y;
    double z4 = l2.p2.z;

    // step 1: the normal vector of l1 and l2
    // n12 = l1 x l2
    struct point n12_vec = cal_normal_vector(l1, l2);
    //cout<<"n12_vec = ("<<n12_vec.x<<","<<n12_vec.y<<","<<n12_vec.z<<")"<<endl;
    struct point n12_unit = normalization(n12_vec);
    //cout<<"n12_unit = ("<<n12_unit.x<<","<<n12_unit.y<<","<<n12_unit.z<<")"<<endl;

    // the direction vector of l1 and l2
    struct point l1_vec;
    l1_vec.x = x2-x1;
    l1_vec.y = y2-y1;
    l1_vec.z = z2-z1;
    struct point l1_unit = normalization(l1_vec);
    //cout<<"l1_unit = ("<<l1_unit.x<<","<<l1_unit.y<<","<<l1_unit.z<<")"<<endl;
    struct point l2_vec;
    l2_vec.x = x4-x3;
    l2_vec.y = y4-y3;
    l2_vec.z = z4-z3;
    //cout<<"l2_vec = ("<<l2_vec.x<<","<<l2_vec.y<<","<<l2_vec.z<<")"<<endl;
    struct point l2_unit = normalization(l2_vec);
    //cout<<"l2_unit = ("<<l2_unit.x<<","<<l2_unit.y<<","<<l2_unit.z<<")"<<endl;

    // the normal vector of l1 and n12
    struct Line l12 = {{0,0,0}, {n12_unit.x,n12_unit.y,n12_unit.z}};
    struct point n12prime_vec = cal_normal_vector(l1, l12);
    struct point n12prime_unit = normalization(n12prime_vec);
    //cout<<"n12prime_unit = ("<<n12prime_unit.x<<","<<n12prime_unit.y<<","<<n12prime_unit.z<<")"<<endl;

    // step 2: find point in l2 in the plane (n12 x l1)
    // define intersection point is p1
    // xp1 = x3 + t1*l2_unit.x
    // yp1 = y3 + t1*l2_unit.y
    // zp1 = z3 + t1*l2_unit.z
    // line (p1,x1) should be perpendicular to n12prime_unit
    // line (p1,x1) dot n12prime_unit = 0

    double t1 = -(x3*n12prime_unit.x+y3*n12prime_unit.y+z3*n12prime_unit.z) / (l2_unit.x*n12prime_unit.x+l2_unit.y*n12prime_unit.y+l2_unit.z*n12prime_unit.z);
    //cout<<"t1 = "<<t1<<endl;
    struct point p1;
    p1.x = x3+t1*l2_unit.x;
    p1.y = y3+t1*l2_unit.y;
    p1.z = z3+t1*l2_unit.z;
    //cout<<"p1 = ("<<p1.x<<","<<p1.y<<","<<p1.z<<")"<<endl;

    // step 3: find point p2 in l1
    // line (p2, p1) is parrallel to n12_unit
    // xp2 = xp1 + t2*n12_unit.x
    // yp2 = yp1 + t2*n12_unit.y
    // zp2 = zp1 + t2*n12_unit.z
    // line (p2,x1) should parrallel to l1
    // note: 
    // line (p2,x1) cross l1 = 0
    double t2 = 0;
    if( (n12_unit.y*l1_unit.z-n12_unit.z*l1_unit.y)!=0 ) t2 = ((p1.z-z1)*l1_unit.y-(p1.y-y1)*l1_unit.z) / (n12_unit.y*l1_unit.z-n12_unit.z*l1_unit.y);
    else if ( (n12_unit.x*l1_unit.y-n12_unit.y*l1_unit.x)!=0 ) t2 = ((p1.y-y1)*l1_unit.x-(p1.x-x1)*l1_unit.y) / (n12_unit.x*l1_unit.y-n12_unit.y*l1_unit.x);
    else if ( (n12_unit.z*l1_unit.x-n12_unit.x*l1_unit.z)!=0 ) t2 = ((p1.x-x1)*l1_unit.z-(p1.z-z1)*l1_unit.x) / (n12_unit.z*l1_unit.x-n12_unit.x*l1_unit.z);
    //cout<<"t2 = "<<t2<<endl;
    struct point p2;
    p2.x = p1.x+t2*n12_unit.x;
    p2.y = p1.y+t2*n12_unit.y;
    p2.z = p1.z+t2*n12_unit.z;
    //cout<<"p2 = ("<<p2.x<<","<<p2.y<<","<<p2.z<<")"<<endl;

    struct Line lper;
    lper.p1 = p1;
    lper.p2 = p2;

    return lper;
}

struct point cal_normal_vector(struct Line l1, struct Line l2){
    double x1 = l1.p1.x;
    double y1 = l1.p1.y;
    double z1 = l1.p1.z;
    double x2 = l1.p2.x;
    double y2 = l1.p2.y;
    double z2 = l1.p2.z;
    double x3 = l2.p1.x;
    double y3 = l2.p1.y;
    double z3 = l2.p1.z;
    double x4 = l2.p2.x;
    double y4 = l2.p2.y;
    double z4 = l2.p2.z;

    double x21 = x2-x1;
    double y21 = y2-y1;
    double z21 = z2-z1;
    double x43 = x4-x3;
    double y43 = y4-y3;
    double z43 = z4-z3;

    double xn = (y21*z43-y43*z21);
    double yn = -(x21*z43-x43*z21);
    double zn = (x21*y43-x43*y21);

    struct point normal_direction;
    normal_direction.x = xn;
    normal_direction.y = yn;
    normal_direction.z = zn;

    return normal_direction;
}

struct point normalization(struct point p1){
    double x1 = p1.x;
    double y1 = p1.y;
    double z1 = p1.z;
    //cout<<"x1 = "<<x1<<" , y1 = "<<y1<<" , z1 = "<<z1<<endl;

    double length = sqrt(x1*x1+y1*y1+z1*z1);
    //cout<<"length = "<<length<<endl;
    struct point unit;
    unit.x = x1/length;
    unit.y = y1/length;
    unit.z = z1/length;

    return unit;
}
