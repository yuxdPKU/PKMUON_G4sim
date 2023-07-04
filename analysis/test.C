#include "dca.h"

int test() {
  double rec_x[4], rec_y[4], rec_z[4];
  
  //no smear
  rec_x[0]=1;rec_y[0]=0;rec_z[0]=0;
  rec_x[1]=0;rec_y[1]=1;rec_z[1]=0;
  rec_x[2]=0;rec_y[2]=0;rec_z[2]=0;
  rec_x[3]=0;rec_y[3]=0;rec_z[3]=-1;

  //smear
//  58.632321 * 58.696420
//  59.917919 * 60.097611
//  161.45163 * 162.36022
//  162.39760 * 163.89755
/*
  rec_x[0]=58.632321;rec_y[0]=58.696420;rec_z[0]=-513.39;
  rec_x[1]=59.917919;rec_y[1]=60.097611;rec_z[1]=-500.05;
  rec_x[2]=161.45163;rec_y[2]=162.36022;rec_z[2]=513.29;
  rec_x[3]=162.39760;rec_y[3]=163.89755;rec_z[3]=526.63;
*/

//  rec_x[0]=-1;rec_y[0]=0;rec_z[0]=0;
//  rec_x[1]=0;rec_y[1]=0;rec_z[1]=0;
//  rec_x[2]=0;rec_y[2]=0;rec_z[2]=0;
//  rec_x[3]=1;rec_y[3]=0;rec_z[3]=0;
  double dca_x, dca_y, dca_z, dca;

  struct Line l1;
  l1.p1.x=rec_x[0];
  l1.p1.y=rec_y[0];
  l1.p1.z=rec_z[0];
  l1.p2.x=rec_x[1];
  l1.p2.y=rec_y[1];
  l1.p2.z=rec_z[1];
  printf("l1: (%f %f %f)\n",l1.p2.x-l1.p1.x,l1.p2.y-l1.p1.y,l1.p2.z-l1.p1.z);
  struct Line l2;
  l2.p1.x=rec_x[2];
  l2.p1.y=rec_y[2];
  l2.p1.z=rec_z[2];
  l2.p2.x=rec_x[3];
  l2.p2.y=rec_y[3];
  l2.p2.z=rec_z[3];
  printf("l2: (%f %f %f)\n",l2.p2.x-l2.p1.x,l2.p2.y-l2.p1.y,l2.p2.z-l2.p1.z);

  line_perpendicular(l1,l2);
  printf("dca=%f p(%f %f %f)\n",pub_dca,pub_dca_x,pub_dca_y,pub_dca_z);



/*
    Line3D line1 = {{1.0, 1.0, 1.0}, {1.0, 0.0, 0.0}};
    Line3D line2 = {{2.0, 3.0, 0.0}, {0.0, 1.0, 0.0}};
    
    Point3D midpoint = calculatePerpendicularMidpoint(&line1, &line2);
    
    printf("Midpoint: (%f, %f, %f)\n", midpoint.x, midpoint.y, midpoint.z);
*/
    
    return 0;
}

