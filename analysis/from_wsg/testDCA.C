#include "dca.h"
void testDCA(){
  double rec_x[4], rec_y[4], rec_z[4];
  
  rec_x[0]=-1;rec_y[0]=0;rec_z[0]=0;
  rec_x[1]=0;rec_y[1]=0;rec_z[1]=0;
  rec_x[2]=0;rec_y[2]=0;rec_z[2]=0;
  rec_x[3]=1;rec_y[3]=0;rec_z[3]=0;
  double dca_x, dca_y, dca_z, dca;

  struct Line l1;
  l1.p1.x=rec_x[0];
  l1.p1.y=rec_y[0];
  l1.p1.z=rec_z[0];
  l1.p2.x=rec_x[1];
  l1.p2.y=rec_y[1];
  l1.p2.z=rec_z[1];
  struct Line l2;
  l2.p1.x=rec_x[2];
  l2.p1.y=rec_y[2];
  l2.p1.z=rec_z[2];
  l2.p2.x=rec_x[3];
  l2.p2.y=rec_y[3];
  l2.p2.z=rec_z[3];

  line_perpendicular(l1,l2);
  printf("dca=%f p(%f %f %f)\n",pub_dca,pub_dca_x,pub_dca_y,pub_dca_z);

}
