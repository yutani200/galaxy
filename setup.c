#include <math.h>
#include "unit.h"
#include "config.h"
#include "struct.h"
#include "setup.h"
#include "run.h"
#include "calc_energy.h"
#include "gaussian.h"

#define STAR 200

void make_spherical_df(int n,double R, double M, double r_v,Istar *star){
  int i,k;
  double W,sigma;
  double dum_x, dum_y, dum_z, r2;
	double rstar,mstar;
  sigma = 1.0;
	rstar = R*PC_CGS/UnitLength;
	mstar = M*Msun_CGS/UnitMass;
	fprintf(stderr,"pre drand\n");
  for(i=0; i<n; i++){
    do{
      dum_x  = 2.0*drand48() - 1.0;
      dum_y  = 2.0*drand48() - 1.0;
      dum_z  = 2.0*drand48() - 1.0;
      r2 = SQR(dum_x) + SQR(dum_y) + SQR(dum_z);
    }while(r2 > 1 || r2 == 0.0);
    star[i].x[0] = rstar*dum_x;
    star[i].x[1] = rstar*dum_y;
    star[i].x[2] = rstar*dum_z;
    star[i].m    = mstar*M/n;
  }
	fprintf(stderr,"pre calc_W\n");
  W = calc_W(n, W, r_v, star);
  sigma = sqrt(2.0*r_v*fabs(W)/3.0);
  for(i=0; i<n; i++){
    for(k=0; k<3; k++){
      star[i].v[k] = sigma*gaussian();
    }
  }
  fprintf(stderr,"v_sigma = %lf\n", sigma);
  return;
}

/*
void twobody(double m[], double x[][3], double v[][3]){
  for(int i=0;i<2;i++){
    m[i] = 500; //Msun
    for(int j=0;j<3;j++){
      x[i][j] = 0;
      v[i][j] = 0;
    }
  }
  m[0] = m[1]*10;
  x[0][0] = -100; //1AU
  x[1][0] = 100;
  v[0][1] = -100*0.210805; //50 km/s to AU/yr
  //v[1][1] = 1*0.210805; //50 km/s to AU/yr
  return;
}
*/
void read_nemo(int NMAX,double Mstar,double Rstar,double eps,Istar *star){
  int i,j,nn,dim;
  double t;
  double mstar,rstar,vstar;
	double Rxy,AccR,cos_p,sin_p;
	double cs;
  FILE *fp;
  if(NMAX == 200)
    fp=fopen("exp4disk200.ascii","r");
  if(NMAX == 1000)
    fp=fopen("pl1k.ascii","r");
  fscanf(fp,"%d%d%lf",&nn,&dim,&t);
  star[i].eps = eps;
  for(i=0;i<nn;i++) fscanf(fp,"%lf",&star[i].m);
  for(i=0;i<nn;i++) fscanf(fp,"%lf %lf %lf",&star[i].x[0],&star[i].x[1],&star[i].x[2]);
  for(i=0;i<nn;i++) fscanf(fp,"%lf %lf %lf",&star[i].v[0],&star[i].v[1],&star[i].v[2]);
  fclose(fp); 
  mstar = Mstar*Msun_CGS/UnitMass/NMAX;
  rstar = Rstar*PC_CGS/UnitLength;
  for(i=0;i<nn;i++){
    star[i].m *= mstar;
    for(int j=0;j<3;j++) star[i].x[j] *= rstar;
  }
  calc_force(nn,star);
	cs = 10;
	for(i=0;i<nn;i++){
    Rxy = sqrt(SQR(star[i].x[0])+SQR(star[i].x[1]));
		if(Rxy > eps){
      cos_p = star[i].x[0]/Rxy;
			sin_p = star[i].x[1]/Rxy;
		}else{
      cos_p = 1e-10;
			sin_p = 1e-10;
		}
		AccR = star[i].a[0]*cos_p + star[i].a[1]*sin_p;
    star[i].v[0] = -sqrt(fabs(AccR*Rxy))*sin_p;
		star[i].v[1] = sqrt(fabs(AccR*Rxy))*cos_p;
		star[i].v[2] = cs*gaussian();
	}
  star[0].m = 1e6*Msun_CGS/UnitMass;
  for(int j=0;j<3;j++) star[0].x[j] = 0.0;
  for(int j=0;j<3;j++) star[0].v[j] = 0.0;
  return;
}
