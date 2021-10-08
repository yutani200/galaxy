#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "struct.h"
#include "setup.h"
#include "run.h"

#define NMAX 200

//AU,Msun,yrで計算する。
int main(void){
  double Tend = 1e4;
  double output_dt = Tend/10000;
  double eps,eps2;
  //static double m[NMAX], x[NMAX][3], v[NMAX][3], a[NMAX][3];
  double Rstar,Mstar;
  Istar star[NMAX];
  
  Rstar = 500; //pc
  Mstar = 1e7; //Msun
  eps = 5;  //pc
  eps2 = SQR(eps);
  read_nemo(NMAX,Mstar,Rstar,eps,star);
  //make_spherical_df(NMAX,Rstar,Mstar,m,x,v,0.1,eps);
	printf("pre run");
  run(NMAX,Tend,output_dt,star);
  return 0;
}
