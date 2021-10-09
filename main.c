#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "struct.h"
#include "setup.h"
#include "run.h"

#define NMAX 200

//pc,1e8Msun,Myrで計算する。
int main(void){
  double Tend = 1e4; //Myr
  double output_dt = Tend/1000;
  double eps,eps2;
  //static double m[NMAX], x[NMAX][3], v[NMAX][3], a[NMAX][3];
  double Rstar,Mstar;
  Istar star[NMAX];
  
  Rstar = 100*PC_CGS; //pc
  Mstar = 1e7*Msun_CGS; //Msun
  eps = 5*PC_CGS;  //pc
  eps2 = SQR(eps);
  //setup
  //read_nemo(NMAX,Mstar,Rstar,eps,star);
  make_spherical_df(NMAX,Rstar,Mstar,0.1,eps,star);
  //run
  fprintf(stderr,"pre run\n");
  //initialize_run(NMAX,)
  run(NMAX,Tend,output_dt,star);
  return 0;
}
