#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "config.h"
#include "struct.h"
#include "setup.h"
#include "run.h"

//pc,1e8Msun,Myrで計算する。
int main(void){
  int NMAX;
  double eps,Rstar,Mstar;
  AP ap, *Pap;
  Pap = &ap;

  Rstar = 100*PC_CGS; //pc
  Mstar = 1e7*Msun_CGS; //Msun
  eps = 5*PC_CGS;  //pc
  NMAX = 200;
  Pap->Tend = 1e4;
  Pap->output_dt = Pap->Tend/100;
  
  ISTAR star[NMAX];
  Pap->NMAX = NMAX;
  //setup
  fprintf(stderr,"begin to setup\n");
  //read_nemo(NMAX,Mstar,Rstar,eps,star);
  make_spherical_df(Pap->NMAX,Rstar,Mstar,0.1,eps,star);
  //run
  fprintf(stderr,"begin to run\n");
  initialize_run(Pap,star);
  run(Pap,star);
  return 0;
}
