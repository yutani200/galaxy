#include <stdio.h>
#include <math.h>
#include "struct.h"
#include "config.h"
#include "unit.h"

double calc_energy(int n,ISTAR *star){
  int i,j,k;
  double e_k, e_w, r2, G;
  // kinetic energy 
  e_k = 0;
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      e_k += star[i].m*(SQR(star[i].v[k]));
    }
  }
  e_k *= 0.5;
  // potential energy
  e_w = 0;
  G = G_CGS*UnitMass*SQR(UnitTime)/CUBE(UnitLength);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      if(i!=j){
        r2 = SQR(star[j].x[0]-star[i].x[0]) + SQR(star[j].x[1]-star[i].x[1]) + SQR(star[j].x[2]-star[i].x[2]);
        e_w = e_w - G*star[i].m*star[j].m/sqrt(r2+SQR(star[j].eps));
      }
    }
  }
  e_w *= 0.5;
  return (e_k + e_w);
}

double calc_W(int n, double W, double r_v,ISTAR *star){
  int i,j;
  double r2,G;
  W = 0.0;
  G = G_CGS*Msun_CGS*SQR(1e6*Year_CGS)/CUBE(PC_CGS);
  for(i=0; i<n-1; i++){
    for(j=i+1; j<n; j++){
      r2 = SQR(star[j].x[0]-star[i].x[0]) + SQR(star[j].x[1]-star[i].x[1]) + SQR(star[j].x[2]-star[i].x[2]);
      W  = W - G*star[i].m*star[j].m/sqrt(r2 + SQR(star[j].eps));
    }
  }
  fprintf(stderr,"W = %lf\n", W);
  return(W);
}
