#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "config.h"
#include "output.c"
#include "unit.h"
#include "calc_energy.h"

double calc_dt(int n,double pre_dt,Istar *star);
double Iterate_dt(int n,double pre_dt,Istar *star);
void leap_frog(int n,double dt,Istar *star);

void calc_force(int n,Istar *star){
  int i,j,k;
  double r[3],r3;
  double G;
  G = G_CGS*UnitMass*SQR(UnitTime)/CUBE(UnitLength);
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      star[i].a[k] = 0.0;
    }   
  }
  for(i=0;i<n-1;i++){
    //calc_r2(n,i,r2,r,x);
    for(j=i+1;j<n;j++){
      for(k=0;k<3;k++){
        r[k] = star[j].x[k] - star[i].x[k];
      }   
      r3 = CUBE(sqrt(SQR(r[0]) + SQR(r[1]) + SQR(r[2])+SQR(star[i].eps)));
      //fprintf(stderr,"r[0] %lf\n",r[0]);
      for(k=0;k<3;k++){
        star[i].a[k] += G*star[j].m*r[k]/r3;
        star[j].a[k] += -G*star[i].m*r[k]/r3;
      }   
    }
  }
  return;
}

double get_dt(int n,double T,double pre_dt,double output_time,Istar *star){
  double vel_tmp,acc_tmp,dt_tmp,T_tmp;
  //fprintf(stderr," %.16lf\n",pre_dt);
  //dt_tmp = Iterate_dt(n,pre_dt,star);
  dt_tmp = calc_dt(n,pre_dt,star);
  T_tmp = T+dt_tmp;
  if(output_time<T_tmp && output_time - T < dt_tmp){
    dt_tmp = output_time - T;
  }
  //fprintf(stderr,"dt %lf max-v %lf\n",dt_tmp,vel_tmp);
  //fprintf(stderr,"%lf %.16lf\n",T,dt_tmp);
  return dt_tmp;
}

double calc_dt(int n,double pre_dt,Istar *star){
  double max=0;
  double acc_tmp,dt;
  for(int i=0;i<n;i++){
    //vel_tmp = sqrt(SQR(v[i][0])+SQR(v[i][1])+SQR(v[i][2]));
		acc_tmp = sqrt(SQR(star[i].a[0])+SQR(star[i].a[1])+SQR(star[i].a[2]));
    //fprintf(stderr,"vel %lf\n",vel_tmp);
    if(max<acc_tmp){
      max = acc_tmp;
    }
  }
  //dt  = (0.01*sqrt(0.1*star[0].eps/acc_tmp) + pre_dt)*0.5;
  dt  = 0.01*sqrt(0.1*star[0].eps/acc_tmp);
  return dt;
}

double Iterate_dt(const int n,double pre_dt,Istar *star){
  int i,j,k;
  double dt;
  Istar star_tmp[200];
  for(i=0;i<n;i++){
    star_tmp[i] = star[i];
  }
  //printf("%p %p\n",star_tmp,star);
  //static double tmp_x[NMAX][3], tmp_v[NMAX][3], tmp_a[NMAX][3];
  dt = calc_dt(n,pre_dt,star_tmp);
  //fprintf(stderr," %.16lf\n",dt);
  for(i=0;i<3;i++){
    for(j=0;j<n;j++){
      star_tmp[j] = star[j];
    }
    leap_frog(n,dt,star_tmp);
    dt = calc_dt(n,pre_dt,star_tmp);
    //printf("i: %d dt %.16lf \n",i,dt);
  }
  //fprintf(stderr," %.16lf\n",pre_dt);
  for(i=0;i<n;i++){
    if(star_tmp[i].x[0] - star[i].x[0] == 0){
      printf("warning\n");
    }
  }
  return dt;
}

void euler(int n,double dt,Istar *star){
  int i,k;
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      star[i].x[k] += star[i].v[k]*dt;
      star[i].v[k] += star[i].a[k]*dt;
    }
  }
  calc_force(n,star);
  return;
}

void sympletic_euler(int n,double dt,Istar *star){
  int i,k;
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      star[i].x[k] += star[i].v[k]*dt;
    }
  }
  calc_force(n,star);
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      star[i].v[k] += star[i].a[k]*dt;
    }
  }
  return;
}

void leap_frog(int n,double dt,Istar *star){
  int i,k;
  double v_half[n][3];
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      v_half[i][k] = star[i].v[k] + 0.5*star[i].a[k]*dt;
      star[i].x[k] += v_half[i][k]*dt;
    }   
  }
  calc_force(n,star);
  for(i=0;i<n;i++){
    for(k=0;k<3;k++){
      star[i].v[k] = v_half[i][k] + 0.5*star[i].a[k]*dt;
    }   
  }
  return;
}

void run(int n,double Tend,double output_dt,Istar *star){
  int nn=1;
  double T,dt,H,H0;
  T = 0.0;
  output(n,T,fabs(1-H/H0),dt,star);
  while(T < Tend){
    dt = get_dt(n,T,dt,output_dt*nn,star);
    //euler(n,dt,star);
    //sympletic_euler(n,m,x,vel,acc,dt,eps2);
    leap_frog(n,dt,star);
    H = calc_energy(n,star);
    if(nn==1.0){
      H0 = H;
    }
    //fprintf(stderr,"%lf %.16lf\n",T,dt);
    if(output_dt*nn==T){
      output(n,T,fabs(1-H/H0),dt,star);
      printf("%lf %.16lf %lf\n",T,fabs(1.0-H/H0),H0);
      nn += 1;
    }
    T += dt;
  }
  return;
}
