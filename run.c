#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "config.h"
#include "output.c"
#include "unit.h"
#include "calc_energy.h"
#include "flags.h"

double calc_dt(int NMAX,ISTAR *star);
double Iterate_dt(int NMAX,ISTAR *star);
void leap_frog(int NMAX,double dt,ISTAR *star);

void calc_force(int NMAX,ISTAR *star){
  int i,j,k;
  double r[3],r3;
  double G;
  G = G_CGS*UnitMass*SQR(UnitTime)/CUBE(UnitLength);
  G = 1;
  for(i=0;i<NMAX;i++){
    for(k=0;k<3;k++){
      star[i].a[k] = 0.0;
    }   
  }
  for(i=0;i<NMAX-1;i++){
    //calc_r2(NMAX,i,r2,r,x);
    //fprintf(stderr,"star[i].x[0] = %f\n",star[i].x[0]);
    for(j=i+1;j<NMAX;j++){
      for(k=0;k<3;k++){
        r[k] = star[j].x[k] - star[i].x[k];
      } 
      r3 = CUBE(sqrt(SQR(r[0]) + SQR(r[1]) + SQR(r[2])+SQR(star[i].eps)));
      //fprintf(stderr,"m = %lf, r3 %lf\n",star[i].m,r[0]);
      for(k=0;k<3;k++){
        star[i].a[k] += G*star[j].m*r[k]/r3;
        star[j].a[k] += -G*star[i].m*r[k]/r3;
      }   
    }
  }
  return;
}

double get_dt(double dt,PAP Pap,ISTAR *star){
  double vel_tmp,acc_tmp,dt_tmp,T_tmp;
  //fprintf(stderr,"111NMAX = %d,T=%lf,output_time=%lf\n",Pap->NMAX,Pap->T,Pap->output_dt*Pap->filenum);
  dt_tmp = calc_dt(Pap->NMAX,star);
  T_tmp = Pap->T+dt_tmp;
  if(Pap->output_dt*Pap->filenum<T_tmp && Pap->output_dt*Pap->filenum - Pap->T < dt_tmp){
    dt_tmp = Pap->output_dt*Pap->filenum - Pap->T;
  }
  //fprintf(stderr,"dt %lf max-v %lf\n",dt_tmp,vel_tmp);
  //fprintf(stderr,"%lf %.16lf\n",T,dt_tmp);
  return dt_tmp;
}

double calc_dt(int NMAX,ISTAR *star){
  double max=0;
  double acc_tmp,vel_tmp,dt;
  for(int i=0;i<NMAX;i++){
    vel_tmp = sqrt(SQR(star[i].v[0])+SQR(star[i].v[1])+SQR(star[i].v[2]));
		//acc_tmp = sqrt(SQR(star[i].a[0])+SQR(star[i].a[1])+SQR(star[i].a[2]));
    //fprintf(stderr,"vel %lf\n",vel_tmp);
    if(max<vel_tmp){
      max = vel_tmp;
    }
    //fprintf(stderr,"222max_vel = %f\n",vel_tmp);
    //fprintf(stderr,"star[%d]:a[1] = %f,a[2] = %f,a[3] = %f\n",i,star[i].a[0],star[i].a[1],star[i].a[2]);
  }
  //dt  = (0.01*sqrt(star[0].eps/vel_tmp) + pre_dt)*0.5;
  dt  = (0.01*star[0].eps/vel_tmp);
  //dt  = (0.01*sqrt(0.1*star[0].eps/acc_tmp) + pre_dt)*0.5;
  //dt  = 0.01*sqrt(star[0].eps/acc_tmp);
  //fprintf(stderr,"max_acc = %f\n",acc_tmp);
  return dt;
}

void euler(int NMAX,double dt,ISTAR *star){
  int i,k;
  for(i=0;i<NMAX;i++){
    for(k=0;k<3;k++){
      star[i].x[k] += star[i].v[k]*dt;
      star[i].v[k] += star[i].a[k]*dt;
    }
  }
  calc_force(NMAX,star);
  return;
}

void sympletic_euler(int NMAX,double dt,ISTAR *star){
  int i,k;
  for(i=0;i<NMAX;i++){
    for(k=0;k<3;k++){
      star[i].x[k] += star[i].v[k]*dt;
    }
  }
  calc_force(NMAX,star);
  for(i=0;i<NMAX;i++){
    for(k=0;k<3;k++){
      star[i].v[k] += star[i].a[k]*dt;
    }
  }
  return;
}

void leap_frog(int NMAX,double dt,ISTAR *star){
  int i,k;
  double v_half[NMAX][3];
  for(i=0;i<NMAX;i++){
    for(k=0;k<3;k++){
      v_half[i][k] = star[i].v[k] + 0.5*star[i].a[k]*dt;
      //fprintf(stderr,"star[i].x[k]=%lf, v[k] = %lf, a[k] = %lf, dt = %lf\n",star[i].x[k],star[i].v[k],star[i].a[k],dt);
      star[i].x[k] += v_half[i][k]*dt;
      //fprintf(stderr,"star[i].x[k]=%lf\n",star[i].x[k]);
    }   
  }
  calc_force(NMAX,star);
  for(i=0;i<NMAX;i++){
    for(k=0;k<3;k++){
      star[i].v[k] = v_half[i][k] + 0.5*star[i].a[k]*dt;
    }   
  }
  return;
}

//初期の加速度を得る。
void initialize_run(PAP Pap,ISTAR* star){
  Pap->nstep = 0;
  Pap->filenum = 0;
  Pap->T = 0.0;
  calc_force(Pap->NMAX,star);
  output(Pap->NMAX,Pap->T,0,0,star);
  return;
}

//Tendまで計算する。
void run(PAP Pap, ISTAR* star){
  double dt,H,H0;
  fprintf(stderr,"Init Output\n");
  while(Pap->T < Pap->Tend){
    dt = get_dt(dt,Pap,star);
#ifdef EULER
    euler(Pap->NMAX,dt,star);
#elif defined(SYMPLETIC_EULER)
    sympletic_euler(Pap->NMAX,dt,star);
#elif defined(LEAP_FROG)
    leap_frog(Pap->NMAX,dt,star);
#endif
    H = calc_energy(Pap->NMAX,star);
    if(Pap->nstep==1.0){
      H0 = H;
    }

    //output
    if(Pap->output_dt*Pap->filenum==Pap->T){
      output(Pap->NMAX,Pap->T,fabs(1-H/H0),dt,star);
      //printf("%lf %.16lf %lf\n",Pap->T,fabs(1.0-H/H0),H0);
      Pap->filenum += 1;
    }
    Pap->T += dt;
    Pap->nstep += 1;
    fprintf(stderr,"step = %d, dt = %f, T = %f\n",Pap->nstep,dt,Pap->T);
  }
  return;
}
