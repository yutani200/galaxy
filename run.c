#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "struct.h"
#include "config.h"
#include "output.c"
#include "unit.h"
#include "calc_energy.h"
#include "flags.h"

double calc_dt(int NMAX,double pre_dt,Istar *star);
double Iterate_dt(int NMAX,double pre_dt,Istar *star);
void leap_frog(int NMAX,double dt,Istar *star);

void calc_force(int NMAX,Istar *star){
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

double get_dt(int NMAX,double T,double pre_dt,double output_time,Istar *star){
  double vel_tmp,acc_tmp,dt_tmp,T_tmp;
  //fprintf(stderr," %.16lf\n",pre_dt);
  //dt_tmp = Iterate_dt(NMAX,pre_dt,star);
  dt_tmp = calc_dt(NMAX,pre_dt,star);
  T_tmp = T+dt_tmp;
  if(output_time<T_tmp && output_time - T < dt_tmp){
    dt_tmp = output_time - T;
  }
  //fprintf(stderr,"dt %lf max-v %lf\n",dt_tmp,vel_tmp);
  //fprintf(stderr,"%lf %.16lf\n",T,dt_tmp);
  return dt_tmp;
}

double calc_dt(int NMAX,double pre_dt,Istar *star){
  double max=0;
  double acc_tmp,vel_tmp,dt;
  for(int i=0;i<NMAX;i++){
    vel_tmp = sqrt(SQR(star[i].v[0])+SQR(star[i].v[1])+SQR(star[i].v[2]));
		//acc_tmp = sqrt(SQR(star[i].a[0])+SQR(star[i].a[1])+SQR(star[i].a[2]));
    //fprintf(stderr,"vel %lf\n",vel_tmp);
    if(max<vel_tmp){
      max = vel_tmp;
    }
    //fprintf(stderr,"star[%d]:a[1] = %f,a[2] = %f,a[3] = %f\n",i,star[i].a[0],star[i].a[1],star[i].a[2]);
  }
  //dt  = (0.01*sqrt(star[0].eps/vel_tmp) + pre_dt)*0.5;
  dt  = (0.01*star[0].eps/vel_tmp);
  fprintf(stderr,"max_vel = %f\n",vel_tmp);
  //dt  = (0.01*sqrt(0.1*star[0].eps/acc_tmp) + pre_dt)*0.5;
  //dt  = 0.01*sqrt(star[0].eps/acc_tmp);
  //fprintf(stderr,"max_acc = %f\n",acc_tmp);
  return dt;
}

double Iterate_dt(const int NMAX,double pre_dt,Istar *star){
  int i,j,k;
  double dt;
  Istar star_tmp[200];
  for(i=0;i<NMAX;i++){
    star_tmp[i] = star[i];
  }
  //printf("%p %p\n",star_tmp,star);
  //static double tmp_x[NMAX][3], tmp_v[NMAX][3], tmp_a[NMAX][3];
  dt = calc_dt(NMAX,pre_dt,star_tmp);
  //fprintf(stderr," %.16lf\n",dt);
  for(i=0;i<3;i++){
    for(j=0;j<NMAX;j++){
      star_tmp[j] = star[j];
    }
    leap_frog(NMAX,dt,star_tmp);
    dt = calc_dt(NMAX,pre_dt,star_tmp);
    //printf("i: %d dt %.16lf \n",i,dt);
  }
  //fprintf(stderr," %.16lf\n",pre_dt);
  for(i=0;i<NMAX;i++){
    if(star_tmp[i].x[0] - star[i].x[0] == 0){
      printf("warning\n");
    }
  }
  return dt;
}

void euler(int NMAX,double dt,Istar *star){
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

void sympletic_euler(int NMAX,double dt,Istar *star){
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

void leap_frog(int NMAX,double dt,Istar *star){
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

/*
//初期の加速度を得る為の関数。
void initialize_run(int NMAX,double T,Istar *star){
  calc_force(NMAX,star);
  output(NMAX,T)
  return;
}
*/

//whileより手前はInitialize_runに入れる。
//後、NMAXやT等の一般的な変数は構造体として纏めて管理する。
void run(int NMAX,double Tend,double output_dt,Istar *star){
  int nstep=1,file_num=1;
  double T,dt,H,H0;
  T = 0.0;
  fprintf(stderr,"Init Output\n");
  output(NMAX,T,fabs(1-H/H0),dt,star);
  calc_force(NMAX,star);
  while(T < Tend){
    dt = get_dt(NMAX,T,dt,output_dt*file_num,star);
#ifdef EULER
    euler(NMAX,dt,star);
#elif defined(SYMPLETIC_EULER)
    sympletic_euler(NMAX,dt,star);
#elif defined(LEAP_FROG)
    leap_frog(NMAX,dt,star);
#endif
    H = calc_energy(NMAX,star);
    if(nstep==1.0){
      H0 = H;
    }
    //fprintf(stderr,"%lf %.16lf\n",T,dt);

    //output
    if(output_dt*file_num==T){
      output(NMAX,T,fabs(1-H/H0),dt,star);
      printf("%lf %.16lf %lf\n",T,fabs(1.0-H/H0),H0);
      file_num += 1;
    }
    fprintf(stderr,"step = %d, dt = %f, T = %f\n",nstep,dt,T);
    T += dt;
    nstep += 1;
  }
  return;
}
