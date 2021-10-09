#ifdef STRUCT_H
#else
#define STRUCT_H
#include <stdio.h>

typedef struct allparticles{
  int NMAX;
  int nstep;
  int filenum;
  double Tend;
  double T;
  double output_dt;
} AP,*PAP;

typedef struct istar{
  double eps;
  double m;
  double x[3];
  double v[3];
  double a[3];
  double dt;
  double pre_dt;
} ISTAR;

#endif
