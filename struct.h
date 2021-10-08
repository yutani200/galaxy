#ifdef STRUCT_H
#else
#define STRUCT_H
#include <stdio.h>

typedef struct{
  double eps;
  double m;
  double x[3];
  double v[3];
  double a[3];
} Istar;

#endif
