#include <math.h>
#include <stdlib.h>

double gaussian(void)
/* Gaussian with mean = 0.0 and dispersion = 0.0 by Box-Muller method */
{
  double x, y, r2;
  double z;

  do{
    x = 2.0*drand48() - 1.0;
    y = 2.0*drand48() - 1.0;
    r2 = x*x + y*y;
  }while(r2 >= 1.0 || r2 == 0.0);
  z = sqrt(-2.0*log(r2)/r2)*x; /* discard another Gaussian */

  return(z);
}
