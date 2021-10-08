#include <math.h>
#include "config.h"
#include "unit.h"

double GetUnitNemoTime(const double Rstar,const double Mstar){
  return sqrt(CUBE(Rstar*PC_CGS)/(G_CGS*Mstar*Msun_CGS))/(1e6*Year_CGS);
}
