#ifdef Unit_H
#define Unit_H
#endif
#include "config.h"

#define UnitLength (PC_CGS)
#define UnitMass (1e8*Msun_CGS)
#define UnitTime (1e6*Year_CGS)

double GetUnitNemoTime(const double,const double);

