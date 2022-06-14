#ifndef DAMPING_H_
#define DAMPING_H_
#include "Vec.h"
#include "World.h"
void Damping_Func(double& , long , long , double& );
void damping_fields(Array3D<double3>& fieldE, Array3D<double3>& fieldB,const Region& domain);
void SetEnergyDampLineNull();
#endif 
