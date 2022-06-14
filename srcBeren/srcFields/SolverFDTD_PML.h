#ifndef SolverFDTD_PML_H_
#define SolverFDTD_PML_H_
#include "Vec.h"
#include "World.h"
void solver_FDTD_PML(Array3D<double3>& fieldE, Array3D<double3>& fieldB,
                     Array3D<double3>& fieldEp, Array3D<double3>& fieldBp, 
                     const Array3D<double3>& fieldJ, const World& world);

#endif 	
