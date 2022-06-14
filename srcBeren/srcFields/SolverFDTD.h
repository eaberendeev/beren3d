#ifndef SolverFDTD_H_
#define SolverFDTD_H_

void solver_FDTD(Array3D<double3>& fieldE, Array3D<double3>& fieldB, const Array3D<double3>& fieldJ, const  World& world);

#endif 	
