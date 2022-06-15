// Author: Evgeny Berendeev
// Budker Institute of Nuclear Physics of Siberian Branch Russian Academy of Sciences
// beren@inp.nsk.su
// (c) 2022, for licensing details see the LICENSE file

#include "util.h"

//Main function simply hands off control to the Simulation class
int main(int argc, char* argv[]) {

	Operator mat(totalSize.prod()*3,totalSize.prod()*3);
	Operator identity(totalSize.prod()*3,totalSize.prod()*3);
	stencil_identity(identity);
	vector<Trip> trips;
	trips.reserve(totalSize.prod()*3);
	stencil_mat_xx(trips);
	stencil_mat_xy(trips);
	mat.setFromTriplets(trips.begin(), trips.end());

	for (int k=0; k<mat.outerSize(); ++k){
	  for (Eigen::SparseMatrix<double>::InnerIterator it(mat,k); it; ++it){
	    std::cout  << it.row() << " " <<  it.col() << " " << it.value()  << "\n";   // col index (here it is equal to k)
	   // it.index(); // inner index, here it is equal to it.row()
	  }
	}
	return(0);
}