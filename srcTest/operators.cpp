// Author: Evgeny Berendeev
// Budker Institute of Nuclear Physics of Siberian Branch Russian Academy of Sciences
// beren@inp.nsk.su
// (c) 2022, for licensing details see the LICENSE file

#include "util.h"

using namespace std;

using Eigen::Vector3d;

typedef Eigen::Triplet<double> Trip;

//create an identity operator
void stencil_identity(Operator& identity) {
	vector<Trip> trips;
	trips.reserve(totalSize.prod()*3);

	//operators should overlap at borders
	//all cells + a few border cells
	for(int x = 0; x < Nx; x++) {
		for(int y = 0; y < Ny; y++) {
			for(int z = 0; z < Nz; z++) {

				//this is a vector operator
				for(int c = 0; c < 1; c++) {
					int cur = vind(x, y, z, c);

					trips.push_back(Trip(cur, cur, 1.0));
				}
			}
		}
	}
	identity.setFromTriplets(trips.begin(), trips.end());
}

void stencil_mat_xx(vector<Trip>& trips) {
	;
	//operators should overlap at borders
	//all cells + a few border cells
	for(int x = 0; x < Nx; x++) {
		for(int y = 0; y < Ny; y++) {
			for(int z = 0; z < Nz; z++) {
				//this is a vector operator
				for(int c = 0; c < 1; c++) {
					int cur = vind(x, y, z, c);
					double valy = -0.25*dt*dt/(dy*dy);
					double valz = - 0.25*dt*dt/(dz*dz);
					trips.push_back(Trip(cur, cur, 1-2.*valy - 2.*valz ));
					if(y!=Ny-1){
						cur = vind(x, y+1, z, c);
						trips.push_back(Trip(cur, cur, valy ));
					}
					if(y!=0) {
						cur = vind(x, y-1, z, c);
						trips.push_back(Trip(cur, cur, valy ));
					}
					if(z!=Nz-1) {
						cur = vind(x, y, z+1, c);
						trips.push_back(Trip(cur, cur, valz ));
					}
					if(z!=0) {
						cur = vind(x, y, z-1, c);
						trips.push_back(Trip(cur, cur, valz ));				}
					}
			}
		}
	}
}
void stencil_mat_xy(vector<Trip>& trips) {
	;
	//operators should overlap at borders
	//all cells + a few border cells
	for(int x = 0; x < Nx; x++) {
		for(int y = 0; y < Ny; y++) {
			for(int z = 0; z < Nz; z++) {
				//this is a vector operator
				for(int c = 0; c < 1; c++) {

					int cur = vind(x, y, z, c);
					int cur1 = cur + Nx*Ny*Nz;
					double valxy = -0.25*dt*dt/(dx*dy);
					trips.push_back(Trip(cur1, cur, -valxy ));
					if(x!=Nx-1){
						cur = vind(x, y+1, z, c);
						trips.push_back(Trip(cur1, cur, valy ));
					}
					if(y!=0) {
						cur = vind(x, y-1, z, c);
						trips.push_back(Trip(cur1, cur, valy ));
					}
					if(z!=Nz-1) {
						cur = vind(x, y, z+1, c);
						trips.push_back(Trip(cur1, cur, valz ));
					}
					if(z!=0) {
						cur = vind(x, y, z-1, c);
						trips.push_back(Trip(cur1, cur, valz ));				}
					}
			}
		}
	}
}