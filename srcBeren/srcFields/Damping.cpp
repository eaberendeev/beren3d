#include "World.h"
#include "Damping.h"

/*void SetEnergyDampLineNull(const Region& domain){
    for( auto i = 0; i < NumCellsR_glob + 2 * SHAPE_SIZE - 1; ++i){
	energy.eDampLineB_Z1[i] = 0.;
	energy.eDampLineE_Z1[i] = 0.;
	energy.eDampLineB_Z2[i] = 0.;
	energy.eDampLineE_Z2[i] = 0.;
    }
    for( auto i = 0; i < NumCellsZ_glob + 2 * SHAPE_SIZE - 1; ++i){
	energy.eDampLineB_R[i] = 0.;
	energy.eDampLineE_R[i] = 0.;
    }
}*/
void Damping_Func(double& source, long i, long maxi, double& energyDamp){
    double koeff=0.8;
    double a, damp;

    a = (1.0-koeff)/(maxi*maxi);
    damp = a * i * i + koeff;

    energyDamp = 0.5*source*source*(1.0 - damp*damp);
    source *= damp;	
} 

void damping_fields(Array3D<double3>& fieldE, Array3D<double3>& fieldB,const Region& domain){
    long i,j,k,i1,j1,k1;
    double energyDamp, DampCells;
    const long cellsYP1 = 0.5 * (domain.numCells.y() - 1.2*PlasmaCellsY_glob);
    const long cellsYP2 = 0.5 * (domain.numCells.y() + 1.2*PlasmaCellsY_glob);
    const long cellsZP1 = 0.5 * (domain.numCells.z() - 1.2*PlasmaCellsZ_glob);
    const long cellsZP2 = 0.5 * (domain.numCells.z() + 1.2*PlasmaCellsZ_glob);    
    long DampCellsX[2];
    long max_indx = domain.numNodes.x()-1;
    long max_indy = domain.numNodes.y()-1;
    long max_indz = domain.numNodes.z()-1;

    DampCellsX[0] = domain.dampCells[0].x() + CELLS_SHIFT + 1;
    DampCellsX[1] = domain.dampCells[1].x() + CELLS_SHIFT + 1;
  	//std::cout << max_indy << " "<<max_indz << "\n";
    if (domain.boundType[0].x() == OPEN) { 
      	for(i = 0; i < DampCellsX[0]; i++){
	  		for (j = 0; j <= max_indy ; j++){
		  		for (k = 0; k <= max_indz ; k++){

					if( ( i >= DampCellsX[0] - domain.dampCellsForce[0].x()) && 
						j > cellsYP1 && j < cellsYP2 && k > cellsZP1 && k < cellsZP2){
					    i1 = i - (DampCellsX[0] - domain.dampCellsForce[0].x());
					    DampCells = domain.dampCellsForce[0].x();
					}
					else{ 
					    i1 = i;
					    DampCells = DampCellsX[0];
					}
					
					Damping_Func(fieldE(i,j,k).x(), i1, DampCells, energyDamp);					
				
					Damping_Func(fieldE(i,j,k).y(), i1, DampCells, energyDamp);
								
					Damping_Func(fieldE(i,j,k).z(), i1, DampCells, energyDamp);
								
					Damping_Func(fieldB(i,j,k).z(), i1, DampCells, energyDamp);

					Damping_Func(fieldB(i,j,k).y(), i1, DampCells, energyDamp);				
					
					Damping_Func(fieldB(i,j,k).x(), i1, DampCells, energyDamp);
		    	}
	    	}
		}
    }

    if (domain.boundType[1].x() == OPEN) { 
      	for(i = max_indx; i > max_indx - DampCellsX[1]; i--){
	  	    for (j = 0; j <= max_indy ; j++){
		  	    for (k = 0; k <= max_indz ; k++){

				    if((i <= max_indx - (DampCellsX[1] - domain.dampCellsForce[1].x())) && 
				    	j > cellsYP1 && j < cellsYP2 && k > cellsZP1 && k < cellsZP2 ){
					    i1 = -i + (max_indx - (DampCellsX[1] - domain.dampCellsForce[1].x()));
					    DampCells = domain.dampCellsForce[1].x();
				    } else{ 
					    i1 = -(i - max_indx);
					    DampCells = DampCellsX[1];
				    }
		    
					
					Damping_Func(fieldE(i,j,k).x(), i1, DampCells, energyDamp);							

					Damping_Func(fieldE(i,j,k).y(), i1, DampCells, energyDamp);
								
					Damping_Func(fieldE(i,j,k).z(), i1, DampCells, energyDamp);
								
					Damping_Func(fieldB(i,j,k).z(), i1, DampCells, energyDamp);

					Damping_Func(fieldB(i,j,k).y(), i1, DampCells, energyDamp);				
					
					Damping_Func(fieldB(i,j,k).x(), i1, DampCells, energyDamp);
		    	}
	    	}
		}
    }

	if (domain.boundType[0].y() == OPEN) { 

	    for( i = 0; i <= max_indx; ++i){			
			for( j = max_indy; j > max_indy - domain.dampCells[1].y(); --j){
				for( k = 0; k <= max_indz; ++k){
								
				    j1 = - j + max_indy;
				    DampCells = domain.dampCells[1].y();

				    Damping_Func(fieldE(i,j,k).x(), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k).y(), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k).z(), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldB(i,j,k).z(), j1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k).y(), j1, DampCells, energyDamp);

						
				    Damping_Func(fieldB(i,j,k).x(), j1, DampCells, energyDamp);

				}
			}
			for( j = 0; j< domain.dampCells[0].y(); ++j){
				for( k = 0; k <= max_indz; ++k){
								
				    j1 = j;//- j + domain.numCells_d2+ADD_NODES;
				    DampCells = domain.dampCells[0].y();

				    Damping_Func(fieldE(i,j,k).x(), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k).y(), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k).z(), j1, DampCells, energyDamp);
									
				    Damping_Func(fieldB(i,j,k).z(), j1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k).y(), j1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k).x(), j1, DampCells, energyDamp);
				}
			}		
	    }
  	}

	if (domain.boundType[0].z() == OPEN) { 

	    for( i = 0; i <= max_indx; ++i){			
			for( j = 0; j <= max_indy; ++j){
				for( k = max_indz; k > max_indz - domain.dampCells[1].z(); --k){
								
				    k1 = - k + max_indz;
				    DampCells = domain.dampCells[1].z();

				    Damping_Func(fieldE(i,j,k).x(), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k).y(), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k).z(), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldB(i,j,k).z(), k1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k).y(), k1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k).x(), k1, DampCells, energyDamp);
				}
			
			for( k = 0; k< domain.dampCells[0].z(); ++k){
								
				    k1 = k;//- j + domain.numCells_d2+ADD_NODES;
				    DampCells = domain.dampCells[0].z();

				    Damping_Func(fieldE(i,j,k).x(), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k).y(), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldE(i,j,k).z(), k1, DampCells, energyDamp);
									
				    Damping_Func(fieldB(i,j,k).z(), k1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k).y(), k1, DampCells, energyDamp);

				    Damping_Func(fieldB(i,j,k).x(), k1, DampCells, energyDamp);
				}
			}		
	    }
  	}
}

