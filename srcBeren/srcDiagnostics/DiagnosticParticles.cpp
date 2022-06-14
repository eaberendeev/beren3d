#include "Diagnostic.h"

void Writer::write_particles3D(long timestep){ 
  char filename[100];

    for( const auto& sp : _species){

      sprintf(filename, (".//Particles//"+sp.name+"//Diag3D//Dens3D%03ld_"+"%03ld").c_str(),_world.MPIconf.rank_line(),timestep / TimeStepDelayDiag2D);
      write_array3D(sp.densityOnGrid, sp.densityOnGrid.size().x() - ADD_NODES, sp.densityOnGrid.size().y(),sp.densityOnGrid.size().z(), filename, _world.MPIconf);
    }
  
}
void Writer::write_particles2D(long timestep){ 
  char filename[100];

    for( const auto& sp : _species){
    	auto size_x = sp.densityOnGrid.size().x() - ADD_NODES;
    	auto size_y = sp.densityOnGrid.size().y();
    	auto size_z = sp.densityOnGrid.size().z();
   		Array2D<double> densityPlaneZ(size_x,size_y);
   		Array2D<double> densityPlaneY(size_x,size_z);
	    
	    long k = size_z/2;
	    for( auto i = 0; i < size_x; i++ ){
	      for( auto j = 0; j < size_y; j++ ){
	            densityPlaneZ(i,j) = sp.densityOnGrid(i,j,k);
	      }
	    }
	   	long j = size_y/2;
	    for( auto i = 0; i < size_x; i++ ){
	        for( auto k = 0; k < size_z; k++ ){
	            densityPlaneY(i,k) = sp.densityOnGrid(i,j,k);
	        }
	    }
      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//DensPlaneZ"+"%03ld").c_str(),timestep / TimeStepDelayDiag2D);
      write_array2D(densityPlaneZ, size_x, size_y, filename, _world.MPIconf);
      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//DensPlaneY"+"%03ld").c_str(),timestep / TimeStepDelayDiag2D);
      write_array2D(densityPlaneY, size_x, size_z, filename, _world.MPIconf);

      sprintf(filename, (".//Particles//"+sp.name+"//Diag2D//Phase2D"+"%03ld").c_str(),timestep / TimeStepDelayDiag2D);
      write_array2D(sp.phaseOnGrid, sp.phaseOnGrid.size().x(), sp.phaseOnGrid.size().y(),filename, _world.MPIconf);

    }
  
}
