#include "SolverFDTD_PML.h"
#include "World.h"
#include "Mesh.h"


inline double pow3(double a){
  return a*a*a;
}

const long pml_x = DampCellsX_glob[0] + CELLS_SHIFT + 1;
const long pml_y = DampCellsY_glob[0] + CELLS_SHIFT + 1;
const long pml_z = DampCellsZ_glob[0] + CELLS_SHIFT + 1;
const double smxmax = -2. * log(1.e-6) / (pml_x * Dx);
const double smymax = -2. * log(1.e-6) / (pml_y * Dy);
const double smzmax = -2. * log(1.e-6) / (pml_z * Dz);

inline double shy(long i){
  if(i >= pml_y) return 0.;  
  
  return double(pml_y - i) / pml_y;
}
inline double shz(long i){
  if(i >= pml_z) return 0.;  
  
  return double(pml_z - i) / pml_z;
}
inline double shx(long i){
  if(i >= pml_x) return 0.; 
  
  return double(pml_x - i) / pml_x;
}

inline bool in_PML_x_left(long i, const Region& region){
	return i <= pml_x && region.boundType[0].x() == OPEN;
}
inline bool in_PML_x_right(long i,const Region& region){
	return i <= pml_x && region.boundType[1].x() == OPEN;
}

inline bool in_PML_y(long i){
	return i <= pml_y;
}
inline bool in_PML_z(long i){
	return i <= pml_z;
}
bool in_PML(long i, long size_x, long j, long size_y,long k, long size_z,  const Region& region){
	return in_PML_x_left(i,region) || in_PML_x_right(size_x-i,region) 
	       || in_PML_y(j) || in_PML_y(size_y-j)
	       || in_PML_y(k) || in_PML_z(size_z-k);
}

double sgm_x(long i, long size_x,const Region& region){
	if( in_PML_x_left(i,region) )
		return pow3( shx(i) ) * smxmax;
	i = size_x - i;
	if(in_PML_x_right(i,region) )
		return pow3( shx(i) ) * smxmax;
	return 0;	
}

double sgm_y(long i, long size_y){
	if( in_PML_y(i) )
		return pow3(shy(i))*smymax;
	i = size_y - i;
	if( in_PML_y(i) )
		return pow3(shy(i))*smymax;
	return 0;	
}
double sgm_z(long i, long size_z){
	if( in_PML_z(i) )
		return pow3(shz(i))*smzmax;
	i = size_z - i;
	if( in_PML_z(i) )
		return pow3(shz(i))*smzmax;
	return 0;	
}


/*
static void update_fieldsB(const Array3D<double3>& fieldE, Array3D<double3>& fieldB, Array3D<double3>& fieldBz, const Region& region){
    long i, j;
    const double dtp = 0.5 * Dt;
    const double rdx = 1. / Dx;
    const double rdy = 1. / Dy;
    const double rdz = 1. / Dz;
    const long size_d1 = fieldB.size_d1();
    const long size_d2 = fieldB.size_d2();
    const long size_d3 = fieldB.size_d3();
    double smx, smy,smz;

    for( i = 0; i < size_x - 1; ++i){
		for( j = 0; j < size_y - 1; ++j){
            for( k = 0; k < size_d3 - 1; ++k ){

		    if(! in_PML(i,size_x,j,size_y,region) ){
    			fieldB(i,j,k).x() -= dtp * (  rdy * (fieldE(i,j+1,k).z() - fieldE(i,j,k).z() ) 
                                              - rdz * (fieldE(i,j,k+1).y() - fieldE(i,j,k).y() )  );
    			fieldB(i,j,k).y() -= dtp * ( rdz * (fieldE(i,j,k+1).x() - fieldE(i,j,k).x() )
                                           - rdx * (fieldE(i+1,j,k).z() - fieldE(i,j,k).z() ) );
    			fieldB(i,j,k).z() -= dtp * ( rdx * (fieldE(i+1,j,k).y() - fieldE(i,j,k).y() )
                                           - rdy * (fieldE(i,j+1,k).x() - fieldE(i,j,k).x() ) );
		    }
		    else{
		    	smx = sgm_x(i,size_x,region);
				smy = sgm_y(j,size_y);
				smz = sgm_z(k,size_z);

				fieldB(i,j).x() = (1.0-dtp*smy)*fieldB(i,j).x() - dtp * rdy * (fieldE(i,j+1).z() - fieldE(i,j).z() );
				fieldB(i,j).y() = (1.0-dtp*smx)*fieldB(i,j).y() + dtp * rdx * (fieldE(i+1,j).z() - fieldE(i,j).z() );		
				
				fieldBz(i,j).x() = (1.0-dtp*smx)*fieldBz(i,j).x() - dtp * rdx * (fieldE(i+1,j).y() - fieldE(i,j).y() );
				fieldBz(i,j).y() = (1.0-dtp*smy)*fieldBz(i,j).y() + dtp * rdy * (fieldE(i,j+1).x() - fieldE(i,j).x() );
				fieldB(i,j).z() = fieldBz(i,j).x() + fieldBz(i,j).y();
		    }
		}
    }

    i = size_x - 1;
    for(j = 0; j < size_y - 1; ++j){

		if(!in_PML(i,size_x,j,size_y,region))
			fieldB(i,j).x() += ( - dtp * rdy * (fieldE(i,j+1).z() - fieldE(i,j).z() ) );
		else{
		    smy =  sgm_y(j,size_y);
		    fieldB(i,j).x() = (1.0-dtp*smy)*fieldB(i,j).x() - dtp * rdy * (fieldE(i,j+1).z() - fieldE(i,j).z() );
		}
    }

	j = size_y - 1;
    for(i = 0; i < size_x - 1; ++i){

		if(!in_PML(i,size_x,j,size_y,region))
			fieldB(i,j).y() += dtp * rdx * (fieldE(i+1,j).z() - fieldE(i,j).z() );
		else{
		    smx = sgm_x(i,size_x,region);
		    fieldB(i,j).y() = (1.0-dtp*smx)*fieldB(i,j).y() + dtp * rdx * (fieldE(i+1,j).z() - fieldE(i,j).z() );
		}
    }
}

*/
void solver_FDTD_PML(Array3D<double3>& fieldE, Array3D<double3>& fieldB,
	Array3D<double3>& fieldEz, Array3D<double3>& fieldBz, const Array3D<double3>& fieldJ, const World& world){
/*    long i, j;
    const double rdx = 1. / Dx;
    const double rdy = 1. / Dy;
    const long size_x = fieldE.size_d1();
    const long size_y = fieldE.size_d2();
    double smx, smy;

    static Array2D<double> Ez0(2,size_y);
    static Array2D<double> Ey0(2,size_y);

    for (j = 0; j < size_y; ++j){
		Ez0(0,j) = fieldE(1,j).z();
		Ez0(1,j) = fieldE(size_x-2,j).z();	
		Ey0(0,j) = fieldE(1,j).y();
		Ey0(1,j) = fieldE(size_x-2,j).y();      
    }

    update_fieldsB(fieldE, fieldB, fieldBz, world.region);
	
    i = 0;
    for(j = 1; j < size_y; ++j){


		if(!in_PML(i,size_x,j,size_y,world.region))
		    fieldE(i,j).x() += ( Dt * rdy * (fieldB(i,j).z() - fieldB(i,j-1).z() ) - Dt * fieldJ(i,j).x());
		else{
				smy = sgm_y(j,size_y);
			    fieldE(i,j).x() = (1.0-Dt*smy)*fieldE(i,j).x() + Dt * rdy * (fieldB(i,j).z() - fieldB(i,j-1).z() ) - Dt*fieldJ(i,j).x();
		}
    }
	
    for(i = 1; i < size_x-1; ++i){
		for(j = 1; j < size_y; ++j){

			if(!in_PML(i,size_x,j,size_y,world.region)){
			    fieldE(i,j).x() += ( Dt * rdy * (fieldB(i,j).z() - fieldB(i,j-1).z() ) - Dt * fieldJ(i,j).x() );
			    fieldE(i,j).y() += (- Dt * rdx * (fieldB(i,j).z() - fieldB(i-1,j).z() ) - Dt * fieldJ(i,j).y() );
			    fieldE(i,j).z() += Dt*( rdx * (fieldB(i,j).y() - fieldB(i-1,j).y() ) - rdy * (fieldB(i,j).x() - fieldB(i,j-1).x() ) - fieldJ(i,j).z() );
			} else{
		    	smx = sgm_x(i,size_x,world.region);
				smy = sgm_y(j,size_y);
			    fieldE(i,j).x() = (1.0-Dt*smy)*fieldE(i,j).x() + Dt * rdy * (fieldB(i,j).z() - fieldB(i,j-1).z() ) - Dt * fieldJ(i,j).x();
			    fieldE(i,j).y() = (1.0-Dt*smx)*fieldE(i,j).y() - Dt * rdx * (fieldB(i,j).z() - fieldB(i-1,j).z() ) - Dt * fieldJ(i,j).y(); 
			    
			    fieldEz(i,j).y() = (1.0-Dt*smy)*fieldEz(i,j).y() - Dt * rdy * (fieldB(i,j).x() - fieldB(i,j-1).x() ) - 0.5 * Dt * fieldJ(i,j).z();
			    fieldEz(i,j).x() = (1.0-Dt*smx)*fieldEz(i,j).x() + Dt * rdx * (fieldB(i,j).y() - fieldB(i-1,j).y() ) - 0.5 * Dt * fieldJ(i,j).z();
			    fieldE(i,j).z() = fieldEz(i,j).x() + fieldEz(i,j).y();
			}
		}
    }
    exchange_fieldsE(fieldE,world.MPIconf);
    exchange_fieldsE(fieldEz,world.MPIconf);
	
    if(world.region.boundType_d1[0] == OPEN){
		for(j = 0; j < size_y; ++j){
		    double Kabc = (Dt - Dx) / (Dt + Dx);
		    fieldE(0,j).z() = Ez0(0,j) + Kabc * (fieldE(1,j).z() - fieldE(0,j).z() );
		    fieldE(0,j).y() = Ey0(0,j) + Kabc * (fieldE(1,j).y() - fieldE(0,j).y() );
	  	}
    }
    if(world.region.boundType_d1[1] == OPEN){
	  for(j = 0; j < size_y; ++j){
	    double Kabc = (Dt - Dx) / (Dt + Dx);
	    fieldE(size_x-1,j).z() = Ez0(1,j) + Kabc * (fieldE(size_x-2,j).z() - fieldE(size_x - 1,j).z() );
	    fieldE(size_x-1,j).y() = Ey0(1,j) + Kabc * (fieldE(size_x-2,j).y() - fieldE(size_x - 1,j).y() );
	  }

    }

    update_fieldsB(fieldE, fieldB, fieldBz, world.region);
		*/	
} 
