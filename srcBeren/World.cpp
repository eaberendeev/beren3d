#include "World.h"

#define RANDOMSTD

#ifdef RANDOMSTD
  #include <random>

  std::mt19937 gen;
  std::uniform_real_distribution<> urd(0, 1); 

  double Uniform01(){
    return urd(gen);
  }

  void SetRandSeed(int val){
    gen.seed(val);
  }
#else
  double Uniform01(){
    return (double)(rand())/RAND_MAX ;
  }
  void SetRandSeed(int val){
      srand(val);
  }
#endif

double Gauss(double sigma){
  double r1 = Uniform01();
  double r2 = Uniform01();
  
  return sigma*sqrt(-2.0*log(r1))*sin(2.0*PI*r2);
}

void MPI_Topology::set_topology(){
  
    MPI_Comm_size(MPI_COMM_WORLD,&_size);
    MPI_Comm_rank(MPI_COMM_WORLD,&_rank);
    
    _sizeDepth = _size / NumAreas;
    _sizeLine = NumAreas;
    
    int colorLine = _rank / NumAreas;
    MPI_Comm_split(MPI_COMM_WORLD, colorLine, _rank, &_commLine);

    int colorDepth = _rank % NumAreas;
    MPI_Comm_split(MPI_COMM_WORLD, colorDepth, _rank, &_commDepth);

    MPI_Comm_rank(_commDepth, &_rankDepth);
    MPI_Comm_rank(_commLine, &_rankLine);
    SetRandSeed(_rank*3+1);

    if( is_master() )
        std::cout << "There are " << NumAreas << " subdomains" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);
}


// Читаем данные из файла параметров, распознаём строки и записываем данные в Params
Region::Region(){

    step_h = double3(Dx,Dy,Dz);
    origin = 0.0;
   
    numCells = long3(NumCellsX_glob, NumCellsY_glob, NumCellsZ_glob);    
    numNodes = long3(numCells.x() + ADD_NODES,
                     numCells.y() + ADD_NODES,
                     numCells.z() + ADD_NODES);
    dampCells[0] = long3(DampCellsX_glob[0],DampCellsY_glob[0],DampCellsZ_glob[0]);
    dampCells[1] = long3(DampCellsX_glob[1],DampCellsY_glob[1],DampCellsZ_glob[1]);

    dampCellsForce[0] = long3(DampCellsXP_glob[0],0.,0.);
    dampCellsForce[1] = long3(DampCellsXP_glob[1],0.,0.);

    boundType[0] = long3(BoundTypeX_glob[0],BoundTypeY_glob[0],BoundTypeZ_glob[0]);
    boundType[1] = long3(BoundTypeX_glob[1],BoundTypeY_glob[1],BoundTypeZ_glob[1]);

}


Region split_region(const Region& regionGlob, int rank, int splitSize){
    Region region = regionGlob;
    
    region.numCells.x() = PlasmaCellsX_glob / splitSize;
    region.origin = region.step_h.x() * rank * region.numCells.x();
    
    if (rank != 0) {
        region.boundType[0].x() = NEIGHBOUR;
        region.dampCells[0].x() = 0;
        region.dampCellsForce[0].x() = 0;
        region.origin += region.step_h.x() * DampCellsX_glob[0];
    } else{
        region.numCells.x() += DampCellsX_glob[0];
    }
    if (rank != splitSize - 1) {
        region.boundType[1].x() = NEIGHBOUR;
        region.dampCells[1].x() = 0;
        region.dampCellsForce[1].x() = 0;

    } else{
        region.numCells.x() += DampCellsX_glob[1];
    }

    region.numNodes.x() = region.numCells.x() + ADD_NODES;

    std::cout << "Current number of cells = "<< region.numCells.x() << ". Start coord =" << region.origin << std::endl;
    return region;
}

long MPI_Topology::accum_sum_line(long value) const{
    long sum = 0;
    int prev, next;
    int id = 1;
    MPI_Status status;
    if (NumAreas == 1) return 0;
    prev = _rankLine - 1;
    next = _rankLine  + 1;
    if( _rankLine == 0){
        prev = MPI_PROC_NULL; 

    }
    else if(_rankLine == _sizeLine - 1){
        next = MPI_PROC_NULL;
    }
    
    MPI_Recv(&sum, 1, MPI_LONG, prev, id, _commLine, &status);
    
    if( _rankLine != 0){
        value += sum; 
    }
    MPI_Send(&value, 1, MPI_LONG, next, id,  _commLine);
    return sum;

}


/*

void World::filter_fields(long timestep){
  long i,j,t;
  long sizeT = long(2 * PI / Dt);
  t = timestep % sizeT;

  if( fieldsAvgLine_y1 > 0){
    
    j = fieldsAvgLine_y1 / Dy;  
    for (i = 0; i < fieldE.size_d1() - ADD_NODES; ++i){
      fieldEAvgLineTop(i,t) = get_fieldE_in_cell(i,j);
      fieldBAvgLineTop(i,t) = get_fieldB_in_cell(i,j);
    }
  }
    if( fieldsAvgLine_y2 > 0){
    
    j = fieldsAvgLine_y2 / Dy;  
    for (i = 0; i < fieldE.size_d1() - ADD_NODES; ++i){
      fieldEAvgLineBottom(i,t) = get_fieldE_in_cell(i,j);
      fieldBAvgLineBottom(i,t) = get_fieldB_in_cell(i,j);
    }
  }
  if( fieldsAvgLine_x1 > 0){
    
    i = (fieldsAvgLine_x1 - region.origin) / Dx;  
      if (!region.in_region(fieldsAvgLine_x1)) return;
    for (j = 0; j < fieldE.size_d2() - ADD_NODES; ++j){
      fieldEAvgLineLeft(j,t) = get_fieldE_in_cell(i,j);
      fieldBAvgLineLeft(j,t) = get_fieldB_in_cell(i,j);
    }
  }

  if( fieldsAvgLine_x2 > 0){
    
  i = (fieldsAvgLine_x2 - region.origin) / Dx;
  if (!region.in_region(fieldsAvgLine_x2)) return;
    for (j = 0; j < fieldE.size_d2() - ADD_NODES; ++j){
      fieldEAvgLineRight(j,t) = get_fieldE_in_cell(i,j);
      fieldBAvgLineRight(j,t) = get_fieldB_in_cell(i,j);
    }
  }

}
double3 World::get_fieldE_in_cell(long i, long j)  const{
  double3 f;
  f.x() = 0.5 * (fieldE(i,j).x() + fieldE(i,j+1).x());
  f.y() = 0.5 * (fieldE(i,j).y() + fieldE(i+1,j).y());
  f.z() = 0.25 * (fieldE(i,j).z() + fieldE(i+1,j).z() + fieldE(i,j+1).z() + fieldE(i+1,j+1).z() );
  return f;
}
double3 World::get_fieldB_in_cell(long i, long j)  const{
  double3 f;
  f.x() = 0.5 *(fieldB(i,j).x() + fieldB(i+1,j).x());
  f.y() = 0.5 * (fieldB(i,j).y() + fieldB(i,j+1).y());
  f.z() = fieldB(i,j).z();
  return f;
}

double3 World::get_fieldE_in_pos(double2 r)  const{
  constexpr auto SMAX = 2*SHAPE_SIZE;
  double sx[SMAX];
  double sy[SMAX];
  double sdx[SMAX];
  double sdy[SMAX];
  double3 ep;
  double xx, yy, arg;
  double snm1, snm2, snm4;
  long xk, yk, indx, indy;
  ep = 0.;
  
  xx = r.x() / Dx;
  yy = r.y() / Dy;
      
  xk = long(xx);
  yk = long(yy);
  
  for(auto n = 0; n < SMAX; ++n){
    arg = -xx + double(xk - CELLS_SHIFT + n);
    sx[n] = Shape(arg);
    sdx[n] = Shape(arg + 0.5);
    arg = -yy + double(yk - CELLS_SHIFT + n);
    sy[n] = Shape(arg);
    sdy[n] = Shape(arg + 0.5);
  }
    
    for( auto n = 0; n < SMAX; ++n){      
	  for( auto m = 0; m < SMAX; ++m){

      
      snm1 = sdx[n] * sy[m];
      snm4 = sx[n] * sy[m];
      snm2 = sx[n] * sdy[m];
      indx = xk + n;
      indy = yk - CELLS_SHIFT + m;
        
        ep += double3(snm1,snm2,snm4) * fieldE(indx,indy);

    }
  }
  
    return ep;
}

void World::laser_source(long timestep){
      /// Laser
  	const auto l_min = region.dampCells_d2;
  	const auto l_max = region.numCells_d2 - region.dampCells_d2 ;
 	for ( const auto& las : lasers){
	    if ( !region.in_region(las.start_d1) ) continue;
	    if ( !las.is_work(timestep) )  continue;
	  	
	  	const long indx = round( (las.start_d1 - region.origin) / Dx) + CELLS_SHIFT + 1;
	    
	    for (auto l = l_min; l <= l_max ; ++l){
	    	double y = (l+0.5)*Dy;
		  		// Надо написать с учётом декомпозиции 
        if(las.type=="Ey"){
    		    fieldE(indx,l).y() = las.get_Ey(las.start_d1, y, timestep);
    		    fieldB(indx,l).z() = sign(las.vg)*las.get_Ey(las.start_d1 + 0.5*Dx,y-0.5*Dy, timestep);
		    } 
        else if(las.type=="Ez"){
            fieldE(indx,l).z() = las.get_Ez(y-0.5*Dy, timestep*Dt-las.delay);
    		    fieldB(indx,l).y() = -sign(las.vg)*las.get_Ez(y-0.5*Dy, Dt*timestep-las.delay);  
	  	  } else{ std::cout << "Critical error! Invalit LASER type!\n"; exit(1);}
      }

    
  	} 
  
}*/
