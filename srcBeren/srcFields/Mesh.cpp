#include "Mesh.h"
#include "World.h"
#include "Shape.h"
#include "SolverFDTD.h"
#include "SolverFDTD_PML.h"
//#include "Damping.h"
//#include "Coil.h"


Mesh::Mesh(const World& world) : Lmat(world.region.total_size()*3,world.region.total_size()*3 ),
          fieldE(world.region.numNodes,3),fieldEn(world.region.numNodes,3),fieldB(world.region.numNodes,3),
          fieldJ(world.region.numNodes,3),fieldB0(world.region.numNodes,3), _world(world){
    
    if (RECOVERY){ 
      read_from_recovery(world.MPIconf);
    } 
    else{
      set_fields();
    }
    _size1 = world.region.numNodes.x();
    _size2 = world.region.numNodes.y();
     _size3 = world.region.numNodes.z();
      _nd = 3;  //  std::vector< std::vector<std::string> > stringParams;
  //  read_params_to_string("Lasers","./LasParams.cfg",stringParams);
  //  for( const auto &params  : stringParams){
  //      lasers.emplace_back(params);
  //  }
}

void Mesh::set_fields(){
	set_uniform_fields();
  //if (BCoil[0] > 0){
  //  set_coils(fieldB,_world);
  //}
  fieldB0 = fieldB;
} 

void Mesh::set_uniform_fields(){
  auto size = fieldE.size();
  fieldE = 0.;
  //auto size_y = fieldE.size_d2();
  for( auto i = 0; i < size.x(); ++i){
    for( auto j = 0; j < size.y(); ++j){
      for( auto k = 0; k < size.z(); ++k){
        fieldB(i,j,k,0) = BUniform[0];
        fieldB(i,j,k,1) = BUniform[1];
        fieldB(i,j,k,2) = BUniform[2];
      }
    }
  }
}

double calc_energy_field(const Field3d& field){
  double potE = 0;
  long i_max = field.size().x() - ADD_NODES;
  long j_max = field.size().y() - ADD_NODES;
  long k_max = field.size().z() - ADD_NODES;
  for(auto i = 0; i < i_max; ++i){
    for(auto j = 0; j < j_max; ++j){
      for(auto k = 0; k < k_max; ++k){
        double3 v = double3(field(i,j,k,0),field(i,j,k,1),field(i,j,k,2));
        potE += dot(v,v );      
      }
    }
  }
  potE *= (0.5 * Dx * Dy *Dz);
  return potE;
}


void Mesh::reduce_current(const MPI_Topology &MPIconf){

  // int left, right;
  // int tag = 1;
  // MPI_Status status;
  // int bufSize = ADD_NODES * fieldJ.size().y()*fieldJ.size().z();
  // long sendIndx; 
  // static double3 *recvBuf = new double3[bufSize];
  
  // MPI_Allreduce ( MPI_IN_PLACE, &fieldJ.data(0), 3*fieldJ.capacity(), MPI_DOUBLE, MPI_SUM, MPIconf.comm_depth() ) ;
  
  // MPI_Barrier(MPI_COMM_WORLD);

  
  // if (MPIconf.rank_line() < MPIconf.size_line() - 1)
  //   right = MPIconf.rank_line() + 1;
  // else
  //   right = 0;
  
  // if (MPIconf.rank_line() > 0)
  //   left = MPIconf.rank_line() - 1;
  // else
  //   left = MPIconf.size_line() - 1;


  // sendIndx = fieldJ.capacity()-bufSize;
  // // fromRight to Left
  // MPI_Sendrecv(&fieldJ.data(sendIndx), 3*bufSize, MPI_DOUBLE_PRECISION, right, tag, 
  //                          recvBuf, 3*bufSize, MPI_DOUBLE_PRECISION, left, tag, 
  //                          MPIconf.comm_line(), &status);
  
  // for (long i = 0; i< bufSize; i++){
  //   fieldJ.data(i) += recvBuf[i];

  // }
  
  // MPI_Sendrecv(&fieldJ.data(0), 3*bufSize, MPI_DOUBLE_PRECISION,left , tag, 
  //                          &fieldJ.data(sendIndx), 3*bufSize, MPI_DOUBLE_PRECISION, right, tag, 
  //                          MPIconf.comm_line(), &status);

  
}
/*
void exchange_fieldsE(Array3D<double3>& fieldE, const MPI_Topology &MPIconf){

  int tag = 1;
  MPI_Status status;
  const auto size = fieldE.size();
  const long bufSize = size.y()*size.z(); 
  long sendIndx, recvIndx; 
  static double3 *recvBuf = new double3[bufSize];
  long shift = ADD_NODES;
  
  auto right = MPIconf.next_line();
  auto left = MPIconf.prev_line();

  sendIndx = (size.x()-shift)*size.y()*size.z();
  // fromRight to Left
  
  MPI_Sendrecv(&fieldE.data(sendIndx), 3*bufSize, MPI_DOUBLE_PRECISION, right, tag, 
                           recvBuf, 3*bufSize, MPI_DOUBLE_PRECISION, left, tag, 
                           MPIconf.comm_line(), &status);
  if( !MPIconf.is_first_line() ){
      for (auto i = 0; i< bufSize; i++){
        fieldE.data(i) = recvBuf[i];
      }
  }
  
  sendIndx = (shift-1)*size.y()*size.z();
  MPI_Sendrecv(&fieldE.data(sendIndx), 3*bufSize, MPI_DOUBLE_PRECISION,left , tag, 
                           recvBuf, 3*bufSize, MPI_DOUBLE_PRECISION, right, tag, 
                           MPIconf.comm_line(), &status);
  
  recvIndx = (size.x()-1)*size.y()*size.z();  
     
  if( !MPIconf.is_last_line() ){

      for (auto i = 0; i< bufSize; i++){
      fieldE.data(recvIndx+i) = recvBuf[i];
      }
  }
  
  
  MPI_Barrier(MPIconf.comm_line() );
  
}
*/

double3 Mesh::get_fieldE_in_cell(long i, long j,long k)  const{
  double3 E;
  E.x() = 0.25 * (fieldE(i,j,  k,0) + fieldE(i,j  ,k+1,0) 
                + fieldE(i,j+1,k,0) + fieldE(i,j+1,k+1,0) );
  
  E.y() = 0.25 * (fieldE(i,  j,k,1) + fieldE(i,  j,k+1,1) + 
                  fieldE(i+1,j,k,1) + fieldE(i+1,j,k+1,1) );
  
  E.z() = 0.125 * (fieldE(i,  j,  k,2) + fieldE(i,  j,  k+1,2) + 
                   fieldE(i,  j+1,k,2) + fieldE(i,  j+1,k+1,2) +
                   fieldE(i+1,j,  k,2) + fieldE(i+1,j,  k+1,2) + 
                   fieldE(i+1,j+1,k,2) + fieldE(i+1,j+1,k+1,2) );
  return E;
}
double3 Mesh::get_fieldB_in_cell(long i, long j,long k)  const{
  double3 B;
  B.x() = 0.25 * (fieldB(i,  j,k,0) + fieldB(i,  j,k+1,0) +
                  fieldB(i+1,j,k,0) + fieldB(i+1,j,k+1,0) );
  
  B.y() = 0.25 * (fieldB(i,j,  k,1) + fieldB(i,j,  k+1,1) +
                  fieldB(i,j+1,k,1) + fieldB(i,j+1,k+1,1) );
  
  B.z() = fieldB(i,j,k,2);
  return B;
}
inline double Shape2(const double& dist){
  double d = fabs(dist);
  if ( d <= 0.5 ) return (0.75 - d * d);
  else 
      if (d < 1.5)  return ((d - 1.5) * (d-1.5) * 0.5);
      else return 0.;
}

/*double3 Mesh::get_fieldE_in_pos(const double2& POS)  const{
  constexpr auto SMAX = 4; //2*SHAPE_SIZE;
  alignas(64) double sx[SMAX];
  alignas(64) double sy[SMAX];
  alignas(64) double sdx[SMAX];
  alignas(64) double sdy[SMAX];
  alignas(64) double3 E;
  double xx, yy, arg;
  double snm1, snm2, snm4;
  long xk, yk, indx, indy;
  
  xx = POS.x() / Dx;
  yy = POS.y() / Dy;
      
  xk = long(xx);
  yk = long(yy);
  
  for(auto n = 0; n < SMAX; ++n){
    arg = -xx + double(xk - CELLS_SHIFT + n);
    sx[n] = Shape2(arg);
    sdx[n] = Shape2(arg + 0.5);
    arg = -yy + double(yk - CELLS_SHIFT + n);
    sy[n] = Shape2(arg);
    sdy[n] = Shape2(arg + 0.5);
  }
    
  for( auto n = 0; n < SMAX; ++n){      
    for( auto m = 0; m < SMAX; ++m){

      snm1 = sdx[n] * sy[m];
      snm4 = sx[n] * sy[m];
      snm2 = sx[n] * sdy[m];
      indx = xk + n;
      indy = yk - CELLS_SHIFT + m;

      E += double3(snm1,snm2,snm4) * fieldE(indx,indy);

    }
  }
  
    return E;
}*/

double3 get_fieldE_in_pos(const Field3d& fieldE,const double3& r) {
    constexpr auto SMAX = 2*SHAPE_SIZE;
    double sx[SMAX], sy[SMAX],sz[SMAX];
    double sdx[SMAX], sdy[SMAX], sdz[SMAX];
    double3 ep;
    double xx, yy, zz,arg;
    double snm1, snm2, snm3;
    long xk, yk, zk, indx, indy, indz;
    ep = 0.;
    
    xx = r.x() / Dx;
    yy = r.y() / Dy;
    zz = r.z() / Dz;
        
    xk = long(xx);
    yk = long(yy);
    zk = long(zz);
  
    for(auto n = 0; n < SMAX; ++n){
      arg = -xx + double(xk - CELLS_SHIFT + n);
      sx[n] = Shape(arg);
      sdx[n] = Shape(arg + 0.5);
      arg = -yy + double(yk - CELLS_SHIFT + n);
      sy[n] = Shape(arg);
      sdy[n] = Shape(arg + 0.5);
      arg = -zz + double(zk - CELLS_SHIFT + n);
      sz[n] = Shape(arg);
      sdz[n] = Shape(arg + 0.5);
    }
    
    for( auto n = 0; n < SMAX; ++n){      
      for( auto m = 0; m < SMAX; ++m){
        for( auto k = 0; k < SMAX; ++k){

          snm1 = sdx[n] * sy[m]  * sz[k];
          snm2 = sx[n]  * sdy[m] * sz[k];
          snm3 = sx[n]  * sy[m]  * sdz[k];

          indx = xk + n;
          indy = yk - CELLS_SHIFT + m;
          indz = zk - CELLS_SHIFT + k; 
            
          ep.x() += snm1 * fieldE(indx,indy,indz,0);
          ep.y() += snm2 * fieldE(indx,indy,indz,1);
          ep.z() += snm3 * fieldE(indx,indy,indz,2);

        }
      }
    }
    return ep;
}

double3 get_fieldE_in_pos_lin(const Field3d& fieldE,const double3& r) {

    double3 E;
    long indx, indy,indz,indx1, indy1,indz1;

    double xx,yy,zz;

    double sx0,sy0,sz0,sdx0,sdy0,sdz0;
    double sx1,sy1,sz1,sdx1,sdy1,sdz1;

    const double rdx = 1. / Dx;
    const double rdy = 1. / Dy;
    const double rdz = 1. / Dz;    
    xx = r.x() * rdx;
    yy = r.y() * rdy;
    zz = r.z() * rdz;
    
    indx = long(xx + 1.);
    indy = long(yy + 1.);
    indz = long(zz + 1.);


    indx1 = long(xx + 0.5);
    indy1 = long(yy + 0.5);
    indz1 = long(zz + 0.5); 
    
    sx1 = (xx - indx + 1.);
    sy1 = (yy - indy + 1.);
    sz1 = (zz - indz + 1.);
    sdx1 = (xx - indx1 + 0.5);
    sdy1 = (yy - indy1 + 0.5);
    sdz1 = (zz - indz1 + 0.5);


    sx0 = 1. - sx1;
    sy0 = 1. - sy1;
    sz0 = 1. - sz1;
    sdx0 = 1. - sdx1;
    sdy0 = 1. - sdy1;
    sdz0 = 1. - sdz1;

    E.x() += sdx0 * ( sy0 * ( sz0 * fieldE(indx1,indy,indz,0) + sz1 * fieldE(indx1,indy,indz+1,0) ) 
                + sy1 * ( sz0 * fieldE(indx1,indy+1,indz,0) + sz1 * fieldE(indx1,indy+1,indz+1,0) ) ) 
        + sdx1 * ( sy0 * ( sz0 * fieldE(indx1+1,indy,indz,0) + sz1 * fieldE(indx1+1,indy,indz+1,0) ) 
                + sy1 * ( sz0 * fieldE(indx1+1,indy+1,indz,0) + sz1 * fieldE(indx1+1,indy+1,indz+1,0) ) );

    E.y() += sx0 * ( sdy0 * ( sz0 * fieldE(indx,indy1,indz,1) + sz1 * fieldE(indx,indy1,indz+1,1) ) 
                + sdy1 * ( sz0 * fieldE(indx,indy1+1,indz,1) + sz1 * fieldE(indx,indy1+1,indz+1,1) ) ) 
        + sx1 * ( sdy0 * ( sz0 * fieldE(indx+1,indy1,indz,1) + sz1 * fieldE(indx+1,indy1,indz+1,1) ) 
                + sdy1 * ( sz0 * fieldE(indx+1,indy1+1,indz,1) + sz1 * fieldE(indx+1,indy1+1,indz+1,1) ) );

    E.z() += sx0 * ( sy0 * ( sdz0 * fieldE(indx,indy,indz1,2) + sdz1 * fieldE(indx,indy,indz1+1,2) ) 
                + sy1 * ( sdz0 * fieldE(indx,indy+1,indz1,2) + sdz1 * fieldE(indx,indy+1,indz1+1,2) ) ) 
        + sx1 * ( sy0 * ( sdz0 * fieldE(indx+1,indy,indz1,2) + sdz1 * fieldE(indx+1,indy,indz1+1,2) ) 
                + sy1 * ( sdz0 * fieldE(indx+1,indy+1,indz1,2) + sdz1 * fieldE(indx+1,indy+1,indz1+1,2) ) );

    return E;
}
double3 get_fieldB_in_pos(const Field3d& fieldB, const double3& r) {

    double3 B;
    long indx, indy,indz,indx1, indy1,indz1;

    double xx,yy,zz;

    double sx0,sy0,sz0,sdx0,sdy0,sdz0;
    double sx1,sy1,sz1,sdx1,sdy1,sdz1;

    const double rdx = 1. / Dx;
    const double rdy = 1. / Dy;
    const double rdz = 1. / Dz;    
    xx = r.x() * rdx;
    yy = r.y() * rdy;
    zz = r.z() * rdz;
    
    indx = long(xx + 1.);
    indy = long(yy + 1.);
    indz = long(zz + 1.);


    indx1 = long(xx + 0.5);
    indy1 = long(yy + 0.5);
    indz1 = long(zz + 0.5); 
    
    sx1 = (xx - indx + 1.);
    sy1 = (yy - indy + 1.);
    sz1 = (zz - indz + 1.);
    sdx1 = (xx - indx1 + 0.5);
    sdy1 = (yy - indy1 + 0.5);
    sdz1 = (zz - indz1 + 0.5);


    sx0 = 1. - sx1;
    sy0 = 1. - sy1;
    sz0 = 1. - sz1;
    sdx0 = 1. - sdx1;
    sdy0 = 1. - sdy1;
    sdz0 = 1. - sdz1;

  B.x() += sx0 * ( sdy0 * ( sdz0 * fieldB(indx,indy1,indz1,0) + sdz1 * fieldB(indx,indy1,indz1+1,0) ) 
               + sdy1 * ( sdz0 * fieldB(indx,indy1+1,indz1,0) + sdz1 * fieldB(indx,indy1+1,indz1+1,0) ) ) 
       + sx1 * ( sdy0 * ( sdz0 * fieldB(indx+1,indy1,indz1,0) + sdz1 * fieldB(indx+1,indy1,indz1+1,0) ) 
               + sdy1 * ( sdz0 * fieldB(indx+1,indy1+1,indz1,0) + sdz1 * fieldB(indx+1,indy1+1,indz1+1,0) ) );

  B.y() += sdx0 * ( sy0 * ( sdz0 * fieldB(indx1,indy,indz1,1) + sdz1 * fieldB(indx1,indy,indz1+1,1) ) 
              + sy1 * ( sdz0 * fieldB(indx1,indy+1,indz1,1) + sdz1 * fieldB(indx1,indy+1,indz1+1,1) ) ) 
      + sdx1 * ( sy0 * ( sdz0 * fieldB(indx1+1,indy,indz1,1) + sdz1 * fieldB(indx1+1,indy,indz1+1,1) ) 
              + sy1 * ( sdz0 * fieldB(indx1+1,indy+1,indz1,1) + sdz1 * fieldB(indx1+1,indy+1,indz1+1,1) ) );

  B.z() += sdx0 * ( sdy0 * ( sz0 * fieldB(indx1,indy1,indz,2) + sz1 * fieldB(indx1,indy1,indz+1,2) ) 
              + sdy1 * ( sz0 * fieldB(indx1,indy1+1,indz,2) + sz1 * fieldB(indx1,indy1+1,indz+1,2) ) ) 
      + sdx1 * ( sdy0 * ( sz0 * fieldB(indx1+1,indy1,indz,2) + sz1 * fieldB(indx1+1,indy1,indz+1,2) ) 
              + sdy1 * ( sz0 * fieldB(indx1+1,indy1+1,indz,2) + sz1 * fieldB(indx1+1,indy1+1,indz+1,2) ) );

    return B;
}


void Mesh::update(long timestep){
    
    fieldB -= fieldB0;
    //laser_source(timestep);
    
    #if DAMP_FIELDS == DAMP
      solver_FDTD(fieldE, fieldB, fieldJ, _world);
    //  damping_fields(fieldE, fieldB,_world.region);

    #elif DAMP_FIELDS == PML
     // solver_FDTD_PML(fieldE, fieldB,fieldEp, fieldBp, fieldJ, _world);
    #endif
    fieldB += fieldB0;


}


void Mesh::laser_source(long timestep){
      /// Laser
   // const auto l_min = _world.region.dampCells[0].y();
   // const auto l_max = _world.region.numCells.y() - _world.region.dampCells[1].y() ;
    /*
    for ( const auto& las : lasers){

      if ( !_world.region.in_region(las.start_d1) ) continue;
          	                 // std::cout <<"dfvd";

      if ( !las.is_work(timestep) )  continue;

      const long indx = round( (las.start_d1 - _world.region.origin) / Dx) + CELLS_SHIFT + 1;
      
      for (auto l = l_min; l <= l_max ; ++l){
        double y = (l+0.5)*Dy;

        if(las.type=="Ey"){
            fieldE(indx,l).y() = las.get_Ey(las.start_d1, y, timestep);
            fieldB(indx,l).z() = sign(las.vg)*las.get_Ey(las.start_d1 + 0.5*Dx,y-0.5*Dy, timestep);
        } 
        else if(las.type=="Ez"){
            fieldE(indx,l).z() = las.get_Ez(y-0.5*Dy, timestep*Dt-las.delay);
            fieldB(indx,l).y() = -sign(las.vg)*las.get_Ez(y-0.5*Dy, Dt*timestep-las.delay);  
        }
      }

    
    } */
  
}

 void Mesh::get_M(){ // !!!!! needs bound condition and if cases!!!!!!
    std::vector<Trip> trips;
    long totalSize = fieldE.capacity();
    trips.reserve(totalSize*3);
    double val;

    const auto size = fieldE.size();
    
    for(long i = 0; i < size.x(); i++){
      for(long j = 0; j < size.x(); j++){
        for(long k = 0; k < size.x(); k++){

          /////  MX  /////
          // Mxx
          val = 1. - 0.25*Dt*Dt*(-2/(Dy*Dy*)-2/(Dz*Dz));
          // Ex i+1/2,j,k
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i,j,k,0),val));
          val = 0.25*Dt*Dt*(Dy*Dy);
          // Ex i+1/2,j+1,k
          if (j != size.y() )
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i,j+1,k,0),-val));
          // Ex i+1/2,j-1,k          
          if (j!=0)
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i,j-1,k,0),-val));
          val = 0.25*Dt*Dt*(Dz*Dz);
          // Ex i+1/2,j,k+1
          if (k != size.z() )
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i,j,k+1,0),-val));
          // Ex i+1/2,j,k-1          
          if (k!=0)
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i,j,k-1,0),-val));

          //Mxy
          val = 0.25*Dt*Dt*(Dx*Dy);
          // Ey i,j+1/2,k
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i,j,k,1),val));
          // Ey i+1,j+1/2,k
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i+1,j,k,1),-val));
          // Ey i,j-1/2,k
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i,j-1,k,1),-val));
          // Ey i+1,j-1/2,k
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i+1,j-1,k,1),val));
          //Mxz
          val = 0.25*Dt*Dt*(Dx*Dz);
          // Ez i,j,k+1/2
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i,j,k,2),val));
          // Ez i+1,j,k+1/2
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i+1,j,k,2),-val));
          // Ez i,j,k-1/2
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i,j,k-1,2),-val));
          // Ez i+1,j,k-1/2
          trips.push_back(Trip(mesh.vind(i,j,k,0),mesh.vind(i+1,j,k-1,2),val));

          /////  MY  /////
          // Myy
          val = 1. - 0.25*Dt*Dt*(-2/(Dx*Dx*)-2/(Dz*Dz));
          // Ey i,j+1/2,k
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i,j,k,1),val));
          val = 0.25*Dt*Dt*(Dx*Dx);
          // Ey i+1,j+1/2,k
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i+1,j,k,1),-val));
          // Ey i-1,j+1/2,k
          val = 0.25*Dt*Dt*(Dz*Dz);
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i-1,j,k,1),-val));
          // Ey i,j+1/2,k+1
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i,j,k+1,1),-val));
          // Ey i,j+1/2,k-1
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i,j,k-1,1),-val));


          // Myx
          val = 0.25*Dt*Dt*(Dy*Dx);
          // Ex i+1/2,j,k
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i,j,k,0),val));     
          // Ex i+1/2,j+1,k     
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i,j+1,k,0),-val));          
          // Ex i-1/2,j,k
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i-1,j,k,0),-val));          
          // Ex i-1/2,j+1,k
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i-1,j+1,k,0),val));          
          
          // Myz
          val = 0.25*Dt*Dt*(Dy*Dz);
          // Ez i,j,k+1/2
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i,j,k,2),val));     
          // Ez i,j+1,k+1/2
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i,j+1,k,2),-val));     
          // Ez i,j,k-1/2
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i,j,k-1,2),-val));               
          // Ez i,j+1,k-1/2
          trips.push_back(Trip(mesh.vind(i,j,k,1),mesh.vind(i,j+1,k-1,2),val));     


          /////  MZ  /////
          // Mzz
          val = 1. - 0.25*Dt*Dt*(-2/(Dx*Dx*)-2/(Dy*Dy));
          // Ez i,j,k+1/2
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i,j,k,2),val));
          // Ez i+1,j,k+1/2
          val = 0.25*Dt*Dt*(Dx*Dx);
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i+1,j,k,2),-val));
          // Ez i-1,j,k+1/2
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i-1,j,k,2),-val));
          val = 0.25*Dt*Dt*(Dy*Dy);
          // Ez i,j+1,k+1/2
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i,j+1,k,2),-val));
          // Ez 1,j-1,k+1/2
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i,j-1,k,2),-val));
          
          // Mzx
          val = 0.25*Dt*Dt*(Dz*Dx);
          // Ex i+1/2,j,k
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i,j,k,0),val));
          // Ex i-1/2,j,k
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i-1,j,k,0),-val));
          // Ex i+1/2,j,k+1
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i,j,k+1,0),-val));
          // Ex i-1/2,j,k+1
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i-1,j,k+1,0),val));

          // Mzy          
          val = 0.25*Dt*Dt*(Dz*Dy);
          // Ey i,j+1/2,k
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i,j,k,1),val));
          // Ey i,j-1/2,k
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i,j-1,k,1),-val));
          // Ey i,j+1/2,k+1
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i,j,k+1,1),-val));
          // Ey i,j-1/2,k+1
          trips.push_back(Trip(mesh.vind(i,j,k,2),mesh.vind(i,j-1,k+1,1),val));

        }
      }
    }
        mesh.Lmat.setFromTriplets(trips.begin(), trips.end());

 }