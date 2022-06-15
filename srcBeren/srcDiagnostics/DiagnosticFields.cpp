#include "Diagnostic.h"
#include <vector>

void Writer::diag_zond(long timestep){
  if (!_world.MPIconf.is_master_depth()) return;
  
  char filename[100];
  static FILE *fZond; 
  ulong n;
  // long i,j,k;
  double xx,yy,zz;  
  double3 E,B,r;


  if (timestep == StartTimeStep){
    sprintf(filename, "./Fields/Zond%03d.dat",_world.MPIconf.rank_line() );
    fZond = fopen(filename, "w");
    fprintf(fZond, "%s ", "## timestep ");
    for (n = 0; n < diagData.params.zondCoords.size(); ++n){
      
      if( ! _world.region.in_region(diagData.params.zondCoords[n].x() ) ) continue;
      xx = diagData.params.zondCoords[n].x();
      yy = diagData.params.zondCoords[n].y();
      zz = diagData.params.zondCoords[n].z();
      fprintf(fZond, "%s%g%s%g%s%g%s %s%g%s%g%s%g%s %s%g%s%g%s%g%s %s%g%s%g%s%g%s %s%g%s%g%s%g%s %s%g%s%g%s%g%s ", 
          "Ex(",xx,",",yy,",",zz,")","Ey(",xx,",",yy,",",zz,")","Ez(",xx,",",yy,",",zz,")",
          "Bx(",xx,",",yy,",",zz,")","By(",xx,",",yy,",",zz,")","Bz(",xx,",",yy,",",zz,")"  );
    }
    fprintf(fZond, "\n");
  }
  
    fprintf(fZond, "%g ",Dt*timestep);
    
    for (n = 0; n < diagData.params.zondCoords.size(); ++n){
      if( ! _world.region.in_region( diagData.params.zondCoords[n].x() ) ) continue;
      r.x() = diagData.params.zondCoords[n].x() -  _world.region.origin;
      r.y() = diagData.params.zondCoords[n].y();
      r.z() = diagData.params.zondCoords[n].z();

      E = get_fieldE_in_pos(_mesh.fieldE,r);
      B = get_fieldB_in_pos(_mesh.fieldB,r);
      fprintf(fZond, "%g %g %g %g %g %g ",  E.x(), E.y(), E.z(), B.x(), B.y(), B.z() );
    }
    
    fprintf(fZond, "\n");
    if(  timestep % TimeStepDelayDiag1D == 0){
      fflush(fZond);
    }

}
/*
static void write_fields_lineX_to_file(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long indy, long indz,
      const MPI_Topology& MPIconf, const Region& domain, MPI_File& fZondLine, long currentStep,long timestep){
    
    MPI_Status status;
    char filename[100];
    auto numCells = (domain.numCells.x() - domain.dampCells[0].x() - domain.dampCells[1].x()) ;
    auto sizeData = 6*numCells;
    long indx;
    
    static Array1D<float>floatData(6*sizeData);
   

    if( timestep == StartTimeStep){
        sprintf(filename, "./Fields/ZondLineX_nodeY%04d_nodeZ%04d.bin",(int)indy,(int)indz);
        MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fZondLine);
        float info = float(numCells*MPIconf.size_line());
    
        if(MPIconf.is_first_line()){   
            MPI_File_write_at(fZondLine, 0, &info, 1, MPI_FLOAT, &status);
        }
    }

    for( auto i = 0; i < numCells; ++i){
        indx = i + domain.dampCells[0].x();
        floatData(6 * i    ) = float(fieldE(indx,indy,indz).x() );
        floatData(6 * i + 1) = float(fieldE(indx,indy,indz).y() );
        floatData(6 * i + 2) = float(fieldE(indx,indy,indz).z() );
        floatData(6 * i + 3) = float(fieldB(indx,indy,indz).x() );
        floatData(6 * i + 4) = float(fieldB(indx,indy,indz).y() );
        floatData(6 * i + 5) = float(fieldB(indx,indy,indz).z() );
    }
  
    int startWrite = MPIconf.rank_line()*sizeData*sizeof(float) + sizeof(float);
  
    int sizeWrite = MPIconf.size_line()*sizeData*sizeof(float);
    
    startWrite += sizeWrite*currentStep;
    
    MPI_File_write_at(fZondLine, startWrite, &floatData.data(0), sizeData, MPI_FLOAT, &status);
      
    //delete[] floatData;
}

void Writer::write_fields_lineX(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long timestep){
  auto numZonds = diagData.params.zondCoordsLineX.size();
  if (!_world.MPIconf.is_master_depth() ) return;
  
  long indy,indz;
  static std::vector<MPI_File> fZondLine;
  static std::vector<long> currentStep; 
  if( timestep == StartTimeStep){
    for (ulong n = 0; n < numZonds; n++ ){
          MPI_File file;
          fZondLine.push_back(file);
          currentStep.push_back(0);
      }
  }
   
  for (ulong n = 0; n < numZonds; n++ ){
    auto coord = diagData.params.zondCoordsLineX[n];
    indy = _mesh.get_node_from_coordY(coord.y() );
    indz = _mesh.get_node_from_coordZ(coord.z() );
    write_fields_lineX_to_file(_mesh.fieldE, _mesh.fieldB, indy, indz,_world.MPIconf, _world.region, fZondLine[n], currentStep[n],timestep);
    currentStep[n]++;
  }
}


static void write_fields_lineY_to_file(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long globIndx, long indz,
      const MPI_Topology& MPIconf, const Region& domain, MPI_File& fZondLine, long currentStep,long timestep){
    
    MPI_Status status;
    char filename[100];
    auto numCells = domain.numNodes.y();
    auto sizeData = 6*numCells;
    long indx = domain.get_index_loc(globIndx);
    long indy;
    static Array1D<float>floatData(6*sizeData);
   

    if( timestep == StartTimeStep){
        sprintf(filename, "./Fields/ZondLineY_nodeX%04d_nodeZ%04d.bin",(int)globIndx,(int)indz);
        MPI_File_open(MPIconf.comm_depth(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fZondLine);
        float info = float(numCells);
    
        if(MPIconf.is_master_depth()){   
            MPI_File_write_at(fZondLine, 0, &info, 1, MPI_FLOAT, &status);
        }
    }

    for( auto i = 0; i < numCells; ++i){
        indy = i;
        floatData(6 * i    ) = float(fieldE(indx,indy,indz).x() );
        floatData(6 * i + 1) = float(fieldE(indx,indy,indz).y() );
        floatData(6 * i + 2) = float(fieldE(indx,indy,indz).z() );
        floatData(6 * i + 3) = float(fieldB(indx,indy,indz).x() );
        floatData(6 * i + 4) = float(fieldB(indx,indy,indz).y() );
        floatData(6 * i + 5) = float(fieldB(indx,indy,indz).z() );
    }
  
    int startWrite =  sizeof(float);
  
    int sizeWrite = sizeData*sizeof(float);
    
    startWrite += sizeWrite*currentStep;
    
    if(MPIconf.is_master_depth()){   
      MPI_File_write_at(fZondLine, startWrite, &floatData.data(0), sizeData, MPI_FLOAT, &status);
    }      
}
static void write_fields_lineZ_to_file(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long globIndx, long indy,
      const MPI_Topology& MPIconf, const Region& domain, MPI_File& fZondLine, long currentStep,long timestep){
    
    MPI_Status status;
    char filename[100];
    auto numCells = domain.numNodes.z();
    auto sizeData = 6*numCells;
    long indx = domain.get_index_loc(globIndx);
    long indz;
    static Array1D<float>floatData(6*sizeData);
   

    if( timestep == StartTimeStep){
        sprintf(filename, "./Fields/ZondLineZ_nodeX%04d_nodeY%04d.bin",(int)globIndx,(int)indy);
        MPI_File_open(MPIconf.comm_depth(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fZondLine);
        float info = float(numCells);
    
        if(MPIconf.is_master_depth()){   
            MPI_File_write_at(fZondLine, 0, &info, 1, MPI_FLOAT, &status);
        }
    }

    for( auto i = 0; i < numCells; ++i){
        indz = i;
        floatData(6 * i    ) = float(fieldE(indx,indy,indz).x() );
        floatData(6 * i + 1) = float(fieldE(indx,indy,indz).y() );
        floatData(6 * i + 2) = float(fieldE(indx,indy,indz).z() );
        floatData(6 * i + 3) = float(fieldB(indx,indy,indz).x() );
        floatData(6 * i + 4) = float(fieldB(indx,indy,indz).y() );
        floatData(6 * i + 5) = float(fieldB(indx,indy,indz).z() );
    }
  
    int startWrite =  sizeof(float);
  
    int sizeWrite = sizeData*sizeof(float);
    
    startWrite += sizeWrite*currentStep;
    
    if(MPIconf.is_master_depth()){   
      MPI_File_write_at(fZondLine, startWrite, &floatData.data(0), sizeData, MPI_FLOAT, &status);
    }      
}


void Writer::write_fields_lineY(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long timestep){
  auto numZonds = diagData.params.zondCoordsLineY.size();
  
  
  long globIndx,indz;
  static std::vector<MPI_File> fZondLine;
  static std::vector<long> currentStep; 
  if( timestep == StartTimeStep){
    for (ulong n = 0; n < numZonds; n++ ){
          MPI_File file;
          fZondLine.push_back(file);
          currentStep.push_back(0);
      }
  }
   
  for (ulong n = 0; n < numZonds; n++ ){
    auto coord = diagData.params.zondCoordsLineY[n];
    if( _world.region.in_region(coord.x() ) ){
      globIndx = _mesh.get_node_from_coordX(coord.x() ); 
      indz = _mesh.get_node_from_coordZ(coord.z() ); 
      write_fields_lineY_to_file(_mesh.fieldE, _mesh.fieldB, globIndx, indz,_world.MPIconf, _world.region, fZondLine[n], currentStep[n],timestep);
      currentStep[n]++;
    }
  }
}


void Writer::write_fields_lineZ(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long timestep){
  auto numZonds = diagData.params.zondCoordsLineY.size();
  
  
  long globIndx,indy;
  static std::vector<MPI_File> fZondLine;
  static std::vector<long> currentStep; 
  if( timestep == StartTimeStep){
    for (const auto& coord : diagData.params.zondCoordsLineZ ){
          MPI_File file;
          fZondLine.push_back(file);
          currentStep.push_back(0);
      }
  }
   
  for (ulong n = 0; n < diagData.params.zondCoordsLineZ.size(); n++ ){
    auto coord = diagData.params.zondCoordsLineZ[n];    
    if( _world.region.in_region(coord.x() ) ){
      globIndx = _mesh.get_node_from_coordX(coord.x() ); 
      indy = _mesh.get_node_from_coordY(coord.y() ); 
      write_fields_lineZ_to_file(_mesh.fieldE, _mesh.fieldB, globIndx, indy,_world.MPIconf, _world.region, fZondLine[n], currentStep[n],timestep);
      currentStep[n]++;
    }
  }
}

static void write_fields_circle_to_file(const Array2D<double3>& circleE, const Array2D<double3>& circleB, long globIndx, long id,
      const MPI_Topology& MPIconf, const Region& domain, MPI_File& fZondLine, long currentStep,long timestep){
    
    MPI_Status status;
    char filename[100];
    auto numCells = circleE.size().y();
    auto sizeData = 6*numCells;
    long indx = domain.get_index_loc(globIndx);
    long indr;
    
    static Array1D<float>floatData(6*sizeData);
   

    if( timestep == StartTimeStep){
        sprintf(filename, "./Fields/ZondCircle_nodeX%04d_radius%04d.bin",(int)globIndx,(int)id);
        MPI_File_open(MPIconf.comm_depth(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fZondLine);
        float info = float(numCells);
    
        if(MPIconf.is_master_depth()){   
            MPI_File_write_at(fZondLine, 0, &info, 1, MPI_FLOAT, &status);
        }
    }

    for( auto i = 0; i < numCells; ++i){
        indr = i;
        floatData(6 * i    ) = float(circleE(indx,indr).x() );
        floatData(6 * i + 1) = float(circleE(indx,indr).y() );
        floatData(6 * i + 2) = float(circleE(indx,indr).z() );
        floatData(6 * i + 3) = float(circleB(indx,indr).x() );
        floatData(6 * i + 4) = float(circleB(indx,indr).y() );
        floatData(6 * i + 5) = float(circleB(indx,indr).z() );
    }
  
    int startWrite =  sizeof(float);
  
    int sizeWrite = sizeData*sizeof(float);
    
    startWrite += sizeWrite*currentStep;
    
    if(MPIconf.is_master_depth()){   
      MPI_File_write_at(fZondLine, startWrite, &floatData.data(0), sizeData, MPI_FLOAT, &status);
    }      
}

void Writer::write_fields_circle( long timestep){
  long globIndx,indr;
  static std::vector<MPI_File> fZondLine;
  static std::vector<long> currentStep; 
  if( timestep == StartTimeStep){
    for (const auto& coord : diagData.params.zondCoordsCircleX ){
        for(const auto& data : diagData.radialDiag){
          MPI_File file;
          fZondLine.push_back(file);
          currentStep.push_back(0);
        }
      }
  }
  for (ulong n = 0; n < diagData.params.zondCoordsCircleX.size(); n++ ){  
    auto coord =  diagData.params.zondCoordsCircleX[n];
    if( _world.region.in_region(coord ) ){
        globIndx = _mesh.get_node_from_coordX(coord ); 
        indr = 0;
        for(const auto& data : diagData.radialDiag){
          indr +=1;
          write_fields_circle_to_file(data.circleE, data.circleB, globIndx, indr,_world.MPIconf, _world.region, fZondLine[n], currentStep[n],timestep);
          currentStep[n]++;
        }

    }
  }
}
*/
/*
static void write_zond_lineX_bin(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long indx,
               const MPI_Topology& MPIconf, const Region& domain, MPI_File& fZondLine, long currentStep,long timestep){
    
    MPI_Status status;
    char filename[50];
    int numCells = (domain.numCells_d2 - domain.dampCells_d2) ;
    auto sizeData = 6*numCells;

    std::cout<<sizeData<<"\n";
    static float* floatData = new float[sizeData];

    //static long currentStep = 0; 

    if( timestep == StartTimeStep){
        sprintf(filename, "./Fields/ZondLineX_node%04d.bin",(int)indx + int(domain.origin/Dx) );

        MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fZondLine);
        float info = float(numCells);
    
        if(MPIconf.is_first_line()){   
            MPI_File_write_at(fZondLine, 0, &info, 1, MPI_FLOAT, &status);
        }
    }

    for( auto j = 0; j < numCells; ++j){
        floatData[6 * j    ] = float(fieldE(indx,j).x() );
        floatData[6 * j + 1] = float(fieldE(indx,j).y() );
        floatData[6 * j + 2] = float(fieldE(indx,j).z() );
        floatData[6 * j + 3] = float(fieldB(indx,j).x() );
        floatData[6 * j + 4] = float(fieldB(indx,j).y() );
        floatData[6 * j + 5] = float(fieldB(indx,j).z() );
    }
  
    int startWrite = MPIconf.rank_line()*sizeData*sizeof(float) + sizeof(float);
  
    int sizeWrite = MPIconf.size_line()*sizeData*sizeof(float);
    
    startWrite += sizeWrite*currentStep;

    MPI_File_write_at(fZondLine, startWrite ,floatData, sizeData, MPI_FLOAT, &status);
   // delete[] floatData;

}

void Writer::diag_zond_lineX_bin(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, long timestep){
 auto numZonds = zondLineX[0];
  if (! _world.MPIconf.is_master_depth() ||  numZonds <= 0) return;
  long indx;
  static MPI_File fZondLine[5];
  static long currentStep[5] = {0,0,0,0,0}; 

  assert (numZonds < 5 && "Number of files less then number of zonds!\n");
   
  for (auto n = 0; n < numZonds; ++n ){
    if( ! _world.region.in_region(zondLineX[n+1]) ) continue;
    indx = long( (zondLineX[n+1] -  _world.region.origin) / Dx);
    write_zond_lineX_bin(_mesh.fieldE, _mesh.fieldB, indx, _world.MPIconf, _world.region, fZondLine[n], currentStep[n],timestep);
    currentStep[n]++;
  }
  
}*/


void Writer::write_fields2D_planeX(const Field3d& fieldE, const Field3d& fieldB, double coordX, const long& timestep){
    if (!_world.region.in_region(coordX ) ) return;
    long globIndex = _mesh.get_node_from_coordX(coordX);
    long index = _world.region.get_index_loc(globIndex);
    
    char filename[100];
    float info;    
    MPI_File fField2D;
    MPI_Status status;
    long indx;

   // long size_x = fieldE.size().x() - ADD_NODES;
    long size_y = fieldE.size().y();
    long size_z = fieldE.size().z();
    long size1 = size_y;
    long size2 = size_z;

    float* floatData[6];
    int sizeData = size1 * size2;

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size1*size2];
    }

    sprintf(filename, ".//Fields//Diag2D//FieldPlaneX_%03ld_%03ld",globIndex,timestep / TimeStepDelayDiag2D); 
    MPI_File_open(_world.MPIconf.comm_depth(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fField2D);
        
     long i = index; 
      for( auto j = 0; j < size_y; j++ ){
          for( auto k = 0; k < size_z; k++ ){
              indx = j*size_z + k;
              floatData[0][indx] = float(fieldE(i,j,k,0) );
              floatData[1][indx] = float(fieldE(i,j,k,1) );
              floatData[2][indx] = float(fieldE(i,j,k,2) );
              floatData[3][indx] = float(fieldB(i,j,k,0) );
              floatData[4][indx] = float(fieldB(i,j,k,1) );
              floatData[5][indx] = float(fieldB(i,j,k,2) );
          }
      }
    
    if(_world.MPIconf.is_master_depth() ){
        info = float(size1);
        MPI_File_write_at(fField2D, 0, &info, 1, MPI_FLOAT, &status);
        info = float(size2);
        MPI_File_write_at(fField2D, sizeof(float), &info, 1, MPI_FLOAT, &status);

        long startWrite =  2 * sizeof(float);
        long sizeWriteField = size1 * size2 * sizeof(float);

        for(auto i = 0; i<6; ++i){
          MPI_File_write_at(fField2D, startWrite + sizeWriteField*i,floatData[i], sizeData, MPI_FLOAT, &status);
        }
    }
    

    
    MPI_File_close(&fField2D);
    
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}



void Writer::write_fields2D_planeZ(const Field3d& fieldE, const Field3d& fieldB, double coordZ, const long& timestep){
    if (!_world.MPIconf.is_master_depth()) return;
    char filename[100];
    float info;    
    MPI_File fField2D;
    MPI_Status status;
    long sizeWriteField;
    long indx;

    long size_x = fieldE.size().x() - ADD_NODES;
    long size_y = fieldE.size().y();
    //long size_z = fieldE.size().z();
    long size1 = size_x;
    long size2 = size_y;


    float* floatData[6]; 
    int sizeData = size1 * size2;

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size1*size2];
    }

    sprintf(filename, ".//Fields//Diag2D//FieldPlaneZ_%03ld_%03ld",_mesh.get_node_from_coordZ(coordZ),timestep / TimeStepDelayDiag2D);    MPI_File_open(_world.MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fField2D);
    
    auto sumSize_x = _world.MPIconf.accum_sum_line(size_x);
    
    long startWrite = sumSize_x * size2 * sizeof(float) + 2 * sizeof(float);

    long k = _mesh.get_node_from_coordZ(coordZ); 
      for( auto i = 0; i < size_x; i++ ){
        for( auto j = 0; j < size_y; j++ ){
            indx = i*size_y + j;
              floatData[0][indx] = float(fieldE(i,j,k,0) );
              floatData[1][indx] = float(fieldE(i,j,k,1) );
              floatData[2][indx] = float(fieldE(i,j,k,2) );
              floatData[3][indx] = float(fieldB(i,j,k,0) );
              floatData[4][indx] = float(fieldB(i,j,k,1) );
              floatData[5][indx] = float(fieldB(i,j,k,2) );
        }
      }
    

    if(_world.MPIconf.is_last_line()){
        info = float(sumSize_x + size1);
        MPI_File_write_at(fField2D, 0, &info, 1, MPI_FLOAT, &status);
        info = float(size2);
        MPI_File_write_at(fField2D, sizeof(float), &info, 1, MPI_FLOAT, &status);
        sizeWriteField = (sumSize_x + size1) * size2 * sizeof(float);
    }
    
    MPI_Bcast(&sizeWriteField, 1, MPI_LONG, _world.MPIconf.last_line(), _world.MPIconf.comm_line() );

    for(auto i = 0; i<6; ++i){
        MPI_File_write_at(fField2D, startWrite + sizeWriteField*i,floatData[i], sizeData, MPI_FLOAT, &status);
    }
    
    MPI_File_close(&fField2D);
    
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}


void Writer::write_fields2D_planeY(const Field3d& fieldE, const Field3d& fieldB, double coordY, const long& timestep){
    if (!_world.MPIconf.is_master_depth()) return;
    char filename[100];
    float info;    
    MPI_File fField2D;
    MPI_Status status;
    long sizeWriteField;
    long indx;

    long size_x = fieldE.size().x() - ADD_NODES;
    //long size_y = fieldE.size().y();
    long size_z = fieldE.size().z();
    long size1 = size_x;
    long size2 = size_z;

    float* floatData[6]; 
    int sizeData = size1 * size2;

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size1*size2];
    }

    sprintf(filename, ".//Fields//Diag2D//FieldPlaneY_%03ld_%03ld",_mesh.get_node_from_coordY(coordY),timestep / TimeStepDelayDiag2D); 
    MPI_File_open(_world.MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fField2D);
    
    auto sumSize_x = _world.MPIconf.accum_sum_line(size_x);
    
    long startWrite = sumSize_x * size2 * sizeof(float) + 2 * sizeof(float);

     long j = _mesh.get_node_from_coordY(coordY); 
      for( auto i = 0; i < size_x; i++ ){
          for( auto k = 0; k < size_z; k++ ){
            indx = i*size_z + k;
              floatData[0][indx] = float(fieldE(i,j,k,0) );
              floatData[1][indx] = float(fieldE(i,j,k,1) );
              floatData[2][indx] = float(fieldE(i,j,k,2) );
              floatData[3][indx] = float(fieldB(i,j,k,0) );
              floatData[4][indx] = float(fieldB(i,j,k,1) );
              floatData[5][indx] = float(fieldB(i,j,k,2) );
          }
      }

    if(_world.MPIconf.is_last_line()){
        info = float(sumSize_x + size1);
        MPI_File_write_at(fField2D, 0, &info, 1, MPI_FLOAT, &status);
        info = float(size2);
        MPI_File_write_at(fField2D, sizeof(float), &info, 1, MPI_FLOAT, &status);
        sizeWriteField = (sumSize_x + size1) * size2 * sizeof(float);
    }
    
    MPI_Bcast(&sizeWriteField, 1, MPI_LONG, _world.MPIconf.last_line(), _world.MPIconf.comm_line() );

    for(auto i = 0; i<6; ++i){
        MPI_File_write_at(fField2D, startWrite + sizeWriteField*i,floatData[i], sizeData, MPI_FLOAT, &status);
    }
    
    MPI_File_close(&fField2D);
    
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}

/*
void Writer::write_fields2D(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, const std::string& axes,long index, const long& timestep){
    if (!_world.MPIconf.is_master_depth()) return;
    char filename[60];
    float info;    
    MPI_File fField2D;
    MPI_Status status;

    if (axes != "xy" && axes !="xz" && axes !="yz") return;

    long size_x = fieldE.size().x() - ADD_NODES;
    long size_y = fieldE.size().y();
    long size_z = fieldE.size().z();
    long size1 = size_x;
    long size2;
    if (axes == "xy"){
      size2 = size_y;
      sprintf(filename, ".//Fields//Diag2D//Field2Dxy_%03ld_%03ld",index,timestep / TimeStepDelayDiag2D);
    }
    else if (axes == "xz"){
      sprintf(filename, ".//Fields//Diag2D//Field2Dxz_%03ld_%03ld",index,timestep / TimeStepDelayDiag2D); 
      size2 = size_z;
    } else if (axes == "yz" || axes == "zy"){
      sprintf(filename, ".//Fields//Diag2D//Field2Dyz_%03ld_%03ld",index,timestep / TimeStepDelayDiag2D); 
      size1 = size_y;
      size2 = size_z;
    } else return;
    int sizeData = size1 * size2;
    long sizeWriteField;
    long indx;
    float* floatData[6];

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size1*size2];
    }

    MPI_File_open(_world.MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fField2D);
    
    auto sumSize_x = _world.MPIconf.accum_sum_line(size_x);
    
    long startWrite = sumSize_x * size2 * sizeof(float) + 2 * sizeof(float);
    
    if (axes == "yz"){
     long i = index; //size_y/2;
      for( auto j = 0; j < size_y; j++ ){
          for( auto k = 0; k < size_z; k++ ){
              indx = j*size_z + k;
              floatData[0][indx] = float(fieldE(i,j,k).x() );
              floatData[1][indx] = float(fieldE(i,j,k).y() );
              floatData[2][indx] = float(fieldE(i,j,k).z() );
              floatData[3][indx] = float(fieldB(i,j,k).x() );
              floatData[4][indx] = float(fieldB(i,j,k).y() );
              floatData[5][indx] = float(fieldB(i,j,k).z() );
          }
      }
    

    } else if (axes == "xy"){
      long k = index; //size_z/2;
      for( auto i = 0; i < size_x; i++ ){
        for( auto j = 0; j < size_y; j++ ){
            indx = i*size_y + j;
              floatData[0][indx] = float(fieldE(i,j,k).x() );
              floatData[1][indx] = float(fieldE(i,j,k).y() );
              floatData[2][indx] = float(fieldE(i,j,k).z() );
              floatData[3][indx] = float(fieldB(i,j,k).x() );
              floatData[4][indx] = float(fieldB(i,j,k).y() );
              floatData[5][indx] = float(fieldB(i,j,k).z() );
        }
      }
    }
    if (axes == "xz"){
     long j = index; //size_y/2;
      for( auto i = 0; i < size_x; i++ ){
          for( auto k = 0; k < size_z; k++ ){
            indx = i*size_z + k;
              floatData[0][indx] = float(fieldE(i,j,k).x() );
              floatData[1][indx] = float(fieldE(i,j,k).y() );
              floatData[2][indx] = float(fieldE(i,j,k).z() );
              floatData[3][indx] = float(fieldB(i,j,k).x() );
              floatData[4][indx] = float(fieldB(i,j,k).y() );
              floatData[5][indx] = float(fieldB(i,j,k).z() );
          }
      }
    }



    if(_world.MPIconf.is_last_line()){
        info = float(sumSize_x + size1);
        MPI_File_write_at(fField2D, 0, &info, 1, MPI_FLOAT, &status);
        info = float(size2);
        MPI_File_write_at(fField2D, sizeof(float), &info, 1, MPI_FLOAT, &status);
        sizeWriteField = (sumSize_x + size1) * size2 * sizeof(float);
    }
    

    MPI_Bcast(&sizeWriteField, 1, MPI_LONG, _world.MPIconf.last_line(), _world.MPIconf.comm_line() );

    for(auto i = 0; i<6; ++i){
        MPI_File_write_at(fField2D, startWrite + sizeWriteField*i,floatData[i], sizeData, MPI_FLOAT, &status);
    }
    
    MPI_File_close(&fField2D);
    
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}*/

void Writer::write_fields2D_circle(const Array2D<double3>& fieldE, const Array2D<double3>& fieldB, long series, const long& timestep){
    if (!_world.MPIconf.is_master_depth()) return;
    char filename[60];
    float info;    
    MPI_File fField2D;
    MPI_Status status;

    long size_x = fieldE.size().x() - ADD_NODES;
    long size_y = fieldE.size().y();
    long size1 = size_x;
    long size2 = size_y;
    sprintf(filename, ".//Fields//Diag2D//Field2Dcircle_%03ld_%03ld",series,timestep / TimeStepDelayDiag2D);
    
    int sizeData = size1 * size2;
    long sizeWriteField;
    long indx;
    float* floatData[6];

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size1*size2];
    }

    MPI_File_open(_world.MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fField2D);

      for( auto i = 0; i < size_x; i++ ){
        for( auto j = 0; j < size_y; j++ ){
            indx = i*size_y + j;
              floatData[0][indx] = float(fieldE(i,j).x() );
              floatData[1][indx] = float(fieldE(i,j).y() );
              floatData[2][indx] = float(fieldE(i,j).z() );
              floatData[3][indx] = float(fieldB(i,j).x() );
              floatData[4][indx] = float(fieldB(i,j).y() );
              floatData[5][indx] = float(fieldB(i,j).z() );
        }
      }
    


    auto sumSize_x = _world.MPIconf.accum_sum_line(size1);
    
    long startWrite = sumSize_x * size2 * sizeof(float) + 2 * sizeof(float);

    if(_world.MPIconf.is_last_line()){
        info = float(sumSize_x + size1);
        MPI_File_write_at(fField2D, 0, &info, 1, MPI_FLOAT, &status);
        info = float(size2);
        MPI_File_write_at(fField2D, sizeof(float), &info, 1, MPI_FLOAT, &status);
        sizeWriteField = (sumSize_x + size1) * size2 * sizeof(float);
    }
    

    MPI_Bcast(&sizeWriteField, 1, MPI_LONG, _world.MPIconf.last_line(), _world.MPIconf.comm_line() );

    for(auto i = 0; i<6; ++i){
        MPI_File_write_at(fField2D, startWrite + sizeWriteField*i,floatData[i], sizeData, MPI_FLOAT, &status);
    }
    
    MPI_File_close(&fField2D);
    
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}

void Writer::write_fields3D(const Field3d& fieldE, const Field3d& fieldB, const long& timestep){
   // if (!_world.MPIconf.is_master_depth()) return;
    
    MPI_File fField3D;
    MPI_Status status;
    
    char filename[100];
    float info;
    long size_x = fieldE.size().x() - ADD_NODES;
    long size_y = fieldE.size().y();
    long size_z = fieldE.size().z();
    int sizeData = size_x *size_y * size_z;
    long sizeWriteField;
    long indx;
    float* floatData[6];

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size_x*size_y*size_z];
    }

    sprintf(filename, ".//Fields//Diag3D//Field3D%03ld_%03ld",_world.MPIconf.rank_line(),timestep / TimeStepDelayDiag2D);

    MPI_File_open(_world.MPIconf.comm_depth(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fField3D);

    for( auto i = 0; i < size_x; i++ ){
      for( auto j = 0; j < size_y; j++ ){
        for( auto k = 0; k < size_z; k++ ){
          indx = i*size_y*size_z + j*size_z + k;
              floatData[0][indx] = float(fieldE(i,j,k,0) );
              floatData[1][indx] = float(fieldE(i,j,k,1) );
              floatData[2][indx] = float(fieldE(i,j,k,2) );
              floatData[3][indx] = float(fieldB(i,j,k,0) );
              floatData[4][indx] = float(fieldB(i,j,k,1) );
              floatData[5][indx] = float(fieldB(i,j,k,2) );
        }
      }
    }
    
    //auto sumSize_x = _world.MPIconf.accum_sum_line(size_x);
    
    long startWrite = 3 * sizeof(float);

    if(_world.MPIconf.is_master_depth()){
        info = float(size_x);
        MPI_File_write_at(fField3D, 0, &info, 1, MPI_FLOAT, &status);
        info = float(size_y);
        MPI_File_write_at(fField3D, sizeof(float), &info, 1, MPI_FLOAT, &status);
        info = float(size_z);
        MPI_File_write_at(fField3D, 2*sizeof(float), &info, 1, MPI_FLOAT, &status);
        sizeWriteField = size_x * size_y * size_z * sizeof(float);
        for(auto i = 0; i<6; ++i){
            MPI_File_write_at(fField3D, startWrite + sizeWriteField*i,floatData[i], sizeData, MPI_FLOAT, &status);
        }
    }
    
    MPI_File_close(&fField3D);
    
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}
/*
void Writer::write_fields3D(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, const long& timestep){
    if (!_world.MPIconf.is_master_depth()) return;
    
    MPI_File fField3D;
    MPI_Status status;
    
    char filename[100];
    float info;
    long size_x = fieldE.size().x() - ADD_NODES;
    long size_y = fieldE.size().y();
    long size_z = fieldE.size().z();
    int sizeData = size_x *size_y * size_z;
    long sizeWriteField;
    long indx;
    float* floatData[6];

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size_x*size_y*size_z];
    }

    sprintf(filename, ".//Fields//Diag3D//Field3D%03ld",timestep / TimeStepDelayDiag2D);

    MPI_File_open(_world.MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fField3D);

    for( auto i = 0; i < size_x; i++ ){
      for( auto j = 0; j < size_y; j++ ){
        for( auto k = 0; k < size_z; k++ ){
          indx = i*size_y*size_z + j*size_z + k;
            floatData[0][indx] = float(fieldE(i,j,k).x() );
            floatData[1][indx] = float(fieldE(i,j,k).y() );
            floatData[2][indx] = float(fieldE(i,j,k).z() );
            floatData[3][indx] = float(fieldB(i,j,k).x() );
            floatData[4][indx] = float(fieldB(i,j,k).y() );
            floatData[5][indx] = float(fieldB(i,j,k).z() );
        }
      }
    }
    
    auto sumSize_x = _world.MPIconf.accum_sum_line(size_x);
    
    long startWrite = sumSize_x * size_y * size_z* sizeof(float) + 3 * sizeof(float);

    if(_world.MPIconf.is_last_line()){
        info = float(sumSize_x + size_x);
        MPI_File_write_at(fField3D, 0, &info, 1, MPI_FLOAT, &status);
        info = float(size_y);
        MPI_File_write_at(fField3D, sizeof(float), &info, 1, MPI_FLOAT, &status);
        info = float(size_z);
        MPI_File_write_at(fField3D, 2*sizeof(float), &info, 1, MPI_FLOAT, &status);
        sizeWriteField = (sumSize_x + size_x) * size_y * size_z * sizeof(float);
    }
    

    MPI_Bcast(&sizeWriteField, 1, MPI_LONG, _world.MPIconf.last_line(), _world.MPIconf.comm_line() );

    //int startWrite = _world.MPIconf.rank_line()*(PlasmaCellsX_glob/_world.MPIconf.size_line())*size_y*sizeof(float) + 2*sizeof(float);
    
    //if(!_world.MPIconf.is_first_line() ){ 
    //  startWrite += DampCellsX_glob[0] * size_y * sizeof(float);
    //}
//    int sizeWrite = (PlasmaCellsX_glob + DampCellsX_glob[0] + DampCellsX_glob[1])*size_y*sizeof(float);
    //int sizeWrite = (PlasmaCellsX_glob + DampCellsX_glob[0] + DampCellsX_glob[1])*size_y*sizeof(float);
    
    for(auto i = 0; i<6; ++i){
        MPI_File_write_at(fField3D, startWrite + sizeWriteField*i,floatData[i], sizeData, MPI_FLOAT, &status);
    }
    
    MPI_File_close(&fField3D);
    
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}*/


/*
void Writer::write_fields2D(const Array3D<double3>& fieldE, const Array3D<double3>& fieldB, const std::string& axes ,const long& timestep){
    if (!_world.MPIconf.is_master_depth()) return;
    
    MPI_File fField3D;
    MPI_Status status;
    
    char filename[50];
    float info;

    long size_x = fieldE.size().x() - ADD_NODES;
    long size_y = fieldE.size().y();
    long size_z = fieldE.size().z();
    int sizeData = size_x *size_y * size_z;
    long sizeWriteField;
    long indx;
    float* floatData[6];

    for(auto i = 0; i<6; i++){
        floatData[i] = new float[size_x*size_y*size_z];
    }

    sprintf(filename, ".//Fields//Diag2D//Field2Dxy%03ld",timestep / TimeStepDelayDiag2D);
    MPI_File_open(_world.MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fField3D);

    for( auto i = 0; i < size_x; i++ ){
      for( auto j = 0; j < size_y; j++ ){
        for( auto k = 0; k < size_z; k++ ){
          indx = i*size_y*size_z + j*size_z + k;
            floatData[0][indx] = float(fieldE(i,j,k).x() );
            floatData[1][indx] = float(fieldE(i,j,k).y() );
            floatData[2][indx] = float(fieldE(i,j,k).z() );
            floatData[3][indx] = float(fieldB(i,j,k).x() );
            floatData[4][indx] = float(fieldB(i,j,k).y() );
            floatData[5][indx] = float(fieldB(i,j,k).z() );
        }
      }
    }
    
    auto sumSize_x = _world.MPIconf.accum_sum_line(size_x);
    
    long startWrite = sumSize_x * size_y * size_z* sizeof(float) + 3 * sizeof(float);

    if(_world.MPIconf.is_last_line()){
        info = float(sumSize_x + size_x);
        MPI_File_write_at(fField3D, 0, &info, 1, MPI_FLOAT, &status);
        info = float(size_y);
        MPI_File_write_at(fField3D, sizeof(float), &info, 1, MPI_FLOAT, &status);
        info = float(size_z);
        MPI_File_write_at(fField3D, 2*sizeof(float), &info, 1, MPI_FLOAT, &status);
        sizeWriteField = (sumSize_x + size_x) * size_y * size_z * sizeof(float);
    }
    

    MPI_Bcast(&sizeWriteField, 1, MPI_LONG, _world.MPIconf.last_line(), _world.MPIconf.comm_line() );

    //int startWrite = _world.MPIconf.rank_line()*(PlasmaCellsX_glob/_world.MPIconf.size_line())*size_y*sizeof(float) + 2*sizeof(float);
    
    //if(!_world.MPIconf.is_first_line() ){ 
    //  startWrite += DampCellsX_glob[0] * size_y * sizeof(float);
    //}
//    int sizeWrite = (PlasmaCellsX_glob + DampCellsX_glob[0] + DampCellsX_glob[1])*size_y*sizeof(float);
    //int sizeWrite = (PlasmaCellsX_glob + DampCellsX_glob[0] + DampCellsX_glob[1])*size_y*sizeof(float);
    
    for(auto i = 0; i<6; ++i){
        MPI_File_write_at(fField3D, startWrite + sizeWriteField*i,floatData[i], sizeData, MPI_FLOAT, &status);
    }
    
    MPI_File_close(&fField3D);
    
    for(auto i = 0; i<6; i++){
        delete[] floatData[i];
    }

}*/