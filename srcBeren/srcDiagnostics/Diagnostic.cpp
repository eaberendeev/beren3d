#include "Diagnostic.h"

// Устанавливаем значения основных параметров в переменной Params в соответствии с содержимым сроки
void DiagData::set_params_from_string(const std::string& line){
    std::vector<std::string> strvec;

    strvec = split(line, ' ');

    if(strvec[0]=="radiationDiagRadiuses") {
      for(ulong i = 1; i < strvec.size();i++ ){
        this->params.radiationDiagRadiuses.push_back( stod(strvec[i] ) );
      }
    }
    if(strvec[0]=="sliceRadiationPlaneX") {
      for(ulong i = 1; i < strvec.size(); i++ ){
        this->params.sliceRadiationPlaneX.push_back( stod(strvec[i] ) );
      }
    }
    if(strvec[0]=="sliceRadiationPlaneY") {
      for(ulong i = 1; i < strvec.size(); i++ ){
        this->params.sliceRadiationPlaneY.push_back( stod(strvec[i] ) );
      }
    }

    if(strvec[0]=="sliceRadiationPlaneZ") {
      for(ulong i = 1; i < strvec.size() ; i++){
        this->params.sliceRadiationPlaneZ.push_back( stod(strvec[i] ) );
      }
    }

    if(strvec[0]=="sliceFieldsPlaneX") {
      for(ulong i = 1; i < strvec.size() ; i++){
        this->params.sliceFieldsPlaneX.push_back( stod(strvec[i] ) );
      }
    }
    if(strvec[0]=="sliceFieldsPlaneY") {
      for(ulong i = 1; i < strvec.size(); i++ ){
        this->params.sliceFieldsPlaneY.push_back( stod(strvec[i] ) );
      }
    }
    if(strvec[0]=="sliceFieldsPlaneZ") {
      for(ulong i = 1; i < strvec.size(); i++){
        this->params.sliceFieldsPlaneZ.push_back( stod(strvec[i] ) );
      }
    }   
    if(strvec[0]=="zondCoords") {
      for(ulong i = 1; i < strvec.size(); i+=3 ){
        auto strx = strvec[i];
        auto stry = strvec[i+1];
        auto strz = strvec[i+2];
        strx.erase(strx.find('('), 1);
        strz.erase(strz.find(')'), 1);
        this->params.zondCoords.push_back( double3(stod(strx ),
                                                    stod(stry ),
                                                    stod(strz ) ) );
      }
    }
    if(strvec[0]=="zondCoordsLineX") {
      for(ulong i = 1; i < strvec.size(); i+=3 ){
        auto strx = strvec[i];
        auto stry = strvec[i+1];
        auto strz = strvec[i+2];
        strx.erase(strx.find('('), 1);
        strz.erase(strz.find(')'), 1);
        this->params.zondCoordsLineX.push_back( double3(stod(strx ),
                                                    stod(stry ),
                                                    stod(strz ) ) );
      }
    }          
    if(strvec[0]=="zondCoordsLineY") {
      for(ulong i = 1; i < strvec.size(); i+=3 ){
        auto strx = strvec[i];
        auto stry = strvec[i+1];
        auto strz = strvec[i+2];
        strx.erase(strx.find('('), 1);
        strz.erase(strz.find(')'), 1);
        this->params.zondCoordsLineY.push_back( double3(stod(strx ),
                                                    stod(stry ),
                                                    stod(strz ) ) );
      }
    }   
    if(strvec[0]=="zondCoordsLineZ") {
      for(ulong i = 1; i < strvec.size(); i+=3 ){
        auto strx = strvec[i];
        auto stry = strvec[i+1];
        auto strz = strvec[i+2];
        strx.erase(strx.find('('), 1);
        strz.erase(strz.find(')'), 1);
        this->params.zondCoordsLineZ.push_back( double3(stod(strx ),
                                                    stod(stry ),
                                                    stod(strz ) ) );
      }
    }       

    if(strvec[0]=="outTime3D") {
      for(ulong i = 1; i < strvec.size(); i++ ){
        this->params.outTime3D.push_back( stol(strvec[i] ) );
      }
    }
    

}


void Writer::output(long timestep){

////// Calculation power of radiation ///////    
/////////////////////////////////////////////
    /// Radiation on plane X
    for (ulong i = 0; i < diagData.powerRadPlaneX.size(); i++){
        diagData.calc_radiation_pointing_planeX(diagData.powerRadPlaneX[i], diagData.params.sliceRadiationPlaneX[i], _mesh,_world.region);
    }
    /// Radiation on plane Y
    for (ulong i = 0; i < diagData.powerRadPlaneY.size(); i++){
        diagData.calc_radiation_pointing_planeY(diagData.powerRadPlaneY[i], diagData.params.sliceRadiationPlaneY[i], _mesh,_world.region);
    }
    /// Radiation on plane Z
    for (ulong i = 0; i < diagData.powerRadPlaneZ.size(); i++){
        diagData.calc_radiation_pointing_planeZ(diagData.powerRadPlaneZ[i], diagData.params.sliceRadiationPlaneZ[i], _mesh,_world.region);
    }
    /// Radiation on circle
    for (auto& radial : diagData.radialDiag){
        radial.calc_radiation_pointing_circle2D(_mesh);
    }

////////////////////////////////////////////////
////////////////////////////////////////////////
//////////////// WRITE DATA ////////////////////
////////////////////////////////////////////////
////////////////////////////////////////////////    


/////// 2D DATA //////////////////////////////
    if (timestep % TimeStepDelayDiag2D == 0){

        write_radiation_circle2D(timestep);
        write_radiation_planes(timestep);

        write_particles2D(timestep);

        for (auto coordX : diagData.params.sliceFieldsPlaneX){
            write_fields2D_planeX(_mesh.fieldE, _mesh.fieldB, coordX, timestep);
        }
        for (auto coordY : diagData.params.sliceFieldsPlaneY){
            write_fields2D_planeY(_mesh.fieldE, _mesh.fieldB, coordY, timestep);
        }     
        for (auto coordZ : diagData.params.sliceFieldsPlaneZ){
            write_fields2D_planeZ(_mesh.fieldE, _mesh.fieldB, coordZ, timestep);
        }
        long series = 0;
        for(const auto& data:diagData.radialDiag){
            write_fields2D_circle(data.circleE, data.circleB, series,timestep);
            series++;
        }
    }
///////////////////////////////////

///// 3D DATA ////////////////////    
    for (const auto& time3D : diagData.params.outTime3D){
        if ( long(Dt*timestep) == time3D ){
            write_particles3D(timestep);
            write_fields3D(_mesh.fieldE, _mesh.fieldB,timestep);
        }    
    }
//////////////////////////////////


////// 1D DATA ////////////////////////////////
    if (timestep % 2 == 0){
      diag_zond(timestep);
      write_fields_lineX(_mesh.fieldE, _mesh.fieldB, timestep);
      write_fields_lineY(_mesh.fieldE, _mesh.fieldB, timestep);
      write_fields_lineZ(_mesh.fieldE, _mesh.fieldB, timestep);
      write_fields_circle(timestep);
    }

    if (timestep % TimeStepDelayDiag1D == 0){

   //     write_radiation_line(timestep);

        diagData.calc_energy(_mesh,_species);
        write_energies(timestep);
        for (auto& radial: diagData.radialDiag){
            radial.powerRadCircle = 0;
        }  
    }

///// RECOVERY DATA ///////////////////
    if (timestep % RecTimeStep == 0 && timestep != StartTimeStep){
        for (const auto& sp: _species){
            sp.write_recovery(_world.MPIconf);
        }
        _mesh.write_recovery(_world.MPIconf);
    }
///////////////////////////////////////
}


void make_folders(){
    mkdir(".//Fields", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Fields//Diag3D", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Fields//Diag2D", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Fields//Diag1D", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Recovery", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Recovery//Fields", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Recovery//Particles", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Anime", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Particles", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir(".//Performance", S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    std::cout << "Create folders for output...\n";
}

Writer::Writer(const World &world, const Mesh &mesh,std::vector<ParticlesArray> &species) : 
    _world(world),_mesh(mesh),_species(species),diagData(world.region) {
  if( !_world.MPIconf.is_master() ) return; 
  
  make_folders(); 
  for( const auto &sp : _species){
    mkdir((".//Particles//" + sp.name).c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir((".//Particles//" + sp.name+"//Diag1D").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir((".//Particles//" + sp.name+"//Diag2D").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    mkdir((".//Particles//" + sp.name+"//Diag3D").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  fDiagEnergies = fopen("Energies.dat", "w");
  
}



void write_array2D(const Array2D<double>& data, long size1, long size2, const char* filename, const MPI_Topology& MPIconf){
    MPI_File fData2D;
    MPI_Status status;

    long sumSize1;
    float info;
    
    if ( !MPIconf.is_master_depth() ) return;

    Array2D<float> floatData(size1,size2);
    long sizeData = floatData.capacity();

    MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fData2D);
    
    for( auto i = 0; i < size1; ++i ){
       for( auto j = 0; j < size2; ++j ){
          floatData(i,j) = float(data(i,j));
      }
    }


    sumSize1 = MPIconf.accum_sum_line(size1);

    if( MPIconf.is_last_line() ){
          info = float(sumSize1 + size1);
          MPI_File_write_at(fData2D, 0,&info, 1,MPI_FLOAT, &status);
          info = float(size2);
          MPI_File_write_at(fData2D, sizeof(float),&info, 1,MPI_FLOAT, &status);
    }
    long startWrite = sumSize1 * size2 * sizeof(float) + 2 * sizeof(float);
        
    MPI_File_write_at(fData2D, startWrite, &floatData.data(0), sizeData, MPI_FLOAT, &status);
    MPI_File_close(&fData2D);
  
}


void write_array3D(const Array3D<double>& data, long size1, long size2, long size3,const char* filename, const MPI_Topology& MPIconf){
    MPI_File fData3D;
    MPI_Status status;
    //long sumSize1,indx;
    float info;
    
    //if ( !MPIconf.is_master_depth() ) return;

    Array3D<float> floatData(size1,size2,size3);
    int sizeData = size1*size2*size3;

    MPI_File_open(MPIconf.comm_depth(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fData3D);

    for( auto i = 0; i < size1; ++i ){
       for( auto j = 0; j < size2; ++j ){
         for( auto k = 0; k < size3; ++k ){
            floatData(i,j,k) = float(data(i,j,k));
        }
      }
    }

    long startWrite = 3 * sizeof(float);

    if( MPIconf.is_master_depth() ){
          info = float(size1);
          MPI_File_write_at(fData3D, 0,&info, 1,MPI_FLOAT, &status);
          info = float(size2);
          MPI_File_write_at(fData3D, sizeof(float),&info, 1,MPI_FLOAT, &status);
          info = float(size3);
          MPI_File_write_at(fData3D, 2*sizeof(float),&info, 1,MPI_FLOAT, &status);

          MPI_File_write_at(fData3D, startWrite, &floatData.data(0), sizeData, MPI_FLOAT, &status);
    }
        
    MPI_File_close(&fData3D);

}

/*
void write_array3D(const Array3D<double>& data, long size1, long size2, long size3,const char* filename, const MPI_Topology& MPIconf){
    MPI_File fData3D;
    MPI_Status status;
    long sumSize1,indx;
    float info;
    
    if ( !MPIconf.is_master_depth() ) return;

    Array3D<float> floatData(size1,size2,size3);
    int sizeData = size1*size2*size3;

    MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fData3D);

    for( auto i = 0; i < size1; ++i ){
       for( auto j = 0; j < size2; ++j ){
         for( auto k = 0; k < size3; ++k ){
            floatData(i,j,k) = float(data(i,j,k));
        }
      }
    }


    sumSize1 = MPIconf.accum_sum_line(size1);

    if( MPIconf.is_last_line() ){
          info = float(sumSize1 + size1);
          MPI_File_write_at(fData3D, 0,&info, 1,MPI_FLOAT, &status);
          info = float(size2);
          MPI_File_write_at(fData3D, sizeof(float),&info, 1,MPI_FLOAT, &status);
          info = float(size3);
          MPI_File_write_at(fData3D, 2*sizeof(float),&info, 1,MPI_FLOAT, &status);
    }
    long startWrite = sumSize1 * size2 * size3 * sizeof(float) + 3 * sizeof(float);
        
    MPI_File_write_at(fData3D, startWrite, &floatData.data(0), sizeData, MPI_FLOAT, &status);
    MPI_File_close(&fData3D);

}*/