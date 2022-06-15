#include "Diagnostic.h"
#include <assert.h>

void DiagData::calc_energy(const Mesh &mesh, const std::vector<ParticlesArray> &species){
  
  for ( auto &sp : species){
    energyParticlesKinetic[sp.name] = sp.get_kinetic_energy();
    energyParticlesInject[sp.name] = sp.get_inject_energy();
  }

  energyFieldE = mesh.get_fieldE_energy();
  energyFieldB = mesh.get_fieldB_energy();

}

void Writer::write_energies(long timestep){
  std::stringstream ss;

  const MPI_Topology &MPIconf = _world.MPIconf;

  if(timestep == 0){
    ss << "Time ";
    for (auto it = diagData.energyParticlesKinetic.begin(); it != diagData.energyParticlesKinetic.end(); ++it){
      ss << "Area_" << it->first << " ";
    }
    for (auto it = diagData.energyParticlesInject.begin(); it != diagData.energyParticlesInject.end(); ++it){
      ss << "Injection_" << it->first << " ";
    }
    ss << "Area_E^2 " << "Area_B^2 " ;
      // << "powerRadXmin "<< "powerRadXmax "<< "powerRadYmin " << "powerRadYmax "<< "powerRadZmin " << "powerRadZmax " ;
       //<< "powerRadAvgXmin "<< "powerRadAvgXmax "<< "powerRadAvgYmin "<< "powerRadAvgYmax "<< "powerRadAvgZmin "<< "powerRadAvgZmax\n";
    for (const auto& data : diagData.radialDiag){
      ss << "powerRadCircle(radius " << data.radiationDiagRadius << ") " ;
    }
    ss << "\n";
  }
  
  ss << timestep*Dt << " ";
    
  for (auto it = diagData.energyParticlesKinetic.begin(); it != diagData.energyParticlesKinetic.end(); ++it){
    double energyParticles = it->second;
    MPI_Allreduce ( MPI_IN_PLACE, &energyParticles, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    ss << energyParticles << " ";
  }
  
  for (auto it = diagData.energyParticlesInject.begin(); it != diagData.energyParticlesInject.end(); ++it){
    double energyInject = it->second;
    MPI_Allreduce ( MPI_IN_PLACE, &energyInject, 1, MPI_DOUBLE, MPI_SUM, MPIconf.comm_line() );
    ss << energyInject << " ";
  }
  std::vector<double> vecEnergy = { diagData.energyFieldE, diagData.energyFieldB,
              //diagData.powerRad.x_min, diagData.powerRad.x_max,
              //diagData.powerRad.y_min, diagData.powerRad.y_max,
              //diagData.powerRad.z_min, diagData.powerRad.z_max,
              //diagData.powerRadAvg.x_min, diagData.powerRadAvg.x_max,
              //diagData.powerRadAvg.y_min, diagData.powerRadAvg.y_max,
              //diagData.powerRadAvg.z_min, diagData.powerRadAvg.z_max, 
              };
  for (const auto& data : diagData.radialDiag){
      vecEnergy.push_back(data.powerRadCircle);
  }
  
  MPI_Allreduce ( MPI_IN_PLACE, &vecEnergy[0], vecEnergy.size(), MPI_DOUBLE, MPI_SUM, MPIconf.comm_line() );

  for(uint i = 0; i< vecEnergy.size(); ++i){
  ss << vecEnergy[i] << " "; 
  }
  ss <<"\n";
  if( MPIconf.is_master() ){
      fprintf(fDiagEnergies, "%s",  ( ss.str() ).c_str() ); 
      std::cout << ss.str();
  }
    if( MPIconf.is_master() && timestep % TimeStepDelayDiag1D == 0){
      fflush(fDiagEnergies);
    }


}
/*
void DiagData::calc_radiation_pointing_box(const Mesh &mesh,const Region &region){
  long i,j,k;
  double3 fieldEAvg;
  double3 fieldBAvg;
  
  powerRad.clear();

  assert(region.numCells.z() <= powerRadLine.x_min.size() && "The analysis of the radiation at the boundary X should be along the axis Z. \n Check the size of the array\n");

  if(region.boundType[0].x() == OPEN){
      i = (zondRadX[0] - region.origin) / Dx;
        for( k = 0; k < region.numCells.z(); k++){
          for( j = 0; j < region.numCells.y(); j++){

            fieldEAvg = mesh.get_fieldE_in_cell(i,j,k);
            fieldBAvg = mesh.get_fieldB_in_cell(i,j,k);
            powerRadLine.x_min(k) += Dt*Dy*Dz*fabs(fieldBAvg.z()*fieldEAvg.y()-fieldEAvg.z()*fieldBAvg.y());
            ///powerRad.x_min += Dt*Dy*Dz*fabs(fieldBAvg.z()*fieldEAvg.y()-fieldEAvg.z()*fieldBAvg.y());
          }
          powerRad.x_min += powerRadLine.x_min(k);
      }
    }
    if(region.boundType[1].x() == OPEN){
        i = (zondRadX[1] - region.origin) / Dx;
        for( k = 0; k < region.numCells.z(); k++){
          for( j = 0; j < region.numCells.y() ; j++){
            fieldEAvg = mesh.get_fieldE_in_cell(i,j,k);
            fieldBAvg = mesh.get_fieldB_in_cell(i,j,k);
            powerRadLine.x_max(k) += Dt*Dy*Dz*fabs(fieldBAvg.z()*fieldEAvg.y()-fieldEAvg.z()*fieldBAvg.y());
            //powerRad.x_max  += Dt*Dy*Dz*fabs(fieldBAvg.z()*fieldEAvg.y()-fieldEAvg.z()*fieldBAvg.y());
          }
          powerRad.x_max += powerRadLine.x_max(k);

        }
    }
    
    j = zondRadY[0] / Dy;
    for (i = 0; i < region.numCells.x() ; i++){
      for (k = 0; k < region.numCells.z() ; k++){
        fieldEAvg = mesh.get_fieldE_in_cell(i,j,k);
        fieldBAvg = mesh.get_fieldB_in_cell(i,j,k);
        //powerRadLine.y_min(i) += Dt*Dx * fabs(fieldEAvg.x()*fieldBAvg.z()-fieldEAvg.z()*fieldBAvg.x());
        powerRadLine.y_min(i) += Dt*Dx*Dz* fabs(fieldEAvg.x()*fieldBAvg.z()-fieldEAvg.z()*fieldBAvg.x());
      }
      powerRad.y_min += powerRadLine.y_min(i);
    }
    j = zondRadY[1] / Dy;
    for (i = 0; i < region.numCells.x() ; i++){
      for (k = 0; k < region.numCells.z() ; k++){

        fieldEAvg = mesh.get_fieldE_in_cell(i,j,k);
        fieldBAvg = mesh.get_fieldB_in_cell(i,j,k);
        powerRadLine.y_max(i) += Dt*Dx * Dz*fabs(fieldEAvg.x()*fieldBAvg.z()-fieldEAvg.z()*fieldBAvg.x());
      }
      powerRad.y_max += powerRadLine.y_max(i);
    }

    k = zondRadZ[0] / Dz;
    for (i = 0; i < region.numCells.x() ; i++){
      for (j = 0; j < region.numCells.y() ; j++){
        fieldEAvg = mesh.get_fieldE_in_cell(i,j,k);
        fieldBAvg = mesh.get_fieldB_in_cell(i,j,k);
        //powerRadLine.y_min(i) += Dt*Dx * fabs(fieldEAvg.x()*fieldBAvg.z()-fieldEAvg.z()*fieldBAvg.x());
        powerRadLine.z_min(i) += Dt*Dx * Dy * fabs(fieldEAvg.x()*fieldBAvg.y()-fieldEAvg.y()*fieldBAvg.x()) ;
      }
        powerRad.z_min += powerRadLine.z_min(i);
    }
    k = zondRadZ[1] / Dz;
    for (i = 0; i < region.numCells.x() ; i++){
      for (j = 0; j < region.numCells.y() ; j++){
        fieldEAvg = mesh.get_fieldE_in_cell(i,j,k);
        fieldBAvg = mesh.get_fieldB_in_cell(i,j,k);
        powerRadLine.z_max(i) += Dt*Dx * Dy* fabs(fieldEAvg.x()*fieldBAvg.y()-fieldEAvg.y()*fieldBAvg.x());
      }
        powerRad.z_max += powerRadLine.z_max(i);
    } 
}
*/


/*
void DiagData::calc_radiation_pointing_planeX(Array2D<double>& dataPlaneX, double coordX, const Mesh &mesh, const Region &region){
  double3 E, B, r;

  if (!region.in_region(coordX) ) return;
    
  r.x() = coordX - region.origin;

  for (auto j = 0; j < dataPlaneX.size().x()-1; ++j ){
    for (auto k = 0; k < dataPlaneX.size().y()-1 ; ++k ){
      r.y() = mesh.get_coord_from_nodeY(j) + 0.5*Dy;
      r.z() = mesh.get_coord_from_nodeZ(k) + 0.5*Dz;

      E = get_fieldE_in_pos(mesh.fieldE,r);
      B = get_fieldB_in_pos(mesh.fieldB,r) - get_fieldB_in_pos(mesh.fieldB0,r);
      dataPlaneX(j,k) += Dt*Dy*Dz* fabs(E.y()*B.z() - E.z()*B.y());
    }
  }

}
void DiagData::calc_radiation_pointing_planeY(Array2D<double>& dataPlaneY, double coordY, const Mesh &mesh, const Region &region){
  double3 E, B, r;
    
  r.y() = coordY;

  for (auto i = 0; i < dataPlaneY.size().x()-1; ++i ){
    for (auto k = 0; k < dataPlaneY.size().y()-1 ; ++k ){
      r.x() = mesh.get_coord_from_nodeX(i) + 0.5*Dx;
      r.z() = mesh.get_coord_from_nodeZ(k) + 0.5*Dz;

      E = get_fieldE_in_pos(mesh.fieldE,r);
      B = get_fieldB_in_pos(mesh.fieldB,r) - get_fieldB_in_pos(mesh.fieldB0,r);
      dataPlaneY(i,k) += Dt*Dx*Dz* fabs(E.x()*B.z() - E.z()*B.x());
    }
  }

}
void DiagData::calc_radiation_pointing_planeZ(Array2D<double>& dataPlaneZ, double coordZ, const Mesh &mesh, const Region &region){
  double3 E, B, r;
    
  r.z() = coordZ;

  for (auto i = 0; i < dataPlaneZ.size().x()-1; ++i ){
    for (auto j = 0; j < dataPlaneZ.size().y()-1 ; ++j ){
      r.x() = mesh.get_coord_from_nodeX(i) + 0.5*Dx;
      r.y() = mesh.get_coord_from_nodeY(j) + 0.5*Dy;

      E = get_fieldE_in_pos(mesh.fieldE,r);
      B = get_fieldB_in_pos(mesh.fieldB,r) - get_fieldB_in_pos(mesh.fieldB0,r);
      dataPlaneZ(i,j) += Dt*Dx*Dy* fabs(E.y()*B.x() - E.x()*B.y());
    }
  }

}



void RadialDiagData::calc_radiation_pointing_circle2D(const Mesh &mesh){
  double angle,Ep,Bp,rad,radX;
  double3 r;
  double3 E,B;
  double center_y = 0.5*Dy*(mesh.fieldE.size().y() - ADD_NODES);
  double center_z = 0.5*Dz*(mesh.fieldE.size().z() - ADD_NODES);

  for (auto i = 0; i < powerRadCircle2D.size().x() - (ADD_NODES - 2); ++i ){
    radX = 0;
    r.x() = Dx * (i - CELLS_SHIFT + 0.5);
    for ( auto j = 0; j < powerRadCircle2D.size().y(); ++j  ){

      angle = 2 * PI * j / powerRadCircle2D.size().y();

      r.y() = center_y + radiationDiagRadius * cos(angle);
      r.z() = center_z + radiationDiagRadius * sin(angle);

      E = get_fieldE_in_pos(mesh.fieldE,r);
      B = get_fieldB_in_pos(mesh.fieldB,r) - get_fieldB_in_pos(mesh.fieldB0,r);
      circleE(i,j) = E;
      circleB(i,j) = B;
      r.y()-= center_y;
      r.z()-= center_z;
      
      Ep = - r.z() / radiationDiagRadius * E.y() + r.y() / radiationDiagRadius * E.z();
      Bp = - r.z() / radiationDiagRadius * B.y() + r.y() / radiationDiagRadius * B.z();
      
        rad = Dt*(2.0*PI*radiationDiagRadius*Dx/powerRadCircle2D.size().y() ) * fabs(Ep*B.x() - E.x()*Bp );
        powerRadCircle2D(i,j) += rad;
        radX += rad;
    }
    powerRadCircleLine(i) += radX;
    powerRadCircle += radX;
  }

}



void Writer::write_radiation_circle2D(long timestep){
    char filename[100];
    long i = 0;
    for (auto& data : diagData.radialDiag){
      sprintf(filename, ".//Fields//Diag2D//RadCircle%03ld_%03ld",i,timestep / TimeStepDelayDiag2D);
      write_array2D(data.powerRadCircle2D, data.powerRadCircle2D.size().x()-ADD_NODES, data.powerRadCircle2D.size().y(), filename, _world.MPIconf);
      data.powerRadCircle2D.clear();
      i++;
    }
}


void Writer::write_radiation_planes(long timestep){
    char filename[100];
    long i = 0;
    for (auto& data : diagData.powerRadPlaneX){
      sprintf(filename, ".//Fields//Diag2D//RadPlaneX%03ld_%03ld",i,timestep / TimeStepDelayDiag2D);
      write_array2D(data, data.size().x()-ADD_NODES, data.size().y(), filename, _world.MPIconf);
      data.clear();
      i++;
    }

    long j = 0;
    for (auto& data : diagData.powerRadPlaneX){
      sprintf(filename, ".//Fields//Diag2D//RadPlaneY%03ld_%03ld",j,timestep / TimeStepDelayDiag2D);
      write_array2D(data, data.size().x()-ADD_NODES, data.size().y(), filename, _world.MPIconf);
      data.clear();
      j++;
    }

    long k = 0;
    for (auto& data : diagData.powerRadPlaneZ){
      sprintf(filename, ".//Fields//Diag2D//RadPlaneZ%03ld_%03ld",k,timestep / TimeStepDelayDiag2D);
      write_array2D(data, data.size().x()-ADD_NODES, data.size().y(), filename, _world.MPIconf);
      data.clear();
      k++;
    }
}


void Writer::write_radiation_line(long timestep){
    const Region &region = _world.region;
    const MPI_Topology &MPIconf = _world.MPIconf;

    if (!MPIconf.is_master_depth()) return;

    static MPI_File fRadX_min,fRadX_max,fRadY_min,fRadY_max,fRadZ_min,fRadZ_max,fRadC;
    MPI_Status status;
    char filename[100];
    int sizeDataX = (region.numCells.x() - region.dampCells[0].x() - region.dampCells[1].x());
   // int sizeDataY = (region.numCells.y() - 2*region.dampCells[0].y());
    int sizeDataZ = (region.numCells.z() - 2*region.dampCells[0].z());
    static Array1D<float> floatDataX(sizeDataX);
    //static float* floatDataY = new float[sizeDataY];
    static Array1D<float> floatDataZ(sizeDataZ);
    static long currentStep = 0; 
    int startWrite, sizeWrite;
    float info;
    long i,indx; 

    if( timestep == StartTimeStep){
        sprintf(filename, ".//Fields//Diag1D//PRadCircle.bin" );
        MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fRadC);
        info = float(MPIconf.size_line()*sizeDataX);
        if(MPIconf.is_first_line()){   
            MPI_File_write_at(fRadC, 0, &info, 1, MPI_FLOAT, &status);
        }

        sprintf(filename, ".//Fields//Diag1D//PRadY_min.bin" );
        MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fRadY_min);
        info = float(MPIconf.size_line()*sizeDataX);
        if(MPIconf.is_first_line()){   
            MPI_File_write_at(fRadY_min, 0, &info, 1, MPI_FLOAT, &status);
        }
        sprintf(filename, ".//Fields//Diag1D//PRadY_max.bin" );
        MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fRadY_max);
        info = float(MPIconf.size_line()*sizeDataX);
        if(MPIconf.is_first_line()){   
            MPI_File_write_at(fRadY_max, 0, &info, 1, MPI_FLOAT, &status);
        }
        sprintf(filename, ".//Fields//Diag1D//PRadZ_min.bin" );
        MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fRadZ_min);
        info = float(MPIconf.size_line()*sizeDataX);
        if(MPIconf.is_first_line()){   
            MPI_File_write_at(fRadZ_min, 0, &info, 1, MPI_FLOAT, &status);
        }
        sprintf(filename, ".//Fields//Diag1D//PRadZ_max.bin" );
        MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fRadZ_max);
        info = float(MPIconf.size_line()*sizeDataX);
        if(MPIconf.is_first_line()){   
            MPI_File_write_at(fRadZ_max, 0, &info, 1, MPI_FLOAT, &status);
        }

        sprintf(filename, ".//Fields//Diag1D//PRadX_min.bin" );
        MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fRadX_min);
        info = float(sizeDataZ);
        if(MPIconf.is_first_line()){   
            MPI_File_write_at(fRadX_min, 0, &info, 1, MPI_FLOAT, &status);
        }

        sprintf(filename, ".//Fields//Diag1D//PRadX_max.bin" );
        MPI_File_open(MPIconf.comm_line(), filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fRadX_max);
        info = float(sizeDataZ);
        if(MPIconf.is_first_line()){   
            MPI_File_write_at(fRadX_max, 0, &info, 1, MPI_FLOAT, &status);
        }
      MPI_Barrier(MPIconf.comm_line());
    }
    
    indx = region.dampCells[0].x();
    for( i = 0; i < sizeDataX ;i++){
      floatDataX(i) = float(diagData.powerRadLine.y_min(i+indx) );
    }
    startWrite = MPIconf.rank_line()*sizeDataX*sizeof(float) + sizeof(float);
    sizeWrite = MPIconf.size_line()*sizeDataX*sizeof(float);
    startWrite += sizeWrite*currentStep;

    for( i = 0; i < sizeDataX ;i++){

      floatDataX(i) = float(diagData.powerRad1D(i+indx) );
    }
    MPI_File_write_at(fRadC, startWrite, &floatDataX.data(0), sizeDataX, MPI_FLOAT, &status);



    MPI_File_write_at(fRadY_min, startWrite, &floatDataX.data(0), sizeDataX, MPI_FLOAT, &status);
    
    for( i = 0; i < sizeDataX ;i++){

      floatDataX(i) = float(diagData.powerRadLine.y_max(i+indx) );
    }
    MPI_File_write_at(fRadY_max, startWrite, &floatDataX.data(0), sizeDataX, MPI_FLOAT, &status);

    for( i = 0; i < sizeDataX ;i++){
      floatDataX(i) = float(diagData.powerRadLine.z_min(i+indx) );
    }
    MPI_File_write_at(fRadZ_min, startWrite, &floatDataX.data(0), sizeDataX, MPI_FLOAT, &status);

    for( i = 0; i < sizeDataX;i++){
      floatDataX(i) = float(diagData.powerRadLine.z_max(i+indx) );
    }
    MPI_File_write_at(fRadZ_max, startWrite, &floatDataX.data(0), sizeDataX, MPI_FLOAT, &status);



    startWrite = sizeDataZ*sizeof(float)*currentStep;
    indx = region.dampCells[0].z();
    if(region.boundType[0].x()==OPEN){
        for(i = region.dampCells[0].z(); i < sizeDataZ; i++){
            floatDataZ(i) = float(diagData.powerRadLine.x_min(i+indx) );
        }
        MPI_File_write_at(fRadX_min, startWrite + sizeof(float), &floatDataZ.data(0), sizeDataZ, MPI_FLOAT, &status);
    }

    if(region.boundType[1].x()==OPEN){
        for(i = region.dampCells[0].z(); i < sizeDataZ; i++){
            floatDataZ(i) = float(diagData.powerRadLine.x_max(i+indx) );
        }
        MPI_File_write_at(fRadX_max, startWrite + sizeof(float), &floatDataZ.data(0), sizeDataZ, MPI_FLOAT, &status);
    }

    currentStep++;
    
}
*/
