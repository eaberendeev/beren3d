#include "Vec.h"
#include "World.h"
#include "Particles.h"
#include "Shape.h"
std::ostream& operator<<(std::ostream& out, const ParticleSimple& particle){
	out << particle.coord.x() << " " << particle.coord.y() << " "<< particle.coord.z() << " " << particle.pulse.x() << " " << particle.pulse.y() << " " << particle.pulse.z();
	return out;
} 
std::ostream& operator<<(std::ostream& out, const ParticleMass& particle){
	out << particle.coord.x() << " " << particle.coord.y() << " " << particle.coord.z() << " "<< particle.pulse.x() << " " << particle.pulse.y() << " " << particle.pulse.z() << " " << particle.mass;
	return out;
} 


std::ostream& operator<<(std::ostream& out,const double3& val){
	  out << " x " << val.x() << " y " << val.y() << " z " << val.z();
	  return out;
} 

long ParticlesArray::counter;
ParticlesArray::ParticlesArray(const std::vector<std::string>& vecStringParams, World& world):
    particlesData(MaxSizeOfParts),
    densityOnGrid(world.region.numNodes),phaseOnGrid(PxMax,PpMax),
    _world(world) 
    {
    for (const auto& line: vecStringParams){
        set_params_from_string(line);
    }
    injectionEnergy = 0.;

    if (RECOVERY){ 
        read_from_recovery(world.MPIconf);
    }else{
        set_distribution();
    }
    density_on_grid_update();
    phase_on_grid_update(); 
      std::cout << particlesData.size() << " " << particlesData.capacity()  <<"\n";

}

double PulseFromKev(double kev, double mass){
  double gama = kev / MC2 + mass;
  return sqrt((gama*gama)- mass);
}


int get_num_of_type_particles(const std::vector<ParticlesArray> &Particles, const std::string& ParticlesType){
  for (int i = 0; i < NumOfPartSpecies; ++i){
    if(Particles[i].name == ParticlesType) return i;
  }
  return -1;
}

void join_density_comm_line(Array3D<double>& dens, const MPI_Topology& MPIconf){

  int tag = 1;
  MPI_Status status;
  auto size = dens.size();
  int bufSize = ADD_NODES*size.y()*size.z();
  long sendIndx; 
  static double *recvBuf = new double[bufSize];
  
  auto right = MPIconf.next_line();
  
  auto left = MPIconf.prev_line();


  sendIndx = (size.x()-ADD_NODES)*size.y()*size.z();
  // fromRight to Left
  MPI_Sendrecv(&dens.data(sendIndx), bufSize, MPI_DOUBLE_PRECISION, right, tag, 
                           recvBuf, bufSize, MPI_DOUBLE_PRECISION, left, tag, 
                           MPIconf.comm_line(), &status);
  
  for (long i = 0; i< bufSize; i++){
    dens.data(i)+= recvBuf[i];
  }
  
  MPI_Sendrecv(&dens.data(0), bufSize, MPI_DOUBLE_PRECISION,left , tag, 
                           &dens.data(sendIndx), bufSize, MPI_DOUBLE_PRECISION, right, tag, 
                           MPIconf.comm_line(), &status);
  
  MPI_Barrier(MPIconf.comm_line());
  
}
void ParticlesArray::density_on_grid_update_pic(){ 
    long yk, xk, zk, i,l,j,k;
    alignas(64) double sx[2];
    alignas(64) double sy[2];
    alignas(64) double sz[2];
    Particle particle;
    long indx, indy,indz;
    long3 size = densityOnGrid.size();

    densityOnGrid.clear();

    for ( j = 0; j < particlesData.size(); ++j){
        particle = particlesData(j);
        xk = particle.coord.x() / Dx;
        yk = particle.coord.y() / Dy;
        zk = particle.coord.z() / Dz;
        indx = long(xk);
        indy = long(yk);
        indz = long(zk);
        sx[0] = (xk - indx);
        sy[0] = (yk - indy);
        sz[0] = (zk - indz);
        sx[1] = 1. - sx[0];
        sy[1] = 1. - sy[0];
        sz[1] = 1. - sz[0];
        for (i = 0; i < 2; ++i){
          for (l = 0; l < 2; ++l){
            for (k = 0; k < 2; ++k){
              densityOnGrid(indx+i+1,indy+l+1,indz+k+1) += mpw(j) * sx[i]*sy[l]*sz[k] / (Dx*Dy*Dz);
            }
          }
        }

    }

    MPI_Allreduce ( MPI_IN_PLACE, &densityOnGrid.data(0), size.x()*size.y()*size.z(), MPI_DOUBLE, MPI_SUM, _world.MPIconf.comm_depth());
    MPI_Barrier(MPI_COMM_WORLD);

    join_density_comm_line(densityOnGrid, _world.MPIconf);
}


void ParticlesArray::density_on_grid_update(){ 
    constexpr auto SMAX = 2*SHAPE_SIZE;
    long yk, xk, zk, k,n,m;
    double sx[SMAX];
    double sy[SMAX];
    double sz[SMAX];
    double arg;
    Particle particle;
    long indx, indy,indz;
    long3 size = densityOnGrid.size();

    
    densityOnGrid.clear();


    for ( auto j = 0; j < particlesData.size(); ++j){
        particle = particlesData(j);
        xk = long(particle.coord.x() / Dx);
        yk = long(particle.coord.y() / Dy);
        zk = long(particle.coord.z() / Dz);

        for(n = 0; n < SMAX; ++n){
            arg = -particle.coord.y() / Dy + double(yk - CELLS_SHIFT + n);
            sy[n] = Shape(arg) ;
            arg = -particle.coord.x() / Dx + double(xk - CELLS_SHIFT + n);
            sx[n] = Shape(arg);
            arg = -particle.coord.z() / Dz + double(zk - CELLS_SHIFT + n);
            sz[n] = Shape(arg);
        }

        for(n = 0; n < SMAX; ++n){
            indx = xk + n;
            for(m = 0; m < SMAX; ++m){
                indy = yk - CELLS_SHIFT + m;
                for(k = 0; k < SMAX; ++k){
                  indz = zk - CELLS_SHIFT + k; 
                  densityOnGrid(indx,indy,indz) += mpw(j) * sx[n] * sy[m] * sz[k] / (Dx*Dy*Dz);
                } 
            }
        }
    }

    MPI_Allreduce ( MPI_IN_PLACE, &densityOnGrid.data(0), size.x()*size.y()*size.z(), MPI_DOUBLE, MPI_SUM, _world.MPIconf.comm_depth());
    MPI_Barrier(MPI_COMM_WORLD);

    join_density_comm_line(densityOnGrid, _world.MPIconf);
}

void ParticlesArray::phase_on_grid_update(){ 
    auto size = phaseOnGrid.size();
    long pk, xk;
    double x, px;
    bool blounder;
    Particle particle;
    double x_min = Dx*(_world.region.dampCells[0].x());
    double x_max = Dx*(_world.region.numCells.x() - _world.region.dampCells[1].x() );
    double pdx = (x_max - x_min) / PxMax;
    double pdp = (phasePXmax -phasePXmin) / PpMax;

    phaseOnGrid.clear();

    for (auto k = 0 ; k < particlesData.size(); k++){
        particle = particlesData(k);   
        x = particle.coord.x();
        px = particle.pulse.x();

        xk = long((x-x_min) / pdx); 
        pk = long((px - phasePXmin) / pdp);

        blounder = (xk<0) || (xk>=PxMax) || (pk<0) || (pk>=PpMax);
        if(!blounder){
            phaseOnGrid(xk,pk) += (mpw(k) / (Dx*Dy*pdx*pdp) );
        }

    }

    MPI_Allreduce ( MPI_IN_PLACE, &phaseOnGrid.data(0), size.x()*size.y(), MPI_DOUBLE, MPI_SUM, _world.MPIconf.comm_depth());
    MPI_Barrier(MPI_COMM_WORLD);

}


void ParticlesArray::set_distribution(){
		set_space_distribution();
		set_pulse_distribution();	
}

// Устанавливаем значения основных параметров в переменной Params в соответствии с содержимым сроки
void ParticlesArray::set_params_from_string(const std::string& line){
    std::vector<std::string> strvec;

    strvec = split(line, ' ');

    if(strvec[0]=="Particles") {
        this->name = strvec[1];
    }

    if(strvec[0]=="Temperature"){
        temperature = stod(strvec[1]);
    }

    if(strvec[0]=="Mass"){
        this->_mass = stod(strvec[1]);
    }

    if(strvec[0]=="Charge"){
       this->charge =  stol(strvec[1]);
    }
    if(strvec[0]=="SmoothMass"){
       option.smoothMass =  stol(strvec[1]);
    }
    if(strvec[0]=="SmoothMassSize"){
       option.smoothMassSize =  stod(strvec[1]);
    }
    if(strvec[0]=="SmoothMassMax"){
       option.smoothMassMax =  stod(strvec[1]);
    }  
    if(strvec[0]=="WidthY"){
        widthY = stod(strvec[1]);
    }
    if(strvec[0]=="WidthZ"){
        widthZ = stod(strvec[1]);
    }
    if(strvec[0]=="Density"){
        this->density = stod(strvec[1]); 
        this->_mpw = density * Dx * Dy * Dz / double(NumPartPerCell);
    }
    if(strvec[0]=="Focus"){
       focus = stod(strvec[1]);
    }
    if(strvec[0]=="Velocity"){
        velocity = stod(strvec[1]);
    }
    if(strvec[0]=="Pot_I"){
        this->pot_I = stod(strvec[1]);
    }
    if(strvec[0]=="Pot_k"){
        this->pot_k = stod(strvec[1]);
    }

    if(strvec[0]=="Px_min"){
        this->phasePXmin = stod(strvec[1]);
    }

    if(strvec[0]=="Px_max"){
         this->phasePXmax = stod(strvec[1]);
    }
    if (strvec[0]=="DistParams"){
        initDist = strvec[1];
        if(strvec[1]=="UniformCosX_dn_k"){
            distParams.push_back(stod(strvec[2]));
            distParams.push_back(stod(strvec[3]));
        }
    }
    if(strvec[0]=="BoundResumption"){
        option.boundResumption = stod(strvec[1]);
    }
    if(strvec[0]=="SourceType"){
        option.sourceType = stod(strvec[1]);
    }
    if(strvec[0]=="sourceAngle"){
        option.sourceAngle = stod(strvec[1]);
    }
}



void ParticlesArray::update(Mesh& mesh,long timestep){
  try {
    inject(timestep);
  }
  catch (std:: string& msg){
    std::cout << "InjectError\n";
  }

  move(mesh,timestep);
  move_virt(mesh,timestep);

  if (timestep % TimeStepDelayDiag2D == 0){
        density_on_grid_update();
        phase_on_grid_update();
  }
      
}
void ParticlesArray::inject(long timestep){
    switch(option.sourceType){
       case SOURCE_UNIFORM:
            source_uniform_from_bound(timestep);
            break;
        case SOURCE_FOCUSED_GAUSS:
            source_focused_gauss(timestep);
            break;
        case SOURCE_NONE:
            break;
        default:
            break;

        }
}

double ParticlesArray::get_kinetic_energy() const{
    double energy = 0;
    double3 pulse;
    double gama;

    for(auto k = 0; k < particlesData.size(); ++k){
        pulse = particlesData(k).pulse;
        gama = sqrt(mass(k)*mass(k) + dot(pulse,pulse) );
        energy += mpw(k) * (gama - mass(k));
    }   

    return energy;
}