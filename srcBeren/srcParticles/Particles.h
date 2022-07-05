#ifndef PARTICLES_H_
#define PARTICLES_H_
#include "World.h"
#include "Vec.h"
#include "Mesh.h"
#include <functional>
#include <assert.h>

typedef Eigen::Triplet<double> Trip;

struct ParticleSimple{
	double3 coord;
    double3 pulse;
    double3 velocity;
    double3 startCoord;

	friend std::ostream& operator<<(std::ostream& out, const ParticleSimple &particle);
    
    void set_global(const Region& domain){
        coord.x() += domain.origin;
    }
    void set_local(const Region& domain){
        coord.x() -= domain.origin;      
    }
    void move(double dt){
        coord += velocity * dt;
    }
};

struct ParticleMPW : ParticleSimple{
    double mpw;
    friend std::ostream& operator<<(std::ostream& out, const ParticleMPW &particle);

};
struct ParticleMass : ParticleSimple{
	double mass;
  	friend std::ostream& operator<<(std::ostream& out, const ParticleMass &particle);
};

struct ParticleMassMPW : ParticleSimple{
    double mass,mpw;
    friend std::ostream& operator<<(std::ostream& out, const ParticleMassMPW &particle);
};

#ifdef PARTICLE_MASS
    #ifdef PARTICLE_MPW
        typedef ParticleMassMPW Particle;
    #else
        typedef ParticleMass Particle;
    #endif
#else
    #ifdef PARTICLE_MPW
        typedef ParticleMPW Particle;
    #else
        typedef ParticleSimple Particle;
    #endif
#endif


struct ParticlesOption{
    long boundResumption;
    long sourceType;
    long smoothMass;
    double smoothMassSize;
    double smoothMassMax;
    double sourceAngle; 
};


class ParticlesArray{

public:
    Array<Particle> particlesData;

    Array3D<double> densityOnGrid;
    Array2D<double> phaseOnGrid;
    
    long charge;
    double density;
    double phasePXmin, phasePXmax;
    double pot_I;
    double pot_k;
    double kineticEnergy;
    double injectionEnergy;
    std::string name;
    double temperature;
    double velocity;
    double widthY,widthZ;
    double focus;
    static long counter;
    ParticlesOption option;
    std::string initDist;
    std::vector<double> distParams;
    void add_particle_scatter(const Particle& particle){
        if(counter % _world.MPIconf.size_depth() == _world.MPIconf.rank_depth() )
            particlesData.push_back(particle);
        ++counter;
    }

    void set_distribution_density(std::function<double(double3 )> set_density);
    void set_smooth_mass();
    ParticlesArray(const std::vector<std::string>& vecStringParams, World& world);
    void set_params_from_string(const std::string& line);
    void density_on_grid_update();
    void density_on_grid_update_pic();
    void phase_on_grid_update();
    void set_distribution();
    void inject(long timestep);
    void update(Mesh& mesh,long timestep);
    double mass() const{
          return _mass;
    }
    double mass(long k) const{
        #ifdef PARTICLE_MASS
          return particlesData(k).mass;
        #else
          return _mass;
        #endif
    }
    double mpw(long k) const{
        #ifdef PARTICLE_MPW
          return particlesData(k).mpw;
        #else
          return _mpw;
        #endif
    }    
    //void read_recovery_particles(const MPI_Topology& MPIconf);
    Particle& operator() (long i) {
        return particlesData(i);
    }

    const Particle& operator() (long i) const{
        return particlesData(i);
    }
    long size() const{
        return particlesData.size();
    }
    double get_kinetic_energy() const;
    double get_inject_energy() const{
        return injectionEnergy;
    }
    void move(Mesh& mesh,long timestep);
    void move(double dt);
    void correctv(const Mesh& mesh);
    void get_velocity(const Mesh& mesh);
    void get_esirkepov_current( Mesh& mesh);
    void get_current_predict( Mesh& mesh);
    void get_L( Mesh& mesh);



    void move_virt(Mesh& mesh,long timestep);
    void bound_resumption(const Particle& particle, const double3& r_new, const double3& p_new);

    void set_space_distribution();
    void set_pulse_distribution();
    void set_uniform_circle(long3 start, long3 end);
    void set_strict_uniform(long3 start, long3 end);
    void set_uniform(long3 start, long3 end);
    void set_strict_cosX(long3 start, long3 end, double delta_n0, double k0);
    void source_uniform_from_bound(long timestep);
    void source_focused_gauss(long timestep);

    void read_from_recovery(const MPI_Topology &MPIconf);
    void write_recovery(const MPI_Topology &MPIconf) const;
protected:
    World &_world;
    double _mass;
    double _mpw; /*macroparticle weight*/
};



double PulseFromKev(double kev, double mass);

int get_num_of_type_particles(const std::vector<ParticlesArray> &Particles, const std::string& ParticlesType);
/// Ionization particles = electron(particles_e) +  ion (particles_i) 
void collision(const Mesh &mesh, const World& world ,std::vector<ParticlesArray> &Particles,long timestep);
//double getFieldEInLaser(double z, double r, const Region& domain, long timestep);

#endif 
