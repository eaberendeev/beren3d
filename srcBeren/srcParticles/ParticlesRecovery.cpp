#include "Diagnostic.h"


void ParticlesArray::read_from_recovery(const MPI_Topology &MPIconf){
  	FILE *fParticles;
	long info;
	Particle particle;
  	char filename[100];
	
  	sprintf(filename, (".//Recovery//Particles//"+name+"%03ld"+".rec").c_str(),MPIconf.rank());
	
	fParticles = fopen(filename, "rb");
	size_t result = fread(&info, sizeof(long), 1, fParticles);
  	assert(result == sizeof(long) && "Number of reading byte is incorrect");

	std::cout << "Raed " << info << " particles "<< name << std::endl;

	while (info>0){	
    result = fread(&particle, sizeof(Particle), 1, fParticles);
    assert(result == sizeof(Particle) && "Number of reading byte is incorrect");
		particlesData.push_back(particle);
		--info;
	}
	fclose(fParticles);

}
void ParticlesArray::write_recovery(const MPI_Topology &MPIconf ) const {
  	FILE *fParticles;
  	char filename[100];
	sprintf(filename, (".//Recovery//Particles//" + name + "%03ld").c_str(),MPIconf.rank());
	
	fParticles = fopen(filename, "wb");
	
	long info = size();
	fwrite(&info, sizeof(long), 1, fParticles);
	fwrite(&particlesData(0), sizeof(Particle)*info, 1, fParticles);
	fclose(fParticles);

}
