#include "World.h"
#include "Particles.h"
#include "Mesh.h"
#include "Diagnostic.h"
#include "Read.h"
#include "Vec.h"



int main(int argc,char **argv){
	Region regionGlob;
	MPI_Topology MPIconf;
	
	MPI_Init(&argc,&argv);	

	MPIconf.set_topology();

	Region region = split_region(regionGlob, MPIconf.rank_line(), MPIconf.size_line() );
/// Make modeling area
	World world(regionGlob, region, MPIconf);	
///	Make mesh
	Mesh mesh(world);	

	std::vector< std::vector<std::string> > stringParams;
///// Make particles
	std::vector<ParticlesArray> species;
	read_params_to_string("Particles","./PartParams.cfg",stringParams);
	for( const auto &params  : stringParams){
	    species.emplace_back(params,world);
	}

	Writer writer(world,mesh,species);

	writer.output(StartTimeStep);
	
	MPI_Barrier(MPI_COMM_WORLD);

	Timer globalTimer("globalFunctions.time");
	for(auto timestep = StartTimeStep + 1; timestep <= MaxTimeStep; ++timestep){

		mesh.prepare();

		for( auto &sp  : species){
			sp.move(0.5*Dt); //  +++ x_n -> x_{n+1/2}
			sp.get_current_predict(mesh); // +++ get J(x_{n+1/2},v_n)_predict 
			sp.get_L(mesh); // --- get Lgg'(x_{n+1/2}) 
		}
		
		mesh.solveE_predict(); // --- solve A*E'_{n+1}=f(E_n, B_n, J(x_{n+1/2})). mesh consist En, En+1_predict
		
		for( auto &sp  : species){

			sp.get_velocity(mesh); // +++ get v'_{n+1} from v_{n} and E'_{n+1}
			sp.move(0.5*Dt); // +++ x_{n+1/2} -> x_{n+1}
			sp.get_esirkepov_current(mesh); // ++++ get J_e(x_n,x_{n+1}) from Esirkepov
		}
		
		mesh.correctE(); // ---- get E_{n+1} from E'_{n+1}, J_e and J_predict. mesh En changed to En+1_final 


		for( auto &sp  : species){
			sp.correctv(mesh); // ++++ get v_{n+1} final from v'_{n+1}
		}
		//mesh.clear();

  		//globalTimer.start("Collision");
		//collision(mesh,world,species,timestep);
		//MPI_Barrier(MPI_COMM_WORLD);
		//globalTimer.finish("Collision");

		// globalTimer.start("particles");
		// for( auto &sp  : species){
		// 	sp.update(mesh,timestep);
		// }

		// globalTimer.finish("particles");
			
		// globalTimer.start("ReduceCurrent");
		// //mesh.reduce_current(MPIconf);
		// globalTimer.finish("ReduceCurrent");

  // 		globalTimer.start("Fields");
		// mesh.update(timestep);
		// globalTimer.finish("Fields");

  // 		globalTimer.start("Output");
		// writer.output(timestep);
		// globalTimer.finish("Output");

		// globalTimer.write(timestep, MPIconf);

	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	
	return 0;
}

