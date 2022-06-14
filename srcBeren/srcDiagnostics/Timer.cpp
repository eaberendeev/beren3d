#include "Timer.h"

void Timer::write(long timestep, const MPI_Topology& MPIconf){
  if( times.empty() || timestep%TimeStepDelayDiag1D !=0) return;
    std::stringstream ss;

    if( firstWrite){
      ss<< "Time ";
        for (auto it = times.begin(); it != times.end(); ++it){
          ss << it->first << " ";
        }
        firstWrite = false;
         
    }
    ss << "\n" << timestep << " ";

    for (auto it = times.begin(); it != times.end(); ++it){
          ss << it->second << " ";
    } 
    
  if( MPIconf.is_master() ) 
    fprintf(fTimes, "%s",  ( ss.str() ).c_str() ); 
    
    reset();

  if( MPIconf.is_master()  && timestep % TimeStepDelayDiag1D == 0 )
    fflush(fTimes);
  
}