#ifndef TIMER_H_
#define TIMER_H_
#include "World.h"

struct Timer{
    FILE *fTimes;
    std::map<std::string,double> times;
    std::map<std::string,double> startT;
    bool firstWrite;
    Timer(const std::string& filename){
        fTimes = fopen( (".//Performance//"+filename).c_str(), "w");
        _reset = true;
        firstWrite = true;
    }

    void start(const std::string &s){
        startT[s] = MPI_Wtime();
    }
    void finish(const std::string &s){

        times[s] += MPI_Wtime() - startT[s];
    }
    void reset(){
        for (auto it = times.begin(); it != times.end(); ++it){
            it->second = 0.;
        }
        _reset = true;
    }
    void write(long timestep, const MPI_Topology& MPIconf); 
protected:
    bool _reset;
};

#endif 	
