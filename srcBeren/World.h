#ifndef WORLD_H_
#define WORLD_H_
#include "Vec.h"
#include "const.h"
#include "defines.h"
#include <mpi.h>
#include <map>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sys/types.h>
#include <sys/stat.h>
#include <assert.h>

double Uniform01();
void SetRandSeed(int val);
double Gauss(double sigma);

template <typename T> 
int sign(T val) {
    return (T(0) < val) - (val < T(0));
}

struct Region {
    double3 step_h;
    double origin;
    long3 numCells;
    long3 numNodes;
    long3 dampCells[2];
    long3 dampCellsForce[2];
    long3 boundType[2];
    bool in_region(double x) const{
        if ( x < origin || x >= origin + numCells.x() * step_h.x() )
            return false;
  
        return true;
    }
    double3 get_coord_loc(const double3 &POS) const{
        double3 POS_loc = POS;
        POS_loc.x() -= origin;
        return POS_loc;
    }
    double3 get_coord_glob(const double3 &POS) const{
        double3 POS_glob = POS;
        POS_glob.x() += origin;
        return POS_glob;    
    }
    long get_index_loc(long indx) const{
        return indx - round(origin / Dx);
    }
    Region();
    long total_size() const {
        return numNodes.x()*numNodes.y()*numNodes.z(); 
    };
};

Region split_region(const Region& regionGlob, int rank, int splitSize);

class MPI_Topology{
public:
    int size() const{ return _size;}
    int rank() const{ return _rank;}
    int size_depth() const{ return _sizeDepth;}
    int size_line() const{ return _sizeLine;}
    int rank_depth() const{ return _rankDepth;}
    int rank_line() const{ return _rankLine;}
    int next_line() const{
        if( _rankLine != last_line() ) 
            return _rankLine + 1;
        else 
            return 0;
    }
    int prev_line() const{
        if( _rankLine != first_line() ) 
            return _rankLine - 1;
        else 
            return last_line();
    }
    int last_line() const {return _sizeLine - 1;}
    int first_line() const {return 0;}

    bool is_master() const {
        return _rank == 0;
    }
    bool is_master_depth() const{
        return _rankDepth == 0;
    }
    bool is_first_line() const{return rank_line() == first_line();}
    bool is_last_line() const{return rank_line() == last_line();}

    MPI_Comm comm_depth() const{ return _commDepth;}
    MPI_Comm comm_line() const{ return _commLine;}
    void set_topology();
    long accum_sum_line(long value) const;


protected:
    int _size, _rank;
    MPI_Comm _commDepth, _commLine;
    int _sizeDepth, _sizeLine, _rankDepth, _rankLine; 
};

struct World{
    World(const Region& regionGlob,const Region& regionSplit, const MPI_Topology& MPIconf): 
                                regionGlob(regionGlob),region(regionSplit),MPIconf(MPIconf){
    };
    Region regionGlob;
    Region region;
    MPI_Topology MPIconf;
    ~World(){
    }

};

#endif 
