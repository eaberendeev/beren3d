#ifndef MESH_H_
#define MESH_H_
#include "World.h"
#include "Read.h"
//#include "Laser.h"
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

double calc_energy_field(const Field3d& field);
double3 get_fieldE_in_pos(const Field3d& fieldE, const double3& r);
double3 get_fieldB_in_pos(const Field3d& fieldB, const double3& r);
void exchange_fieldsE(Field3d& field,const MPI_Topology &MPIconf);

struct Mesh{
    Mesh(const World& world);

    Operator Lmat;

    //Sources and fields on the grid
    Field3d chargeDensity;

    Field3d fieldE;
    Field3d fieldEn;
    Field3d fieldB;
    Field3d fieldJ;

    Field3d fieldB0;
    int vind(long i, long j, long k, long d) const {
        return d + _nd*(i * _size2 * _size3 + j * _size3 + k);
    };
    void set_fields();

    void set_uniform_fields();
    void read_from_recovery(const MPI_Topology &MPIconf);
    void write_recovery(const MPI_Topology &MPIconf) const;
    void update(long timestep);

    double3 get_fieldE_in_cell(long i, long j, long k)  const;
    double3 get_fieldB_in_cell(long i, long j, long k)  const;
    double get_fieldE_energy() const{
        return calc_energy_field(fieldE);
    };
    void laser_source(long timestep);
    double get_fieldB_energy() const{
        return calc_energy_field(fieldB)-calc_energy_field(fieldB0);
    };
    void reduce_current(const MPI_Topology &MPIconf);
    ~Mesh(){
    }
    void clear(){
        fieldJ = 0.;
    }
    double get_coord_from_nodeX(long index) const{
        return (index - CELLS_SHIFT) * Dx;
    }
    double get_coord_from_nodeY(long index) const{
        return (index - CELLS_SHIFT) * Dy;
    }
    double get_coord_from_nodeZ(long index) const{
        return (index - CELLS_SHIFT) * Dz;
    }
    long get_node_from_coordX(double coord) const{
        return long(coord / Dx + CELLS_SHIFT); 
    }
    long get_node_from_coordY(double coord) const{
        return long(coord / Dy + CELLS_SHIFT); 
    }
    long get_node_from_coordZ(double coord) const{
        return long(coord / Dz + CELLS_SHIFT); 
    }
protected:
    const World &_world;
    long _size1, _size2, _size3, _nd;
};


#endif 
