#!/bin/bash
#
#PBS -l walltime=96:00:00
#PBS -l nodes=0:ppn=32
#PBS -N 
#

MPI_NP=

# cd to work dir
cd $PBS_O_WORKDIR

/share/apps/bin/mpiexec ./

#cat $PBS_NODEFILE
