#!/bin/bash
#PBS -l walltime=2:00:00
#PBS -l select=1:ncpus=1
#PBS -N taurus
#

MPI_NP=40
EXE=speed_sts

# cd to work dir
cd $PBS_O_WORKDIR

# start mpi job interconnect=infiniband
#mpirun -IBV -prot -rdma -np $MPI_NP ./$EXE
#mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./Nerpa
./VideoPlot.sh
# Checking all setups:
#echo "mpi-selector set to:"
#mpi-selector --query

#echo "mpirun:"
#which mpirun

#echo "PBS_NODEFILE:"
#cat $PBS_NODEFILE
