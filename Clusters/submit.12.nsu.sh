if [ $1 -lt 12 ]
then
   nd=1
else
   nd=$[($1+11)/12]
fi
WD=`pwd`
echo '#!/bin/bash
#
#PBS -l walltime=120:00:00
#PBS -l select='$nd':ncpus=12:mpiprocs=12:mem=24492MB
#PBS -N '$2'
#

MPI_NP='$1'
EXE=speed_sts

# cd to work dir
cd $PBS_O_WORKDIR

# start mpi job interconnect=infiniband
#mpirun -IBV -prot -rdma -np $MPI_NP ./$EXE
mpirun -hostfile $PBS_NODEFILE -np $MPI_NP ./'$2'
# Checking all setups:
#echo "mpi-selector set to:"
#mpi-selector --query

#echo "mpirun:"
#which mpirun

#echo "PBS_NODEFILE:"
#cat $PBS_NODEFILE'>&submit.sh
chmod 744 submit.sh
