if [ $1 -lt 24 ]
then
   nd=1
else
   nd=$[($1+23)/24]
fi
WD=`pwd`
echo '#!/bin/bash
#
#PBS -l walltime=120:00:00
#PBS -l select='$nd':ncpus=24:mpiprocs=24:mem=180000m
#PBS -N '$2'
#PBS -q xl230g9q
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
