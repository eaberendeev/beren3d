if [ $1 -lt 36 ]
then
   nd=1
else
   nd=$[($1+35)/36]
fi
WD=`pwd`
echo '#!/bin/bash
#
#PBS -l walltime=96:00:00
#PBS -l nodes='$nd':ppn=36
#PBS -N '$2'
#

MPI_NP='$1'

# cd to work dir
cd $PBS_O_WORKDIR

/share/apps/bin/mpiexec ./'$2'

#cat $PBS_NODEFILE'>&submit.sh
chmod 744 submit.sh
