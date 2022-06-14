#!/bin/bash
if [ $1 -lt 32 ]
then
   nd=1
else
   nd=$[($1+31)/32]
fi


echo '#!/bin/bash
# set the number of nodes
#SBATCH --nodes='$nd'
## hyperthreading off
#SBATCH --threads-per-core=1
# set max wallclock time
#SBATCH --time=9999
# set name of job
#SBATCH --job-name=Nerpa
# set queue name
#SBATCH -p broadwell
#SBATCH --reservation=bdw
##SBATCH --ntasks-per-node=1
# run the application
mpirun -n='$1' ./'$2'' > start.sh
