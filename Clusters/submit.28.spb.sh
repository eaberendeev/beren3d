if [ $1 -lt 28 ]
then
   nd=1
else
   nd=$[($1+27)/28]
fi
WD=`pwd`
echo "#!/bin/bash
#
#SBATCH --nodes=$nd
#SBATCH --tasks-per-node=28
###SBATCH --cpus-per-task=1
#SBATCH -p tornado
#SBATCH -t 10-00:00:00
#SBATCH -J $2
#SBATCH -o $2-%j.out
#SBATCH -e $2-%j.err
if [ -f /etc/profile.d/modules-basis.sh ]; then
source /etc/profile.d/modules-basis.sh
fi 

mpirun -np $1 ./$2 " > submit.sh
chmod 744 submit.sh
