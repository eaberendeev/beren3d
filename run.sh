#!/bin/bash
rm *.tmp
rm *.cfg
python set_params.py

echo "......Reading of auxilary variables......"
name=Beren3D
WorkDir=$(cat workdir.tmp)      # WORK DIRECTORY
proc=$(cat proc.tmp)      # Number of procs
queue=$(cat queue.tmp)      # queue type
clu=$(cat cluster.tmp)      # cluster

if [[ -z $WorkDir ]]
then
  exit
fi


echo "......Copy files to work directory......"   
   rm *.tmp
   rm -r $WorkDir
   mkdir $WorkDir
   cp -r PlotScripts ./$WorkDir
   cp -r srcBeren ./$WorkDir
   cp -r Scripts ./$WorkDir
   cp -r Clusters ./$WorkDir
   cp -r Recovery ./$WorkDir
   mv *.h ./$WorkDir/srcBeren
   mv *.par ./$WorkDir/srcBeren
   cp *.py ./$WorkDir
   mv *.cfg ./$WorkDir
   cp *.sh ./$WorkDir
   cd ./$WorkDir

echo "......Compile......"
cd srcBeren

if [[ $clu = nks1p ]]
then
module purge
module load intel/2017.4.196 parallel/mpi.intel.broadwell/2017.4.196 compilers/intel/2017.4.196
fi 

make -f Makefile_cpu clean
make -f Makefile_cpu
echo "......End compile......"
ls | grep $name
cp $name ../
cd ..

if [[ $queue = home ]]
then

#export OMP_NUM_THREADS=$nt
mpirun -np $proc ./$name

else

./Clusters/submit."$queue"."$clu".sh  $proc $name
  if [[ $clu = nks1p ]]
  then
    sbatch ./start.sh
  fi
#qsub submit.sh
fi

