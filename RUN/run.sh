#PBS -l nodes=np_nodes:ppn=np_threads
#PBS -N expname
#PBS -o expname.out
#PBS -e expname.err

#NPROCS=`wc -l < $PBS_NODEFILE`

source $MODULESHOME/init/sh
module purge
module load intel14/comp intel14/hdf5-1.8.12 intel14/openmpi-1.8 intel14/netcdf-4.3.2 intel14/petsc-3.6.2 intel14/pnetcdf-1.3.1

cd $PBS_O_WORKDIR

export EXPHDR_tmp='expname ../../DATA/expname'

date > runtime

mpirun ./vvm < INPUT | tee OUTPUT

date >> runtime


