source /opt/intel/oneapi/mkl/latest/env/vars.sh
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc/25.5
export CUDA_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/cuda
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
rm -rf build
mkdir build
cd build
cmake -DOPENMP=OFF -DMKL=ON -DACC=ON -DCMAKE_INSTALL_PREFIX=/home/shuo/bin/gronor-25.06_gn2 ..
make 
make install
