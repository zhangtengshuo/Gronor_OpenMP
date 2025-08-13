# HDF5
export HDF5_ROOT=/home/shuo/bin/hdf5_1.14.6
# Intel MKL
source /opt/intel/oneapi/mkl/latest/env/vars.sh
# NVHPC
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc/25.5
export CUDA_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/cuda
export CC=mpicc
export CXX=mpicxx
export FC=mpif90

cd /home/shuo/develop/gronor.26.06.dev/gronor-25.06_openmp_shuo
rm -rf build && mkdir build && cd build
cmake -DOPENMP=OFF -DMKL=ON -DACC=ON -CMAKE_BUILD_TYPE=Debug -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_INSTALL_PREFIX=/home/shuo/bin/gronor-25.06_HDF5 ..
make -j8
make install
