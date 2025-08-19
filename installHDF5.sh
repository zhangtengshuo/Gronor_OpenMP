# Intel MKL
source /opt/intel/oneapi/mkl/latest/env/vars.sh
# NVHPC
module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc/25.5
export CUDA_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/25.5/cuda
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
# disable conda
if command -v conda &> /dev/null; then
    if [ -n "$CONDA_DEFAULT_ENV" ]; then 
    conda deactivate
    fi
fi

rm -rf build && mkdir build && cd build
cmake -DOPENMP=OFF -DMKL=ON -DACC=ON -DCMAKE_BUILD_TYPE=Debug \
      -DHDF5_ROOT=/home/shuo/bin/hdf5_1.14.6 \
      -DCMAKE_INSTALL_PREFIX=/home/shuo/bin/gronor-25.06_HDF5 ..
make 
make install