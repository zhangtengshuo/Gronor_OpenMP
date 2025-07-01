#!/bin/bash

module load /opt/nvidia/hpc_sdk/modulefiles/nvhpc/25.5
export CUDA_HOME=/opt/nvidia/hpc_sdk/Linux_x86_64/2025/cuda
export CC=nvc
export CXX=nvc++
export FC=nvfortran