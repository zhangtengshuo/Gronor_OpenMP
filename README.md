# GronOR Non-Orthogonal Configuration Interaction

This is a modified version of [GronOR](https://gitlab.com/gronor/gronor).  
The original version was primarily optimized for supercomputing clusters and did not fully utilize GPU capabilities. The current modifications adapt the parallel framework for smaller-scale computations (which might seem somewhat unconventional).

## Current Modifications and Future Goals

1. **GPU subroutine optimization**:  
   - Already implemented optimization for `cofac1`. Achieved 10x speedup with OpenACC (excluding SVD/EVD operations).  
   - `gntwo` appears to have limited optimization potential under OpenACC - the current batch approach may not be optimal.  

2. **SVD/EVD stays on CPU**:  
   - The computation scale is too small to benefit significantly from GPU acceleration.  
   - All GPU solvers have been removed - these operations will now run on CPU.  

3. **OpenMP integration**:  
   - **Multithreading now handles work rank tasks.**  
   - Yes, OpenMP has been reintroduced (after being removed previously).  
   - Generates sufficient SVD/EVD computation data to keep the GPU busy with `gntwo` calculations.  
   - Addresses CPU/GPU idle time due to mutual waiting.  

4. **Parallel framework**:  
   - I'm being stubborn here - using a hybrid `MPI` + `OpenMP` + `OpenACC` approach (determined to make it work before refactoring).  
   - Future optimization direction: `MPI` + `OpenMP` + `OpenMP target offload`.

## Authors 
     T. P. Straatsma, Oak Ridge National Laboratory, Oak Ridge, TN  
	 C. de Graaf, University Rovira i Virgili, Tarragona, Spain  
	 A. Sanchez, University Rovira i Virgili, Tarragona, Spain  
	 R. K. Kathir, University of Groningen, Groningen, Netherlands  

## Modifier
	 Teng-Shuo Zhang, Zhejiang University of Technology, HangZhou, China  

## Reference
    T. P. Straatsma, R. Broer, A. Sanchez-Mansilla, C. Sousa, and C. de Graaf,   
	“GronOR: Scalable and Accelerated Non-Orthogonal Configuration Interaction for Molecular fragment Wave Functions”,   
	Journal of Chemical Theory and Computation, 18, 3549-3565 (2022).  

## Install
### Download the origin GronOR from the GitLab repository:  
```bash
git clone --recursive git@gitlab.com:gronor/gronor.git
#or
git clone --recursive https://www.gitlab.com/gronor/gronor.git
```
The resulting master branch is the most recent release of gronor. To use an earlier release (e.g. version 23.08) use  
```bash
git checkout tags/23.08
```

### Downloading the modified version from the Github repository:  
```bash
git clone https://github.com/zhangtengshuo/Gronor_OpenMP.git
```

The initial directory structure is as follows:
```plain
gronor 
 ├─ src (source directory with sub-directories)
 ├─ aux (auxiliary programs)
 ├─ include (with a few include files)
 ├─ examples (example input files for full OpenMolcas/GronOR runs)
 ├─ scripts (scripts directory)
 ├─ molcas_interface (auxiliary programs to interface with OpenMolcas)
 ├─ CMakeLists.txt (the cmake build file)
 └─ CTestConfig.cmake (cmake script to setup automated testing)
```

To build, do the following within the gronor directory;
```bash
mkdir build
cd build
cmake [flags] ../
make -j 10
```
This will expand the directory structure as follows:
```plain
gronor 
 ├─ src (source directory with sub-directories)
 ├─ aux (auxiliary programs)
 ├─ include (with a few include files)
 ├─ examples (example input files for full OpenMolcas/GronOR runs)
 ├─ scripts (scripts directory)
 ├─ molcas_interface (auxiliary programs to interface with OpenMolcas)
 ├─ CMakeLists.txt (the cmake build file)
 ├─ CTestConfig.cmake (cmake script to setup automated testing)
 └─ build 
      ├─ bin     (directory with the gronor binary)
      ├─ lib   	  (directory with the gronor libraries)
      └─ CMakeFiles  (cmake files created during build)	  
```

The following rules need to be followed for the src subdirectories:

1. programs in a single source file in a subdirectory with the same name
2. library files in a subdirectory will be in a single library


Build flags are the following:

`-DMPI=ON` is default as GronOR requires a minimum of 2 MPI ranks to run  
`-DACC=ON` will activate OpenACC directives to be interpreted for PGI compilers   
`-DOMP=ON` will activate OpenMP directives to be interpreted for most compilers   
`-DOMPTGT=ON` will activate OpenMP target directives to be interpreted for most compilers  
`-DCUDA=ON` will activate basic access required for CUSOLVER and/or CUSOLVERJ   
`-DCUSOLVER=ON` will activate the CUSOLVER library, and requires `-DOPENACC=ON`  
`-DCUSOLVERJ=ON` will activate the CUSOLVER library including the iterative Jacobi solvers, and requires `-DOPENACC=ON`  
`-DMKL=ON` will link against MKL libraries for compatible compilers  
`-DPROFILING=ON` will activate compiler flags enabling profiling tools (e.g. symbol tables)  

GronOR is interfaced with OpenMolcas for integrals and CASSCF orbital coefficients.
Some minor changes in the OpenMolcas source code are required for GronOR to properly 
function per the instructions in the file OpenMolcas-GronOR.pdf.

GronOR is made available as open source software under the Apache License Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0) 
and any use of the software has to be in compliance with this license. Unless required by applicable law or agreed to in writing, 
software distributed under the license is distributed on an ‘as is’ bases, without warranties or conditions of any kind, either 
express or implied. 
See the license for the specific language governing permissions and limitations under the license.
