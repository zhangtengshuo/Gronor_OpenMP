# GronOR Non-Orthogonal Configuration Interaction

## Authors 
     T. P. Straatsma, Oak Ridge National Laboratory, Oak Ridge, TN  
	 C. de Graaf, University Rovira i Virgili, Tarragona, Spain  
	 A. Sanchez, University Rovira i Virgili, Tarragona, Spain  
	 R. K. Kathir, University of Groningen, Groningen, Netherlands  

## Reference
    T. P. Straatsma, R. Broer, A. Sanchez-Mansilla, C. Sousa, and C. de Graaf,   
	“GronOR: Scalable and Accelerated Non-Orthogonal Configuration Interaction for Molecular fragment Wave Functions”,   
	Journal of Chemical Theory and Computation, 18, 3549-3565 (2022).  

## Install
Downloading from the GitLab repository:  
```bash
git clone --recursive git@gitlab.com:gronor/gronor.git
#or
git clone --recursive https://www.gitlab.com/gronor/gronor.git
```
The resulting master branch is the most recent release of gronor. To use an earlier release (e.g. version 23.08) use  
```bash
git checkout tags/23.08
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
`-DCUSOLVERJ=ON` will activate the CUSOLVER library including the iterative Jacobi solvers, and requires `-`DOPENACC=ON`  
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
