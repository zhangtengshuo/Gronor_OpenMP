
# Global CMakeLists.txt for building GronOR using cmake

# Set minimum required for cmake
#    Note: As of 3.13 changes in the behavior of option() under policy CMP0077

cmake_minimum_required (VERSION 3.11.00)

# Define project name with associated languages

message ("Programming Environment:")
message (STATUS "CMake Version: ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}")

project (GronOR VERSION 24.06 LANGUAGES Fortran C CXX)

# CMake optional parameters with defaults
#    Note: MPI is not optional and removed from options list

set(TARGET "" CACHE STRING "Target computer system")

option (ACC "Enable OPENACC off-loading." OFF)
option (OPENACC "Enable OPENACC off-loading." OFF)
option (OPENMP "Enable multi-threading." OFF)
option (OMPTGT "Enable OPENMP off-loading." OFF)
option (CUSOLVER "Enable CUSOLVER library." OFF)
option (CUSOLVERJ "Enable CUSOLVER library with Jacobi." OFF)
option (HIPSOLVER "Enable HIPSOLVER library." OFF)
option (HIPSOLVERJ "Enable HIPSOLVER library with Jacobi." OFF)
option (ROCSOLVER "Enable ROCSOLVER library." OFF)
option (MKL "Enable MKL library." OFF)
option (LAPACK "Enable LAPACK library." OFF)
option (MAGMA "Enable MAGMA library." OFF)
option (CUDA "Enable CUDA Functions." OFF)
option (CUBLAS "Enable CUDA BLAS library." OFF)
option (NVBLAS "Enable NVidia BLAS library." OFF)
option (CTEST "Enable CTest." OFF)
option (SINGLEP "Enable single precision." OFF)
option (PROFILING "Enable profiling." OFF)
option (DEBUG "Debugging." OFF)
option (AVX512 "Compiler option for AVX512." OFF)
option (GPU "GPU Vendor" NVIDIA)
option (GUI "Graphical User Interface" OFF)

# Set global session identifier

# Detect system information

message ("Detecting system information:")

if(DEFINED ENV{HOSTNAME})
  message(STATUS "Found hostname: $ENV{HOSTNAME}")
  if($ENV{HOSTNAME} MATCHES ".juwels$")
    set(HOST "JUWELS")
  endif()
  if($ENV{HOSTNAME} MATCHES ".vega.izum.si$")
	  set(HOST "VEGA")
  endif()
  if($ENV{HOSTNAME} MATCHES ".crusher.olcf.ornl.gov$")
	  set(HOST "CRUSHER")
  endif()  
  if($ENV{HOSTNAME} MATCHES ".frontier.olcf.ornl.gov$")
	  set(HOST "FRONTIER")
  endif() 
  if($ENV{HOSTNAME} MATCHES ".leonardo.local$")
    set(HOST "LEONARDO")
  endif() 
else()
  message(STATUS "Not found: hostname")
endif()

if(NOT "${TARGET}" STREQUAL "")
  set(HOST "${TARGET}")
  message(STATUS "Target host is ${HOST}")
endif()

# Operating system

set (OS ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})

message (STATUS "Operating System: ${OS}")

# Address mode

if (${CMAKE_SIZEOF_VOID_P} EQUAL 8)
  set (ADDRMODE 64)
#  add_definitions (-D_I8_)
else ()
  set (ADDRMODE 32)
endif ()

message (STATUS "Address Mode: ${ADDRMODE}")

# Platform settings

if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
  add_definitions (-D_LINUX_)
  if ((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86") OR
      (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "i686"))
    set (PLATFORM "LINUX")
  elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "x86_64")
    set (PLATFORM "LINUX64")
  elseif (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ia64")
    set (PLATFORM "LINUX64_IA")
  elseif ((${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc_64") OR
      (${CMAKE_SYSTEM_PROCESSOR} STREQUAL "ppc64le"))
    set (PLATFORM "PPC64")
  endif ()
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
  add_definitions (-D_LINUX_ -D_DARWIN_)
  set (CMAKE_MACOSX_RPATH 0)
  set (PLATFORM "MacOS")
elseif (${CMAKE_SYSTEM_NAME} STREQUAL "CYGWIN")
  add_definitions (-D_LINUX_ -D_CYGWIN_)
  if (ADDRMODE EQUAL 32)
    set (PLATFORM "WIN")
  elseif (ADDRMODE EQUAL 64)
    set (PLATFORM "WIN64")
  endif ()
endif ()

if (NOT PLATFORM)
  message (FATAL_ERROR "unsupported platform")
else ()
  message (STATUS "Platform: ${PLATFORM}")
endif ()

# General session configuration

set (ENV{CMAKE_SESSION} "GronOR")

message ("Configuring session:")
message (STATUS "Project: ${PROJECT_NAME}")
message (STATUS "Version: ${PROJECT_VERSION}")
message (STATUS "Set environment CMAKE_SESSION: $ENV{CMAKE_SESSION}")

# Enforce out-of-source builds only

set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)

get_filename_component (sourcedir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component (binarydir "${CMAKE_BINARY_DIR}" REALPATH)

if ("${sourcedir}" STREQUAL "${binarydir}")
  message (FATAL_ERROR "In-source builds in ${CMAKE_BINARY_DIR} are not "
    "allowed, please remove ./CMakeCache.txt and ./CMakeFiles/, create a "
    "separate build directory and run cmake from there.")
endif ()

message (STATUS "Source directory: ${CMAKE_SOURCE_DIR}")
message (STATUS "Binary directory: ${CMAKE_BINARY_DIR}")
set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/CMakeFiles/mod")

# Define GronOR compilation

set (LMODCOMP "$ENV{LMOD_FAMILY_COMPILER}/$ENV{LMOD_FAMILY_COMPILER_VERSION}")
set (LMODMPI "$ENV{LMOD_FAMILY_MPI}/$ENV{LMOD_FAMILY_MPI_VERSION}")

# The build type, specifies the different optimization types.

if (NOT DEFINED CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
  set (CMAKE_BUILD_TYPE "Release" CACHE STRING
    "Type of possible builds are: None (CFLAGS/FFLAGS can be used), Debug, Garble, RelWithDebInfo, Release, Fast."
    FORCE)
endif ()

message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS None Debug Garble RelWithDebInfo Release Fast)

# The place where GronOR is installed with make install.

if (NOT DEFINED CMAKE_INSTALL_PREFIX)
  set (CMAKE_INSTALL_PREFIX "/opt/gronor" CACHE STRING
    "Location where GronOR will be installed.")
else()
#  message(STATUS "GRONORVERSION ${PROJECT_NAME}/${PROJECT_VERSION}")
  set (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/gronor/${PROJECT_VERSION}")
endif ()

message (STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
message (STATUS "Install  rpath: ${CMAKE_INSTALL_RPATH}")

message ("Configuring MPI parallellization:")
find_package (MPI REQUIRED)

include_directories (${MPI_Fortran_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})

link_directories(${MPI_Fortran_LIBRARIES})

message (STATUS "MPI_C_INCLUDE_PATH: ${MPI_C_INCLUDE_PATH}")
message (STATUS "MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")
message (STATUS "MPI_Fortran_INCLUDE_PATH: ${MPI_Fortran_INCLUDE_PATH}")
message (STATUS "MPI_Fortran_MODULE_DIR: ${MPI_Fortran_MODULE_DIR}")
message (STATUS "MPI_Fortran_LIBRARIES: ${MPI_Fortran_LIBRARIES}")
message (STATUS "MPIEXEC: ${MPIEXEC}")


message (STATUS "HOST: ${HOST}")
if(HOST STREQUAL "VEGA")
  set (MPI_IMPLEMENTATION "openmpi")
elseif(HOST STREQUAL "JUWELS")
  set (MPI_IMPLEMENTATION "openmpi")
elseif(HOST STREQUAL "LEONARDO")
  set (MPI_IMPLEMENTATION "openmpi")
else ()
execute_process (
  COMMAND ${MPIEXEC} "--version"
  OUTPUT_VARIABLE MPIEXEC_VERSION
  ERROR_VARIABLE MPIEXEC_VERSION
  )
if (MPIEXEC_VERSION MATCHES "OpenRTE")
  set (MPI_IMPLEMENTATION "openmpi")
elseif (MPIEXEC_VERSION MATCHES "Intel\\(R\\) MPI")
  set (MPI_IMPLEMENTATION "impi")
elseif (MPIEXEC_VERSION MATCHES "HYDRA")
  set (MPI_IMPLEMENTATION "mpich")
elseif (MPIEXEC_VERSION MATCHES "mpiexec \\(IBM Spectrum MPI\\)")
  set (MPI_IMPLEMENTATION "openmpi")
elseif (MPIEXEC_VERSION MATCHES "mpirun \\(Open MPI\\)")
  set (MPI_IMPLEMENTATION "openmpi")
elseif (MPIEXEC_VERSION MATCHES "mpiexec \\(Open MPI\\)")
  set (MPI_IMPLEMENTATION "openmpi")
elseif (MPIEXEC_VERSION MATCHES "mpiexec \\(OpenRTE\\)")
  set (MPI_IMPLEMENTATION "openmpi")
elseif(HOST STREQUAL "JUWELS")
  set (MPI_IMPLEMENTATION "openmpi")
elseif(HOST STREQUAL "LEONARDO")
  set (MPI_IMPLEMENTATION "openmpi")
elseif(HOST STREQUAL "VEGA")
  set (MPI_IMPLEMENTATION "openmpi")
elseif(HOST STREQUAL "CRUSHER")
  set (MPI_IMPLEMENTATION "cray-mpich")
elseif(HOST STREQUAL "FRONTIER")
  set (MPI_IMPLEMENTATION "cray-mpich")
else ()
  message (FATAL_ERROR "Unknown or unsupported MPI implementation:\n${MPIEXEC_VERSION}")
endif ()
endif ()
message (STATUS "MPI_IMPLEMENTATION: ${MPI_IMPLEMENTATION}")


# Find external packages, programs and libraries

message ("External packages, programs and libraries:")

# Set CUDA_HOME to ${CUDATOOLKIT_HOME} (Needed on e.g. Piz Daint)

if("$ENV{CUDA_HOME}" STREQUAL "")
  if(NOT "$ENV{CUDATOOLKIT_HOME}" STREQUAL "")
    set(ENV{CUDA_HOME} $ENV{CUDATOOLKIT_HOME})
  endif()
  if(NOT "$ENV{CUDA_DIR}" STREQUAL "")
    set(ENV{CUDA_HOME} $ENV{CUDA_DIR})
  endif()
endif()
if(NOT DEFINED ENV{CUDA_HOME})
  message (STATUS "Not found: cuda")
else()
  message (STATUS "Found cuda: $ENV{CUDA_HOME}")
endif()

find_program (GIT "git")
message (STATUS "Found git: ${GIT}")

set(MOLCAS "ON")

if(NOT DEFINED ENV{MOLCAS})
  set(MOLCAS "")
  message (STATUS "Not found: molcas")
else()
#  string (REPLACE "/build" "" molcas_root $ENV{MOLCAS})
  set(molcas_root $ENV{MOLCAS_BUILD})
  message (STATUS "Found Molcas: ${molcas_root}")
endif()

if(NOT DEFINED ENV{MKLROOT})
  set(MOLCAS "")
  message (STATUS "Not found: mkl")
else()
  message (STATUS "Found mkl: $ENV{MKLROOT}")
endif()

if(NOT DEFINED ENV{OLCF_MAGMA_ROOT})
  message (STATUS "Not found: magma")
else()
  message (STATUS "Found magma: $ENV{OLCF_MAGMA_ROOT}")
  set(MAGMALIB "$ENV{OLCF_MAGMA_ROOT}/lib")
  set(MAGMAINC "$ENV{OLCF_MAGMA_ROOT}/include")
endif()

if(HOST STREQUAL "FRONTIER")
  if(NOT DEFINED ENV{CRAY_LIBSCI_PREFIX})
	  message (STATUS "Not found: libsci LAPACK")
  else()
	  message (STATUS "Found LAPACK: $ENV{CRAY_LIBSCI_PREFIX}")
	  set(LAPACKROOT $ENV{CRAY_LIBSCI_PREFIX})
  endif()
endif()

if(NOT DEFINED ENV{GA})
  set(MOLCAS "")
  message (STATUS "Not found: ga")
else()
  message (STATUS "Found ga: $ENV{GA}")
endif()

if(NOT DEFINED ENV{UTILS})
  message (STATUS "Not found: utils")
  find_package(BLAS)
  find_package(LAPACK)
  #  find_package(SCALAPACK)
else()
  message (STATUS "Found utils: $ENV{UTILS}")
  set(UTILS "$ENV{UTILS}")
endif()

if (MOLCAS STREQUAL "")
  unset(MOLCAS)
endif()
  
if(DEFINED MOLCAS)
  message (STATUS "Molcas dependent auxiliary programs will be compiled")
  include_directories ("${molcas_root}/src/Include")
  message (STATUS "Added include directory: ${molcas_root}/src/Include")
else()
  message (STATUS "Molcas dependent auxiliary programs will not be compiled")
endif()

#     Note: build and debug options

option (DEBUGPRINT "Enable excessive print statements for debugging." OFF)
option (DEBUGTRACE "Enable tracing print statements for debugging." OFF)

option (BUILD_SHARED_LIBS "Build dynamically-linked GronOR library." OFF)
option (BUILD_STATIC_LIB "Build statically-linked GronOR library too." OFF)

if(GUI)
  message ("Configuring Java Programming Environment")
  find_package(Java REQUIRED)
  include(UseJava)
  if(DEFINED ENV{JAVA_CLASSES})
    set(CLASSDIR $ENV{JAVA_CLASSES})
    file (GLOB dir_jars "${CLASSDIR}" "${CLASSDIR}/*")
    list (LENGTH dir_jars num_dir_jars)
    foreach (jardir ${dir_jars})
      file (GLOB jars "${jardir}/*.jar")
      foreach (jarfile ${jars})
        if(${jarfile} MATCHES "cformat.jar$")
          list (APPEND JAVA_EXTERNAL_JARS ${jarfile})
        endif()
        if(${jarfile} MATCHES "plotapplication.jar$")
          list (APPEND JAVA_EXTERNAL_JARS ${jarfile})
        endif()
      endforeach(jarfile)
    endforeach(jardir)
  endif()
  list (LENGTH JAVA_EXTERNAL_JARS numjars)
  message(STATUS "JAVA_EXTERNAL_JARS ${numjars} ${JAVA_EXTERNAL_JARS}")
  set (JAR_NAME gronor_gui)
  set (JAVA_SOURCE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if (DEFINED BUILD_SHARED_LIBS)
  mark_as_advanced (CLEAR BUILD_STATIC_LIB)
else ()
  mark_as_advanced (FORCE BUILD_STATIC_LIB)
endif ()

#      Note: Avoid using CXX compiler for linking

set (CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES FALSE)



message ("Configuring base directories")
  
  # Set cmake module search paths

message (STATUS "CMAKE_CURRENT_LIST DIR: ${CMAKE_CURRENT_LIST_DIR}")

set (GRONOR_DIR ${CMAKE_CURRENT_LIST_DIR})

if (NOT ${PROJECT_SOURCE_DIR} STREQUAL ${GRONOR_DIR})
  set (basedirs ${PROJECT_SOURCE_DIR})
endif ()
list (APPEND basedirs ${GRONOR_DIR})

message(STATUS "Base directories ${basedirs}")

foreach (BASE_DIR ${basedirs})
  file (GLOB source_roots_tmp RELATIVE ${PROJECT_SOURCE_DIR} "${BASE_DIR}/src*")
  list (SORT source_roots_tmp)
  list (APPEND source_roots ${source_roots_tmp})
  file (GLOB test_roots_tmp RELATIVE ${PROJECT_SOURCE_DIR} "${BASE_DIR}/test*")
  list (SORT test_roots_tmp)
  list (APPEND test_roots ${test_roots_tmp})
endforeach ()

if(DEFINED MOLCAS)
  message (STATUS "Adding molcas dependent auxiliary programs")
  foreach (BASE_DIR ${basedirs})
    file (GLOB source_roots_tmp RELATIVE ${PROJECT_SOURCE_DIR} "${BASE_DIR}/molcas*")
    list (SORT source_roots_tmp)
    list (APPEND source_roots ${source_roots_tmp})
  endforeach ()
endif()

message (STATUS "Adding auxiliary programs")
foreach (BASE_DIR ${basedirs})
  file (GLOB source_roots_tmp RELATIVE ${PROJECT_SOURCE_DIR} "${BASE_DIR}/aux*")
  list (SORT source_roots_tmp)
  list (APPEND source_roots ${source_roots_tmp})
endforeach ()

list (REVERSE basedirs)

message (STATUS "Source root directories: ${source_roots}")
message (STATUS "Test root directories: ${test_roots}")

set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set (MAIN_INCLUDE_DIR ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY})

# Copy the gronor_config.h to the include directory after modifying the proper version number
configure_file (${GRONOR_DIR}/include/gronor_config.h.in ${MAIN_INCLUDE_DIR}/gronor_config.fh)

# Create the compiler_mpi.h to the include directory
file(WRITE ${MAIN_INCLUDE_DIR}/gronor_compiler.fh "#define _LMODCOMPILER_ \"${LMODCOMP}\"\n")
file(APPEND ${MAIN_INCLUDE_DIR}/gronor_compiler.fh "#define _LMODMPI_ \"${LMODMPI}\"\n")
file(APPEND ${MAIN_INCLUDE_DIR}/gronor_compiler.fh "#define _CMAKE_ \"${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}\"\n")

STRING(TIMESTAMP COMPTIME "%y/%m/%d %H:%M:%S")

file(WRITE ${MAIN_INCLUDE_DIR}/gronor_compile_time.fh "#define _COMPILE_TIME_ \"${COMPTIME}\"\n")

# Add the include directories from each source base directory to compiler include path
foreach (BASE_DIR ${basedirs})
  if (EXISTS "${BASE_DIR}/src/include")
    include_directories (${BASE_DIR}/src/include)
  endif ()
endforeach ()

# Add the main include directoriy to compiler include path
include_directories (${MAIN_INCLUDE_DIR})

# Define a function to find source sub-directories
function (find_source name)
  foreach (src ${source_roots})
    if (EXISTS ${PROJECT_SOURCE_DIR}/${src}/${name})
      file (GLOB tmplist RELATIVE ${PROJECT_SOURCE_DIR}/${src}/${name} ${PROJECT_SOURCE_DIR}/${src}/${name}/*)
      # ignore possible leftover files from previous configure+make
      list (REMOVE_ITEM tmplist Makefile 00dependencies 00sources)
      if (tmplist)
        set (found "${src}/${name}")
        break ()
      endif ()
    endif ()
  endforeach ()
  # do not overwrite existing definitions
  if (DEFINED ${name}_src)
    set (${name}_src "${${name}_src}" PARENT_SCOPE)
  else ()
    if (DEFINED found)
      set (${name}_src "${found}" PARENT_SCOPE)
    endif ()
  endif ()
endfunction ()

# Define a function to insert list
function (insert_before _list _item _new)
  list (GET ${_item} 0 _first)
  list (FIND ${_list} ${_first} _index)
  if (_index GREATER -1)
    list (INSERT ${_list} ${_index} ${_new})
  else ()
    list (APPEND ${_list} ${_new})
  endif ()
  set (${_list} ${${_list}} PARENT_SCOPE)
endfunction ()

message(STATUS "GPU Vendor is ${GPU}")

# Generate a list of supported c, c++ and fortran compilers with definitions of compiler flags

# GNU Compiler Collection (GCC)

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "GNU")
set (CXXFLAGS_GNU_BASIC "")
set (CXXFLAGS_GNU_OPENMP "-fopenmp")
set (CXXFLAGS_GNU_PROFILING "")
# C++ build targets
set (CXXFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (CXXFLAGS_GNU_GARBLE "-O2 -g -Wall")
set (CXXFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (CXXFLAGS_GNU_RELEASE "-O2")
set (CXXFLAGS_GNU_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "GNU")
set (CFLAGS_GNU_BASIC "-std=gnu99")
set (CFLAGS_GNU_OPENMP "-fopenmp")
set (CFLAGS_GNU_PROFILING "")
# C build targets
set (CFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (CFLAGS_GNU_GARBLE "-O2 -g -Wall")
set (CFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (CFLAGS_GNU_RELEASE "-O2")
set (CFLAGS_GNU_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "GNU")
set (FFLAGS_GNU_BASIC "-DGNU -fallow-argument-mismatch")
set (FFLAGS_GNU_PREPROCESS "-cpp")
set (FFLAGS_GNU_OPENMP "-fopenmp")
set (FFLAGS_GNU_PROFILING "")
set (FFLAGS_GNU_ILP64 "-fdefault-integer-8")
# Fortran build targets
set (FFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (FFLAGS_GNU_GARBLE "-O2 -g -Wall -finit-real=snan -finit-integer=730432726 -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans")
set (FFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (FFLAGS_GNU_RELEASE "-O2")
set (FFLAGS_GNU_FAST "-O3")
# Fix for GNU Fortran compiler version >= 4.8 to prevent wrong optimization
# version check for Fortran doesn't work yet, assume it is the same as the C compiler GCC version
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" AND
    CMAKE_Fortran_COMPILER_ID STREQUAL "GNU" AND
    (NOT CMAKE_C_COMPILER_VERSION VERSION_LESS "4.8"))
  set (FFLAGS_GNU_BASIC "${FFLAGS_GNU_BASIC} -fno-aggressive-loop-optimizations")
  set (FFLAGS_GNU_BOUNDS "-fsanitize=address -fno-omit-frame-pointer")
endif ()
# Add runtime checks
# (No boundary checks because of the Work arrays)
if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "6")
    set (FFLAGS_GNU_GARBLE "${FFLAGS_GNU_GARBLE} -fcheck=array-temps,do,mem,pointer,recursion")
  else ()
    set (FFLAGS_GNU_GARBLE "${FFLAGS_GNU_GARBLE} -fcheck=all,no-bounds")
  endif ()
endif ()

# LLVM C frontend (Clang)

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Clang")
set (CXXFLAGS_Clang_BASIC "-std=gnu99")
set (CXXFLAGS_Clang_OPENMP "-fopenmp")
set (CXXFLAGS_Clang_PROFILING "")
# C++ build targets
set (CXXFLAGS_Clang_DEBUG "-O0 -g -Wall")
set (CXXFLAGS_Clang_GARBLE "-O2 -g -Wall")
set (CXXFLAGS_Clang_RELWITHDEBINFO "-O2 -g -Wall")
set (CXXFLAGS_Clang_RELEASE "-O2")
set (CXXFLAGS_Clang_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Clang")
set (CFLAGS_Clang_BASIC "-std=gnu99")
set (CFLAGS_Clang_OPENMP "-fopenmp")
set (CFLAGS_Clang_PROFILING "")
# C build targets
set (CFLAGS_Clang_DEBUG "-O0 -g -Wall")
set (CFLAGS_Clang_GARBLE "-O2 -g -Wall")
set (CFLAGS_Clang_RELWITHDEBINFO "-O2 -g -Wall")
set (CFLAGS_Clang_RELEASE "-O2")
set (CFLAGS_Clang_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "Flang")
set (FFLAGS_Flang_BASIC "-DFLANG -fbackslash")
set (FFLAGS_Flang_PREPROCESS "")
#set (FFLAGS_Flang_OPENMP "-mp")
set (FFLAGS_Flang_PROFILING "")
set (FFLAGS_Flang_ILP64 "-fdefault-integer-8")
# Fortran build targets
set (FFLAGS_Flang_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_Flang_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_Flang_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_Flang_RELEASE "-O2")
set (FFLAGS_Flang_FAST "-fast -fastsse")

if(HOST STREQUAL "CRUSHER")
set (FFLAGS_Flang_OPENMP "-DAMD -fopenmp -DGPUAMD")
else()
if(HOST STREQUAL "FRONTIER")
set (FFLAGS_Flang_OPENMP "-DAMD -fopenmp -DGPUAMD -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a")
else()
set (FFLAGS_Flang_OPENMP "-mp -h")
endif()
endif()

# Intel C++/C/Fortran Compilers

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Intel")
set (CXXFLAGS_Intel_BASIC "")
set (CXXFLAGS_Intel_OPENMP "-qopenmp")
set (CXXFLAGS_Intel_PROFILING "-g -parallel-source-info=2")
# C++ build targets
set (CXXFLAGS_Intel_DEBUG "-debug -w3")
set (CXXFLAGS_Intel_GARBLE "-O2 -debug -w3")
set (CXXFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -w3")
set (CXXFLAGS_Intel_RELEASE "-O2")
set (CXXFLAGS_Intel_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Intel")
set (CFLAGS_Intel_BASIC "-std=gnu99")
set (CFLAGS_Intel_OPENMP "-qopenmp")
set (CFLAGS_Intel_PROFILING "-g -parallel-source-info=2")
# C build targets
set (CFLAGS_Intel_DEBUG "-debug -w3")
set (CFLAGS_Intel_GARBLE "-O2 -debug -w3")
set (CFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -w3")
set (CFLAGS_Intel_RELEASE "-O2")
set (CFLAGS_Intel_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "Intel")
set (FFLAGS_Intel_BASIC "-DINTEL")
set (FFLAGS_Intel_PREPROCESS "-fpp")
set (FFLAGS_Intel_OPENMP "-qopenmp -DOMP")
set (FFLAGS_Intel_PROFILING "-g -parallel-source-info=2")
set (FFLAGS_Intel_MKL "-DMKL")
set (FFLAGS_Intel_ILP64 "-i8 -r8 -heap-arrays")
#set (FFLAGS_Intel_ILP64 "-r8 -heap-arrays")
# Fortran build targets
set (FFLAGS_Intel_DEBUG "-init=snan -debug -traceback -warn all,nodeclarations -check bounds")
set (FFLAGS_Intel_GARBLE "-O2 -debug -traceback -warn all,nodeclarations -check all,nobounds,noarg_temp_created")
set (FFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -traceback -warn all,nodeclarations")
set (FFLAGS_Intel_RELEASE "-O2 -traceback")
set (FFLAGS_Intel_FAST "-xHOST")
set (FFLAGS_Intel_AVX512 "-xCOMMON-AVX512")
# Intel versions prior to 15 used -openmp
if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
  set (CXXFLAGS_Intel_OPENMP "-openmp")
endif ()
if (CMAKE_C_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_C_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
  set (CFLAGS_Intel_OPENMP "-openmp")
endif ()
if (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
  set (FFLAGS_Intel_OPENMP "-openmp")
endif ()

# IntelLLVM C++/C/Fortran Compilers

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "IntelLLVM")
set (CXXFLAGS_IntelLLVM_BASIC "")
set (CXXFLAGS_IntelLLVM_OPENMP "-qopenmp")
set (CXXFLAGS_IntelLLVM_PROFILING "-g -parallel-source-info=2")
# C++ build targets
set (CXXFLAGS_IntelLLVM_DEBUG "-debug -w3")
set (CXXFLAGS_IntelLLVM_GARBLE "-O2 -debug -w3")
set (CXXFLAGS_IntelLLVM_RELWITHDEBINFO "-O2 -debug -w3")
set (CXXFLAGS_IntelLLVM_RELEASE "-O2")
set (CXXFLAGS_IntelLLVM_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "IntelLLVM")
set (CFLAGS_IntelLLVM_BASIC "-std=gnu99")
set (CFLAGS_IntelLLVM_OPENMP "-qopenmp")
set (CFLAGS_IntelLLVM_PROFILING "-g -parallel-source-info=2")
# C build targets
set (CFLAGS_IntelLLVM_DEBUG "-debug -w3")
set (CFLAGS_IntelLLVM_GARBLE "-O2 -debug -w3")
set (CFLAGS_IntelLLVM_RELWITHDEBINFO "-O2 -debug -w3")
set (CFLAGS_IntelLLVM_RELEASE "-O2")
set (CFLAGS_IntelLLVM_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "IntelLLVM")
set (FFLAGS_IntelLLVM_BASIC "-DINTEL")
set (FFLAGS_IntelLLVM_PREPROCESS "-fpp")
set (FFLAGS_IntelLLVM_OPENMP "-qopenmp -DOMP")
set (FFLAGS_IntelLLVM_PROFILING "-g -parallel-source-info=2")
set (FFLAGS_IntelLLVM_MKL "-DMKL")
set (FFLAGS_IntelLLVM_ILP64 "-i8 -r8 -heap-arrays")
#set (FFLAGS_IntelLLVM_ILP64 "-r8 -heap-arrays")
# Fortran build targets
set (FFLAGS_IntelLLVM_DEBUG "-init=snan -debug -traceback -warn all,nodeclarations -check bounds")
set (FFLAGS_IntelLLVM_GARBLE "-O2 -debug -traceback -warn all,nodeclarations -check all,nobounds,noarg_temp_created")
set (FFLAGS_IntelLLVM_RELWITHDEBINFO "-O2 -debug -traceback -warn all,nodeclarations")
set (FFLAGS_IntelLLVM_RELEASE "-O2 -traceback")
set (FFLAGS_IntelLLVM_FAST "-xHOST")
set (FFLAGS_IntelLLVM_AVX512 "-xCOMMON-AVX512")
# IntelLLVM versions prior to 15 used -openmp
if (CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM" AND
    CMAKE_CXX_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
  set (CXXFLAGS_IntelLLVM_OPENMP "-openmp")
endif ()
if (CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM" AND
    CMAKE_C_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
  set (CFLAGS_IntelLLVM_OPENMP "-openmp")
endif ()
if (CMAKE_Fortran_COMPILER_ID STREQUAL "IntelLLVM" AND
    CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "15.0.0.20140528")
  set (FFLAGS_IntelLLVM_OPENMP "-openmp")
endif ()

# Portland Group (PGI) Compilers

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "PGI")
set (CXXFLAGS_PGI_BASIC "-c99 -pgf90libs")
set (CXXFLAGS_PGI_OPENMP "-mp")
set (CXXFLAGS_PGI_PROFILING "")
# C++ build targets
set (CXXFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_PGI_RELEASE "-O2")
set (CXXFLAGS_PGI_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "PGI")
set (CFLAGS_PGI_BASIC "-c99 -pgf90libs")
set (CFLAGS_PGI_OPENMP "-mp")
set (CFLAGS_PGI_PROFILING "")
# C build targets
set (CFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_PGI_RELEASE "-O2")
set (CFLAGS_PGI_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "PGI")
set (FFLAGS_PGI_BASIC "-DPGI -Mbackslash")
set (FFLAGS_PGI_PREPROCESS "-Mpreprocess")
set (FFLAGS_PGI_OPENMP "-mp")
set (FFLAGS_PGI_PROFILING "")
set (FFLAGS_PGI_ILP64 "-i8")
set (FFLAGS_PGI_CUSOLVER "-DCUSOLVER -lcusolver")
set (FFLAGS_PGI_CUSOLVERJ "-DCUSOLVER -DCUSOLVERJ -lcusolver")
# Fortran build targets
set (FFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_PGI_RELEASE "-O2")
set (FFLAGS_PGI_FAST "-fast -fastsse")

# NVIDIA Compilers

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "NVIDIA")
set (CXXFLAGS_NVIDIA_BASIC "-c99 -pgf90libs")
set (CXXFLAGS_NVIDIA_OPENMP "-mp")
set (CXXFLAGS_NVIDIA_PROFILING "")
# C++ build targets
set (CXXFLAGS_NVIDIA_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_NVIDIA_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_NVIDIA_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_NVIDIA_RELEASE "-O2")
set (CXXFLAGS_NVIDIA_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "NVIDIA")
set (CFLAGS_NVIDIA_BASIC "-c99 -pgf90libs")
set (CFLAGS_NVIDIA_OPENMP "-mp")
set (CFLAGS_NVIDIA_PROFILING "")
# C build targets
set (CFLAGS_NVIDIA_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_NVIDIA_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_NVIDIA_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_NVIDIA_RELEASE "-O2")
set (CFLAGS_NVIDIA_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "NVIDIA")
set (FFLAGS_NVIDIA_BASIC "-DNVIDIA -Mbackslash")
set (FFLAGS_NVIDIA_PREPROCESS "-Mpreprocess")
set (FFLAGS_NVIDIA_OPENMP "-mp")
set (FFLAGS_NVIDIA_PROFILING "-g")
set (FFLAGS_NVIDIA_MKL "-DMKL")
set (FFLAGS_NVIDIA_NVTX "-g -lnvhpcwrapnvtx")
set (FFLAGS_NVIDIA_ILP64 "-i8")
set (FFLAGS_NVIDIA_ACC "-acc -Minfo=accel -gpu=cc80")
set (FFLAGS_NVIDIA_CUSOLVER "-DCUSOLVER -lcusolver")
set (FFLAGS_NVIDIA_CUSOLVERJ "-DCUSOLVER -DCUSOLVERJ -lcusolver")
# Fortran build targets
set (FFLAGS_NVIDIA_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_NVIDIA_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_NVIDIA_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_NVIDIA_RELEASE "-O2")
set (FFLAGS_NVIDIA_FAST "-fast -fastsse")

# NVHPC Compilers

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "NVHPC")
set (CXXFLAGS_NVHPC_BASIC "-c99 -pgf90libs")
if(OMPTGT)
set (CXXFLAGS_NVHPC_OPENMP "-mp=gpu")
else()
set (CXXFLAGS_NVHPC_OPENMP "-mp")
endif()
set (CXXFLAGS_NVHPC_PROFILING "")
# C++ build targets
set (CXXFLAGS_NVHPC_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_NVHPC_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_NVHPC_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_NVHPC_RELEASE "-O2")
set (CXXFLAGS_NVHPC_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "NVHPC")
set (CFLAGS_NVHPC_BASIC "-c99 -pgf90libs")
set (CFLAGS_NVHPC_OPENMP "-mp")
set (CFLAGS_NVHPC_PROFILING "")
# C build targets
set (CFLAGS_NVHPC_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_NVHPC_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_NVHPC_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_NVHPC_RELEASE "-O2")
set (CFLAGS_NVHPC_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "NVHPC")
set (FFLAGS_NVHPC_BASIC "-DNVHPC -Mbackslash -Mlist")
set (FFLAGS_NVHPC_PREPROCESS "-Mpreprocess")
if(OMPTGT)
set (FFLAGS_NVHPC_OPENMP "-mp=gpu")
else()
set (FFLAGS_NVHPC_OPENMP "-mp")
endif()
set (FFLAGS_NVHPC_PROFILING "-g")
set (FFLAGS_NVHPC_MKL "-DMKL")
set (FFLAGS_NVHPC_NVTX "-g -lnvhpcwrapnvtx")
set (FFLAGS_NVHPC_ILP64 "-i8")
if(${GPU} STREQUAL "AMD")
   set (FFLAGS_NVHPC_ACC "-acc -DGPUAMD")
else()
   set (FFLAGS_NVHPC_ACC "-acc -Minfo=accel -gpu=cc80 -DGPUNVIDIA")
endif()
set (FFLAGS_NVHPC_CUSOLVER "-DCUSOLVER -L/autofs/nccs-svm1_sw/summit/nvhpc_sdk/Linux_ppc64le/23.9/math_libs/12.2/targets/ppc64le-linux/lib -lcusolver")
set (FFLAGS_NVHPC_CUSOLVERJ "-DCUSOLVER -DCUSOLVERJ -L/autofs/nccs-svm1_sw/summit/nvhpc_sdk/Linux_ppc64le/23.9/math_libs/12.2/targets/ppc64le-linux/lib -lcusolver")
set (FFLAGS_NVHPC_MAGMA "-DMAGMA -lmagma")
# Fortran build targets
set (FFLAGS_NVHPC_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_NVHPC_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_NVHPC_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_NVHPC_RELEASE "-O2")
set (FFLAGS_NVHPC_FAST "-fast -fastsse")

# Cray Compilers

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Cray")
set (CXXFLAGS_Cray_BASIC "-c99 -pgf90libs")
set (CXXFLAGS_Cray_OPENMP "-h omp")
set (CXXFLAGS_Cray_PROFILING "")
# C++ build targets
set (CXXFLAGS_Cray_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_Cray_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_Cray_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_Cray_RELEASE "-O2")
set (CXXFLAGS_Cray_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Cray")
set (CFLAGS_Cray_BASIC "-c99 -pgf90libs")
set (CFLAGS_Cray_OPENMP "-h omp")
set (CFLAGS_Cray_PROFILING "")
# C build targets
set (CFLAGS_Cray_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_Cray_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_Cray_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_Cray_RELEASE "-O2")
set (CFLAGS_Cray_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "Cray")
#set (FFLAGS_Cray_BASIC "-Mbackslash")
set (FFLAGS_Cray_BASIC "-DCRAY")
set (FFLAGS_Cray_LAPACK "-DLAPACK")
set (FFLAGS_Cray_MAGMA "-DMAGMA -L${MAGMALIB} -lmagma")
#set (FFLAGS_Cray_PREPROCESS "-Mpreprocess")

if(HOST STREQUAL "CRUSHER")
set (FFLAGS_Cray_OPENMP "-DCRAY -DOMP5 -h omp -h list=a -DGPUAMD")
else()
if(HOST STREQUAL "FRONTIER")
set (FFLAGS_Cray_OPENMP "-DCRAY -DOMP5 -h omp -h list=a -DGPUAMD")
else()
set (FFLAGS_Cray_OPENMP "-DCRAY -h omp -h list=a")
endif()
endif()

if(HOST STREQUAL "CRUSHER")
set (FFLAGS_Cray_OPENMP "-DCRAY -DOMP5 -h omp -h list=a -DGPUAMD")
else()
if(HOST STREQUAL "FRONTIER")
set (FFLAGS_Cray_OPENMP "-DCRAY -DOMP5 -h omp -h list=a -DGPUAMD")
else()
set (FFLAGS_Cray_OPENMP "-DCRAY -h omp -h list=a")
endif()
endif()

set (FFLAGS_Cray_PROFILING "")
set (FFLAGS_Cray_ILP64 "-s integer64")
if(${GPU} STREQUAL "AMD")
   set (FFLAGS_Cray_ACC "-h acc -DGPUAMD -h list=a")
else()
   set (FFLAGS_Cray_ACC "-h acc -h list=a -DGPUNVIDIA")
endif()
set (FFLAGS_Cray_CUSOLVER "-DCUSOLVER -lcusolver")
set (FFLAGS_Cray_CUSOLVERJ "-DCUSOLVER -DCUSOLVERJ -lcusolver")
set (FFLAGS_Cray_HIPSOLVER "-M878 -DHIPSOLVER -I$ENV{OLCF_HIPFORT_ROOT}/include/hipfort/amdgcn -L/opt/rocm-5.3.0/hipsolver/lib -lhipsolver")
set (FFLAGS_Cray_HIPSOLVERJ "-M878 -DHIPSOLVER -I$ENV{OLCF_HIPFORT_ROOT}/include/hipfort/amdgcn -L/opt/rocm-5.3.0/hipsolver/lib -DHIPSOLVERJ -lhipsolver")
set (FFLAGS_Cray_ROCSOLVER "-M878 -DROCSOLVER -I$ENV{OLCF_HIPFORT_ROOT}/include/hipfort/amdgcn -L$ENV{OLCF_HIPFORT_ROOT}/lib -lhipfort-amdgcn -L/opt/rocm-5.3.0/lib -lrocsolver -lrocblas -lhipblas -lamdhip64")
set (FFLAGS_Cray_ROCSOLVERJ "-M878 -DROCSOLVER -I$ENV{OLCF_HIPFORT_ROOT}/include/hipfort/amdgcn -L$ENV{OLCF_HIPFORT_ROOT}/lib -lhipfort-amdgcn -L/opt/rocm-5.3.0/rocsolver/lib -DHIPSOLVERJ -lhipsolver")
#set (FFLAGS_Cray_ROCSOLVER "-M878 -DROCSOLVER -I$ENV{OLCF_HIPFORT_ROOT}/include/hipfort/amdgcn -L$ENV{OLCF_HIPFORT_ROOT}/lib -lhipfort-amdgcn -L$ENV{ROCM_PATH}/lib -lhipblas -lamdhip64")
# Fortran build targets
set (FFLAGS_Cray_DEBUG "-O0 -g")
set (FFLAGS_Cray_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_Cray_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_Cray_RELEASE "-O2")
set (FFLAGS_Cray_FAST "-fast -fastsse")

# IBM XL Compilers

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "IBM")
set (CXXFLAGS_IBM_BASIC "-c99 -pgf90libs")
set (CXXFLAGS_IBM_OPENMP "-mp")
set (CXXFLAGS_IBM_PROFILING "")
# C++ build targets
set (CXXFLAGS_IBM_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_IBM_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_IBM_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_IBM_RELEASE "-O2")
set (CXXFLAGS_IBM_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "IBM")
set (CFLAGS_IBM_BASIC "-c99 -pgf90libs -qxlf77")
set (CFLAGS_IBM_OPENMP "-mp")
# C build targets
set (CFLAGS_IBM_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_IBM_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_IBM_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_IBM_RELEASE "-O2")
set (CFLAGS_IBM_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "IBM")
set (FFLAGS_IBM_BASIC "-DIBM -Mbackslash -qlanglvl=77std -qintsize=8")
set (FFLAGS_IBM_PREPROCESS "-Mpreprocess")
set (FFLAGS_IBM_OPENMP "-mp -qoffload -qsmp")
set (FFLAGS_IBM_PROFILING "")
set (FFLAGS_IBM_ILP64 "-i8")
set (FFLAGS_IBM_ACC "-acc -Minfo=accel -Mcuda=rdc")
set (FFLAGS_IBM_CUSOLVER "-DCUSOLVER -lcusolver")
set (FFLAGS_IBM_CUSOLVERJ "-DCUSOLVER -DCUSOLVERJ -lcusolver")
set (FFLAGS_IBM_MAGMA "-DMAGMA -lmagma")
# Fortran build targets
set (FFLAGS_IBM_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_IBM_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_IBM_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_IBM_RELEASE "-O2")
set (FFLAGS_IBM_FAST "-fast -fastsse")

# IBM XL Compilers

# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "XL")
set (CXXFLAGS_XL_BASIC "")
set (CXXFLAGS_XL_OPENMP "-qsmp=omp")
set (CXXFLAGS_XL_PROFILING "")
# C++ build targets
set (CXXFLAGS_XL_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_XL_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_XL_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_XL_RELEASE "-O2")
set (CXXFLAGS_XL_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "XL")
set (CFLAGS_XL_BASIC "")
set (CFLAGS_XL_OPENMP "-qsmp=omp")
# C build targets
set (CFLAGS_XL_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_XL_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_XL_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_XL_RELEASE "-O2")
set (CFLAGS_XL_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "XL")
set (FFLAGS_XL_BASIC "-DIBM -qintsize=8")
set (FFLAGS_XL_PREPROCESS "")
set (FFLAGS_XL_OPENMP "-qsmp=omp -qoffload")
set (FFLAGS_XL_PROFILING "")
set (FFLAGS_XL_ILP64 "")
set (FFLAGS_XL_ACC "-acc -Minfo=accel -Mcuda=rdc")
set (FFLAGS_XL_CUSOLVER "-DCUSOLVER -lcusolver")
set (FFLAGS_XL_CUSOLVERJ "-DCUSOLVER -DCUSOLVERJ -lcusolver")
set (FFLAGS_XL_MAGMA "-DMAGMA -lmagma")
# Fortran build targets
set (FFLAGS_XL_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_XL_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_XL_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_XL_RELEASE "-O2")
set (FFLAGS_XL_FAST "-fast -fastsse")

# Test if a supported c compiler is used
foreach (C_COMPILER ${SUPPORTED_C_COMPILERS})
  if (C_COMPILER STREQUAL "${CMAKE_C_COMPILER_ID}")
    set (CCID ${CMAKE_C_COMPILER_ID})
  endif ()
endforeach ()
if (NOT CCID)
  message (FATAL_ERROR "Unsupported C compiler: ${CMAKE_C_COMPILER_ID}")
else()
  message (STATUS "C compiler used is supported: ${CMAKE_C_COMPILER_ID}") 
endif ()

# Test if a supported c++ compiler is used
foreach (CXX_COMPILER ${SUPPORTED_CXX_COMPILERS})
  if (CXX_COMPILER STREQUAL "${CMAKE_CXX_COMPILER_ID}")
    set (CXXID ${CMAKE_CXX_COMPILER_ID})
  endif ()
endforeach ()
if (NOT CXXID)
  message (FATAL_ERROR "Unsupported C++ compiler: ${CMAKE_CXX_COMPILER_ID}")
else()
  message (STATUS "C++ compiler used is supported: ${CMAKE_CXX_COMPILER_ID}") 
endif ()

# Test if a supported Fortran compiler is used
foreach (Fortran_COMPILER ${SUPPORTED_Fortran_COMPILERS})
  if (Fortran_COMPILER STREQUAL "${CMAKE_Fortran_COMPILER_ID}")
    set (FCID ${CMAKE_Fortran_COMPILER_ID})
    #	  message (STATUS "Fortran compiler ID is ${FCID}")
  endif ()
endforeach ()
if (NOT FCID)
  message (FATAL_ERROR "Unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
else()
  message (STATUS "Fortran compiler used is supported: ${CMAKE_Fortran_COMPILER_ID}")
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Cray")
   set (ENV{CRAY} "Y")
#   message(STATUS "Fortran compiler is ${CMAKE_Fortran_COMPILER_ID}")
   else()
   unset(ENV{CRAY})
#   message(STATUS "Fortran compiler is ${CMAKE_Fortran_COMPILER_ID}")
  endif()
endif ()

# Change selected compiler defaults set by cmake

# overwrite any CMake defaults (some are strange, so we set our own)
if (NOT CMAKE_CXX_FLAGS_USER)
  set (CMAKE_CXX_FLAGS_DEFAULT "$ENV{CXXFLAGS} ${CXXFLAGS_${CXXCID}_BASIC}"
    CACHE STRING "C++ compiler flags." FORCE)
  set (CMAKE_CXX_FLAGS_DEBUG ${CXXFLAGS_${CXXCID}_DEBUG}
    CACHE STRING "C++ compiler flags used for debug builds." FORCE)
  set (CMAKE_CXX_FLAGS_GARBLE ${CXXFLAGS_${CXXCID}_GARBLE}
    CACHE STRING "C++ compiler flags used for garble builds." FORCE)
  set (CMAKE_CXX_FLAGS_RELWITHDEBINFO ${CXXFLAGS_${CXXCID}_RELWITHDEBINFO}
    CACHE STRING "C++ compiler flags used for release builds with debug info." FORCE)
  set (CMAKE_CXX_FLAGS_RELEASE ${CXXFLAGS_${CXXCID}_RELEASE}
    CACHE STRING "C++ compiler flags used for release builds." FORCE)
  set (CMAKE_CXX_FLAGS_FAST ${CXXFLAGS_${CXXCID}_FAST}
    CACHE STRING "C++ compiler flags used for fast release builds." FORCE)
  set (CMAKE_CXX_FLAGS_USER "TRUE"
    CACHE INTERNAL "set to FALSE to reset to default C++ flags" FORCE)
endif ()
if (EXPERT)
  mark_as_advanced(CLEAR
    CMAKE_CXX_FLAGS_DEFAULT
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_GARBLE
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_FAST
    )
else ()
  mark_as_advanced(FORCE
    CMAKE_CXX_FLAGS_DEFAULT
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_GARBLE
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_FAST
    )
endif ()
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEFAULT}")

# overwrite any CMake defaults (some are strange, so we set our own)
if (NOT CMAKE_C_FLAGS_USER)
  set (CMAKE_C_FLAGS_DEFAULT "$ENV{CFLAGS} ${CFLAGS_${CCID}_BASIC}"
    CACHE STRING "C compiler flags." FORCE)
  set (CMAKE_C_FLAGS_DEBUG ${CFLAGS_${CCID}_DEBUG}
    CACHE STRING "C compiler flags used for debug builds." FORCE)
  set (CMAKE_C_FLAGS_GARBLE ${CFLAGS_${CCID}_GARBLE}
    CACHE STRING "C compiler flags used for garble builds." FORCE)
  set (CMAKE_C_FLAGS_RELWITHDEBINFO ${CFLAGS_${CCID}_RELWITHDEBINFO}
    CACHE STRING "C compiler flags used for release builds with debug info." FORCE)
  set (CMAKE_C_FLAGS_RELEASE ${CFLAGS_${CCID}_RELEASE}
    CACHE STRING "C compiler flags used for release builds." FORCE)
  set (CMAKE_C_FLAGS_FAST ${CFLAGS_${CCID}_FAST}
    CACHE STRING "C compiler flags used for fast release builds." FORCE)
  set (CMAKE_C_FLAGS_USER "TRUE"
    CACHE INTERNAL "set to FALSE to reset to default C flags" FORCE)
endif ()
if (EXPERT)
  mark_as_advanced(CLEAR
    CMAKE_C_FLAGS_DEFAULT
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_GARBLE
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_FAST
    )
else ()
  mark_as_advanced(FORCE
    CMAKE_C_FLAGS_DEFAULT
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_GARBLE
    CMAKE_C_FLAGS_RELWITHDEBINFO
    CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_FAST
    )
endif ()
set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS_DEFAULT}")


# overwrite any CMake defaults (some are strange, so we set our own)
if (NOT CMAKE_Fortran_FLAGS_USER)
  set (CMAKE_Fortran_FLAGS_DEFAULT "$ENV{FFLAGS} ${FFLAGS_${FCID}_BASIC}"
    CACHE STRING "Fortran compiler flags." FORCE)
  set (CMAKE_Fortran_FLAGS_DEBUG ${FFLAGS_${FCID}_DEBUG}
    CACHE STRING "Fortran compiler flags used for debug builds." FORCE)
  set (CMAKE_Fortran_FLAGS_GARBLE ${FFLAGS_${FCID}_GARBLE}
    CACHE STRING "Fortran compiler flags used for garble builds." FORCE)
  set (CMAKE_Fortran_FLAGS_RELWITHDEBINFO ${FFLAGS_${FCID}_RELWITHDEBINFO}
    CACHE STRING "Fortran compiler flags used for release builds with debug info." FORCE)
  set (CMAKE_Fortran_FLAGS_RELEASE ${FFLAGS_${FCID}_RELEASE}
    CACHE STRING "Fortran compiler flags used for release builds." FORCE)
  set (CMAKE_Fortran_FLAGS_FAST ${FFLAGS_${FCID}_FAST}
    CACHE STRING "Fortran compiler flags used for fast release builds." FORCE)
  set (CMAKE_Fortran_FLAGS_USER "TRUE"
    CACHE INTERNAL "set to FALSE to reset to default Fortran flags" FORCE)
endif ()
if (EXPERT)
  mark_as_advanced(CLEAR
    CMAKE_Fortran_FLAGS_DEFAULT
    CMAKE_Fortran_FLAGS_DEBUG
    CMAKE_Fortran_FLAGS_GARBLE
    CMAKE_Fortran_FLAGS_RELWITHDEBINFO
    CMAKE_Fortran_FLAGS_RELEASE
    CMAKE_Fortran_FLAGS_FAST
    )
else ()
  mark_as_advanced(FORCE
    CMAKE_Fortran_FLAGS_DEFAULT
    CMAKE_Fortran_FLAGS_DEBUG
    CMAKE_Fortran_FLAGS_GARBLE
    CMAKE_Fortran_FLAGS_RELWITHDEBINFO
    CMAKE_Fortran_FLAGS_RELEASE
    CMAKE_Fortran_FLAGS_FAST
    )
endif ()
set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS_DEFAULT}")


# Fortran accelerator flags
if("$ENV{HOST}" MATCHES "daint")
  set (FFLAGS_NVHPC_ACC "-acc -Minfo=accel -ta=tesla,cc60 -DCUDA")
  set (FFLAGS_NVIDIA_ACC "-acc -Minfo=accel -ta=tesla,cc60 -DCUDA")
  set (FFLAGS_PGI_ACC "-acc -Minfo=accel -ta=tesla,cc60 -DCUDA")
  message (STATUS "Custom OpenACC compiler flags for DAINT: ${FFLAGS_NVHPC_ACC}")
elseif(HOST STREQUAL "JUWELS")
  set (FFLAGS_NVHPC_ACC "-acc -Minfo=accel -Mcuda=rdc -DCUDA")
  set (FFLAGS_NVIDIA_ACC "-acc -Minfo=accel -Mcuda=rdc -DCUDA")
  set (FFLAGS_PGI_ACC "-acc -Minfo=accel -Mcuda=rdc -DCUDA")
  message (STATUS "Custom OpenACC compiler flags for JUWELS: ${FFLAGS_NVHPC_ACC}")
elseif(HOST STREQUAL "LEONARDO")
  set (FFLAGS_NVHPC_ACC "-acc -Minfo=accel -gpu=cc80 -DCUDA")
  set (FFLAGS_NVIDIA_ACC "-acc -Minfo=accel -Mcuda=rdc -DCUDA")
  set (FFLAGS_PGI_ACC "-acc -Minfo=accel -Mcuda=rdc -DCUDA")
  message (STATUS "Custom OpenACC compiler flags for LEONARDO: ${FFLAGS_NVHPC_ACC}")
else()
   if(${GPU} STREQUAL "AMD")
     set (FFLAGS_Cray_ACC "-h acc -h list=a -DAMD")
   else() 
     set (FFLAGS_NVHPC_ACC "-acc -Minfo=accel -Mcuda=rdc -Mlist -DCUDA -DGPUNVIDIA")
#     set (FFLAGS_NVHPC_ACC "-acc -Minfo=accel -cuda -gpu=cc70 -Mlist -DCUDA -DGPUNVIDIA")
     set (FFLAGS_NVIDIA_ACC "-acc -Minfo=accel -Mcuda=rdc -Mlist -DCUDA -DGPUNVIDIA")
     set (FFLAGS_PGI_ACC "-acc -Minfo=accel -Mcuda=rdc -Mlist -DCUDA -DGPUNVIDIA")
     set (FFLAGS_GNU_ACC "-fopenacc -lcuda -fopt-info-optimized-omp -DGPUNVIDIA")
     set (FFLAGS_GCC_ACC "-fopenacc -lcuda -fopt-info-optimized-omp -DGPUNVIDIA")
   endif()
#  message (STATUS "Generic OpenACC compiler flags: ${FFLAGS_NVHPC_ACC}")
endif()

# molcas dependent auxiliary programs need the molcas libraries
if(DEFINED MOLCAS)
  list (APPEND EXTERNAL_LIBRARIES "-L$ENV{MOLCAS}/lib -lmolcas")
  list (APPEND EXTERNAL_LIBRARIES "-L$ENV{GA}/lib -lga -larmci")
  list (APPEND EXTERNAL_LIBRARIES "-L$ENV{MKLLIB} -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lmkl_blacs_intelmpi_ilp64 -lmkl_scalapack_ilp64")
endif()

# Define the appropriate compiler flags

if(DEFINED UTILS)
  message ("TESTING MKL FOR NON_UTILS")
  if(${MKL} MATCHES "OFF")
   list (APPEND EXTERNAL_LIBRARIES "-L$ENV{UTILS}/lib  -llapack -lblas -ltimer -lchem -linp -lsolvers -lcml -ltps")
  else()
   list (APPEND EXTERNAL_LIBRARIES "-L$ENV{UTILS}/lib  -ltimer -lchem -linp -lsolvers -lcml -ltps")
  endif()
  include_directories ("$ENV{UTILS}/include")
endif()

message(STATUS "Host: ${HOST} $ENV{CRAY}")

if(HOST STREQUAL "CRUSHER")
  if($ENV{CRAY} MATCHES "Y")
#   list (APPEND EXTERNAL_LIBRARIES "-ltcmalloc")
#   message(STATUS "Adding tcmalloc")
  endif()
endif()

if(HOST STREQUAL "FRONTIER")
  if($ENV{CRAY} MATCHES "Y")
#   list (APPEND EXTERNAL_LIBRARIES "-ltcmalloc")
#   message(STATUS "Adding tcmalloc")
  endif()
endif()

message(STATUS "External libraries defined: ${EXTERNAL_LIBRARIES}")

# Define the appropriate compiler flags

if (NOT CMAKE_GENERATOR STREQUAL "Ninja")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_PREPROCESS}")
endif ()

if (ADDRMODE EQUAL 64)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_ILP64}")
endif ()

if (OPENMP)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_OPENMP}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_OPENMP}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_OPENMP}")
endif ()

if (OMP)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_OPENMP}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_OPENMP}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_OPENMP}")
endif ()

if (PROFILING)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_PROFILING}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_PROFILING}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_PROFILING}")
endif ()

if (DEBUG)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_DEBUG}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_DEBUG}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_DEBUG}")
endif ()

if (AVX512)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_AVX512}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_AVX512}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_AVX512}")
endif ()

if (FAST)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_FAST}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_FAST}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_FAST}")
endif ()

if (BOUNDS)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_BOUNDS}")
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_BOUNDS}")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_BOUNDS}")
  add_definitions (-D_MALLOC_INTERCEPT_)
endif ()

if(ACC)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_ACC} -DACC")
endif()

if(OPENACC)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_ACC} -DACC")
endif()

if(OMPTGT)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_OPENMP} -DOMP -DOMPTGT")
endif()

if(OPENOMP)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_OPENMP} -DOMP")
endif()

if(OPENOMP)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_OPENMP} -DOMP")
endif()

if(CUDA)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCUDA")
endif()

if(MAGMA)
  set (ENV{MAGMA} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_MAGMA}")
endif()

if(CUSOLVER)
  set (ENV{CUSOLVER} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_CUSOLVER}")
endif()

if(CUSOLVERJ)
  set (ENV{CUSOLVER} "Y")
  set (ENV{CUSOLVERJ} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_CUSOLVERJ}")
endif()

if(HIPSOLVER)
  set (ENV{HIPSOLVER} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_HIPSOLVER}")
endif()

if(HIPSOLVERJ)
  set (ENV{HIPSOLVER} "Y")
  set (ENV{HIPSOLVERJ} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_HIPSOLVERJ}")
endif()

if(ROCSOLVER)
  set (ENV{ROCSOLVER} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_ROCSOLVER}")
endif()

if(SINGLEP)
  set (ENV{SINGLEP} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DSINGLEP")
endif()

if(MKL)
  set (ENV{MKL} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_MKL}")
endif()

if(LAPACK)
  set (ENV{LAPACK} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_LAPACK}")
endif()

if(${PLATFORM} STREQUAL "PPC64")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DLINUX")
endif()

if(${PLATFORM} STREQUAL "LINUX64")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DLINUX")
endif()

if(${PLATFORM} STREQUAL "LINUX")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DLINUX")
endif()

# Set target computer system flag

if(NOT "${TARGET}" STREQUAL "")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -D${TARGET}")
endif()

# Allow user-specified modifications

if (EXTRA_CXXFLAGS)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXXFLAGS}")
endif ()
if (EXTRA_CFLAGS)
  set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${EXTRA_CFLAGS}")
endif ()
if (EXTRA_FFLAGS)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${EXTRA_FFLAGS}")
endif ()

# Install directory

if (NOT IS_ABSOLUTE "${CMAKE_INSTALL_PREFIX}")
  set (CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}")
endif ()
get_filename_component (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)

# Special print statements

if (DEBUGPRINT)
  add_definitions (-D_DEBUG_)
endif ()
if (DEBUGTRACE)
  add_definitions (-D_DEBUG_TRACE_)
endif ()

# Parallel settings

list (APPEND EXTERNAL_LIBRARIES ${MPI_Fortran_LIBRARIES})


# BLAS/LAPACK settings
#=====================

message ("Configuring linear algebra libraries:")

# If one uses BLAS/LAPACK, you can use the CMake FindBLAS/FindLAPACK modules to
# set everything, which is determined by setting the BLA_VENDOR environment
# variable. Unfortunately, currently, the modules that comes with CMake doesn't
# seem to support any ilp64 version of MKL, as it only has e.g. Intel10_64lp
# For now, I'm using some non-portable thing.

set (MKLROOT "" CACHE PATH "MKL root directory.")
if (MKL)
  message (STATUS "Using Intel Math Kernel Library (MKL)")

  # quick bypass for 32bit, don't bother with it (for now?)
  if (NOT ADDRMODE EQUAL 64)
    message (FATAL_ERROR "Use of MKL for 32bit installation not implemented")
  endif ()

  mark_as_advanced (CLEAR MKLROOT)
  if (MKLROOT STREQUAL "")
    set (MKLROOT $ENV{MKLROOT} CACHE PATH "MKL root directory." FORCE)
    if (NOT MKLROOT)
      message (FATAL_ERROR
        "You must set environment variable MKLROOT, "
        "or specify -DMKLROOT=/path/to/mkl_root_dir "
        "when running cmake."
        )
    endif ()
  endif ()
  # at this point, MKLROOT should be defined and not empty
  message (STATUS "MKLROOT = ${MKLROOT}")

  # here we check if MKLROOT has changed after it was last used
  # succesfully, if so then we need to reprocess the checks here.
  if (NOT MKLROOT STREQUAL MKLROOT_LAST)
    # reset MKL paths
    set (MKL_INCLUDE_PATH "${MKLROOT}/include" CACHE PATH
      "location of MKL include files." FORCE)
    if (ADDRMODE EQUAL 64)
      set (libpath "${MKLROOT}/lib/intel64")
    elseif (ADDRMODE EQUAL 32)
      set (libpath "${MKLROOT}/lib/ia32")
    endif ()
    set (MKL_LIBRARY_PATH ${libpath} CACHE PATH
      "location of MKL libraries." FORCE)
      message(STATUS "MKL_LIBRARY_PATH is ${MKL_LIBRARY_PATH}")
    # uncache variables
    unset (LIBMKL_CORE CACHE)
    unset (LIBMKL_INTERFACE CACHE)
    unset (LIBMKL_SEQUENTIAL CACHE)
    unset (LIBMKL_THREADING CACHE)
    unset (LIBMKL_SCALAPACK CACHE)
    unset (LIBMKL_BLACS CACHE)
    # cache last used MKLROOT
    set (MKLROOT_LAST ${MKLROOT} CACHE INTERNAL "last value." FORCE)
  endif ()
  message (STATUS "MKL_INCLUDE_PATH = ${MKL_INCLUDE_PATH}")
  message (STATUS "MKL_LIBRARY_PATH = ${MKL_LIBRARY_PATH}")

  # core library
  find_library (LIBMKL_CORE NAMES "mkl_core"
    PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  # compiler-specific library interface
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    find_library (LIBMKL_INTERFACE NAMES "mkl_gf_ilp64"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    find_library (LIBMKL_INTERFACE NAMES "mkl_intel_ilp64"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI") # apparently PGI uses this too
    find_library (LIBMKL_INTERFACE NAMES "mkl_intel_ilp64"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA") # apparently NVIDIA uses this too
    find_library (LIBMKL_INTERFACE NAMES "mkl_intel_ilp64"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC") # apparently NVIDIA uses this too
    find_library (LIBMKL_INTERFACE NAMES "mkl_intel_ilp64"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  endif ()
  # sequential/compiler-specific threading interface
  find_library (LIBMKL_SEQUENTIAL NAMES "mkl_sequential"
    PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set (LIBMKL_OMP_LINK_FLAGS "${FFLAGS_GNU_OPENMP}")
    find_library (LIBMKL_THREADING NAMES "mkl_gnu_thread"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set (LIBMKL_OMP_LINK_FLAGS "${FFLAGS_Intel_OPENMP}")
    find_library (LIBMKL_THREADING NAMES "mkl_intel_thread"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVIDIA")
    set (LIBMKL_OMP_LINK_FLAGS "${FFLAGS_NVIDIA_OPENMP}")
    find_library (LIBMKL_THREADING NAMES "mkl_intel_thread"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    set (LIBMKL_OMP_LINK_FLAGS "${FFLAGS_NVHPC_OPENMP}")
    find_library (LIBMKL_THREADING NAMES "mkl_intel_thread"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  endif ()
  # find scalapack/blacs for parallel lapack support
  find_library (LIBMKL_SCALAPACK NAMES "mkl_scalapack_ilp64"
    PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  if (MPI_IMPLEMENTATION STREQUAL "openmpi")
    find_library (LIBMKL_BLACS NAMES "mkl_blacs_openmpi_ilp64"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (MPI_IMPLEMENTATION STREQUAL "impi")
    find_library (LIBMKL_BLACS NAMES "mkl_blacs_intelmpi_ilp64"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  elseif (MPI_IMPLEMENTATION STREQUAL "mpich")
    find_library (LIBMKL_BLACS NAMES "mkl_blacs_intelmpi_ilp64"
      PATHS ${MKL_LIBRARY_PATH} NO_DEFAULT_PATH)
  endif ()

  # generate actual library list with paths
  if (MPI)
    add_definitions (-DSCALAPACK)
    list (APPEND MKL_LIBRARIES ${LIBMKL_SCALAPACK})
  endif ()
  list (APPEND MKL_LIBRARIES ${LIBMKL_INTERFACE})
  list (APPEND MKL_LIBRARIES ${LIBMKL_CORE})
  if (OPENMP)
    set (CMAKE_EXE_LINKER_FLAGS ${LIBMKL_OMP_LINK_FLAGS})
    list (APPEND MKL_LIBRARIES ${LIBMKL_THREADING})
  else ()
    list (APPEND MKL_LIBRARIES ${LIBMKL_SEQUENTIAL})
  endif ()
  if (MPI)
    list (APPEND MKL_LIBRARIES ${LIBMKL_BLACS})
  endif ()

  add_definitions (-D_MKL_)
  include_directories (${MKL_INCLUDE_PATH})
  set (LINALG_LIBRARIES ${MKL_LIBRARIES})
else ()
  mark_as_advanced(FORCE MKLROOT)
  message(STATUS "No MKL libraries found")
endif ()

if (LAPACK)
  message (STATUS "Using External LAPACK Library")
endif()

if(${BLAS_FOUND} MATCHES "TRUE")
 if(${FCID} MATCHES "XL")
 else()
   list (APPEND LINALG_LIBRARIES ${BLAS_LIBRARIES})
 endif()
endif()

if(${LAPACK_FOUND} MATCHES "TRUE")
 list (APPEND LINALG_LIBRARIES ${LAPACK_LIBRARIES})
endif()

list (APPEND EXTERNAL_LIBRARIES ${LINALG_LIBRARIES})

message (STATUS "Linear Algebra Libraries: ${LINALG_LIBRARIES}")
message (STATUS "External Libraries: ${EXTERNAL_LIBRARIES}")

# global module directory

set (MAIN_MOD_DIR ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)
if (SINGLE_MOD_DIR)
  include_directories (${MAIN_MOD_DIR}/_single)
else ()
  include_directories (${MAIN_MOD_DIR})
endif ()

# Print summary

message ("Build Summary:")
message (STATUS "Build type: ${CMAKE_BUILD_TYPE}")
string (TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TARGET)
if (NOT BUILD_TARGET STREQUAL "NONE")
  message (STATUS "C compiler: ${CMAKE_C_COMPILER}")
  message (STATUS "C compiler flags: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${BUILD_TARGET}}")
  message (STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
  message (STATUS "C++ compiler flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${BUILD_TARGET}}")
  message (STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
  message (STATUS "Fortran compiler flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${BUILD_TARGET}}")
else ()
  message (STATUS "C compiler: ${CMAKE_C_COMPILER}")
  message (STATUS "C compiler flags: ${CMAKE_C_FLAGS}")
  message (STATUS "C++ compiler: ${CMAKE_CXX_COMPILER}")
  message (STATUS "C++ compiler flags: ${CMAKE_CXX_FLAGS}")
  message (STATUS "Fortran compiler: ${CMAKE_Fortran_COMPILER}")
  message (STATUS "Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
endif ()
get_directory_property(defs COMPILE_DEFINITIONS)
message (STATUS "Definitions: ${defs}")
message (STATUS "Debug definitions: ${DEBUG_DEFS}")

if(MKL)
elseif(UTILS)
elseif(LAPACK)
else()
  if(${BLAS_FOUND} MATCHES "TRUE")
    if(${FCID} MATCHES "XL")
      list(APPEND dir_utils "blas")
    endif()
  else()
    list(APPEND dir_utils "blas")
  endif()
  if(${LAPACK_FOUND} MATCHES "TRUE")
  else()
    list(APPEND dir_utils "lapack")
  endif()
endif()
if(UTILS)
else()
  list(APPEND dir_utils "utils")
  list(APPEND dir_utils "inp")
endif()

# Libraries and Programs

# Make a list of source directories

message ("Configuring libraries list:")

foreach (dir ${source_roots})
  file (GLOB dir_src RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${PROJECT_SOURCE_DIR}/${dir}/[^.]*")
  list(REMOVE_ITEM dir_src "config") 
  list(REMOVE_ITEM dir_src "include")
  if("${AUXILIARY}" STREQUAL "")
    list(REMOVE_ITEM dir_src "aux")
  endif()
  if(MKL)
    list(REMOVE_ITEM dir_src "blas") 
    list(REMOVE_ITEM dir_src "lapack") 
  endif()
  if(LAPACK)
    list(REMOVE_ITEM dir_src "blas") 
    list(REMOVE_ITEM dir_src "lapack") 
  endif()
  if(UTILS)
    list(REMOVE_ITEM dir_src "blas") 
    list(REMOVE_ITEM dir_src "lapack")  
    list(REMOVE_ITEM dir_src "inp")  
    list(REMOVE_ITEM dir_src "utils")
  endif()

# If installed BLAS and LAPACK are found they do not need to be generated from source here

  if(${BLAS_FOUND} MATCHES "TRUE")
    if(${FCID} MATCHES "XL")
    else()
      list(REMOVE_ITEM dir_src "blas")
    endif()
  endif()
  if(${LAPACK_FOUND} MATCHES "TRUE")
   list(REMOVE_ITEM dir_src "lapack") 
  endif()

  foreach (srcdir ${dir_src})
    set (fildir ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
    file (GLOB_RECURSE srcfiles RELATIVE "${PROJECT_SOURCE_DIR}" "${fildir}/*.F" "${fildir}/*.f" "${fildir}/*.F90" "${fildir}/*.f90" "${fildir}/*.c")
    file (GLOB_RECURSE prgfiles RELATIVE "${PROJECT_SOURCE_DIR}" "${fildir}/${srcdir}.F" "${fildir}/${srcdir}.f" "${fildir}/${srcdir}.F90" "${fildir}/${srcdir}.f90" "${fildir}/${srcdir}.c")
    foreach (prg ${prgfiles})
      list (REMOVE_ITEM srcfiles ${prg})
    endforeach(prg)
    set (SRCDIRPATH ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
    list (LENGTH srcfiles numsrc)
    list (LENGTH prgfiles numprg)
    if(${numsrc} GREATER 0)
      message (STATUS "Create library for ${SRCDIRPATH}: ${numsrc} library files")
    endif()
    if(numsrc GREATER 0)
      add_library(${srcdir} STATIC ${srcfiles})
      list(APPEND libraries ${srcdir})
      list(APPEND dir_libs ${srcdir})
    endif()
  endforeach(srcdir)
endforeach(dir)

message (STATUS "Libraries constructed: ${libraries}")

foreach (dir ${dir_utils})
  list(REMOVE_ITEM dir_libs ${dir})
endforeach(dir)

list(LENGTH dir_utils numutils)
if(numutils GREATER 0)
  foreach (dir ${dir_libs})
    add_dependencies(${dir} ${dir_utils})
    target_link_libraries (${dir} ${dir_utils})
    message (STATUS "Dependencies ${dir}: ${dir_utils}")
  endforeach(dir)
endif()

message ("Configuring program list:")

foreach (dir ${source_roots})
  file (GLOB dir_src RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${PROJECT_SOURCE_DIR}/${dir}/[^.]*")
  list(REMOVE_ITEM dir_src "config") 
  list(REMOVE_ITEM dir_src "include") 
  foreach (srcdir ${dir_src})
    set (fildir ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
    file (GLOB_RECURSE srcfiles RELATIVE "${PROJECT_SOURCE_DIR}" "${fildir}/*.F" "${fildir}/*.f" "${fildir}/*.f90" "${fildir}/*.c")
    file (GLOB_RECURSE prgfiles RELATIVE "${PROJECT_SOURCE_DIR}" "${fildir}/${srcdir}.F" "${fildir}/${srcdir}.f" "${fildir}/${srcdir}.f90" "${fildir}/${srcdir}.c")
    foreach (prg ${prgfiles})
      list (REMOVE_ITEM srcfiles ${prg})
    endforeach(prg)
    set (SRCDIRPATH ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
    list (LENGTH srcfiles numsrc)
    list (LENGTH prgfiles numprg)
    if(numprg GREATER 0)
      add_executable(${srcdir} ${prgfiles})
      message (STATUS "Create program executable ${srcdir} for ${prgfiles}")
      add_dependencies (${srcdir} ${libraries})
      target_link_libraries (${srcdir} ${libraries} ${EXTERNAL_LIBRARIES})
      list(APPEND programs ${srcdir})
      list(APPEND dir_prgs ${srcdir})
    endif()
  endforeach(srcdir)
endforeach(dir)

message (STATUS "Programs constructed: ${programs}")

message ("Configuring java list:")

foreach (dir ${source_roots})
  file (GLOB dir_src RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${PROJECT_SOURCE_DIR}/${dir}/[^.]*")
  list(REMOVE_ITEM dir_src "config") 
  list(REMOVE_ITEM dir_src "include") 
  foreach (srcdir ${dir_src})
    set (fildir ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
    file (GLOB_RECURSE javafiles RELATIVE "${PROJECT_SOURCE_DIR}" "${fildir}/*.java")
    set (SRCDIRPATH ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
    list (LENGTH javafiles numjava)
    if(numjava GREATER 0)
      message (STATUS "Create java classes for ${SRCDIRPATH}: ${numjava} java files ${javafiles}")
      list (APPEND JAVA_SOURCE_FILES ${javafiles})
    endif()
  endforeach(srcdir)
endforeach(dir)

message ("Configuring install list:")

install(TARGETS ${libraries})
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
install(TARGETS ${programs})
set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/bin")

enable_testing()

include(CTest)

message ("Configuring testing:")
foreach (dir ${test_roots})
  file (GLOB dir_test RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${PROJECT_SOURCE_DIR}/${dir}/[^.]*")
  foreach (testdir ${dir_test})
    set (fildir ${PROJECT_SOURCE_DIR}/${dir}/${testdir})
    file (GLOB_RECURSE testfiles RELATIVE "${fildir}" "${fildir}/*.run")
    foreach (test ${testfiles})
      string(REPLACE ".run" "" testroot ${test})
      message (STATUS "Test included:" ${fildir} "/" ${testroot})
      add_test (NAME ${testroot} WORKING_DIRECTORY ${fildir} COMMAND python ${PROJECT_SOURCE_DIR}/scripts/test.py ${testroot})
    endforeach(test)
  endforeach(testdir)
endforeach(dir)

message ("Finalization:")
set (LINK_DIRECTORIES ${LIBRARY_OUTPUT_PATH})
message (STATUS "Link directories: ${LINK_DIRECTORIES}")


if(GUI)
  if(${numjars} GREATER 1)
    add_jar(${JAR_NAME} SOURCES ${JAVA_SOURCE_FILES} INCLUDE_JARS ${JAVA_EXTERNAL_JARS} ENTRY_POINT gronor.gronor_Main)
    install_jar(${JAR_NAME} DESTINATION "${CMAKE_INSTALL_PREFIX}/jar")
  endif()
endif()
