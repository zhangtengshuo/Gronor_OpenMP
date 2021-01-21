
# Global CMakeLists.txt for building GronOR using cmake

cmake_minimum_required (VERSION 2.8.11)

# Enforce out-of-source builds only

set (CMAKE_DISABLE_SOURCE_CHANGES ON)
set (CMAKE_DISABLE_IN_SOURCE_BUILD ON)

get_filename_component (sourcedir "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component (binarydir "${CMAKE_BINARY_DIR}" REALPATH)

if ("${sourcedir}" STREQUAL "${binarydir}")
        message(FATAL_ERROR "In-source builds in ${CMAKE_BINARY_DIR} are not "
                "allowed, please remove ./CMakeCache.txt and ./CMakeFiles/, create a "
                "separate build directory and run cmake from there.")
endif ()

# Set global session identifier

set (ENV{CMAKE_SESSION} "GronOR")

message ("Session: $ENV{CMAKE_SESSION}")

message ("-- Source: ${CMAKE_SOURCE_DIR}")
#message ("Binary dir is ${CMAKE_BINARY_DIR}")

# Find required tools

find_program (GIT "git")

message ("-- Found git: ${GIT}")

# Define GronOR compilation options

# The build type, specifies the different optimization types.

if (NOT DEFINED CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "")
   set (CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Type of build, options are: None (CFLAGS/FFLAGS can be used), Debug, Garble, RelWithDebInfo, Release, Fast."
      FORCE)
endif ()

message ("-- Build type: ${CMAKE_BUILD_TYPE}")

set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS None Debug Garble RelWithDebInfo Release Fast)

# The place where GronOR is installed with make install.

if (NOT DEFINED CMAKE_INSTALL_PREFIX)
        set (CMAKE_INSTALL_PREFIX "/opt/gronor" CACHE STRING
                "Location where GronOR will be installed.")
endif ()

message("-- Install prefix: ${CMAKE_INSTALL_PREFIX}")

option (DEBUGPRINT "Enable excessive print statements for debugging." OFF)
option (DEBUGTRACE "Enable tracing print statements for debugging." OFF)

set (DEBUG_DEFS_DESC "Additional debugging definitions.")
if (DEFINED DEBUG_DEFS)
        separate_arguments(DEBUG_DEFS)
        set (DEBUG_DEFS "${DEBUG_DEFS}" CACHE STRING "${DEBUG_DEFS_DESC}" FORCE)
else ()
        set (DEBUG_DEFS "" CACHE STRING "${DEBUG_DEFS_DESC}")
endif ()

option (BUILD_SHARED_LIBS "Build dynamically-linked GronOR library." OFF)
option (BUILD_STATIC_LIB "Build statically-linked GronOR library too." OFF)
if (BUILD_SHARED_LIBS)
        mark_as_advanced (CLEAR BUILD_STATIC_LIB)
else ()
        mark_as_advanced (FORCE BUILD_STATIC_LIB)
      endif ()
      
option (MPI "Enable MPI parallellization." ON)
option (OPENACC "Enable OPENACC off-loading." OFF)
option (_OL_ACC_ "Enable OPENACC off-loading." OFF)
option (OPENMP "Enable multi-threading." OFF)
option (_OL_OMP_ "Enable OPENMP off-loading." OFF)
option (CUSOLVER "Enable CUSOLVER library." OFF)
option (CUSOLVERJ "Enable CUSOLVER library with Jacobi." OFF)
option (MKL "Enable MKL library." OFF)
option (CUDA "Enable CUDA Functions." OFF)
option (CUBLAS "Enable CUDA BLAS library." OFF)
option (NVBLAS "Enable NVidia BLAS library." OFF)
option (CTEST "Enable CTest." OFF)
option (SINGLEP "Enable single precision." OFF)


# Define project name with associated languages

project (GronOR Fortran C CXX)

# Piz Daint specific stuff

if("$ENV{CUDA_HOME}" STREQUAL "")
  if("$ENV{CUDATOOLKIT_HOME}" STREQUAL "")
    message ("-- CUDA_HOME not set")
  else()
    set(ENV{CUDA_HOME} $ENV{CUDATOOLKIT_HOME})
    message ("-- CUDA_HOME set to CUDATOOLKIT_HOME: $ENV{CUDA_HOME}")
  endif()
else()
  message ("-- CUDA_HOME found: $ENV{CUDA_HOME}")
endif()

# Avoid using CXX compiler for linking

set (CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES FALSE)

# Set cmake module search paths

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH})

#message ("Project source dir is ${PROJECT_SOURCE_DIR}")
#message ("Project binary dir is ${PROJECT_BINARY_DIR}")
#message ("CMake current list dir is ${CMAKE_CURRENT_LIST_DIR}")

set (GRONOR_DIR ${CMAKE_CURRENT_LIST_DIR})

if (NOT ${PROJECT_SOURCE_DIR} STREQUAL ${GRONOR_DIR})
        set (basedirs ${PROJECT_SOURCE_DIR})
endif ()

list (APPEND basedirs ${GRONOR_DIR})


foreach (BASE_DIR ${basedirs})
  file (GLOB source_roots_tmp RELATIVE ${PROJECT_SOURCE_DIR} "${BASE_DIR}/src*")
  list (SORT source_roots_tmp)
  list (APPEND source_roots ${source_roots_tmp})
  file (GLOB test_roots_tmp RELATIVE ${PROJECT_SOURCE_DIR} "${BASE_DIR}/test*")
  list (SORT test_roots_tmp)
  list (APPEND test_roots ${test_roots_tmp})
endforeach ()
      
list (REVERSE basedirs)

#message ("basedirs is ${basedirs}")
#message ("source_roots is ${source_roots}")
#message ("test_roots is ${test_roots}")

set (LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set (EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

#message ("library output path is ${LIBRARY_OUTPUT_PATH}")
#message ("executable output path is ${EXECUTABLE_OUTPUT_PATH}")

set (MAIN_INCLUDE_DIR ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/include)

configure_file (${GRONOR_DIR}/include/gronor_version.h.in ${MAIN_INCLUDE_DIR}/gronor_version.h)

foreach (BASE_DIR ${basedirs})
        if (EXISTS "${BASE_DIR}/src/include")
                include_directories (${BASE_DIR}/src/include)
        endif ()
endforeach ()
include_directories (${MAIN_INCLUDE_DIR})

#message ("main include dir is ${MAIN_INCLUDE_DIR}")

# set function to find source directories

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
                #else ()
                #        message (WARNING "${name} not found, expect trouble")
                endif ()
        endif ()
endfunction ()

# set function to insert list

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

# Detect system information

message ("Detecting system info:")

# Operating system

set (OS ${CMAKE_SYSTEM_NAME}-${CMAKE_SYSTEM_PROCESSOR})

message ("-- OS: ${OS}")

# Address mode

if (${CMAKE_SIZEOF_VOID_P} EQUAL 8)
        set (ADDRMODE 64)
        add_definitions (-D_I8_)
else ()
        set (ADDRMODE 32)
endif ()

message ("-- ADDRMODE: ${ADDRMODE}")

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
        message ("-- PLATFORM: ${PLATFORM}")
endif ()

# GNU Compiler Collection (GCC)

# defines __GNUC__
# C++ compiler

list (APPEND SUPPORTED_CXX_COMPILERS "GNU")
set (CXXFLAGS_GNU_BASIC "")
set (CXXFLAGS_GNU_OPENMP "-fopenmp")
# build targets
set (CXXFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (CXXFLAGS_GNU_GARBLE "-O2 -g -Wall")
set (CXXFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (CXXFLAGS_GNU_RELEASE "-O2")
set (CXXFLAGS_GNU_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "GNU")
set (CFLAGS_GNU_BASIC "-std=gnu99")
set (CFLAGS_GNU_OPENMP "-fopenmp")
# build targets
set (CFLAGS_GNU_DEBUG "-O0 -g -Wall")
set (CFLAGS_GNU_GARBLE "-O2 -g -Wall")
set (CFLAGS_GNU_RELWITHDEBINFO "-O2 -g -Wall")
set (CFLAGS_GNU_RELEASE "-O2")
set (CFLAGS_GNU_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "GNU")
set (FFLAGS_GNU_BASIC "")
set (FFLAGS_GNU_PREPROCESS "-cpp")
set (FFLAGS_GNU_OPENMP "-fopenmp")
set (FFLAGS_GNU_ILP64 "-fdefault-integer-8")
# build targets
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

# defines __clang__
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Clang")
set (CXXFLAGS_Clang_BASIC "-std=gnu99")
set (CXXFLAGS_Clang_OPENMP "-fopenmp")
# build targets
set (CXXFLAGS_Clang_DEBUG "-O0 -g -Wall")
set (CXXFLAGS_Clang_GARBLE "-O2 -g -Wall")
set (CXXFLAGS_Clang_RELWITHDEBINFO "-O2 -g -Wall")
set (CXXFLAGS_Clang_RELEASE "-O2")
set (CXXFLAGS_Clang_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Clang")
set (CFLAGS_Clang_BASIC "-std=gnu99")
set (CFLAGS_Clang_OPENMP "-fopenmp")
# build targets
set (CFLAGS_Clang_DEBUG "-O0 -g -Wall")
set (CFLAGS_Clang_GARBLE "-O2 -g -Wall")
set (CFLAGS_Clang_RELWITHDEBINFO "-O2 -g -Wall")
set (CFLAGS_Clang_RELEASE "-O2")
set (CFLAGS_Clang_FAST "-O2")

# Intel C++/C/Fortran Compilers

# defines __INTEL_COMPILER
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Intel")
set (CXXFLAGS_Intel_BASIC "")
set (CXXFLAGS_Intel_OPENMP "-qopenmp")
# build targets
set (CXXFLAGS_Intel_DEBUG "-debug -w3")
set (CXXFLAGS_Intel_GARBLE "-O2 -debug -w3")
set (CXXFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -w3")
set (CXXFLAGS_Intel_RELEASE "-O2")
set (CXXFLAGS_Intel_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Intel")
set (CFLAGS_Intel_BASIC "-std=gnu99")
set (CFLAGS_Intel_OPENMP "-qopenmp")
# build targets
set (CFLAGS_Intel_DEBUG "-debug -w3")
set (CFLAGS_Intel_GARBLE "-O2 -debug -w3")
set (CFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -w3")
set (CFLAGS_Intel_RELEASE "-O2")
set (CFLAGS_Intel_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "Intel")
set (FFLAGS_Intel_BASIC "-xHost")
set (FFLAGS_Intel_PREPROCESS "-fpp")
set (FFLAGS_Intel_OPENMP "-qopenmp")
set (FFLAGS_Intel_MKL "-DMKL")
set (FFLAGS_Intel_ILP64 "-i8 -r8 -heap-arrays")
# build targets
set (FFLAGS_Intel_DEBUG "-init=snan -debug -traceback -warn all,nodeclarations")
set (FFLAGS_Intel_GARBLE "-O2 -debug -traceback -warn all,nodeclarations -check all,nobounds,noarg_temp_created")
set (FFLAGS_Intel_RELWITHDEBINFO "-O2 -debug -traceback -warn all,nodeclarations")
set (FFLAGS_Intel_RELEASE "-O2 -traceback")
set (FFLAGS_Intel_FAST "-fast")

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

# Portland Group (PGI) Compilers

# defines __PGI
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "PGI")
set (CXXFLAGS_PGI_BASIC "-c99 -pgf90libs")
set (CXXFLAGS_PGI_OPENMP "-mp")
# build targets
set (CXXFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_PGI_RELEASE "-O2")
set (CXXFLAGS_PGI_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "PGI")
set (CFLAGS_PGI_BASIC "-c99 -pgf90libs")
set (CFLAGS_PGI_OPENMP "-mp")
# build targets
set (CFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_PGI_RELEASE "-O2")
set (CFLAGS_PGI_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "PGI")
set (FFLAGS_PGI_BASIC "-Mbackslash")
set (FFLAGS_PGI_PREPROCESS "-Mpreprocess")
set (FFLAGS_PGI_OPENMP "-mp")
set (FFLAGS_PGI_ILP64 "-i8")
if("$ENV{HOST}" MATCHES "daint")
set (FFLAGS_PGI_ACC "-acc -Minfo=accel -ta=tesla,cc60")
message ("-- Custom OpenACC compiler flags: ${FFLAGS_PGI_ACC}")
else()
set (FFLAGS_PGI_ACC "-acc -Minfo=accel -Mcuda=rdc")
message ("-- Generic OpenACC compiler flags: ${FFLAGS_PGI_ACC}")
endif()
set (FFLAGS_PGI_CUSOLVER "-DCUSOLVER -lcusolver")
set (FFLAGS_PGI_CUSOLVERJ "-DCUSOLVER -DCUSOLVERJ -lcusolver")
# build targets
set (FFLAGS_PGI_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_PGI_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_PGI_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_PGI_RELEASE "-O2")
set (FFLAGS_PGI_FAST "-fast -fastsse")

# NVIDIA Compilers

# defines __NVIDIA
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "NVIDIA")
set (CXXFLAGS_NVIDIA_BASIC "-c99 -pgf90libs")
set (CXXFLAGS_NVIDIA_OPENMP "-mp")
# build targets
set (CXXFLAGS_NVIDIA_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_NVIDIA_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_NVIDIA_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_NVIDIA_RELEASE "-O2")
set (CXXFLAGS_NVIDIA_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "NVIDIA")
set (CFLAGS_NVIDIA_BASIC "-c99 -pgf90libs")
set (CFLAGS_NVIDIA_OPENMP "-mp")
# build targets
set (CFLAGS_NVIDIA_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_NVIDIA_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_NVIDIA_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_NVIDIA_RELEASE "-O2")
set (CFLAGS_NVIDIA_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "NVIDIA")
set (FFLAGS_NVIDIA_BASIC "-Mbackslash")
set (FFLAGS_NVIDIA_PREPROCESS "-Mpreprocess")
set (FFLAGS_NVIDIA_OPENMP "-mp")
set (FFLAGS_NVIDIA_ILP64 "-i8")
set (FFLAGS_NVIDIA_ACC "-acc -Minfo=accel -Mcuda=rdc")
set (FFLAGS_NVIDIA_CUSOLVER "-DCUSOLVER -lcusolver")
set (FFLAGS_NVIDIA_CUSOLVERJ "-DCUSOLVER -DCUSOLVERJ -lcusolver")
# build targets
set (FFLAGS_NVIDIA_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_NVIDIA_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_NVIDIA_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_NVIDIA_RELEASE "-O2")
set (FFLAGS_NVIDIA_FAST "-fast -fastsse")

# Cray Compilers

# defines __Cray
# C++ compiler
list (APPEND SUPPORTED_CXX_COMPILERS "Cray")
set (CXXFLAGS_Cray_BASIC "-c99 -pgf90libs")
set (CXXFLAGS_Cray_OPENMP "-mp")
# build targets
set (CXXFLAGS_Cray_DEBUG "-O0 -g -Minform=warn")
set (CXXFLAGS_Cray_GARBLE "-O2 -gopt -Minform=warn")
set (CXXFLAGS_Cray_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CXXFLAGS_Cray_RELEASE "-O2")
set (CXXFLAGS_Cray_FAST "-O2")
# C compiler
list (APPEND SUPPORTED_C_COMPILERS "Cray")
set (CFLAGS_Cray_BASIC "-c99 -pgf90libs")
set (CFLAGS_Cray_OPENMP "-mp")
# build targets
set (CFLAGS_Cray_DEBUG "-O0 -g -Minform=warn")
set (CFLAGS_Cray_GARBLE "-O2 -gopt -Minform=warn")
set (CFLAGS_Cray_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (CFLAGS_Cray_RELEASE "-O2")
set (CFLAGS_Cray_FAST "-O2")
# Fortran compiler
list (APPEND SUPPORTED_Fortran_COMPILERS "Cray")
set (FFLAGS_Cray_BASIC "-Mbackslash")
set (FFLAGS_Cray_PREPROCESS "-Mpreprocess")
set (FFLAGS_Cray_OPENMP "-mp")
set (FFLAGS_Cray_ILP64 "-i8")
set (FFLAGS_Cray_ACC "-acc -Minfo=accel -Mcuda=rdc")
set (FFLAGS_Cray_CUSOLVER "-DCUSOLVER -lcusolver")
# build targets
set (FFLAGS_Cray_DEBUG "-O0 -g -Minform=warn")
set (FFLAGS_Cray_GARBLE "-O2 -gopt -Minform=warn")
set (FFLAGS_Cray_RELWITHDEBINFO "-O2 -gopt -Minform=warn")
set (FFLAGS_Cray_RELEASE "-O2")
set (FFLAGS_Cray_FAST "-fast -fastsse")


# C compiler
foreach (C_COMPILER ${SUPPORTED_C_COMPILERS})
        if (C_COMPILER STREQUAL "${CMAKE_C_COMPILER_ID}")
                set (CCID ${CMAKE_C_COMPILER_ID})
        endif ()
endforeach ()

if (NOT CCID)
        message (FATAL_ERROR "unsupported C compiler: ${CMAKE_C_COMPILER_ID}")
endif ()

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

# Fortran compiler
foreach (Fortran_COMPILER ${SUPPORTED_Fortran_COMPILERS})
        if (Fortran_COMPILER STREQUAL "${CMAKE_Fortran_COMPILER_ID}")
          set (FCID ${CMAKE_Fortran_COMPILER_ID})
#	  message ("-- Fortran compiler ID is ${FCID}")
        endif ()
endforeach ()

if (NOT FCID)
        message (FATAL_ERROR "unsupported Fortran compiler: ${CMAKE_Fortran_COMPILER_ID}")
endif ()

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

if (BOUNDS)
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CXXFLAGS_${CXXCID}_BOUNDS}")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${CFLAGS_${CCID}_BOUNDS}")
        set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_BOUNDS}")
        add_definitions (-D_MALLOC_INTERCEPT_)
endif ()

if(OPENACC)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_ACC} -D_OL_ACC_")
endif()

if(_OL_ACC_)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_ACC} -D_OL_ACC_")
endif()

if(_OL_OMP_)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_OPENMP} -D_OL_OMP_")
endif()

if(OPENOMP)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_OPENMP}")
endif()

if(CUDA)
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DCUDA")
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
      
if(SINGLEP)
  set (ENV{SINGLEP} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DSINGLEP")
endif()
      
if(MKL)
  set (ENV{MKL} "Y")
  set (CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${FFLAGS_${FCID}_MKL}")
endif()
      
if(${PLATFORM} STREQUAL "PPC64")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DLINUX")
   message("-- Added LINUX compile flag")
endif()
if(${PLATFORM} STREQUAL "LINUX64")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DLINUX")
   message("-- Added LINUX compile flag")
endif()
if(${PLATFORM} STREQUAL "LINUX")
   set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -DLINUX")
   message("-- Added LINUX compile flag")
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

# Process options                                                              #

# Install directory

if (NOT IS_ABSOLUTE "${CMAKE_INSTALL_PREFIX}")
        set (CMAKE_INSTALL_PREFIX "${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_PREFIX}")
endif ()
get_filename_component (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)

#message ("CMake install prefix is ${CMAKE_INSTALL_PREFIX}")

# Special print statements

if (DEBUGPRINT)
        add_definitions (-D_DEBUG_)
endif ()
if (DEBUGTRACE)
        add_definitions (-D_DEBUG_TRACE_)
endif ()

# Parallel settings

if (MPI)
        message ("Configuring with MPI parallellization:")

        find_package (MPI REQUIRED)

        add_definitions (-D_GRONOR_MPP_)

        include_directories (${MPI_Fortran_INCLUDE_PATH} ${MPI_C_INCLUDE_PATH})

        link_directories(${MPI_Fortran_LIBRARIES})

        message ("-- MPI_C_INCLUDE_PATH: ${MPI_C_INCLUDE_PATH}")
        message ("-- MPI_Fortran_INCLUDE_PATH: ${MPI_Fortran_INCLUDE_PATH}")
        message ("-- MPI_C_LIBRARIES: ${MPI_C_LIBRARIES}")
        message ("-- MPI_Fortran_LIBRARIES: ${MPI_Fortran_LIBRARIES}")
        message ("-- MPIEXEC: ${MPIEXEC}")

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
        elseif (MPIEXEC_VERSION MATCHES "slurm 20.02.6-1.20201123git821b81e9af.el8")
                set (MPI_IMPLEMENTATION "openmpi")
        else ()
                message (FATAL_ERROR "Unknown or unsupported MPI implementation:\n${MPIEXEC_VERSION}")
        endif ()
        message ("-- MPI_IMPLEMENTATION: ${MPI_IMPLEMENTATION}")

        list (APPEND EXTERNAL_LIBRARIES ${MPI_Fortran_LIBRARIES})
endif ()

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
        message ("-- Using Intel Math Kernel Library (MKL)")

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
        message ("-- MKLROOT = ${MKLROOT}")

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
        message ("-- MKL_INCLUDE_PATH = ${MKL_INCLUDE_PATH}")
        message ("-- MKL_LIBRARY_PATH = ${MKL_LIBRARY_PATH}")

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
endif ()

list (APPEND EXTERNAL_LIBRARIES ${LINALG_LIBRARIES})

message ("-- LINALG_LIBRARIES: ${LINALG_LIBRARIES}")

# global module directory

set (MAIN_MOD_DIR ${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mod)
if (SINGLE_MOD_DIR)
        include_directories (${MAIN_MOD_DIR}/_single)
else ()
        include_directories (${MAIN_MOD_DIR})
      endif ()
      
# launcher settings

# The runtime can be controlled by the user through two variables:
# SER_LAUNCHER: a command used to launch a serial program
# MPI_LAUNCHER: a command used to launch an MPI program
# Since these are advanced options, the user would typically modify
# them in the CMakeCache.txt. However, they can also be activated
# by other options for internal use, which is why their presence
# is detected instead of just caching a default string.

# Scripts are tricky, we don't want to set automatic wrappers

set (RUNSCRIPT "$program $input" CACHE STRING "Command to launch scripts.")

# Specify the launch command for a serial run

if (NOT SER_LAUNCHER)
        set (SER_LAUNCHER ""                CACHE STRING "Command to launch serial programs.")
else ()
        set (SER_LAUNCHER "${SER_LAUNCHER}" CACHE STRING "Command to launch serial programs.")
endif ()
set (RUNBINARYSER "${SER_LAUNCHER} $program")

# Specify the launch command for an MPI run

if (MPI)
        if (NOT MPI_LAUNCHER)
                set (MPI_LAUNCHER "${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} $GRONOR_NPROCS ${MPIEXEC_PREFLAGS}" CACHE STRING "Command to launch MPI programs.")
        endif ()
        set (RUNBINARY "${MPI_LAUNCHER} $program")
else ()
        set (RUNBINARY "${SER_LAUNCHER} $program")
endif ()

string (STRIP ${RUNSCRIPT} RUNSCRIPT)
string (STRIP ${RUNBINARY} RUNBINARY)
string (STRIP ${RUNBINARYSER} RUNBINARYSER)

message ("-- RUNSCRIPT:    ${RUNSCRIPT}")
message ("-- RUNBINARY:    ${RUNBINARY}")
message ("-- RUNBINARYSER: ${RUNBINARYSER}")

# Print summary

#message ("Build type: ${CMAKE_BUILD_TYPE}")
#string (TOUPPER ${CMAKE_BUILD_TYPE} BUILD_TARGET)
if (NOT BUILD_TARGET STREQUAL "NONE")
        message ("-- C compiler: ${CMAKE_C_COMPILER}")
        message ("-- C compiler flags: ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${BUILD_TARGET}}")
        message ("-- Fortran compiler: ${CMAKE_Fortran_COMPILER}")
        message ("-- Fortran compiler flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${BUILD_TARGET}}")
else ()
        message ("-- C compiler: ${CMAKE_C_COMPILER}")
        message ("-- C compiler flags: ${CMAKE_C_FLAGS}")
        message ("-- Fortran compiler: ${CMAKE_Fortran_COMPILER}")
        message ("-- Fortran compiler flags: ${CMAKE_Fortran_FLAGS}")
endif ()
get_directory_property(defs COMPILE_DEFINITIONS)
message ("-- Definitions: ${defs}")
message ("-- Debug definitions: ${DEBUG_DEFS}")

if(MKL)
else()
  list(APPEND dir_utils "blas")
  list(APPEND dir_utils "lapack")
endif()
list(APPEND dir_utils "utils")
list(APPEND dir_utils "inp")

# Libraries and Programs

# Make a list of source directories

#message("Source_roots are ${source_roots}")

foreach (dir ${source_roots})
  file (GLOB dir_src RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${PROJECT_SOURCE_DIR}/${dir}/[^.]*")
  list(REMOVE_ITEM dir_src "config") 
  list(REMOVE_ITEM dir_src "include") 
  if(MKL)
    list(REMOVE_ITEM dir_src "blas") 
    list(REMOVE_ITEM dir_src "lapack") 
  endif()
  foreach (srcdir ${dir_src})
#    message("Checking files in ${srcdir}")
    set (fildir ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
    file (GLOB_RECURSE srcfiles RELATIVE "${PROJECT_SOURCE_DIR}" "${fildir}/*.F" "${fildir}/*.f" "${fildir}/*.c")
    file (GLOB_RECURSE prgfiles RELATIVE "${PROJECT_SOURCE_DIR}" "${fildir}/${srcdir}.F" "${fildir}/${srcdir}.f" "${fildir}/${srcdir}.c")
    foreach (prg ${prgfiles})
      list (REMOVE_ITEM srcfiles ${prg})
    endforeach(prg)
#    message("Files found   : ${srcfiles}")
#    message("Programs found: ${prgfiles}")
#    message("Trying ${PROJECT_SOURCE_DIR}/${dir}/${srcdir}")
    set (SRCDIRPATH ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
#    message("Trying files in ${SRCDIRPATH}")
    list (LENGTH srcfiles numsrc)
    list (LENGTH prgfiles numprg)
    message ("-- Source directory ${SRCDIRPATH}: ${numsrc} library files")
    if(numsrc GREATER 0)
#      message ("Add library ${srcdir} : ${srcfiles}")
      add_library(${srcdir} STATIC ${srcfiles})
      list(APPEND libraries ${srcdir})
      list(APPEND dir_libs ${srcdir})
    endif()
  endforeach(srcdir)
#  list(APPEND source_dirs ${dir_src})
endforeach(dir)

message("-- Libraries constructed: ${libraries}")

foreach (dir ${dir_utils})
  list(REMOVE_ITEM dir_libs ${dir})
endforeach(dir)
foreach (dir ${dir_libs})
  add_dependencies(${dir} ${dir_utils})
  target_link_libraries (${dir} ${dir_utils})
  message("-- Dependencies ${dir}: ${dir_utils}")
endforeach(dir)
#message("programs  constructed: ${programs}")
#message("dir_utils: ${dir_utils}")
#message("dir_libs : ${dir_libs}")



foreach (dir ${source_roots})
  file (GLOB dir_src RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${PROJECT_SOURCE_DIR}/${dir}/[^.]*")
  list(REMOVE_ITEM dir_src "config") 
  list(REMOVE_ITEM dir_src "include") 
  foreach (srcdir ${dir_src})
#    message("Checking files in ${srcdir}")
    set (fildir ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
    file (GLOB_RECURSE srcfiles RELATIVE "${PROJECT_SOURCE_DIR}" "${fildir}/*.F" "${fildir}/*.f" "${fildir}/*.c")
    file (GLOB_RECURSE prgfiles RELATIVE "${PROJECT_SOURCE_DIR}" "${fildir}/${srcdir}.F" "${fildir}/${srcdir}.f" "${fildir}/${srcdir}.c")
    foreach (prg ${prgfiles})
      list (REMOVE_ITEM srcfiles ${prg})
    endforeach(prg)
#    message("Files found   : ${srcfiles}")
#    message("Programs found: ${prgfiles}")
#    message("Trying ${PROJECT_SOURCE_DIR}/${dir}/${srcdir}")
    set (SRCDIRPATH ${PROJECT_SOURCE_DIR}/${dir}/${srcdir})
#    message("Trying files in ${SRCDIRPATH}")
    list (LENGTH srcfiles numsrc)
    list (LENGTH prgfiles numprg)
    message ("-- Source directory ${SRCDIRPATH}: ${numprg} program files")
    if(numprg GREATER 0)
#      message ("Add library ${srcdir} : ${srcfiles}")
#      add_dependencies(${scrdir} ${libraries})
      add_executable(${srcdir} ${prgfiles})
      message("-- Added executable ${srcdir} to ${prgfiles}")
      add_dependencies (${srcdir} ${libraries})
      target_link_libraries (${srcdir} ${libraries} ${EXTERNAL_LIBRARIES})
      list(APPEND programs ${srcdir})
      list(APPEND dir_prgs ${srcdir})
    endif()
  endforeach(srcdir)
#  list(APPEND source_dirs ${dir_src})
endforeach(dir)

message("-- Programs constructed: ${programs}")

enable_testing()

include(CTest)

message("Configuring testing:")
foreach (dir ${test_roots})
  file (GLOB dir_test RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${PROJECT_SOURCE_DIR}/${dir}/[^.]*")
  foreach (testdir ${dir_test})
#    message("Test Directory " ${testdir})
    set (fildir ${PROJECT_SOURCE_DIR}/${dir}/${testdir})
    file (GLOB_RECURSE testfiles RELATIVE "${fildir}" "${fildir}/*.run")
    foreach (test ${testfiles})
      string(REPLACE ".run" "" testroot ${test})
#      message("PROJECT_SOURCE_DIR " ${PROJECT_SOURCE_DIR} )
#      message("add_test NAME " ${testroot} )
#      message("WORKING_DIRECTORY " ${fildir})
      message("-- add_test " ${fildir} "/" ${testroot})
      add_test (NAME ${testroot} WORKING_DIRECTORY ${fildir} COMMAND python ${PROJECT_SOURCE_DIR}/scripts/test.py ${testroot})
    endforeach(test)
  endforeach(testdir)
endforeach(dir)

#foreach (prgm ${programs})
#  add_dependencies (${prgm} ${libraries})
#  message("-- Dependencies ${prgm}: ${libraries}")
#  target_link_libraries (${prgm} ${libraries})
#endforeach (prgm)

message ("Finalization:")
set (LINK_DIRECTORIES ${LIBRARY_OUTPUT_PATH})
message ("-- Link directories: ${LINK_DIRECTORIES}")
#message("Source directories are ${source_dirs}")

#foreach (dir ${source_dirs})
#  set (srcdir ${PROJECT_SOURCE_DIR}/src/
#  file(GLOB sources RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${dir}/[^.]*")
#endforeach(dir)
#HIERO

# Define all directories that contain utilities

#list(APPEND dir_utils "blas")
#list(APPEND dir_utils "lapack")
#list(APPEND dir_utils "utils")

# Autodiscover program and utility directories

#foreach (dir ${source_roots})
#        # find all files not starting with a dot
#        file (GLOB dir_progs RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${dir}/[^.]*")
##        file (GLOB dir_utils RELATIVE "${PROJECT_SOURCE_DIR}/${dir}" "${dir}/[^.]*_util")

#	message ("Testing ${PROJECT_SOURCE_DIR}/${dir}")
	
#        # programs are those that do not end in "_util"
#        if (dir_utils)
#                list (REMOVE_ITEM dir_progs ${dir_utils})
#        endif ()
#        # and consider only directories not beginning with uppercase
#        foreach (item ${dir_progs})
#                if (NOT IS_DIRECTORY "${PROJECT_SOURCE_DIR}/${dir}/${item}")
#                        list (REMOVE_ITEM dir_progs ${item})
#                endif ()
#                if (item MATCHES "^[A-Z]")
#                        list (REMOVE_ITEM dir_progs ${item})
#                endif ()
#        endforeach ()
##        foreach (item ${dir_utils})
##                if (NOT IS_DIRECTORY "${PROJECT_SOURCE_DIR}/${dir}/${item}")
##                        list (REMOVE_ITEM dir_utils ${item})
##                endif ()
##                if (item MATCHES "^[A-Z]")
##                        list (REMOVE_ITEM dir_utils ${item})
##                endif ()
##        endforeach ()
#list (APPEND utils ${dir_utils})
#list (REMOVE_ITEM dir_progs "include")
#list (APPEND progs ${dir_progs})
#list (REMOVE_ITEM progs "include")
#	message ("dir_progs is ${dir_progs}")
#	message ("dir_utils is ${dir_utils}")
#endforeach ()
##list (REMOVE_DUPLICATES utils)
#list (REMOVE_DUPLICATES progs)
#list (REMOVE_ITEM dir_progs "include")
#	message ("dir_progs is ${dir_progs}")


# set the location of the sources for the utilities

#foreach (util ${utils})
#        find_source(${util})
#        if (NOT DEFINED ${util}_src)
#                list (REMOVE_ITEM utils ${util})
#                message (WARNING "\"${util}\" not found in source directories, skipping")
#        endif ()
#endforeach ()

# finally, sort the list

#list (SORT utils)

# set the location of the sources for the programs

#foreach (prog ${progs})
#        find_source(${prog})
#        if (NOT DEFINED ${prog}_src)
#                list (REMOVE_ITEM progs ${prog})
#                message (WARNING "\"${prog}\" not found in source directories, skipping")
#        endif ()
#endforeach ()

#message ("progs is ${progs}")

# finally, sort the list

#list (SORT progs)

#message ("utils is ${utils}")
#message ("progs is ${progs}")

# Disable some warnings that are buggy or on code we don't want to modify

#if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
#        set (netlib_nowarn_flags " -w")
#        set (ext_nowarn_flags " -w")
#elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
#        set (netlib_nowarn_flags " -w")
#        set (ext_nowarn_flags " -w")
#elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "SunPro")
#        set (netlib_nowarn_flags " -erroff=%all")
#        set (ext_nowarn_flags " -erroff=%all")
#elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
#        set (netlib_nowarn_flags " -Kieee -w")
#        set (ext_nowarn_flags " -w")
#elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "NAG")
#        set (netlib_nowarn_flags " -w")
#        set (ext_nowarn_flags " -w")
#endif ()

# Blacklisted sources

#if (DEFINED BLACKLIST)
#        foreach (module ${BLACKLIST})
#                file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR}
#                        ${${module}_src}/*.f
#                        ${${module}_src}/*.f90
#                        )
#                set_property (
#                        SOURCE ${source}
#                        APPEND_STRING
#                        PROPERTY
#                        COMPILE_FLAGS " ${BLACKFLAGS}"
#                        )
#        endforeach ()
#        message ("Blacklisting enables, someone doesn't like certain modules...")
#        message ("-- modules: ${BLACKLIST}")
#        message ("-- *additional* flags: ${BLACKFLAGS}")
#endif ()

# Build the utilities and the GronOR library

#set (util "Modules")
#file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${GRONOR_DIR}/src/${util}/*.f ${GRONOR_DIR}/src/${util}/*.f90)
#message ("TEST SOURCE is ${source}")
#add_library (mods_obj OBJECT ${source})
#list (APPEND ALL_SOURCES ${source})
#if (DEFINED ${util}_defs)
#        set_target_properties (mods_obj PROPERTIES COMPILE_DEFINITIONS "${${util}_defs}")
#endif ()
#
#     if (BUILD_SHARED_LIBS)
#        set_target_properties (mods_obj PROPERTIES
#                POSITION_INDEPENDENT_CODE "TRUE"
#        )
#endif ()
#if (SINGLE_MOD_DIR)
#        set_target_properties (mods_obj PROPERTIES
#                Fortran_MODULE_DIRECTORY ${MAIN_MOD_DIR}/_single)
#else ()
#        set_target_properties (mods_obj PROPERTIES
#                Fortran_MODULE_DIRECTORY ${MAIN_MOD_DIR})
#endif ()
#target_include_directories (mods_obj PRIVATE
#        ${GRONOR_DIR}/src/${util}
#        )
#list (APPEND utils_obj mods_obj)
#list (APPEND objlibs "$<TARGET_OBJECTS:mods_obj>")


# now go through utility list

#foreach (util ${utils})
#  message("util is ${util}")
#  message("Search ${PROJECT_SOURCE_DIR}/src/${util}")
#  set (UTILS_DIR ${PROJECT_SOURCE_DIR}/src/${util})
#  message("PROJECT_SOURCE_DIR/src is ${PROJECT_SOURCE_DIR}/${source_roots}/${util}")
#  file (GLOB_RECURSE util_sources RELATIVE ${PROJECT_SOURCE_DIR} ${UTILS_DIR}/*.f ${UTILS_DIR}/*.F ${UTILS_DIR}/*.c)
#  message("srcs is ${util_sources}")
#  add_library(${util} SHARED ${util_sources})
#endforeach(util)

#foreach (util ${utils})
#        file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${${util}_src}/*.f ${${util}_src}/*.f90 ${${util}_src}/*.c)
#        list (APPEND ALL_SOURCES ${source})
#        # ignore empty directories
#        if (NOT "${source}" STREQUAL "")
#                add_library (${util}_obj OBJECT ${source})
#                if (SINGLE_MOD_DIR)
#                        set (mod_dir ${MAIN_MOD_DIR}/_single)
#                else ()
#                        set (mod_dir ${MAIN_MOD_DIR}/${util})
#                endif ()
#                set_target_properties (${util}_obj PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})
#
#                if (BUILD_SHARED_LIBS)
#                        set_target_properties (${util}_obj PROPERTIES
#                                POSITION_INDEPENDENT_CODE "TRUE")
#                endif ()
#
#                # utility-specific compile definitions
#                if (DEFINED ${util}_defs)
#                        set_target_properties (${util}_obj PROPERTIES COMPILE_DEFINITIONS "${${util}_defs}")
#                endif ()
#                # utility-specific include directories
#                list (APPEND ${util}_incs ${PROJECT_SOURCE_DIR}/${${util}_src} ${mod_dir})
#                foreach (inc ${${util}_incs})
#                        target_include_directories (${util}_obj PRIVATE ${inc})
#                endforeach ()
#
#                add_dependencies (${util}_obj mods_obj)
#
#                list (APPEND utils_obj ${util}_obj)
#                list (APPEND objlibs "$<TARGET_OBJECTS:${util}_obj>")
#                #add_library (${util} $<TARGET_OBJECTS:${util}_obj>)
#
#                # if (EXISTS ${PROJECT_SOURCE_DIR}/${${util}_src}/custom.cmake)
#                #         include(${PROJECT_SOURCE_DIR}/${${util}_src}/custom.cmake)
#                # endif ()
#        endif ()
#endforeach (util)
#
## build an actual GronOR library
#
#add_library (gronor ${objlibs})
#set (GRONOR_LIBRARIES "gronor")
#if (BUILD_SHARED_LIBS)
#        # flatten namespace on Mac OS X
#        if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
#                set_target_properties (molcas PROPERTIES
#                        LINK_FLAGS "-Wl,-flat_namespace")
#        endif ()
#        # also create a static library!
#        if (BUILD_STATIC_LIB)
#                add_library (gronor_static STATIC ${objlibs})
#                list (APPEND GRONOR_LIBRARIES "gronor_static")
#        endif ()
#endif ()
#
# if necessary, build BLAS/LAPACK
#
#if (LINALG STREQUAL "Internal")
#        set (util "LinAlg_internal")
#        if (NOT DEFINED ${util}_src)
#                message (FATAL_ERROR "\"${util}\" not found in source directories")
#        endif ()
#        foreach (lib "lapack" "blas")
#                set (source ${${util}_src}/${lib}.f)
#                add_library (${lib} ${source})
#                list (APPEND ALL_SOURCES ${source})
#                list (APPEND ${util}_incs ${internal_linalg})
#                foreach (inc ${${util}_incs})
#                        target_include_directories (${lib} PRIVATE ${inc})
#                endforeach ()
#                if (BUILD_SHARED_LIBS)
#                        set_target_properties (${lib} PROPERTIES
#                                POSITION_INDEPENDENT_CODE "TRUE")
#                endif ()
#        endforeach ()
#endif ()

# Build the program modules


#add_custom_target(only_objs DEPENDS ${utils_obj})

#foreach (prog ${progs})
#        file (GLOB source RELATIVE ${PROJECT_SOURCE_DIR} ${${prog}_src}/*.f ${${prog}_src}/*.f90 ${${prog}_src}/*.F90 ${${prog}_src}/*.c)
#        list (APPEND ALL_SOURCES ${source})
#        # ignore empty programs
#        if ("${source}" STREQUAL "")
#                list (REMOVE_ITEM progs ${prog})
#        else ()
#                # FIXME: the checks below remove any main.* file, which is needed in
#                # order to be able to make a library out of the program sources. This
#                # is required for supermodules.
#                #------------------------------------------------------------------
#                foreach (lext "f" "f90" "c")
#                        if (EXISTS ${PROJECT_SOURCE_DIR}/${${prog}_src}/main.${lext})
#                                list (REMOVE_ITEM source ${${prog}_src}/main.${lext})
#                                set (${prog}_main ${${prog}_src}/main.${lext})
#                        endif ()
#                endforeach ()

#                if (NOT "${source}" STREQUAL "")
#                        add_library (${prog} OBJECT ${source})
#                        if (SINGLE_MOD_DIR)
#                                set (mod_dir ${MAIN_MOD_DIR}/_single)
#                        else ()
#                                set (mod_dir ${MAIN_MOD_DIR}/${prog})
#                        endif ()
#                        set_target_properties (${prog} PROPERTIES Fortran_MODULE_DIRECTORY ${mod_dir})

#                        # programs in "proglib" have their libraries compiled as utils
#                        if (BUILD_SHARED_LIBS AND ";${proglib};" MATCHES ";${prog};")
#                                set_target_properties (${prog} PROPERTIES
#                                        POSITION_INDEPENDENT_CODE "TRUE")
#                        endif ()

#                        # program-specific compile definitions
#                        if (DEFINED ${prog}_defs)
#                                set_target_properties (${prog} PROPERTIES COMPILE_DEFINITIONS "${${prog}_defs}")
#                        endif ()
#                        # program-specific include directories
#                        list (APPEND ${prog}_incs ${PROJECT_SOURCE_DIR}/${${prog}_src} ${mod_dir})
#                        foreach (inc ${${prog}_incs})
#                                target_include_directories (${prog} PRIVATE ${inc})
#                        endforeach ()

#                        add_dependencies (${prog} ${utils_obj})
#                        add_dependencies (only_objs ${prog})
#                endif ()
#        endif ()
#endforeach ()

## set file-specific flags
#if (BOUNDS AND (${FFLAGS_GNU_BOUNDS} MATCHES "-fsanitize=address"))
#        set_property(SOURCE ${rassi_src}/prprop.f APPEND_STRING PROPERTY COMPILE_FLAGS " --param=max-vartrack-size=120000000")
#endif ()

## set up proper RPATH for executables, use RUNPATH if on Linux
#set (CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
#set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#foreach (prog ${progs})
#        # FIXME: Actual program target only depends on main file now, but if
#        # supermodules are gone, the dependency should be on ${source} instead.
#        #----------------------------------------------------------------------
#        # programs in "proglib" have their library included in libmolcas
#        if (TARGET ${prog} AND NOT ";${proglib};" MATCHES ";${prog};")
#                list (APPEND ${prog}_deps $<TARGET_OBJECTS:${prog}>)
#        endif ()
#        add_executable (${prog}.exe ${${prog}_main} ${${prog}_deps})
#        list (APPEND progs_exe ${prog}.exe)
#        if (SINGLE_MOD_DIR)
#                set_target_properties (${prog}.exe PROPERTIES
#                        Fortran_MODULE_DIRECTORY ${MAIN_MOD_DIR}/_single)
#        else ()
#                set_target_properties (${prog}.exe PROPERTIES
#                        Fortran_MODULE_DIRECTORY ${MAIN_MOD_DIR}/${prog})
#        endif ()
#        foreach (inc ${${prog}_incs})
#                target_include_directories (${prog}.exe PRIVATE ${inc})
#        endforeach ()

#        target_link_libraries (${prog}.exe molcas)
#        target_link_libraries (${prog}.exe ${EXTERNAL_LIBRARIES})
#        if (EXTERNAL_PROJECTS)
#                add_dependencies (${prog}.exe ${EXTERNAL_PROJECTS})
#        endif ()

#        if (CUBLAS)
#                set_target_properties (${prog}.exe PROPERTIES LINKER_LANGUAGE "CXX")
#        endif ()

#        list (APPEND PROGRAM_EXECUTABLES ${prog}.exe)
#endforeach (prog)
