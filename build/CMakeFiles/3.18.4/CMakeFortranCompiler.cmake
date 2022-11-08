set(CMAKE_Fortran_COMPILER "/home/qq/coen/bin/nv/bin/mpifort")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "PGI")
set(CMAKE_Fortran_COMPILER_VERSION "20.9.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "Linux")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "")
set(CMAKE_RANLIB "/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "")
set(CMAKE_COMPILER_IS_GNUG77 )
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "mpifort;mpi;nvf;nvomp;dl;nvhpcatm;atomic;pthread;nvcpumath;nsnvc;nvc;rt;pthread;m;gcc;c;gcc;gcc_s")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/home/qq/coen/bin/nv/lib;/opt/nvidia/hpc_sdk/Linux_x86_64/20.9/compilers/lib;/usr/lib64;/usr/lib/gcc/x86_64-redhat-linux/4.8.5")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
