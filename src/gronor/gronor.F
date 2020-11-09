!> GronOR Main Program
      
      program gronor
 
      use mpi
      use cidist

      implicit none
!tps#include 'mpif.h'

      call gronor_environment()

      call gronor_main()

      end
