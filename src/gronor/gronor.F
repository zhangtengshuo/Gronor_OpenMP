!> GronOR Main Program
 
      program gronor

      use mpi
      use cidist

      implicit none

      call gronor_environment()

      call gronor_main()

      end
