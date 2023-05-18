!> GronOR Main Program
 
      program gronor

      use mpi
      use cidist

      implicit none

      external :: gronor_environment,gronor_main
      
      call gronor_environment()

      call gronor_main()

      end
