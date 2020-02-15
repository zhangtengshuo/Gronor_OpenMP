      program gronor
 
      use mpi
      use cidist

      implicit none
!tps#include 'mpif.h'
      integer (kind=4) :: ierror,ierr

      call gronor_environment()

      call gronor_main()

      if(me.eq.master) call mpi_abort(MPI_COMM_WORLD,ierror,ierr)

!      call mpi_finalize(ierr)

      end
