!
!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ bases, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!> @brief
!! Routine that provides all possible calls to Eigensolver library routines
!!
!! @author  T. P. Straatsma, ORNL
!! @date    2025
!!

subroutine gronor_evd()

  !> Routine that provides all possible calls to Eigensolver library routines
  !! including routines executed on the CPU or on GPU accelerators
  !!
  !! For CPU-executed ranks the vectors, matrices, and workspaces are expected in CPU memory
  !! For GPU-accelerated ranks the vectors, matrices, and workspaces are expected in GPU memory
  !! GPU-accelerated ranks can use CPU-based routines, in which case this routine will copy the
  !! required input data from the GPU to the CPU and return result data from CPU to GPU memory
  !!
  !! The library routine that will be used is determined by variable ev_solver which is set in
  !! routine gronor_solver_init based on the input provide in gronor_input. Currently implemented
  !! options for ev_solver are:
  !!
  !! SOLVER_EISPACK      tred2 and tql2 as provided in the source code will run on the CPU
  !! SOLVER_MKL          dsyev from external an Intel MKL library will run on the CPU
  !! SOLVER_MKLD         dsyevd from external an Intel MKL library will run on the CPU
  !!
  !! SOLVER_MKL and SOLVER LAPACK cannot be available in the same executable because of the name conflict
  !!
  !! The boolean flag levcpu is set in gronor_solver_init and indicates the solver routine runs
  !! on the CPU (levcpu=.true.) or the GPU (levcpu=.false.)
  !! The integer iamacc specifies is set in gronor_main and indicates if rank
  !! is CPU resident (iamacc=0) or GPU resident (iamacc=1)
  !!
  !! The integer workspace_i and double workspace_d workspaces are allocated with total sizes
  !! len_work_int and len_work_dbl, respectively. These workspace are shared between all possible
  !! solver routines, including the eigensolver in gronor_svd. Separate copies exist in CPU memory
  !! and in GPU memory. The maximum required workspace sizes are set in gronor_solver_init
  
  use cidef
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_solvers
  use iso_c_binding

  ! library specific modules

#ifdef MKL
  use mkl_solver
#endif

  ! variable declarations

  implicit none

  external :: tred2,tql2
#ifdef MKL
  external :: dsyevd
#endif
  
  integer :: i,j
  integer :: ierr

  ! library specific declarations

  if(iamacc.eq.1) then
     if(levcpu) then
#ifdef ACC
!$acc update host (a)
#endif
      do i=1,nelecs
        sdiag(i)=0.0d0
      enddo

    endif
#ifdef ACC
!$acc kernels present(sdiag)
#endif

      do i=1,nelecs
        sdiag(i)=0.0d0
      enddo
#ifdef ACC
!$acc end kernels
#endif

   else

      do i=1,nelecs
        sdiag(i)=0.0d0
      enddo

   endif

  ! ========== EISPACK =========

  if(ev_solver.eq.SOLVER_EISPACK) then
    call tred2(nelecs,nelecs,a,nelecs,nelecs,diag, &
        nelecs,sdiag,nelecs,a,nelecs,nelecs)
    call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs, &
        a,nelecs,nelecs,ierr)
  endif
  
! ============ MKL ===========
#ifdef MKL
  if(ev_solver.eq.SOLVER_MKL.or.ev_solver.eq.SOLVER_MKLD) then
    ndimm=nelecs
    if(ev_solver.eq.SOLVER_MKL) then
      call dsyev('N','L',ndimm,a,nelecs,diag,workspace_d,len_work_dbl,ierr)
    endif
    if(ev_solver.eq.SOLVER_MKLD) then
      call dsyevd('N','L',ndimm,a,nelecs,diag,workspace_d,len_work_dbl, &
          workspace_i,len_work_int,ierr)
    endif
 endif
#endif 

  if(iamacc.eq.1.and.levcpu) then
#ifdef ACC
!$acc update device (diag)
#endif
  endif

  return
end subroutine gronor_evd
