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

!>    Driver routine for worker ranks
!!    @brief Driver for calculation Hamiltonian matrix elements on worker ranks
!!    @author T. P. Straatsma (ORNL)subroutine gronor_solver_init()

subroutine gronor_solver_init(ntemp)
  
  use mpi
  use cidef
  use cidist
  use gnome_integrals
  use gnome_data
  use gnome_parameters
  use gnome_solvers
  use iso_c_binding
  
#ifdef CUSOLVER
  use cusolverDn
  use cuda_cusolver
  use cudafor
#endif

#ifdef MKL
  use mkl_solver
#endif
  
  implicit none
  
#ifdef MKL
  external :: dgesvd,dsyevd
  integer (kind=4) :: ierr
#endif

  integer :: ntemp
  character(len=255) :: string
  
  integer (kind=8) :: lworki,lwork1m,lwork2m
  integer (kind=4) :: lwork1,lwork2
  
  real(kind=8) :: worksize(2),worksize2(2)
  integer (kind=4) :: iworksize(2)

  nelecs=ntemp

  len_work_int=0
  len_work_dbl=0
  len_work2_dbl=0

! Cusolver initialization for the svd
  
  if(idbg.gt.50) then
    call swatch(date,time)
    write(lfndbg,'(a,1x,a,a,2i4)') date(1:8),time(1:8)," Solver init for ",sv_solver,ev_solver
    flush(lfndbg)
  endif

  lsvcpu=.false.
  levcpu=.false.
  lsvtrns=.true.
  
  if(sv_solver.eq.SOLVER_EISPACK) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_MKL) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_MKLD) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_MKLJ) lsvcpu=.true.

  if(sv_solver.eq.SOLVER_EISPACK) lsvtrns=.false.
  
  if(ev_solver.eq.SOLVER_EISPACK) levcpu=.true.
  if(ev_solver.eq.SOLVER_MKL) levcpu=.true.
  if(ev_solver.eq.SOLVER_MKLD) levcpu=.true.
  if(ev_solver.eq.SOLVER_MKLJ) levcpu=.true.

  ! MKL initialization
  
#ifdef MKL
    ndimm=nelecs
    mdimm=mbasel
    lwork1m=-1
    lwork2m=-1
    lworki=-1
    if(sv_solver.eq.SOLVER_MKL) then
      call dgesvd('All','All',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,worksize,lwork1m,ierr)
      lwork1m=int(worksize(1))
    endif
    if(sv_solver.eq.SOLVER_MKLD) then
      call dgesdd('All',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,worksize,lwork1m,iworksize,ierr)
      lwork1m=int(worksize(1))
      lworki=max(1000,8*nelecs)
    endif
    if(sv_solver.eq.SOLVER_MKLJ) then
!      call dgesvj('L','U','V',ndimm,ndimm,a,ndimm,ev,ndimm,w,ndimm,workspace_d,lwork1m,ierr)
      lwork1m=max(int(worksize(1)),6,2*nelecs,lwork1m)
      lworki=max(int(iworksize(1)),lworki)
    endif
    if(ev_solver.eq.SOLVER_MKL) then
      call dsyev('N','L',ndimm,a,ndimm,w,worksize,lwork2m,ierr)
      lwork2m=int(worksize(1))
    endif
    if(ev_solver.eq.SOLVER_MKLD) then
      call dsyevd('N','L',ndimm,a,ndimm,w,worksize,lwork2m,iworksize,lworki,ierr)
      lwork2m=int(worksize(1))
      lworki=max(1000,int(iworksize(1)),lworki)
    endif
    lwork1m=max(1,lwork1m,lwork2m)
    lworki=max(1,lworki)
    len_work_dbl=max(len_work_dbl,lwork1m)
    len_work_int=max(len_work_int,lworki)
#endif
    
    len_work_dbl=max(1,len_work_dbl)
    len_work_int=max(1,len_work_int)
    len_work2_dbl=max(1,len_work2_dbl)

    allocate(workspace_d(len_work_dbl))
    allocate(workspace2_d(len_work2_dbl))
    allocate(workspace_i(len_work_int))
    allocate(workspace_i4(len_work_int))

    return
  end subroutine gronor_solver_init

subroutine gronor_solver_finalize()

  use mpi
  use inp
  use cidef
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals
  use iso_c_binding
  use iso_fortran_env
  use gnome_solvers

  
  return
end subroutine gronor_solver_finalize

subroutine gronor_solver_create_handle()
  
  use mpi
  use inp
  use cidef
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals
  use iso_c_binding
  use iso_fortran_env
  use gnome_solvers
  
  ! Only accelerated ranks need to define cusolver handles
  
  return
end subroutine gronor_solver_create_handle
