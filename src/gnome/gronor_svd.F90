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
!! Routine that provides all possible calls to Singular Value Decomposition library routines
!! @author  T. P. Straatsma, ORNL
!! @date    2025
!!

subroutine gronor_svd()

  !> Routine that provides all possible calls to Singular Value Decomposition library routines
  !! including routines executed on the CPU or on GPU accelerators 
  !!
  !! For CPU-executed ranks the vectors, matrices, and workspaces are expected in CPU memory
  !! For GPU-accelerated ranks the vectors, matrices, and workspaces are expected in GPU memory
  !! GPU-accelerated ranks can use CPU-based routines, in which case this routine will copy the
  !! required input data from the GPU to the CPU and return result data from CPU to GPU memory
  !!
  !! The library routine that will be used is determined by variable sv_solver which is set in
  !! routine gronor_solver_init based on the input provide in gronor_input. Currently implemented
  !! options for sv_solver are:
  !!
  !! SOLVER_EISPACK      svd as provided in the source code will run on the CPU
  !! SOLVER_MKL          dgesvd from external an Intel MKL library will run on the CPU
  !! SOLVER_MKLD         dgesdd from external an Intel MKL library will run on the CPU
  !! SOLVER_MKLJ         dgesvj from external an Intel MKL library will run on the CPU
  !!
  !! SOLVER_MKL and SOLVER LAPACK cannot be available in the same executable because of the name conflict
  !!
  !! The boolean flag lsvcpu is set in gronor_solver_init and indicates the solver routine runs
  !! on the CPU (lsvcpu=.true.) or the GPU (lsvcpu=.false.)
  !! The integer iamacc specifies is set in gronor_main and indicates if rank
  !! is CPU resident (iamacc=0) or GPU resident (iamacc=1)
  !!
  !! The integer workspace_i and double workspace_d workspaces are allocated with total sizes
  !! len_work_int and len_work_dbl, respectively. These workspace are shared between all possible
  !! solver routines, including the eigensolver in gronor_evd. Separate copies exist in CPU memory
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

  external :: svd,dgesvd
  
  integer :: i,j
  integer :: ierr
  integer (kind=4) :: istat

  ! library specific declarations

  lwork4=int(len_work_dbl,kind=4)
  
  if(iamacc.eq.1.and.lsvcpu) then
#ifdef ACC
!$acc update host (a)
#endif
  endif

  ! ========== EISPACK =========

  if(sv_solver.eq.SOLVER_EISPACK) then
    call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,ev,nelecs,.true., &
        u,nelecs,nelecs,.true.,w,nelecs,nelecs,ierr,sdiag,nelecs)
  endif
  
  ! ============ MKL ===========
#ifdef MKL
  if(sv_solver.eq.SOLVER_MKL.or.sv_solver.eq.SOLVER_MKLD.or.sv_solver.eq.SOLVER_MKLJ) then
    ndimm=nelecs
    if(sv_solver.eq.SOLVER_MKL) then
      call dgesvd('All','All',ndimm,ndimm,a,ndimm,ev,u,ndimm,wt,ndimm, &
          workspace_d,len_work_dbl,ierr)
    endif
    if(sv_solver.eq.SOLVER_MKLD) then
      call dgesdd('All',ndimm,ndimm,a,ndimm,ev,u,ndimm,wt,ndimm, &
          workspace_d,len_work_dbl,workspace_i,ierr)
    endif
    if(sv_solver.eq.SOLVER_MKLJ) then
      call dgesvj('L','U','V',ndimm,ndimm,a,ndimm,ev,ndimm,wt,ndimm, &
          workspace_d,len_work_dbl,ierr)
    endif
    !lsvtrns
  endif
#endif 
   
!! Evaluate right had matrix from its transpose if solver provided the transpose
  if(lsvtrns) then
    if(lsvcpu) then

    elseif(iamacc.eq.1) then

#ifdef ACC
!$acc kernels present(w,wt)
#endif
  do i=1,nelecs
    do j=1,nelecs
      w(i,j)=wt(j,i)
    enddo
  enddo
#ifdef ACC
!$acc end kernels
#endif
      
    endif
  endif

!! Update device if SVD was performed on the host for accelerated ranks
  if(iamacc.eq.1.and.lsvcpu) then
#ifdef ACC
!$acc update device (ev,u,w)
#endif
  endif

  return  
end subroutine gronor_svd
