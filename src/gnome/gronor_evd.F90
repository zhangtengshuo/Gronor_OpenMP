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
!! Cofactor matrix evaluation and factorization
!!
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!

subroutine gronor_evd()
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_solvers

  ! library specific modules

#ifdef MKL
  use mkl_solver
#endif

#ifdef LAPACK  
  use lapack_solver
#else
#ifdef MAGMA  
  use magma_solver
#endif
#endif
  
#ifdef CUSOLVER
  use cusolverDn
  use cuda_cusolver
  use cudafor
#endif

  ! variable declarations

  implicit none

  external :: tred2,tql2
#ifdef MKL
  external :: dsyevd
#endif
  
  integer :: i
  integer :: ierr

  ! library specific declarations

#ifdef CUSOLVER
  character (len=1), target :: jobu, jobvt
#endif
  
  ! ========== EISPACK =========

  if(ev_solver.eq.SOLVER_EISPACK) then
#ifdef ACC
!$acc kernels present(sdiag)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop
#else
!$omp target teams distribute parallel do
#endif
#endif
      do i=1,nelecs
        sdiag(i)=0.0d0
      enddo
#ifdef ACC
!$acc end kernels
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
    if(iamacc.eq.1) then
#ifdef ACC
!$acc update host (a)
#endif
#ifdef OMPTGT
!$omp target update from(a)
#endif    
    endif
    call tred2(nelecs,nelecs,a,nelecs,nelecs,diag, &
        nelecs,sdiag,nelecs,a,nelecs,nelecs)
    call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs, &
        a,nelecs,nelecs,ierr)
    if(iamacc.eq.1) then
#ifdef ACC
!$acc update device (ev,u,w)
#endif
#ifdef OMPTGT
!$omp target update to(ev,u,w)
#endif
    endif
  endif
  
! ============ MKL ===========
#ifdef MKL
  if(ev_solver.eq.SOLVER_MKL.or.ev_solver.eq.SOLVER_MKLD) then
    if(iamacc.eq.1) then
#ifdef ACC
!$acc update host (a)
#endif
#ifdef OMPTGT
!$omp target update from(a)
#endif    
    endif
    ndimm=nelecs
    if(ev_solver.eq.SOLVER_MKL) then
      call dsyevd('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m,workspace_i,lworki,ierr)
    endif
    if(ev_solver.eq.SOLVER_MKLD) then
      call dsyevd('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m,workspace_i,lworki,ierr)
    endif
    if(iamacc.eq.1) then
#ifdef ACC
!$acc update device (ev,u,w)
#endif
#ifdef OMPTGT
!$omp target update to(ev,u,w)
#endif
    endif
  endif
#endif 
  
! ============ LAPACK ===========
#ifdef LAPACK
  if(ev_solver.eq.SOLVER_LAPACK.or.ev_solver.eq.SOLVER_LAPACKD) then
    if(iamacc.eq.1) then
#ifdef ACC
!$acc update host (a)
#endif
#ifdef OMPTGT
!$omp target update from(a)
#endif    
    endif
    ndimm=nelecs
    if(ev_solver.eq.SOLVER_LAPACK) then
      call dsyev('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m,ierr)
    elseif(ev_solver.eq.SOLVER_LAPACKD) then
      call dsyevd('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m, &
          workspace_i,lworki,ierr)
    endif
    if(iamacc.eq.1) then
#ifdef ACC
!$acc update device (ev,u,w)
#endif
#ifdef OMPTGT
!$omp target update to(ev,u,w)
#endif
    endif
  endif
#endif 
  
! ============ MAGMA ===========
#ifdef MAGMA
  if(ev_solver.eq.SOLVER_MAGMA) then
    if(iamacc.eq.1) then
#ifdef ACC
!$acc data create(workspace_d,workspace_i)
!$acc wait
!$acc host_data use_device(a,diag,workspace_d,workspace_i)
#endif
#ifdef OMPTGT
!!!!!$omp target data use_device_addr(a,ev,u,w,dev_info_d,workspace_d)
#endif    
      ndimm=nelecs
      call magma_dsyevd_gpu('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m,workspace_i,lworki,ierr)
#ifdef ACC
!$acc end host_data
!$acc wait
!$acc end data   
#endif
#ifdef OMPTGT
!!!!!$omp end target data
#endif
    else        
      ndimm=nelecs
      call magma_dsyevd('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m,workspace_i,lworki,ierr)
    endif
  endif
#endif 
  
  ! ========= CUSOLVER =========

#ifdef CUSOLVER
  if(ev_solver.eq.SOLVER_CUSOLVER) then
    ndim=nelecs
    mdim=mbasel
    jobu = 'A'  ! all m columns of U
    jobvt= 'A'  ! all m columns of VT
#ifdef ACC
!$acc data copy(dev_info_d) create(workspace_d)
!$acc wait
!$acc host_data use_device(a,diag,dev_info_d,workspace_d)
#endif
#ifdef OMPTGT
!$omp target data use_device_addr(a,diag,dev_info_d,workspace_d)
#endif
    cusolver_status = cusolverDnDsyevd(cusolver_handle, &
        CUSOLVER_EIG_MODE_NOVECTOR,CUBLAS_FILL_MODE_LOWER, &
        ndim,a,ndim,diag,workspace_d,lwork2,dev_info_d)
#ifdef ACC
!$acc end host_data
!$acc wait
!$acc end data   
#endif
#ifdef OMPTGT
!$omp end target data
#endif
    if(cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
        write(*,*) 'cusolverDnDsyevd failed',cusolver_status
  endif
#endif
  
  ! ======== CUSOLVERJ =========

#ifdef CUSOLVERJ  
  if(ev_solver.eq.SOLVER_CUSOLVERJ) then
    ndim=nelecs
    mdim=mbasel
    jobz = CUSOLVER_EIG_MODE_NOVECTOR
    uplo = CUBLAS_FILL_MODE_LOWER

#ifdef ACC
!$acc data copy(dev_info_d,syevj_params) create(workspace_d)
!$acc host_data use_device(a,diag,dev_info_d,workspace_d)
#endif
#ifdef OMPTGT
!$omp target data use_device_addr(a,ev,u,w,dev_info_d,workspace_d)
#endif
    cusolver_status = cusolverDnDsyevj &
        (cusolver_handle, jobz, uplo, ndim,a,ndim,diag, &
        workspace_d,lwork2,dev_info_d,syevj_params)
#ifdef ACC
!$acc end host_data
!$acc end data   
#endif
#ifdef OMPTGT
!$omp end target data
#endif
    cusolver_status=cudaDeviceSynchronize()
    if(cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
        write(*,*) 'cusolverDnDsyevj failed',cusolver_status
  endif
#endif
  
  ! ======== HIPSOLVER =========

#ifdef HIPSOLVER
  if(ev_solver.eq.SOLVER_HIPSOLVER) then
    ndim=nelecs
    mdim=mbasel
  endif
#endif

  ! ======= HIPSOLVERJ =========
  
#ifdef HIPSOLVERJ
  if(ev_solver.eq.SOLVER_HIPSOLVERJ) then
    ndim=nelecs
    mdim=mbasel
  endif
#endif

  ! ======== ROCSOLVER =========
  
#ifdef ROCSOLVER
  if(ev_solver.eq.SOLVER_ROCSOLVER) then
    ndim=nelecs
    mdim=mbasel
    istatus=rocsolver_dsyev(rocsolver_handle,evect,uplo, &
        m,c_loc(at),m,c_loc(dt),c_loc(et),rocinfo)
    call hipCheck(hipDeviceSynchronize())
  endif
#endif

  ! ======= ROCSOLVERJ =========
  
#ifdef ROCSOLVERJ
  if(ev_solver.eq.SOLVER_ROCSOLVERJ) then
    ndim=nelecs
    mdim=mbasel
    
  endif
#endif
    

  return
end subroutine gronor_evd



subroutine gronor_evd_omp()
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_solvers

  ! library specific modules

#ifdef MKL
  use mkl_solver
#endif
#ifdef LAPACK
  use lapack_solver
#endif

  ! variable declarations

  implicit none

  external :: tred2,tql2
#ifdef MKL
  external :: dsyevd
#endif
#ifdef LAPACK
  external :: dsyevd
#endif
  
  integer :: i
  integer :: ierr

  ! ========== EISPACK =========

  if(ev_solver.eq.SOLVER_EISPACK) then
#ifdef OMP
!$omp parallel do shared(sdiag)
#endif
    do i=1,nelecs
      sdiag(i)=0.0d0
    enddo
#ifdef OMP
!$omp end parallel do
#endif
    call tred2(nelecs,nelecs,a,nelecs,nelecs,diag, &
        nelecs,sdiag,nelecs,a,nelecs,nelecs)
    call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,a,nelecs, &
        nelecs,ierr)
  endif
  
  ! ============ MKL ===========
  
#ifdef MKL
  if(ev_solver.eq.SOLVER_MKL) then
    ndimm=nelecs
    call dsyev('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m,ierr)
  elseif(ev_solver.eq.SOLVER_MKLD) then
    ndimm=nelecs
    call dsyevd('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m, &
        workspace_i,lworki,ierr)
  endif
#endif 
  
  ! ============ LAPACK ===========
  
#ifdef LAPACK
  if(ev_solver.eq.SOLVER_LAPACK) then
    ndimm=nelecs
    call dsyev('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m,ierr)
  elseif(ev_solver.eq.SOLVER_LAPACKD) then
    ndimm=nelecs
    call dsyevd('N','L',ndimm,a,nelecs,diag,workspace_d,lwork1m, &
        workspace_i,lworki,ierr)
  endif
#endif 

  return
end subroutine gronor_evd_omp
