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

subroutine gronor_svd()
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

#ifdef CUSOLVER
  use cusolverDn
  use cuda_cusolver
  use cudafor
#endif
#ifdef HIPSOLVER
  use hipvars
  use hipsolver_enums
  use amd_hipsolver
#endif
#ifdef ROCSOLVER
  use rocvars
  use hipfort
  use hipfort_hipmalloc
  use hipfort_rocblas_enums
  use hipfort_rocblas
  use hipfort_rocsolver_enums
  use hipfort_rocsolver
#endif

  ! variable declarations

  implicit none

  external :: svd,dgesvd
  
  integer :: i,j
  integer :: ierr

  ! library specific declarations

#ifdef CUSOLVER
  character (len=1), target :: jobu, jobvt
#endif

#ifdef HIPSOLVER
  integer (kind=1) :: jobu, jobvt
#endif

#ifdef ROCSOLVER
  integer :: jobu, jobvt, jobz
#endif

#ifdef LAPACK
  integer (kind=4) :: lapack_info
#endif


  if(iamacc.eq.1) then
     if(lsvcpu) then
#ifdef ACC
!$acc update host (a)
#endif
#ifdef OMPTGT
!$omp target update from(a)
#endif 
#ifdef OMP
!$omp parallel do shared(sdiag)
#endif
      do i=1,nelecs
        sdiag(i)=0.0d0
      enddo
#ifdef OMP
!$omp end parallel do
#endif   
     endif
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
   else
#ifdef OMP
!$omp parallel do shared(sdiag)
#endif
      do i=1,nelecs
        sdiag(i)=0.0d0
      enddo
#ifdef OMP
!$omp end parallel do
#endif
   endif
   
  ! ========== EISPACK =========

  if(sv_solver.eq.SOLVER_EISPACK) then
    call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,                         &
        &       ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,           &
        &       ierr,sdiag,nelecs)
  endif
  
  ! ============ MKL ===========
#ifdef MKL
  if(sv_solver.eq.SOLVER_MKL.or.sv_solver.eq.SOLVER_MKLD.or.sv_solver.eq.SOLVER_MKLJ) then
    ndimm=nelecs
    if(sv_solver.eq.SOLVER_MKL) then
      call dgesvd('All','All',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm, &
          workspace_d,len_work_dbl,ierr)
    endif
    if(sv_solver.eq.SOLVER_MKLD) then
      call dgesdd('All',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm, &
          workspace_d,len_work_dbl,workspace_i,ierr)
    endif
    if(sv_solver.eq.SOLVER_MKLJ) then
      call dgesvj('L','U','V',ndimm,ndimm,a,ndimm,ev,ndimm,w,ndimm, &
          workspace_d,len_work_dbl,ierr)
    endif
  endif
#endif 
  
  ! ============ LAPACK ===========
#ifdef LAPACK
  if(sv_solver.eq.SOLVER_LAPACK.or.sv_solver.eq.SOLVER_LAPACKD) then
    ndimm=nelecs
    if(sv_solver.eq.SOLVER_LAPACK) then
      call dgesvd('A','A',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm, &
          workspace_d,len_work_dbl,ierr)
    endif
    if(sv_solver.eq.SOLVER_LAPACKD) then
      call dgesdd('All',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm, &
          workspace_d,len_work_dbl,workspace_i,ierr)
    endif
  endif
#endif 
  
  ! ========= CUSOLVER =========

#ifdef CUSOLVER
  if(sv_solver.eq.SOLVER_CUSOLVER) then
    ndim=nelecs
    mdim=mbasel
    jobu = 'A'  ! all m columns of U
    jobvt= 'A'  ! all m columns of VT
#ifdef ACC
!$acc data copy(dev_info_d) create(workspace_d)
!$acc host_data use_device(a,ev,u,w,dev_info_d,workspace_d,rwork)
#endif
#ifdef OMPTGT
!$omp target data use_device_addr(a,ev,u,w,dev_info_d,workspace_d,rwork)
#endif
    cusolver_status=cusolverDnDgesvd(cusolver_handle,jobu,jobvt, &
        ndim,ndim,a,ndim,ev,u,ndim,w,ndim,workspace_d,           &
        int(len_work_dbl,kind=4),rwork,dev_info_d)
#ifdef ACC
!$acc end host_data
!$acc wait
!$acc end data
#endif
#ifdef OMPTGT
!$omp end target data
#endif
    if(cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
        write(*,*) 'cusolverDnDgesvd failed',cusolver_status
#ifdef ACC
!$acc kernels present(temp,u,w,ta)
!$acc loop collapse(2)
#endif
    do i=1,nelecs
      do j=1,nelecs
        temp(i,j)=w(j,i)
      enddo
    enddo
#ifdef ACC
!$acc loop collapse(2)
#endif
    do i=1,nelecs
      do j=1,nelecs
        w(i,j)=temp(i,j)
      enddo
    enddo
#ifdef ACC
!$acc end kernels
#endif

    
  endif
#endif
  
  ! ======== CUSOLVERJ =========

#ifdef CUSOLVERJ  
  if(sv_solver.eq.SOLVER_CUSOLVERJ) then
    ndim=nelecs
    mdim=mbasel
    jobz=CUSOLVER_EIG_MODE_VECTOR
#ifdef ACC
!$acc data copy(dev_info_d,gesvdj_params) create(workspace_d)
!$acc host_data use_device(a,ev,u,w,dev_info_d,workspace_d)
#endif
#ifdef OMPTGT
!$omp target data use_device_addr(a,ev,u,w,dev_info_d,workspace_d,rwork)
#endif
    cusolver_status=cusolverDnDgesvdj(cusolver_handle,jobz,econ, &
               ndim,ndim,a,ndim,ev,u,ndim,w,ndim,workspace_d,    &
               int(len_work_dbl,kind=4),dev_info_d,gesvdj_params)
#ifdef ACC
!$acc end host_data
!$acc end data
#endif
#ifdef OMPTGT
!$omp end target data
#endif
    cusolver_status=cudaDeviceSynchronize()
    if(cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
        write(*,*) 'cusolverDnDgesvdj failed',cusolver_status
    cusolver_status = cusolverDnXgesvdjGetSweeps &
        (cusolver_handle, gesvdj_params, exec_sweeps)
    cusolver_status = cusolverDnXgesvdjGetResidual &
        (cusolver_handle, gesvdj_params, residual)
  endif
#endif
  
  ! ======== HIPSOLVER =========

#ifdef HIPSOLVER
  if(sv_solver.eq.SOLVER_HIPSOLVER) then
    ndim=nelecs
    mdim=mbasel
  endif
#endif

  ! ======= HIPSOLVERJ =========
  
#ifdef HIPSOLVERJ
  if(sv_solver.eq.SOLVER_HIPSOLVERJ) then
    ndim=nelecs
    mdim=mbasel
  endif
#endif

  ! ======== ROCSOLVER =========
  
#ifdef ROCSOLVER
  if(sv_solver.eq.SOLVER_ROCSOLVER) then
    ndim=nelecs
    mdim=mbasel
    call hipcheck(rocsolver_dgesvd(rocsolver_handle, &
        ROCBLAS_SVECT_ALL,ROCBLAS_SVECT_ALL,m,n,c_loc(at),n, &
        c_loc(st),c_loc(ut),n,c_loc(vtt),n,c_loc(work), &
        ROCBLAS_OUTOFPLACE,rocinfo))
    call hipCheck(hipDeviceSynchronize())
  endif
#endif

  ! ======= ROCSOLVERJ =========
  
#ifdef ROCSOLVERJ
  if(sv_solver.eq.SOLVER_ROCSOLVERJ) then
    ndim=nelecs
    mdim=mbasel
  endif
#endif

  

  if(iamacc.eq.1) then
     if(lsvcpu) then
#ifdef ACC
!$acc update device (ev,u,w)
#endif
#ifdef OMPTGT
!$omp target update to(ev,u,w)
#endif
     endif
     if(lsvtrns) then
#ifdef ACC
!$acc kernels present(temp,u,w,ta)
!$acc loop collapse(2)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop
#else
!$omp target teams distribute parallel do
#endif
#endif
        do i=1,nelecs
           do j=1,nelecs
              temp(i,j)=w(i,j)
           enddo
        enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc loop collapse(2)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop
#else
!$omp target teams distribute parallel do
#endif
#endif
        do i=1,nelecs
           do j=1,nelecs
              w(j,i)=temp(i,j)
           enddo
        enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
     endif
  else
     if(lsvtrns) then
#ifdef OMP
!$omp parallel shared(temp,w,nelecs)
!$omp do collapse(2)
#endif
        do i=1,nelecs
           do j=1,nelecs
              temp(i,j)=w(i,j)
           enddo
        enddo
#ifdef OMP
!$omp end do
!$omp do collapse(2)
#endif
        do i=1,nelecs
           do j=1,nelecs
              w(j,i)=temp(i,j)
           enddo
        enddo
#ifdef OMP
!$omp end do
!$omp end parallel
#endif
     endif
  endif

  return
end subroutine gronor_svd
