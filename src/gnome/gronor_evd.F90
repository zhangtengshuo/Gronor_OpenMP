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
  !! SOLVER_LAPACK       dsyevd from an external LAPACK library will run on the CPU
  !! SOLVER_CUSOLVER     cusolverDnDsyevd from the NVIDIA CUSOLVER library wil run on NVDIA GPUs
  !! SOLVER_CUSOLVERJ    cusolverDnDsyesvj from the NVIDIA CUSOLVER library wil run on NVDIA GPUs
  !! SOLVER_ROCSOLVER    rocsolver_dsyevd from the AMD ROCSOLVER library will run on the GPU
  !! SOLVER_ROCSOLVERD   same as SOLVER ROCSOLVER
  !! SOLVER_ROCSOLVERX   same as SOLVER_ROCSOLVER
  !! SOLVER_HIPSOLVER    planned
  !! SOLVER_HIPSOLVER    planned
  !! SOLVER_MAGMA        planned
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

#ifdef LAPACK  
  use lapack_solver
#else
#ifdef MAGMA
  use magma
  use magma_dfortran
  use magma_solver
#endif
#endif
  
#ifdef CUSOLVER
  use cusolverDn
  use cuda_cusolver
  use cudafor
#endif

#ifdef ROCSOLVER
  use iso_fortran_env
  use rocvars
  use rocsolver_interfaces_enums
  use rocsolver_interfaces
#endif

#ifdef HIPSOLVER
  use hipvars
  use iso_fortran_env
  use hipfort
  use hipfort_check
  use hipfort_rocblas_enums
  use hipfort_rocblas
  use hipfort_rocsolver_enums
  use hipfort_rocsolver
#endif

  ! variable declarations

  implicit none

  external :: tred2,tql2
#ifdef MKL
  external :: dsyevd
#endif
#ifdef CUSOLVER
  external cusolverdndsyevj
#endif
#ifdef MAGMA
  integer (kind=4) :: magma_info
#endif
  
  integer :: i,j
  integer :: ierr

  ! library specific declarations

#ifdef CUSOLVER
  character (len=1), target :: jobu, jobvt
#endif

#ifdef ROCSOLVER
  real(kind=8), allocatable, target :: at(:,:),dt(:),et(:,:)
#endif

  if(iamacc.eq.1) then
     if(levcpu) then
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
  
! ============ LAPACK ===========
#ifdef LAPACK
  if(ev_solver.eq.SOLVER_LAPACK.or.ev_solver.eq.SOLVER_LAPACKD) then
    ndimm=nelecs
    if(ev_solver.eq.SOLVER_LAPACK) then
      call dsyev('N','L',ndimm,a,nelecs,diag,workspace_d,len_work_dbl,ierr)
    elseif(ev_solver.eq.SOLVER_LAPACKD) then
      call dsyevd('N','L',ndimm,a,nelecs,diag,workspace_d,len_work_dbl, &
          workspace_i,len_work_int,ierr)
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
!$acc host_data use_device(a,diag,workspace_d,workspace_i4,workspace2_d)
#endif
#ifdef OMPTGT
!$omp target data use_device_addr(a,ev,u,w,dev_info_d,workspace_d,workspace_i4,workspace2_d)
#endif    
      ndimm=nelecs
      ndim4=nelecs
      lwork4=len_work_dbl
      liwork4=len_work_int
      call magmaf_dsyevd_gpu('N','L',ndim4,c_loc(a),ndim4,diag,workspace2_d,ndim4, &
          workspace_d,lwork4,workspace_i4,liwork4,magma_info)
#ifdef ACC
!$acc end host_data
!$acc wait
!$acc end data   
#endif
#ifdef OMPTGT
!$omp end target data
#endif
    else        
      ndimm=nelecs
      ndim4=nelecs
      lwork4=len_work_dbl
      liwork4=len_work_int
      call magmaf_dsyevd('N','L',ndim4,a,ndim4,diag, &
      workspace_d,lwork4,workspace_i4,liwork4,magma_info)
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
        ndim,a,ndim,diag,workspace_d,int(len_work_dbl,kind=4),dev_info_d)
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

  if(ev_solver.eq.SOLVER_CUSOLVERJ) then
    ndim=nelecs
    mdim=mbasel
    jobz = CUSOLVER_EIG_MODE_NOVECTOR
    uplo = CUBLAS_FILL_MODE_LOWER

#ifdef ACC
!$acc data copy(dev_info_d,syevj_params) create(workspace_d)
!$acc host_data use_device(a,diag,dev_info_d,workspace_d,syevj_params)
#endif
#ifdef OMPTGT
!$omp target data use_device_addr(a,ev,u,w,dev_info_d,workspace_d,syevj_params)
#endif
    cusolver_status = cusolverDnDsyevj &
        (cusolver_handle, jobz, uplo, ndim,a,ndim,diag, &
        workspace_d,int(len_work_dbl,kind=4),dev_info_d,syevj_params)
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
!$omp target data use_device_addr(a,diag,workspace_d,rocinfo)
    rocsolver_status=rocsolver_dsyevd(rocsolver_handle,ROCBLAS_EVECT_NONE,ROCBLAS_FILL_LOWER, &
        ndim,c_loc(a),ndim,c_loc(diag),c_loc(workspace_d),rocinfo)
!$omp end target data
!     rocsolver_status=hipDeviceSynchronize()
  endif

  if(ev_solver.eq.SOLVER_ROCSOLVERD) then
    ndim=nelecs
    mdim=mbasel
    rocsolver_status=rocsolver_dsyevd(rocsolver_handle,ROCBLAS_EVECT_NONE,ROCBLAS_FILL_LOWER, &
        ndim,c_loc(a),ndim,c_loc(diag),c_loc(workspace_d),rocinfo)
!    rocsolver_status=hipDeviceSynchronize()
  endif

  if(ev_solver.eq.SOLVER_ROCSOLVERX) then
    ndim=nelecs
    mdim=mbasel
    rocsolver_status=rocsolver_dsyevd(rocsolver_handle,ROCBLAS_EVECT_NONE,ROCBLAS_FILL_LOWER, &
        ndim,c_loc(a),ndim,c_loc(diag),c_loc(workspace_d),rocinfo)
!    rocsolverstatus=hipDeviceSynchronize()
  endif
#endif

! ======== CRAYLIBSCI =========

#ifdef CRAYLIBSCI
  if(ev_solver.eq.SOLVER_CRAYLIBSCID_CPU) then
    ndim=nelecs
    mdim=mbasel
    call dsyevd_acc(jobz, uplo, ndim, a, ndim, diag, workspace_d, lwork, workspace_i, liwork, info)
  endif
  if(ev_solver.eq.SOLVER_CRAYLIBSCID_ACC) then
    ndim=nelecs
    mdim=mbasel
!$omp target enter data map(to:a)
!$omp target data use_device_addr(a,diag,workspace_d,workspace_i,info)
    call dsyevd_acc(jobz, uplo, ndim, a, ndim, diag, workspace_d, lwork, workspace_i, liwork, info)
!$omp end target data
  endif
#endif

  if(iamacc.eq.1.and.levcpu) then
#ifdef ACC
!$acc update device (diag)
#endif
#ifdef OMPTGT
!$omp target update to(diag)
#endif
  endif

  return
end subroutine gronor_evd
