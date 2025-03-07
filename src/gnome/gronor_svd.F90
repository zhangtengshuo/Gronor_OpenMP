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
  !! SOLVER_LAPACK       dgesvd from an external LAPACK library will run on the CPU
  !! SOLVER_CUSOLVER     cusolverDnDgesvd from the NVIDIA CUSOLVER library wil run on NVDIA GPUs
  !! SOLVER_CUSOLVERJ    cusolverDnDgesvdj from the NVIDIA CUSOLVER library wil run on NVDIA GPUs
  !! SOLVER_ROCSOLVER    rocsolver_dgesvd from the AMD ROCSOLVER library will run on the GPU
  !! SOLVER_ROCSOLVERD   same as SOLVER ROCSOLVER
  !! SOLVER_ROCSOLVERX   same as SOLVER_ROCSOLVER
  !! SOLVER_HIPSOLVER    planned
  !! SOLVER_HIPSOLVER    planned
  !! SOLVER_MAGMA        planned
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
#ifdef HIPSOLVER
  use hipvars
  use hipsolver_enums
  use amd_hipsolver
#endif
#ifdef ROCSOLVER
  use iso_fortran_env
  use rocvars
  use rocsolver_interfaces_enums
  use rocsolver_interfaces
#endif

  ! variable declarations

  implicit none

  external :: svd,dgesvd
  
  integer :: i,j
  integer :: ierr
  integer (kind=4) :: istat

  ! library specific declarations

#ifdef CUSOLVER
  character (len=1), target :: jobu, jobvt
  external :: cusolverdndgesvdj,cusolverdnxgesvdjgetsweeps,cusolverdnxgesvdjgetresidual
#endif
#ifdef MAGMA
  integer (kind=4) :: magma_info
#endif

#ifdef HIPSOLVER
  integer (kind=1) :: jobu, jobvt
#endif

#ifdef ROCSOLVER
  integer :: jobu, jobvt
#endif

#ifdef LAPACK
  integer (kind=4) :: lapack_info
#endif

  lwork4=int(len_work_dbl,kind=4)
  
  if(iamacc.eq.1.and.lsvcpu) then
#ifdef ACC
!$acc update host (a)
#endif
#ifdef OMPTGT
!$omp target update from(a)
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
  
  ! ============ LAPACK ===========
#ifdef LAPACK
  if(sv_solver.eq.SOLVER_LAPACK.or.sv_solver.eq.SOLVER_LAPACKD) then
    ndimm=nelecs
    ndim=nelecs
    if(sv_solver.eq.SOLVER_LAPACK) then
      call dgesvd('A','A',ndim,ndim,a,ndim,ev,u,ndim,wt,ndim, &
          workspace_d,lwork4,ierr)
    endif
    if(sv_solver.eq.SOLVER_LAPACKD) then
      call dgesdd('All',ndim,ndim,a,ndim,ev,u,ndim,wt,ndim, &
          workspace_d,lwork4,workspace_i,ierr)
    endif
    !lsvtrns
  endif
#endif 
  
! ============ MAGMA ===========
  
#ifdef MAGMA
  
  ndimm=nelecs
  mdimm=mbasel
  ndim=nelecs
  ndim4=nelecs

  if(sv_solver.eq.SOLVER_MAGMA) then
    if(iamacc.eq.1) then
      call magmaf_dgesdd('A',ndim4,ndim4,a,ndim4,ev,u,ndim4,w,ndim4, &
          workspace_d,lwork4,workspace_i4,magma_info) 
    else
      call magmaf_dgesdd('A',ndim4,ndim4,a,ndim4,ev,u,ndim4,w,ndim4, &
          workspace_d,lwork4,workspace_i4,magma_info) 
    endif
  endif
  if(sv_solver.eq.SOLVER_MAGMAD) then
    if(iamacc.eq.1) then
      call magmaf_dgesvd('A','A',ndim4,ndim4,a,ndim4,ev,u,ndim4,w,ndim4, &
          workspace_d,lwork4,magma_info)
    else
      call magmaf_dgesvd('A','A',ndim4,ndim4,a,ndim4,ev,u,ndim4,w,ndim4, &
          workspace_d,lwork4,magma_info)
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
!$acc host_data use_device(a,ev,u,wt,dev_info_d,workspace_d,rwork)
#endif
#ifdef OMPTGT
!$omp target data use_device_addr(a,ev,u,wt,dev_info_d,workspace_d,rwork)
#endif
    cusolver_status=cusolverDnDgesvd(cusolver_handle,jobu,jobvt, &
        ndim,ndim,a,ndim,ev,u,ndim,wt,ndim,workspace_d, &
        lwork4,rwork,dev_info_d)
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
    !lsvtrns
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
        lwork4,dev_info_d,gesvdj_params)
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
!    istatus=rocsolver_dgesvd(rocsolver_handle, &
!        ROCBLAS_SVECT_ALL,ROCBLAS_SVECT_ALL,ndim,ndim,c_loc(a),ndim, &
!        c_loc(ev),c_loc(u),ndim,c_loc(wt),ndim,c_loc(work), &
!        ROCBLAS_OUTOFPLACE,rocinfo)
!    istatus=hipDeviceSynchronize()
    
!$omp target data use_device_addr(a,ev,u,wt,workspace_d,rocinfo)
    istat=rocsolver_dgesvd(rocsolver_handle, &
        ROCBLAS_SVECT_ALL,ROCBLAS_SVECT_ALL,ndim,ndim,c_loc(a),ndim, &
        c_loc(ev),c_loc(u),ndim,c_loc(wt),ndim,c_loc(workspace_d), &
        ROCBLAS_INPLACE,rocinfo)
!$omp end target data
!    call hipcheck(hipDeviceSynchronize())
  endif
  if(sv_solver.eq.SOLVER_ROCSOLVERX) then
    ndim=nelecs
    mdim=mbasel
!    istatus=rocsolver_dgesvd(rocsolver_handle, &
!        ROCBLAS_SVECT_ALL,ROCBLAS_SVECT_ALL,ndim,ndim,c_loc(a),ndim, &
!        c_loc(ev),c_loc(u),ndim,c_loc(wt),ndim,c_loc(work), &
!        ROCBLAS_OUTOFPLACE,rocinfo)
!    istatus=hipDeviceSynchronize()
    istat=rocsolver_dgesvd(rocsolver_handle, &
        ROCBLAS_SVECT_ALL,ROCBLAS_SVECT_ALL,ndim,ndim,c_loc(a),ndim, &
        c_loc(ev),c_loc(u),ndim,c_loc(wt),ndim,c_loc(work), &
        ROCBLAS_OUTOFPLACE,rocinfo)
!    call hipcheck(hipDeviceSynchronize())
  endif
  !lsvtrns
#endif

!! Evaluate right had matrix from its transpose if solver provided the transpose
  if(lsvtrns) then
    if(lsvcpu) then

#ifdef OMP
!$omp parallel shared(w,wt,nelecs)
!$omp do collapse(2)
#endif
      do i=1,nelecs
        do j=1,nelecs
          w(i,j)=wt(j,i)
        enddo
      enddo
#ifdef OMP
!$omp end do
!$omp end parallel
#endif
    elseif(iamacc.eq.1) then

!======================================
      
!#ifdef ACC
!!$acc update host(wt)
!#endif
!#ifdef OMPTGT
!!$omp target update from(wt) 
!#endif
!   
!#ifdef OMP
!!$omp parallel shared(w,wt,nelecs)
!!$omp do collapse(2)
!#endif
!      do i=1,nelecs
!        do j=1,nelecs
!          w(i,j)=wt(j,i)
!        enddo
!      enddo
!#ifdef OMP
!!$omp end do
!!$omp end parallel
!#endif
!      
!#ifdef ACC
!!$acc update device(w)
!#endif
!#ifdef OMPTGT
!!$omp target update to(w) 
!#endif  

!======================================
      
#ifdef ACC
!$acc kernels present(w,wt)
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
      w(i,j)=wt(j,i)
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

!=========================================
      
    endif
  endif

!! Update device if SVD was performed on the host for accelerated ranks
  if(iamacc.eq.1.and.lsvcpu) then
#ifdef ACC
!$acc update device (ev,u,w)
#endif
#ifdef OMPTGT
!$omp target update to(ev,u,w) 
#endif  
  endif

  return  
end subroutine gronor_svd
