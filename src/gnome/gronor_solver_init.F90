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
#ifdef CUSOLVER
  use cusolverDn
  use cuda_cusolver
  use cudafor
#endif
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
  
  implicit none
  
#ifdef MKL
  external :: dgesvd,dsyevd
  integer (kind=4) :: ierr
#endif
#ifdef LAPACK
  external :: dgesvd,dsyevd
  integer (kind=4) :: lapack_info,ierr
#else  
#ifdef MAGMA
  external :: magma_dsyevd,magma_dsyevd_gpu
  integer (kind=4) :: lapack_info,ierr
#endif
#endif  

  integer :: ntemp
  character(len=255) :: string
  
  real(kind=8) :: worksize(2)
  integer (kind=8) :: iworksize(2)

  nelecs=ntemp

! Cusolver initialization for the svd
  
  if(idbg.gt.50) then
    call swatch(date,time)
    write(lfndbg,'(a,1x,a,a,2i4)') date(1:8),time(1:8)," Solver init for ",sv_solver,ev_solver
    flush(lfndbg)
  endif
  
  if(iamacc.ne.0) then
#ifdef CUSOLVER

    ndim=nelecs
    mdim=mbasel
    lwork1=0
    lwork2=0

#ifdef CUSOLVERJ
    if(sv_solver.eq.SOLVER_CUSOLVER) then
#endif
#ifdef ACC
!$acc data copyin(w,ta) create(dev_info_d)
!$acc host_data use_device(ta)
#endif
      cusolver_status = cusolverDnDgesvd_bufferSize(cusolver_handle,ndim,ndim,lwork1)
#ifdef ACC
!$acc end host_data
#endif
      if (cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
          write(*,*) 'cusolverDnDgesvd_bufferSize failed'

#ifdef ACC
!$acc end data
#endif
#ifdef CUSOLVERJ
    elseif(sv_solver.eq.SOLVER_CUSOLVERJ) then

      ndim=nelecs
      mdim=mbasel
      tol = tolsvj
      max_sweeps = iswsvj

      jobz=CUSOLVER_EIG_MODE_VECTOR

#ifdef ACC
!$acc data create(dev_info_d,u,w,ev,ta)
#endif
      cusolver_status = cudaStreamCreateWithFlags(stream,cudaStreamNonBlocking)

      cusolver_status = cusolverDnSetStream(cusolver_handle,stream)

      cusolver_status = cusolverDnCreateGesvdjInfo(gesvdj_params)

      cusolver_status = cusolverDnXgesvdjSetTolerance(gesvdj_params,tol)

      cusolver_status = cusolverDnXgesvdjSetMaxSweeps(gesvdj_params,max_sweeps)

#ifdef ACC
!$acc host_data use_device(ta,ev,u,w)
#endif
      cusolver_status = cusolverDnDgesvdj_bufferSize(cusolver_handle,jobz,econ, &
          ndim,ndim,ta,mdim,ev,u,ndim,w,ndim,lwork1,gesvdj_params)

#ifdef ACC
!$acc end host_data
#endif
      if (cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
          print *,"cusolverDnDgesvdj_bufferSize failed",cusolver_status

#ifdef ACC
!$acc end data
#endif
    endif
#endif

! Cusolver initialization for the syevd

#ifdef CUSOLVERJ
    if(ev_solver.eq.SOLVER_CUSOLVER) then
#endif
#ifdef ACC
!$acc data copyin(w,ta) create(dev_info_d)
!$acc host_data use_device(ta,w)
#endif
      cusolver_status = cusolverDnDsyevd_bufferSize(cusolver_handle,CUSOLVER_EIG_MODE_NOVECTOR, &
          CUBLAS_FILL_MODE_LOWER,ndim,ta,mdim,w,lwork2)
#ifdef ACC
!$acc end host_data
#endif
      if (cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
                   print *,"cusolverDnDsyevd_bufferSize failed",cusolver_status

#ifdef ACC
!$acc end data
#endif
#ifdef CUSOLVERJ
    elseif(ev_solver.eq.SOLVER_CUSOLVERJ) then

      ! Jacobi EVD

      ndim=nelecs
      mdim=mbasel
      tol = tolevj
      max_sweeps = iswevj

      jobz = CUSOLVER_EIG_MODE_NOVECTOR
      uplo = CUBLAS_FILL_MODE_LOWER

#ifdef ACC
!$acc data create(dev_info_d,w,ta,ev)
#endif
      cusolver_status = cudaStreamCreateWithFlags(stream,cudaStreamNonBlocking)

      cusolver_status = cusolverDnSetStream(cusolver_handle,stream)

      cusolver_status = cusolverDnCreateSyevjInfo(syevj_params)

      cusolver_status = cusolverDnXsyevjSetTolerance(syevj_params,tol)

      cusolver_status = cusolverDnXsyevjSetMaxSweeps(syevj_params,max_sweeps)

#ifdef ACC
!$acc host_data use_device(ta,w)
#endif
      cusolver_status = cusolverDnDsyevj_bufferSize(cusolver_handle,jobz,uplo, &
          ndim,ta,mdim,ev,lwork2,syevj_params)

#ifdef ACC
!$acc end host_data
#endif
      if (cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
          print *,"cusolverDnDsyevj_bufferSize failed",cusolver_status

#ifdef ACC
!$acc end data
#endif
    endif
#endif

    lwork1=max(lwork1,lwork2)

    call gronor_update_device_info()

    if(memavail.gt.0.and.8*lwork1.gt.memavail) then
      write(string,'(a,i10,a,i10)') "Available ",memavail," device memory insufficient for", &
          8*lwork1," needed as workspace for CUSOLVER solvers"
      call gronor_abort(500,string)
    endif
    
    allocate(workspace_d(lwork1))

#endif
  endif

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
      lworki=8*ndimm
    endif
    if(ev_solver.eq.SOLVER_MKL) then
      call dsyevd('N','L',ndimm,a,ndimm,w,worksize,lwork2m,iworksize,lworki,ierr)
      lwork2m=int(worksize(1))
      lworki=int(iworksize(1))
    endif
    if(ev_solver.eq.SOLVER_MKLD) then
      call dsyevd('N','L',ndimm,a,ndimm,w,worksize,lwork2m,iworksize,lworki,ierr)
      lwork2m=int(worksize(1))
      lworki=int(iworksize(1))
    endif
    lwork1m=max(8*ndimm,lwork1m,lwork2m)
    lworki=max(8*ndimm,lworki)
    allocate(workspace_d(lwork1m))
    allocate(workspace_i(lworki))
#endif

#ifdef LAPACK
    ndimm=nelecs
    mdimm=mbasel
    lwork1m=-1
    lwork2m=-1
    lworki=-1
  
    if(sv_solver.eq.SOLVER_LAPACK) then
      call dgesvd('A','A',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,worksize,lwork1m,lapack_info)
      lwork1m=max(int(worksize(1)),1+3*nelecs,5*nelecs)+1024*nelecs
    endif
    if(sv_solver.eq.SOLVER_LAPACKD) then
      call dgesdd('A',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,worksize,lwork1m,iworksize, &
          lapack_info)
      lwork1m=max(int(worksize(1)),7*nelecs+4*nelecs*nelecs)+1024*nelecs
      lworki=8*nelecs
    endif
    if(ev_solver.eq.SOLVER_LAPACK) then
      call dsyev('V','L',ndimm,a,ndimm,w,worksize,lwork2m,lapack_info)
      lwork2m=int(worksize(1))+1024*nelecs
    endif
    if(ev_solver.eq.SOLVER_LAPACKD) then
      call dsyevd('V','L',ndimm,a,ndimm,w,worksize,lwork2m,iworksize,lworki,lapack_info)
      lwork2m=max(int(worksize(1)),1+6*nelecs+2*nelecs*nelecs)
      lworki=max(int(iworksize(1)),3+5*nelecs)
    endif
    lwork1m=max(8,lwork1m,lwork2m)
    lworki=max(8,lworki)
    allocate(workspace_d(lwork1m))
    allocate(workspace_i(lworki))
#endif  

#ifdef MAGMA
    ndimm=nelecs
    mdimm=mbasel
    lwork1m=-1
    lwork2m=-1
    lworki=-1
!    if(sv_solver.eq.SOLVER_LAPACK) then
!      call dgesvd('A','A',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,worksize,lwork1m,lapack_info)
!      lwork1m=int(worksize(1))+1024*nelecs
!    endif
    if(ev_solver.eq.SOLVER_MAGMA) then
      if(iamacc.eq.0) then
        call magma_dsyevd('V','L',ndimm,a,ndimm,w,worksize,lwork2m,iworksize,lworki,lapack_info)
      else
        call magma_dsyevd_gpu('V','L',ndimm,a,ndimm,w,worksize,lwork2m,iworksize,lworki,lapack_info)
      endif
      lwork2m=int(worksize(1))+1024*nelecs
      lworki=int(iworksize(1))+1024*nelecs
    endif
    lwork1m=max(8,lwork1m,lwork2m)
    lworki=max(8,lworki)
    allocate(workspace_d(lwork1m))
    allocate(workspace_i(lworki))
#endif
    return
  end subroutine gronor_solver_init

subroutine gronor_solver_final()

  if(iamacc.gt.0) then
#ifdef CUSOLVER
    cusolver_status = cusolverDnDestroy(cusolver_handle)
    if (cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
        write(*,*) 'cusolver_handle destruction failed'
#endif
#ifdef HIPSOLVER
    hipsolver_status = hipsolverDestroy(hipsolver_handle)
    if (hipsolver_status /= HIPSOLVER_STATUS_SUCCESS) &
        write(*,*) 'hipsolver_handle destruction failed'
#endif
  endif

  return
end subroutine gronor_solver_final

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
  
#ifdef CUSOLVER
  use cusolverDn
  use cuda_cusolver
#endif
  
  ! Only accelerated ranks need to define cusolver handles
  
  if(iamacc.gt.0) then
    
#ifdef CUSOLVER
    cusolver_status=cusolverDnCreate(cusolver_handle)
    if(cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
        print *, me,'cusolver_handle creation failed'
    if(numdev.gt.1) then
!      cpfre=c_loc(memfre)
!      cptot=c_loc(memtot)
      ! istat=cudaMemGetInfo(cpfre,cptot)
    endif
#endif
    
  endif
        
  return
end subroutine gronor_solver_create_handle
