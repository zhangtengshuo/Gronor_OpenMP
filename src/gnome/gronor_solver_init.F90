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
#ifdef ROCSOLVER
  use rocvars
  use rocsolver_interfaces_enums
  use rocsolver_interfaces
!  use hipfort
!  use hipfort_check
!  use hipfort_rocblas_enums
!  use hipfort_rocblas
!  use hipfort_rocsolver_enums
!  use hipfort_rocsolver
#endif

#ifdef HIPSOLVER
  use hipfort
  use hipfort_check
  use hipfort_rocblas_enums
  use hipfort_rocblas
  use hipfort_rocsolver_enums
  use hipfort_rocsolver
#endif
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
  external :: magmaf_dsyevd,magmaf_dsyevd_gpu
  integer (kind=4) :: magma_info,ierr
#endif
#endif
#ifdef CUSOLVER
  external :: cusolverdncreategesvdjinfo,cusolverdnxgesvdjsettolerance
  external :: cusolverdnxgesvdjsetmaxsweeps,cusolverdndgesvdj_buffersize
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
  lsvtrns=.false.
  
  if(sv_solver.eq.SOLVER_EISPACK) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_LAPACK) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_LAPACKD) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_LAPACKQ) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_LAPACKJ) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_LAPACKJH) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_MKL) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_MKLD) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_MKLJ) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_CRAYLIBSCID_CPU) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_MAGMA) lsvcpu=.true.
  if(sv_solver.eq.SOLVER_MAGMAD) lsvcpu=.true.

  if(sv_solver.eq.SOLVER_LAPACK) lsvtrns=.true.
  if(sv_solver.eq.SOLVER_LAPACKD) lsvtrns=.true.
  if(sv_solver.eq.SOLVER_LAPACKQ) lsvtrns=.true.
  if(sv_solver.eq.SOLVER_LAPACKJ) lsvtrns=.true.
  if(sv_solver.eq.SOLVER_LAPACKJH) lsvtrns=.true.
  if(sv_solver.eq.SOLVER_MKL) lsvtrns=.true.
  if(sv_solver.eq.SOLVER_MKLD) lsvtrns=.true.
  if(sv_solver.eq.SOLVER_MKLJ) lsvtrns=.true.
  if(sv_solver.eq.SOLVER_CRAYLIBSCID_CPU) lsvtrns=.true.
  if(sv_solver.eq.SOLVER_CRAYLIBSCID_ACC) lsvtrns=.true.
  
  if(ev_solver.eq.SOLVER_EISPACK) levcpu=.true.
  if(ev_solver.eq.SOLVER_LAPACK) levcpu=.true.
  if(ev_solver.eq.SOLVER_LAPACKD) levcpu=.true.
  if(ev_solver.eq.SOLVER_LAPACKQ) levcpu=.true.
  if(ev_solver.eq.SOLVER_LAPACKJ) levcpu=.true.
  if(ev_solver.eq.SOLVER_LAPACKJH) levcpu=.true.
  if(ev_solver.eq.SOLVER_MKL) levcpu=.true.
  if(ev_solver.eq.SOLVER_MKLD) levcpu=.true.
  if(ev_solver.eq.SOLVER_MKLJ) levcpu=.true.
  if(ev_solver.eq.SOLVER_CRAYLIBSCID_CPU) levcpu=.true.
  
  if(iamacc.ne.0) then
#ifdef CUSOLVER

    ndim=nelecs
    mdim=mbasel
    lwork1=0
    lwork2=0

    if(sv_solver.eq.SOLVER_CUSOLVER) then

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
    len_work_dbl=max(len_work_dbl,lwork1)

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
    len_work_dbl=max(len_work_dbl,lwork1)
    endif

! Cusolver initialization for the syevd

    if(ev_solver.eq.SOLVER_CUSOLVER) then

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

    call gronor_update_device_info()

    if(memavail.gt.0.and.8*lwork1.gt.memavail) then
      write(string,'(a,i10,a,i10)') "Available ",memavail," device memory insufficient for", &
          8*lwork1," needed as workspace for CUSOLVER solvers"
      call gronor_abort(500,string)
    endif

    len_work_dbl=max(len_work_dbl,lwork1,lwork2)

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
      lworki=max(8*nelecs,lworki)
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
      lworki=max(int(iworksize(1)),lworki)
    endif
    lwork1m=max(0,lwork1m,lwork2m)
    lworki=max(0,lworki)
    len_work_dbl=max(len_work_dbl,lwork1m)
    len_work_int=max(len_work_int,lworki)
#endif

#ifdef LAPACK
    ndimm=nelecs
    mdimm=mbasel
    lwork1m=-1
    lwork2m=-1
    lworki=-1
    
    if(sv_solver.eq.SOLVER_LAPACK) then
      call dgesvd('A','A',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,worksize,lwork1m,lapack_info)
      lwork1m=int(worksize(1))
    endif
    if(sv_solver.eq.SOLVER_LAPACKD) then
      call dgesdd('A',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,worksize,lwork1m,iworksize, &
          lapack_info)
      lwork1m=int(worksize(1))
      lworki=8*nelecs
    endif
    if(ev_solver.eq.SOLVER_LAPACK) then
      call dsyev('V','L',ndimm,a,ndimm,w,worksize,lwork2m,lapack_info)
      lwork2m=int(worksize(1))
    endif
    if(ev_solver.eq.SOLVER_LAPACKD) then
      call dsyevd('V','L',ndimm,a,ndimm,w,worksize,lwork2m,iworksize,lworki,lapack_info)
      lwork2m=int(worksize(1))
      lworki=max(int(iworksize(1)),lworki)
    endif
    lwork1m=max(0,lwork1m,lwork2m)
    lworki=max(0,lworki)
    len_work_dbl=max(len_work_dbl,lwork1m)
    len_work_int=max(len_work_int,lworki)
#endif

#ifdef ROCSOLVER
    lwork1m=max(0,5*nelecs*nelecs)
    lworki=max(0,lworki)
    len_work_dbl=max(len_work_dbl,lwork1m)
    len_work_int=max(len_work_int,lworki)
!!!!! !$omp target enter data map(rocinfo)
#endif

#ifdef CRAYLIBSCI
    ndimm=nelecs
    mdimm=mbasel
    lwork1m=-1
    lwork2m=-1
    lworki=-1
    
    if(sv_solver.eq.SOLVER_CRAYLIBSCID_CPU) then
!      call dgesvd('A','A',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,worksize,lwork1m,lapack_info)
!      lwork1m=int(worksize(1))
    endif
    if(sv_solver.eq.SOLVER_CRAYLIBSCID_ACC) then
!      call dgesdd('A',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,worksize,lwork1m,iworksize, &
!          lapack_info)
!      lwork1m=int(worksize(1))
!      lworki=8*nelecs
    endif
    if(ev_solver.eq.SOLVER_CRAYLIBSCID_CPU) then
      call dsyevd_cpu('V','L',ndimm,a,ndimm,w,worksize,lwork2m,info)
      lwork2m=int(worksize(1))
    endif
    if(ev_solver.eq.SOLVER_CRAYLIBSCID_ACC) then
!$omp target enter data map(to:a) map(alloc:w,worksize,iworksize,info)
!$omp target data use_device_addr(a,w,worksize,iworksize,info)
      call dsyevd_acc('V','L',ndimm,a,ndimm,w,worksize,lwork2m,iworksize,lworki,info)
!$omp end target data
      lwork2m=int(worksize(1))
      lworki=max(int(iworksize(1)),lworki)
    endif
    lwork1m=max(0,lwork1m,lwork2m)
    lworki=max(0,lworki)
    len_work_dbl=max(len_work_dbl,lwork1m)
    len_work_int=max(len_work_int,lworki)
!$omp target exit data map(delete:worksize,iworksize)
#endif

#ifdef MAGMA
    ndimm=nelecs
    mdimm=mbasel
    ndim=nelecs
    lwork1m=-1
    lwork2m=-1
    lworki=-1
    ndim4=nelecs
    lwork4=-1
    liwork4=-1
    ndimm=nelecs
    ndim4=nelecs
    lwork4=-1
    liwork4=-1
    
    
    if(sv_solver.eq.SOLVER_MAGMA) then
      if(iamacc.eq.1) then
        call magmaf_dgesvd('A','A',ndim4,ndim4,a,ndim4,ev,u,ndim4,w,ndim4, &
            worksize,lwork4,magma_info)
      else
        call magmaf_dgesvd('A','A',ndim4,ndim4,a,ndim4,ev,u,ndim4,w,ndim4, &
            worksize,lwork4,magma_info)
      endif
      lwork1m=int(worksize(1))
      lworki=8*nelecs
    endif
      
    if(sv_solver.eq.SOLVER_MAGMAD) then
      if(iamacc.eq.1) then
        call magmaf_dgesdd('A',ndim4,ndim4,a,ndim4,ev,u,ndim4,w,ndim4, &
            worksize,lwork4,iworksize,magma_info)
      else
        call magmaf_dgesdd('A',ndim4,ndim4,a,ndim4,ev,u,ndim4,w,ndim4, &
            worksize,lwork4,iworksize,magma_info)
      endif
      lwork1m=int(worksize(1))
      lworki=8*nelecs
    endif
    
    if(ev_solver.eq.SOLVER_MAGMA) then
      if(iamacc.eq.1) then
#ifdef ACC
!$acc data create(workspace_d,workspace_i,workspace2_d)
!$acc wait
!!!!$acc host_data use_device(a,diag,workspace_d,workspace_i4,workspace2_d)
#endif
#ifdef OMPTGT
!$omp target data use_device_addr(a,ev,u,w,dev_info_d,workspace_d,workspace_i4,workspace2_d)
#endif     
        ndimm=nelecs
        ndim4=nelecs
        lwork4=-1
        liwork4=-1
        call magmaf_dsyevd_gpu('N','L',ndim4,c_loc(a),ndim4,diag,worksize2,ndim4, &
            worksize,lwork4,iworksize,liwork4,magma_info)
#ifdef ACC
!!!!$acc end host_data
!$acc wait
!$acc end data   
#endif
#ifdef OMPTGT
!$omp end target data
#endif
      else       
        ndimm=nelecs
        ndim4=nelecs
        lwork4=-1
        liwork4=-1
        call magmaf_dsyevd('N','L',ndim4,a,ndim4,diag, &
            worksize,lwork4,iworksize,liwork4,magma_info)
      endif
      lwork2m=int(worksize(1))
      lworki=int(iworksize(1))
    endif
    lwork1m=max(0,lwork1m,lwork2m)
    lworki=max(0,lworki)
    len_work_dbl=max(len_work_dbl,lwork1m)
    len_work_int=max(len_work_int,lworki)
    len_work2_dbl=nelecs*nelecs
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

subroutine gronor_solver_final()

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

#ifdef ROCSOLVER
  use rocvars
  use rocsolver_interfaces_enums
  use rocsolver_interfaces
!  use hipfort
!  use hipfort_check
!  use hipfort_rocblas_enums
!  use hipfort_rocblas
!  use hipfort_rocsolver_enums
!  use hipfort_rocsolver
#endif

#ifdef HIPSOLVER
  use hipfort
  use hipfort_check
  use hipfort_rocblas_enums
  use hipfort_rocblas
  use hipfort_rocsolver_enums
  use hipfort_rocsolver
#endif

  
  if(iamacc.gt.0) then
#ifdef CUSOLVER
    cusolver_status = cusolverDnDestroy(cusolver_handle)
    if (cusolver_status /= CUSOLVER_STATUS_SUCCESS) &
        write(*,*) 'cusolver_handle destruction failed'
#endif
#ifdef HIPSOLVER
    if (hipsolverDestroy(hipsolver_handle) /= HIPSOLVER_STATUS_SUCCESS) &
        write(*,*) 'hipsolver_handle destruction failed'
#endif
#ifdef ROCSOLVER
    if (rocsolver_destroy_handle(rocsolver_handle) /= 0) &
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

#ifdef ROCSOLVER
  use rocvars
  use rocsolver_interfaces_enums
  use rocsolver_interfaces
!  use hipfort
!  use hipfort_check
!  use hipfort_rocblas_enums
!  use hipfort_rocblas
!  use hipfort_rocsolver_enums
!  use hipfort_rocsolver
#endif

#ifdef HIPSOLVER
  use hipfort
  use hipfort_check
  use hipfort_rocblas_enums
  use hipfort_rocblas
  use hipfort_rocsolver_enums
  use hipfort_rocsolver
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
!      istat=cudaMemGetInfo(cpfre,cptot)
    endif
#endif

#ifdef ROCSOLVER
    istatus=rocsolver_create_handle(rocsolver_handle)
!    call hipcheck(rocblas_create_handle(rocsolver_handle))
#endif

#ifdef HIPSOLVER
    call hipsolvercreate(hipsolver_handle)
#endif
    
  endif
        
  return
end subroutine gronor_solver_create_handle
