      module const
      integer :: nr(20,3)
      data nr/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1,
     &        0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1,
     &        0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1/
      real (kind=8) :: pi,twopi,dnorm,piroot,piterm,enorm
      data pi,twopi,dnorm/3.1415926535898d0,6.2831853071796d0,32.0d0/
      data piroot,piterm,enorm/0.88622692545276d0,34.986836655250d0,
     & 1.0d-12/
      end module const

      module cidist

      implicit none

!     MPI distribution

!     me : global process id
!     meg : group process id
!     np : global number of processes
!     ng : number of process groups
!     npg : number of processes per group
!     comm_group : communicator of group current process belongs to
!     comm_first : communicator of all first processes in each group
!     mpibuf : size of the MPI buffer for integrals in reals

      integer (kind=4) :: me, np, master
      integer (kind=4) :: numdev, mydev
      integer (kind=4) :: group_batch, group_heads, group_world
      integer (kind=4), allocatable :: grp_comm(:)
      integer (kind=4), allocatable :: ranks_list(:), ranks_heads(:)

      integer :: nnodes,nrsets,nranks,ncycls,nrnsets

      integer :: numacc, numnon, numgrp
      integer :: ngr,meg,npg,comm_group,comm_first,mgr,nalive
      integer :: mygroup,myhead,ime,mpibuf
      integer :: iamhead, iamacc, iamactive
      integer (kind=8), allocatable  :: map1(:,:),map2(:,:)
      integer, allocatable :: thisgroup(:), allgroups(:,:), allheads(:)
      integer, allocatable :: numrecs(:)

!     OpenMP distribution

!     num_treads : number of threads set by environment OMP_NUM_THREADS

      integer :: num_threads

!     Accelerator data

      integer (kind=8), target :: memfre, memtot, memavail

      end module cidist

      module cidef

      implicit none

      integer :: nmol,mstates,nbase,nspin

      character (len=255) :: root,filinp,filout,filpro,fildbg,filday
      character (len=255) :: filsym,filciv,filvec,filint,filtst
      character (len=255) :: filf74,filf75,filone,filtwo,fildat
      integer :: lfninp,lfnout,lfnsym,lfnciv,lfnvec,lfnint
      integer :: lfnpro,lfndbg,lfnf74,lfnf75,lfnone,lfntwo,lfndat
      integer :: lfnday,lfntst

      character (len=255) :: user,host,date,time,cwd,command
      integer :: nprocs

      real (kind=8), allocatable :: civm(:,:), vecsm(:,:,:)
      integer, allocatable :: ioccm(:,:,:), nbasm(:), nactm(:)
      integer, allocatable :: inactm(:), idetm(:)
      real (kind=8), allocatable :: civb(:,:), vecsb(:,:,:)
      integer, allocatable :: ioccb(:,:)
      integer, allocatable ::  ncombv(:,:),inactb(:), nactb(:), idetb(:)
      real (kind=8), allocatable :: hbase(:,:), sbase(:,:), tbase(:,:)
      real (kind=8), allocatable :: hev(:)
      integer, allocatable :: nsing(:,:,:)
      integer :: iocch(32)

      end module cidef

      module gnome_parameters

      implicit none

      integer :: iocc(2),iord(2),indx(2),icalc,ipr,ipro,ipvec,idbg,itim
      integer :: iload,ibatch,itest,ifault,isolver,jsolver,idevel
      integer :: iswsvj,iswevj
      integer :: naccel,nacc0,nacc1,inpcib,intfil,ncols
      integer :: iday,idipole,itp4,ione,intadd,nummps,numomp,ixpert
      integer :: ins2
      integer :: ibase0,jbase0,idet0,jdet0
      integer :: ndeti,ndetj,nacti,nactj,inacti,inactj
      real (kind=8) :: thresh, thresh_CI, thresh_SIN
      real (kind=8) :: tolsvj,tolevj

      logical :: corres

      end module gnome_parameters

      module gnome_data

      integer, parameter :: intsrnk=250000000

      integer, parameter :: mncom=400,nqv=780,nsmf=590,nsrep=20
      integer, parameter :: ncf=10000,nkern=300,nprm=1182

      integer :: nnucl,nbas,numint


      integer :: nelec(2), ntcl(2), ntop(2)
      integer :: nalfa, nveca,nvecb, ntcla,ntclb, ntopa,ntopb
      integer :: nclose(2)
      integer :: nopen(2),  iocopen(100,2)                  

      character (len=80) :: text, name(2), namint(2)
      real (kind=8) :: potnuc
      character (len=4), allocatable :: centn(:)
      real (kind=8), allocatable :: xcord(:), ycord(:), zcord(:)
      real (kind=8), allocatable :: znuc(:)


      integer, allocatable :: ic(:), it(:), ictr(:)
      integer, allocatable :: icentn(:), itypen(:), icount(:), istand(:)

      integer :: ninact(2), nact(2), mnact, ittr(3)

      integer :: nbasis, mvec, nvec(2)

      integer, allocatable :: ioccup(:,:),itemp(:), ioccn(:,:)
      real (kind=8), allocatable :: vec(:,:,:),vtemp(:,:,:)
      real (kind=8), allocatable :: vna(:),vnb(:)
      real (kind=8) :: bias,deta,smat,hmat
      real (kind=8) :: etot,e1,e2,e2c
      real (kind=8) :: hh, ss
      real (kind=8) :: e2buff,e2summ
      integer :: ttest

      integer :: ntesta, ntestb
      integer :: nelecs,n1bas,nstdim, mbasel, ijend
      real (kind=8), allocatable :: va(:,:),vb(:,:),tb(:,:),a(:,:)
      real (kind=8), allocatable :: u(:,:),w(:,:),ev(:),sdiag(:),diag(:)
      real (kind=8), allocatable :: bsdiag(:),bdiag(:)
      real (kind=8), allocatable :: csdiag(:),cdiag(:)
      real (kind=8), allocatable :: st(:)
      real (kind=8), allocatable :: veca(:), vecb(:)
      real (kind=8), allocatable :: s12d(:,:)
      real (kind=8), allocatable :: w1(:),w2(:,:)
      real (kind=8), allocatable :: temp(:,:)

#ifdef SINGLEP
      real (kind=8), allocatable :: ta(:,:)
      real (kind=4), allocatable :: taa(:,:),aaa(:,:),tat(:,:)
      real (kind=4), allocatable :: aat(:,:),tt(:,:)
      real (kind=4), allocatable :: sm(:,:)
#else
      real (kind=8), allocatable :: taa(:,:),aaa(:,:),ta(:,:),tat(:,:)
      real (kind=8), allocatable :: aat(:,:),tt(:,:),sm(:,:)
#endif

      integer :: ising

      integer,external :: iloop
      real (kind=8), allocatable :: result(:,:), resultt(:,:)

      real (kind=8), allocatable :: work(:)
      integer (kind=8) :: lwrk, lwork, info

      real (kind=8) :: buffer(8)

      integer (kind=8) :: numdet
      integer (kind=8), allocatable :: ndxdet(:,:), lstdet(:,:)
      
      end module gnome_data

      module gnome_integrals

      integer :: intone,int1,int2,mint2,nt(3),myints(7)
      integer :: igfr,igto
      real (kind=8), allocatable :: s(:,:),t(:),v(:)
#ifdef SINGLEP
      real (kind=4), allocatable :: g(:),gg(:)
#else
      real (kind=8), allocatable :: g(:),gg(:)
#endif
      integer (kind=2), allocatable :: lab(:,:)
      integer, allocatable :: ig(:)
      integer (kind=4), allocatable :: ndxtv(:)
      integer :: nqcp,mbuf
      integer (kind=4), allocatable :: new_comm(:)
      integer (kind=4) :: new_group,grp_group
      integer, allocatable :: ntarget(:)

      end module gnome_integrals

#ifdef _OPENACC
      module cuda_functions
       use iso_c_binding
       interface
        integer (c_int) function cudaMemGetInfo(fre,tot)
     &                  bind(C,name="cudaMemGetInfo")
         use iso_c_binding
         implicit none
         type(c_ptr),value :: fre
         type(c_ptr),value :: tot
        end function cudaMemGetInfo

        integer (c_int) function cudaMalloc ( buff, size )
     &                  bind (C, name="cudaMalloc" )
         use iso_c_binding
         implicit none
         type (c_ptr)  :: buff
         integer (c_size_t), value :: size
        end function cudaMalloc

        integer (c_int) function cudaMemcpy ( dst, src, count, kind )
     &                  bind (C, name="cudaMemcpy" )
!       note: cudaMemcpyHostToDevice=1
!       note: cudaMemcpyDeviceToHost=2
         use iso_c_binding
         type (c_ptr), value :: dst, src
         integer (c_size_t), value :: count, kind
        end function cudaMemcpy

        integer (c_int) function cudaFree(buff)
     &                  bind(C, name="cudaFree")
         use iso_c_binding
         implicit none
         type (c_ptr), value :: buff
        end function cudaFree
       end interface
      end module cuda_functions
#endif

#ifdef CUSOLVER
      module cuda_cusolver
      
      use cusolverDN
      use cudafor
      
      type gesvdjInfo
      type(c_ptr) :: svInfo
      end type gesvdjInfo
      
      type syevjInfo
      type(c_ptr) :: svInfo
      end type syevjInfo
      
      interface
        integer(c_int) function cusolverDnCreateGesvdjInfo(info)
     &       bind(C,name='cusolverDnCreateGesvdjInfo')
        import gesvdjInfo
        type(gesvdjInfo) :: info
        end function cusolverDnCreateGesvdjInfo
      end interface
      
      interface
        integer(c_int) function cusolverDnXgesvdjSetTolerance
     &       (info, tolerance)
     &       bind(C,name='cusolverDnXgesvdjSetTolerance')
        use iso_c_binding
        import gesvdjInfo
        type(gesvdjInfo)       :: info
        real(c_double), value  :: tolerance
        end function cusolverDnXgesvdjSetTolerance
      end interface
      
      interface
        integer(c_int) function cusolverDnXgesvdjSetMaxSweeps
     &       (info, max_sweeps)
     &       bind(C,name='cusolverDnXgesvdjSetMaxSweeps')
        use iso_c_binding
        import gesvdjInfo
        type(gesvdjInfo)       :: info
        integer(c_int),value   :: max_sweeps
        end function cusolverDnXgesvdjSetMaxSweeps
      end interface
      
      interface
        integer(c_int) function cusolverDnXgesvdjGetSweeps
     &       (cusolver_Hndl, info, executed_sweeps)
     &       bind(C,name='cusolverDnXgesvdjGetSweeps')
        use iso_c_binding
        import cusolverDnHandle
        import gesvdjInfo
        type(cusolverDnHandle), value :: cusolver_Hndl
        type(gesvdjInfo)              :: info
        integer(c_int)                :: executed_sweeps
        end function cusolverDnXgesvdjGetSweeps
      end interface
      
      interface
        integer(c_int) function cusolverDnXgesvdjGetResidual
     &       (cusolver_Hndl, info, residual)
     &       bind(C,name='cusolverDnXgesvdjGetResidual')
        use iso_c_binding
        import cusolverDnHandle
        import gesvdjInfo
        type(cusolverDnHandle), value :: cusolver_Hndl
        type(gesvdjInfo)              :: info
        real(c_double)                :: residual
        end function cusolverDnXgesvdjGetResidual
      end interface
      
      interface
        integer(c_int) function cusolverDnDgesvdj_bufferSize
     &       (cusolver_Hndl, jobz, econ,
     &       m, n, a, lda, s, u, ldu, v, ldv, lwork, info)
     &       bind(C,name='cusolverDnDgesvdj_bufferSize')
        use iso_c_binding
        import cusolverDnHandle
        import gesvdjInfo
        type(cusolverDnHandle), value :: cusolver_Hndl
        integer(c_int), value  :: jobz
        integer(c_int), value  :: econ
        integer(c_int), value  :: m, n, lda, ldu, ldv
        integer(c_int)         :: lwork
        real(c_double), device :: a(*), s(*), u(*), v(*)
        type(gesvdjInfo)       :: info
        end function cusolverDnDgesvdj_bufferSize
      end interface
      
      interface
        integer(c_int) function cusolverDnDgesvdj
     &       (cusolver_Hndl, jobz, econ,
     &       m, n, a, lda, s, u, ldu, v, ldv, work, lwork,
     &       devinfo, info )
     &       bind(C,name='cusolverDnDgesvdj')
        use iso_c_binding
        import cusolverDnHandle
        import gesvdjInfo
        type(cusolverDnHandle), value :: cusolver_Hndl
        integer(c_int), value  :: jobz
        integer(c_int), value  :: econ, m, n, lda, ldu, ldv
        real(c_double), device :: a(*), s(*), u(*), v(*), work(*)
        integer(c_int)         :: lwork
        integer(c_int)         :: devinfo
        type(gesvdjInfo)       :: info
        end function cusolverDnDgesvdj
      end interface
      
      interface
        integer(c_int) function cusolverDnCreateSyevjInfo(info)
     &       bind(C,name='cusolverDnCreateSyevjInfo')
        import syevjInfo
        type(syevjInfo) :: info
        end function cusolverDnCreateSyevjInfo
      end interface
      
      interface
        integer(c_int) function cusolverDnXsyevjSetTolerance
     &       (info, tolerance)
     &       bind(C,name='cusolverDnXsyevjSetTolerance')
        use iso_c_binding
        import syevjInfo
        type(syevjInfo)       :: info
        real(c_double), value  :: tolerance
        end function cusolverDnXsyevjSetTolerance
      end interface
      
      interface
        integer(c_int) function cusolverDnXsyevjSetMaxSweeps
     &       (info, max_sweeps)
     &       bind(C,name='cusolverDnXsyevjSetMaxSweeps')
        use iso_c_binding
        import syevjInfo
        type(syevjInfo)       :: info
        integer(c_int),value   :: max_sweeps
        end function cusolverDnXsyevjSetMaxSweeps
      end interface
      
      interface
        integer(c_int) function cusolverDnXsyevjGetSweeps
     &       (cusolver_Hndl, info, executed_sweeps)
     &       bind(C,name='cusolverDnXsyevjGetSweeps')
        use iso_c_binding
        import cusolverDnHandle
        import syevjInfo
        type(cusolverDnHandle), value :: cusolver_Hndl
        type(syevjInfo)              :: info
        integer(c_int)                :: executed_sweeps
        end function cusolverDnXsyevjGetSweeps
      end interface
      
      interface
        integer(c_int) function cusolverDnXsyevjGetResidual
     &       (cusolver_Hndl, info, residual)
     &       bind(C,name='cusolverDnXsyevjGetResidual')
        use iso_c_binding
        import cusolverDnHandle
        import syevjInfo
        type(cusolverDnHandle), value :: cusolver_Hndl
        type(syevjInfo)              :: info
        real(c_double)                :: residual
        end function cusolverDnXsyevjGetResidual
      end interface
      
      interface
        integer(c_int) function cusolverDnDsyevj_bufferSize
     &       (cusolver_Hndl, jobz, uplo,
     &       n, a, lda, v, lwork, info)
     &       bind(C,name='cusolverDnDsyevj_bufferSize')
        use iso_c_binding
        import cusolverDnHandle
        import syevjInfo
        type(cusolverDnHandle), value :: cusolver_Hndl
        integer(c_int), value  :: jobz
        integer(c_int), value  :: uplo
        integer(c_int), value  :: n, lda
        integer(c_int)         :: lwork
        real(c_double), device :: a(*), v(*)
        type(syevjInfo)       :: info
        end function cusolverDnDsyevj_bufferSize
      end interface
      
      interface
        integer(c_int) function cusolverDnDsyevj
     &       (cusolver_Hndl, jobz, uplo,
     &       n, a, lda, v, work, lwork,
     &       devinfo, info )
     &       bind(C,name='cusolverDnDsyevj')
        use iso_c_binding
        import cusolverDnHandle
        import syevjInfo
        type(cusolverDnHandle), value :: cusolver_Hndl
        integer(c_int), value  :: jobz
        integer(c_int), value  :: uplo
        integer(c_int), value  :: n, lda
        real(c_double), device :: a(*), v(*), work(*)
        integer(c_int)         :: lwork
        integer(c_int)         :: devinfo
        type(syevjInfo)       :: info
        end function cusolverDnDsyevj
      end interface
      
      type(cusolverDnHandle) :: cusolver_handle
      type(gesvdjInfo)       :: gesvdj_params
      type(syevjInfo)        :: syevj_params
      integer (kind=4)       :: cusolver_status
      integer (kind=4)       :: worksize,worksize2,ndim,mdim
      real (kind=8)          :: tol, residual
      integer (kind=4)       :: max_sweeps, exec_sweeps
      integer (kind=4), parameter :: econ=0
      real (kind=8), allocatable :: rwork(:)
      real (kind=8), allocatable :: workspace_d(:)
      integer (kind=4) :: dev_info_d
      integer(kind=cuda_stream_kind) :: stream
      integer (kind=4) :: jobz
      integer (kind=4) :: uplo=0
      
      end module cuda_cusolver
#endif
