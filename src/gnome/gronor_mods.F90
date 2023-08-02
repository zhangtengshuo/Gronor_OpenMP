!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ basis, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!> @brief
!! GronOR module definitions
!!
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!

module const
  integer :: nr(20,3)
  data nr/0,1,0,0,2,0,0,1,1,0,3,0,0,2,2,1,0,1,0,1, &
      0,0,1,0,0,2,0,1,0,1,0,3,0,1,0,2,2,0,1,1, &
      0,0,0,1,0,0,2,0,1,1,0,0,3,0,1,0,1,2,2,1/
  real (kind=8) :: pi,twopi,dnorm,piroot,piterm,enorm
  data pi,twopi,dnorm/3.1415926535898d0,6.2831853071796d0,32.0d0/
  data piroot,piterm,enorm/0.88622692545276d0,34.986836655250d0,1.0d-12/
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

  integer (kind=4) :: me,np,mstr
  integer (kind=4) :: numdev,mydev
  integer (kind=4) :: group_batch,group_heads,group_world
  integer (kind=4), allocatable :: ranks_list(:),ranks_heads(:)

  integer :: nnodes,nrsets,nranks,ncycls,nrnsets,ngpus

  integer :: numacc,numnon,numgrp,numgrn,numidle
  integer :: ngr,meg,npg,comm_group,comm_first,mgr,nalive,nabort
  integer :: mynode,mygroup,myhead,ime,mpibuf
  integer :: iamhead,iamacc,iamactive
  integer (kind=4), allocatable  :: map1(:,:),map2(:,:)
  integer, allocatable :: thisgroup(:),allgroups(:,:),allheads(:)
  integer, allocatable :: numrecs(:)
  integer (kind=8), allocatable :: ipbuf(:,:),itbuf(:,:)

  !     OpenMP distribution

  !     num_treads : number of threads set by environment OMP_NUM_THREADS

  integer :: num_threads

  !     Accelerator data

  integer (kind=8), target :: memfre,memtot,memavail
  integer (kind=4), allocatable :: igrn(:)

  integer (kind=8) :: managers,numwrk,maxbuf,numbuf
  integer (kind=8), allocatable :: mgrbuf(:,:)
  integer (kind=8), allocatable :: mgrwrk(:,:)
  integer (kind=8), allocatable :: mipbuf(:,:)
  
  character (len=8) :: machine

  integer (kind=4), parameter :: master=1
  integer (kind=4), parameter :: manager=2
  integer (kind=4), parameter :: worker=3
  integer (kind=4), parameter :: idle=4

  character (len=1), parameter :: crole(4)=(/"M","m","w","i"/)

  integer (kind=4) :: role
  integer (kind=8) :: nperman,numman

end module cidist

module makebasedata
  implicit none
  integer                        :: spinFrag
  integer,allocatable            :: micro_ndets(:)
  real(kind=8),allocatable       :: coef(:),micro_coef(:)
  real(kind=8),allocatable       :: coef_new(:)
  character(len=255),allocatable :: occ(:),micro_occ(:)
  character(len=255),allocatable :: occ_new(:)
  character(len=8)               :: spin_mult(10)
  data spin_mult/'singlet','doublet','triplet','quartet','quintet','sextet',    &
                 'septet','octet','nonet','decaplet'/
end module makebasedata

module cidef

  implicit none

  integer :: nmol,mstates,nbase,nspin

  character (len=255) :: root,filinp,filout,filpro,fildbg,filday
  character (len=255) :: filsys,filciv,filvec,filint,filtst,fildet
  character (len=255) :: filone,filtwo,fildat,fillog,filtim
  character (len=255) :: filcpr,filarx,filrnk,filcml,filxrx
  integer :: lfninp,lfnout,lfnsys,lfnciv,lfnvec,lfnint,lfndet
  integer :: lfnpro,lfndbg,lfnone,lfntwo,lfndat,lfntim,lfnabt
  integer :: lfnday,lfntst,lfnlog,lfncpr,lfnarx,lfnrnk,lfncml
  integer :: lfnwrn,lfnxrx

  character (len=255) :: user,host,date,time,cwd,command
  integer :: nprocs
  integer :: nciaux

  integer :: maxci,maxcib,maxnact,mactb,maxvec,nci,numact

  real (kind=8), allocatable :: civm(:,:),vecsm(:,:,:)
  integer, allocatable :: ioccm(:,:,:),nbasm(:),nactm(:)
  integer, allocatable :: inactm(:),idetm(:),spinm(:),nFrozen(:)
  integer, allocatable :: nElectrons(:)
  integer, allocatable :: inter_couplings(:,:)
  real (kind=8), allocatable :: civb(:,:),vecsb(:,:,:)
  real (kind=8), allocatable :: ecasscf(:),ecaspt2(:)
  integer (kind=1), allocatable :: iocc(:,:,:)
  integer (kind=8), allocatable :: ncombv(:,:),inactb(:)
  integer (kind=8), allocatable :: nactb(:),idetb(:),alldets(:)
  real (kind=8), allocatable :: hbase(:,:),sbase(:,:),tbase(:,:)
  real (kind=8), allocatable :: dqbase(:,:,:)
  real (kind=8), allocatable :: hev(:)
  real (kind=8), allocatable :: nociwf(:,:)
  logical, allocatable :: bpdone(:,:)
  integer, allocatable :: nsing(:,:,:)
  integer :: iocch(32)
  character(len=255), allocatable :: occm_string(:,:)

  logical :: vecdet,lcpr

  integer :: ncorr,nwt,lablen
  real(kind=8), allocatable :: ecorr(:)           ! the shifts
  real(kind=8), allocatable :: hcorr(:,:) ! the shifted Hamiltonian
  real(kind=8)              :: mnuc(9),maxcoef
  character(len=16),allocatable  :: fragLabel(:)
  character(len=128),allocatable :: mebfLabel(:)
  character(len=128)             :: header,key
  logical                        :: mebfLabels
  integer, allocatable  :: nCntr(:,:)

end module cidef

module gnome_parameters

  implicit none

  integer :: icalc,ipr,ipro,ipvec,idbg,itim
  integer :: itest,ifault,isolver,jsolver,idevel,idist,labmax
  integer :: ntask,ntaska,nbatch,nbatcha
  integer :: ndbg,mdbg,load,loada
  integer :: iaslvr,jaslvr,inslvr,jnslvr
  integer :: iswsvj,iswevj
  integer :: naccel,nacc0,nacc1,inpcib,intfil,ncols
  integer :: iday,idipole,itp4,nummps,numgpu,ixpert
  integer :: ins2
  integer :: ibase0,jbase0,idet0,jdet0
  integer :: ndeti,ndetj,nacti,nactj,inacti,inactj
  real (kind=8) :: tau_MO,tau_CI,tau_SIN,thresh,thresh_SIN
  real (kind=8) :: tolsvj,tolevj
  integer :: numtasks
  real (kind=8) :: tmax

  integer (kind=4) :: itreq
  integer (kind=8) :: irbuf(4)
  logical :: odbg,odupl,oterm,otreq,otimeout
  logical :: corres

  character (len=128) :: mebfroot,combas
  character (len=16), allocatable :: fragname(:)
  character (len=16), allocatable :: fragstate(:,:)
  character (len=128), allocatable :: fragfile(:)
  character (len=128), allocatable :: vecfile(:)
  character (len=128), allocatable :: detfile(:)

end module gnome_parameters

module gnome_data

  integer, parameter :: intsrnk=250000000

  integer, parameter :: mncom=400,nqv=780,nsmf=590,nsrep=20
  integer, parameter :: ncf=10000,nkern=300,nprm=1182

  integer :: nnucl,nbas,numint,numfiles

  integer :: nelec(2),ntcl(2),ntop(2)
  integer :: nalfa,nveca,nvecb,ntcla,ntclb,ntopa,ntopb
  integer :: nclose(2)
  integer :: nopen(2),iocopen(100,2)

  character (len=80) :: text,name(2),namint(2)
  real (kind=8) :: potnuc,zNucTot
  character (len=4), allocatable :: centn(:)
  real (kind=8), allocatable :: xcord(:),ycord(:),zcord(:)
  real (kind=8), allocatable :: znuc(:)
  real (kind=8)              :: com(3)


  integer, allocatable :: ic(:),it(:),ictr(:)
  integer, allocatable :: icentn(:),itypen(:),icount(:),istand(:)

  integer :: ninact(2),nact(2),mnact,ittr(3)

  integer :: nbasis,mvec,nvec(2)

  integer :: nb0,nb1

  integer, allocatable :: ioccup(:,:),itemp(:),ioccn(:,:)
  real (kind=8), allocatable :: vec(:,:,:),vtemp(:,:,:)
  real (kind=8) :: bias,deta,smat,hmat
  real (kind=8) :: etot,e1,e2,e2c,etotb,fac,fctr
  real (kind=8) :: mpoles(9)
  real (kind=8) :: e1tot,e2tot,sstot
  real (kind=8) :: hh,ss
  real (kind=8) :: e2buff,e2summ
  integer :: ttest

  integer :: ntesta,ntestb
  integer :: nelecs,n1bas,nstdim,mbasel,ijend
  real (kind=8), allocatable :: va(:,:),vb(:,:),tb(:,:)
#ifdef HIPSOLVER
  real (kind=8), allocatable, target :: a(:,:)
  real (kind=8), allocatable, target :: u(:,:),w(:,:),ev(:)
#else
#ifdef ROCSOLVER
  real (kind=8), allocatable, target :: a(:,:)
  real (kind=8), allocatable, target :: u(:,:),w(:,:),ev(:)
#else
  real (kind=8), allocatable :: a(:,:)
  real (kind=8), allocatable :: u(:,:),w(:,:),ev(:)
#endif
#endif
  real (kind=8), allocatable :: sdiag(:)
  real (kind=8), allocatable, target :: diag(:)
  real (kind=8), allocatable :: bsdiag(:),bdiag(:)
  real (kind=8), allocatable :: csdiag(:),cdiag(:)
  real (kind=8), allocatable :: st(:)
  real (kind=8), allocatable :: veca(:),vecb(:)
  real (kind=8), allocatable :: s12d(:,:)
  real (kind=8), allocatable :: w1(:),w2(:,:)
  real (kind=8), allocatable :: temp(:,:)

  real (kind=8), allocatable :: taa(:,:)
#ifdef HIPSOLVER
  real (kind=8), allocatable, target :: ta(:,:)
#else
#ifdef ROCSOLVER
  real (kind=8), allocatable, target :: ta(:,:)
#else
  real (kind=8), allocatable :: ta(:,:)
#endif
#endif
  real (kind=8), allocatable :: sm(:,:)
  real (kind=8), allocatable :: aaa(:,:),aat(:,:)
  real (kind=8), allocatable :: tt(:,:)
#ifdef SINGLEP
  real (kind=8), allocatable :: diagl(:,:),bdiagl(:,:)
  real (kind=8), allocatable :: csdiagl(:,:),bsdiagl(:,:)
  real (kind=4), allocatable :: sml(:,:,:),prefac(:)
  real (kind=4), allocatable :: aaal(:,:,:),aatl(:,:,:)
  real (kind=4), allocatable :: ttl(:,:,:),tatl(:,:,:)
  real (kind=4), allocatable :: tal(:,:,:)
#else
  real (kind=8), allocatable :: diagl(:,:),bdiagl(:,:)
  real (kind=8), allocatable :: csdiagl(:,:),bsdiagl(:,:)
  real (kind=8), allocatable :: sml(:,:,:),prefac(:)
  real (kind=8), allocatable :: aaal(:,:,:),aatl(:,:,:)
  real (kind=8), allocatable :: ttl(:,:,:),tatl(:,:,:)
  real (kind=8), allocatable :: tal(:,:,:)
#endif
#ifdef SINGLEP
  real (kind=4), allocatable :: sm0(:,:,:),prefac0(:)
  real (kind=4), allocatable :: aaa0(:,:,:),aat0(:,:,:)
  real (kind=4), allocatable :: tt0(:,:,:)
  real (kind=4), allocatable :: ta0(:,:,:)
#else
  real (kind=8), allocatable :: sm0(:,:,:),prefac0(:)
  real (kind=8), allocatable :: aaa0(:,:,:),aat0(:,:,:)
  real (kind=8), allocatable :: tt0(:,:,:)
  real (kind=8), allocatable :: ta0(:,:,:)
#endif
#ifdef SINGLEP
  real (kind=8), allocatable :: diag1(:,:),bdiag1(:,:)
  real (kind=8), allocatable :: csdiag1(:,:),bsdiag1(:,:)
  real (kind=4), allocatable :: sm1(:,:,:),prefac1(:)
  real (kind=4), allocatable :: aaa1(:,:,:),aat1(:,:,:)
  real (kind=4), allocatable :: tt1(:,:,:)
  real (kind=4), allocatable :: ta1(:,:,:)
#else
  real (kind=8), allocatable :: diag1(:,:),bdiag1(:,:)
  real (kind=8), allocatable :: csdiag1(:,:),bsdiag1(:,:)
  real (kind=8), allocatable :: sm1(:,:,:),prefac1(:)
  real (kind=8), allocatable :: aaa1(:,:,:),aat1(:,:,:)
  real (kind=8), allocatable :: tt1(:,:,:)
  real (kind=8), allocatable :: ta1(:,:,:)
#endif

  integer :: ising

  real (kind=8), allocatable :: result(:,:),resultt(:,:)

  real (kind=8), allocatable :: work(:)
  integer (kind=8) :: lwrk,lwork,info

  real (kind=8) :: buffer(17)

  integer (kind=8) :: numdet,melen,memax,icur,jcur
  integer (kind=4), allocatable :: melist(:,:)
  integer (kind=8), allocatable :: ndxdet(:,:)

  real (kind=8) :: gbmelist

  real (kind=8), allocatable :: c2sum(:,:)
  integer (kind=8), allocatable :: nbdet(:,:)
#ifdef ROCSOLVER 
  real (kind=8), allocatable, target :: rwork(:)
#else
#ifdef HIPSOLVER 
  real (kind=8), allocatable, target :: rwork(:)
#else
  real (kind=8), allocatable :: rwork(:)
#endif
#endif

end module gnome_data

module gnome_integrals

  integer :: intone,int1,int2,mint2,nt(3),myints(7),intndx,jntndx
  integer :: igfr,igto
  real (kind=8), allocatable :: s(:,:),t(:),v(:),dqm(:,:)
#ifdef SINGLEP
  real (kind=4), allocatable :: g(:),gg(:)
#else
  real (kind=8), allocatable :: g(:),gg(:)
#endif
  integer (kind=2), allocatable :: lab(:,:)
  integer (kind=8), allocatable :: ndx(:)
  integer (kind=8), allocatable :: ndxk(:)
  integer, allocatable :: ig(:)
  integer (kind=4), allocatable :: ndxtv(:)
  integer :: mbuf,mclab
  integer (kind=4), allocatable :: new_comm(:)
  integer (kind=4), allocatable :: new_comm1(:)
  integer (kind=4), allocatable :: new_comm2(:,:)
  integer (kind=4) :: new_group,grp_group
  integer, allocatable :: ntarget(:)

  integer (kind=8) :: mlab,icomm1

  logical :: lcomm1

end module gnome_integrals

#ifdef CUDA
module cuda_functions
  use iso_c_binding
  interface
    integer (c_int) function cudaMemGetInfo(fre,tot) bind(C,name="cudaMemGetInfo")
      use iso_c_binding
      implicit none
      type(c_ptr),value :: fre
      type(c_ptr),value :: tot
    end function cudaMemGetInfo

    integer (c_int) function cudaMalloc ( buff,size ) bind (C,name="cudaMalloc" )
      use iso_c_binding
      implicit none
      type (c_ptr)  :: buff
      integer (c_size_t), value :: size
    end function cudaMalloc

    integer (c_int) function cudaMemcpy ( dst,src,count,kind ) bind (C,name="cudaMemcpy" )
      !       note: cudaMemcpyHostToDevice=1
      !       note: cudaMemcpyDeviceToHost=2
      use iso_c_binding
      type (c_ptr), value :: dst,src
      integer (c_size_t), value :: count,kind
    end function cudaMemcpy

    integer (c_int) function cudaFree(buff) bind(C,name="cudaFree")
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

#ifdef CUSOLVERJ
  type gesvdjInfo
    type(c_ptr) :: svInfo
  end type gesvdjInfo

  type syevjInfo
    type(c_ptr) :: svInfo
  end type syevjInfo

  interface
    integer(c_int) function cusolverDnCreateGesvdjInfo(info) &
        bind(C,name='cusolverDnCreateGesvdjInfo')
      import gesvdjInfo
      type(gesvdjInfo) :: info
    end function cusolverDnCreateGesvdjInfo
  end interface

  interface
    integer(c_int) function cusolverDnXgesvdjSetTolerance(info,tolerance) &
        bind(C,name='cusolverDnXgesvdjSetTolerance')
      use iso_c_binding
      import gesvdjInfo
      type(gesvdjInfo)       :: info
      real(c_double), value  :: tolerance
    end function cusolverDnXgesvdjSetTolerance
  end interface

  interface
    integer(c_int) function cusolverDnXgesvdjSetMaxSweeps(info,max_sweeps) &
        bind(C,name='cusolverDnXgesvdjSetMaxSweeps')
      use iso_c_binding
      import gesvdjInfo
      type(gesvdjInfo)       :: info
      integer(c_int),value   :: max_sweeps
    end function cusolverDnXgesvdjSetMaxSweeps
  end interface

  interface
    integer(c_int) function cusolverDnXgesvdjGetSweeps(cusolver_Hndl,info,executed_sweeps) &
        bind(C,name='cusolverDnXgesvdjGetSweeps')
      use iso_c_binding
      import cusolverDnHandle
      import gesvdjInfo
      type(cusolverDnHandle), value :: cusolver_Hndl
      type(gesvdjInfo)              :: info
      integer(c_int)                :: executed_sweeps
    end function cusolverDnXgesvdjGetSweeps
  end interface

  interface
    integer(c_int) function cusolverDnXgesvdjGetResidual(cusolver_Hndl,info,residual) &
        bind(C,name='cusolverDnXgesvdjGetResidual')
      use iso_c_binding
      import cusolverDnHandle
      import gesvdjInfo
      type(cusolverDnHandle), value :: cusolver_Hndl
      type(gesvdjInfo)              :: info
      real(c_double)                :: residual
    end function cusolverDnXgesvdjGetResidual
  end interface

  interface
    integer(c_int) function cusolverDnDgesvdj_bufferSize(cusolver_Hndl,jobz,econ, &
        m,n,a,lda,s,u,ldu,v,ldv,lwork,info) bind(C,name='cusolverDnDgesvdj_bufferSize')
      use iso_c_binding
      import cusolverDnHandle
      import gesvdjInfo
      type(cusolverDnHandle), value :: cusolver_Hndl
      integer(c_int), value    :: jobz
      integer(c_int), value    :: econ
      integer(c_int), value    :: m,n,lda,ldu,ldv
      integer(c_int)           :: lwork
      real(c_double), device   :: a(:,:),s(:),u(:,:),v(:,:)
      type(gesvdjInfo)         :: info
    end function cusolverDnDgesvdj_bufferSize
  end interface

  interface
    integer(c_int) function cusolverDnDgesvdj(cusolver_Hndl,jobz,econ,m,n,a,lda,s,u,ldu,v,ldv, &
        work,lwork,devinfo,info ) bind(C,name='cusolverDnDgesvdj')
      use iso_c_binding
      import cusolverDnHandle
      import gesvdjInfo
      type(cusolverDnHandle), value :: cusolver_Hndl
      integer(c_int), value   :: jobz
      integer(c_int), value   :: econ,m,n,lda,ldu,ldv
      real(c_double), device  :: a(:,:),s(:),u(:,:),v(:,:),work(:)
      integer(c_int)          :: lwork
      integer(c_int)          :: devinfo
      type(gesvdjInfo)        :: info
    end function cusolverDnDgesvdj
  end interface

  interface
    integer(c_int) function cusolverDnCreateSyevjInfo(info) bind(C,name='cusolverDnCreateSyevjInfo')
      import syevjInfo
      type(syevjInfo) :: info
    end function cusolverDnCreateSyevjInfo
  end interface

  interface
    integer(c_int) function cusolverDnXsyevjSetTolerance(info,tolerance) &
        bind(C,name='cusolverDnXsyevjSetTolerance')
      use iso_c_binding
      import syevjInfo
      type(syevjInfo)        :: info
      real(c_double), value  :: tolerance
    end function cusolverDnXsyevjSetTolerance
  end interface

  interface
    integer(c_int) function cusolverDnXsyevjSetMaxSweeps(info,max_sweeps) &
        bind(C,name='cusolverDnXsyevjSetMaxSweeps')
      use iso_c_binding
      import syevjInfo
      type(syevjInfo)       :: info
      integer(c_int),value  :: max_sweeps
    end function cusolverDnXsyevjSetMaxSweeps
  end interface

  interface
    integer(c_int) function cusolverDnXsyevjGetSweeps(cusolver_Hndl,info,executed_sweeps) &
        bind(C,name='cusolverDnXsyevjGetSweeps')
      use iso_c_binding
      import cusolverDnHandle
      import syevjInfo
      type(cusolverDnHandle), value :: cusolver_Hndl
      type(syevjInfo)               :: info
      integer(c_int)                :: executed_sweeps
    end function cusolverDnXsyevjGetSweeps
  end interface

  interface
    integer(c_int) function cusolverDnXsyevjGetResidual(cusolver_Hndl,info,residual) &
        bind(C,name='cusolverDnXsyevjGetResidual')
      use iso_c_binding
      import cusolverDnHandle
      import syevjInfo
      type(cusolverDnHandle), value :: cusolver_Hndl
      type(syevjInfo)               :: info
      real(c_double)                :: residual
    end function cusolverDnXsyevjGetResidual
  end interface

  interface
    integer(c_int) function cusolverDnDsyevj_bufferSize(cusolver_Hndl,jobz,uplo, &
        n,a,lda,v,lwork,info) bind(C,name='cusolverDnDsyevj_bufferSize')
      use iso_c_binding
      import cusolverDnHandle
      import syevjInfo
      type(cusolverDnHandle), value :: cusolver_Hndl
      integer(c_int), value         :: jobz
      integer(c_int), value         :: uplo
      integer(c_int), value         :: n,lda
      integer(c_int)                :: lwork
      real(c_double) , device       :: a(:,:),v(:)
      type(syevjInfo)               :: info
    end function cusolverDnDsyevj_bufferSize
  end interface

  interface
    integer(c_int) function cusolverDnDsyevj(cusolver_Hndl,jobz,uplo,n,a,lda,v,work,lwork, &
        devinfo,info) bind(C,name='cusolverDnDsyevj')
      use iso_c_binding
      import cusolverDnHandle
      import syevjInfo
      type(cusolverDnHandle), value :: cusolver_Hndl
      integer(c_int), value   :: jobz
      integer(c_int), value   :: uplo
      integer(c_int), value   :: n,lda
      real(c_double), device  :: a(:,:),v(:),work(:)
      integer(c_int)          :: lwork
      integer(c_int)          :: devinfo
      type(syevjInfo)         :: info
    end function cusolverDnDsyevj
  end interface
#endif

  type(cusolverDnHandle) :: cusolver_handle
#ifdef CUSOLVERJ
  type(gesvdjInfo)       :: gesvdj_params
  type(syevjInfo)        :: syevj_params
#endif
  integer (kind=4)       :: cusolver_status
  integer (kind=4)       :: lwork1,lwork2,ndim,mdim
  real (kind=8)          :: tol,residual
  integer (kind=4)       :: max_sweeps,exec_sweeps
  integer (kind=4), parameter :: econ=0
  real (kind=8), allocatable :: workspace_d(:)
  integer (kind=4) :: dev_info_d
  integer(kind=cuda_stream_kind) :: stream
  integer (kind=4) :: jobz
  integer (kind=4) :: uplo=0

end module cuda_cusolver
#endif

#ifdef HIPSOLVER
module hipsolver_enums
  use iso_c_binding

  !---------------------!
  !   hipSOLVER types   !
  !---------------------!

  enum, bind(c)
    enumerator :: HIPSOLVER_OP_N = 111
    enumerator :: HIPSOLVER_OP_T = 112
    enumerator :: HIPSOLVER_OP_C = 113
  end enum

  enum, bind(c)
    enumerator :: HIPSOLVER_FILL_MODE_UPPER = 121
    enumerator :: HIPSOLVER_FILL_MODE_LOWER = 122
  end enum

  enum, bind(c)
    enumerator :: HIPSOLVER_SIDE_LEFT  = 141
    enumerator :: HIPSOLVER_SIDE_RIGHT = 142
  end enum

  enum, bind(c)
    enumerator :: HIPSOLVER_EIG_MODE_NOVECTOR = 201
    enumerator :: HIPSOLVER_EIG_MODE_VECTOR   = 202
  end enum

  enum, bind(c)
    enumerator :: HIPSOLVER_EIG_TYPE_1 = 211
    enumerator :: HIPSOLVER_EIG_TYPE_2 = 212
    enumerator :: HIPSOLVER_EIG_TYPE_3 = 213
  end enum

  enum, bind(c)
    enumerator :: HIPSOLVER_EIG_RANGE_ALL = 221
    enumerator :: HIPSOLVER_EIG_RANGE_V   = 222
    enumerator :: HIPSOLVER_EIG_RANGE_I   = 223
  end enum

  enum, bind(c)
    enumerator :: HIPSOLVER_STATUS_SUCCESS           = 0
    enumerator :: HIPSOLVER_STATUS_NOT_INITIALIZED   = 1
    enumerator :: HIPSOLVER_STATUS_ALLOC_FAILED      = 2
    enumerator :: HIPSOLVER_STATUS_INVALID_VALUE     = 3
    enumerator :: HIPSOLVER_STATUS_MAPPING_ERROR     = 4
    enumerator :: HIPSOLVER_STATUS_EXECUTION_FAILED  = 5
    enumerator :: HIPSOLVER_STATUS_INTERNAL_ERROR    = 6
    enumerator :: HIPSOLVER_STATUS_NOT_SUPPORTED     = 7
    enumerator :: HIPSOLVER_STATUS_ARCH_MISMATCH     = 8
    enumerator :: HIPSOLVER_STATUS_HANDLE_IS_NULLPTR = 9
    enumerator :: HIPSOLVER_STATUS_INVALID_ENUM      = 10
    enumerator :: HIPSOLVER_STATUS_UNKNOWN           = 11
  end enum

end module hipsolver_enums

module hipvars
  use iso_c_binding
  type(c_ptr)            :: hipsolver_handle
  integer (kind=4)       :: hipsolver_status
  integer (kind=4)       :: lwork1,mdim
  integer (kind=4)       :: lwork2
  integer (kind=4)       :: ndim
  real (kind=8)          :: tol,residual
  integer (kind=4)       :: max_sweeps,exec_sweeps
  integer (kind=4), parameter :: econ=0
  real (kind=8), allocatable, target :: workspace_d(:)
  integer (kind=4) :: dev_info_d
  !     integer(kind=cuda_stream_kind) :: stream
  integer (kind=8) :: jobz
  integer (kind=8) :: uplo=0      
end module hipvars

module amd_hipsolver
  use iso_c_binding
  use hipsolver_enums
  !      use hipsolverDn     
  !      type(hipsolverDnHandle) :: hipsolver_handle

  !       type(c_ptr)            :: hipsolver_handle
  !       integer (kind=4)       :: hipsolver_status
  !       integer (kind=4)       :: lwork1,lwork2,ndim,mdim
  !       real (kind=8)          :: tol,residual
  !       integer (kind=4)       :: max_sweeps,exec_sweeps
  !       integer (kind=4), parameter :: econ=0
  !       real (kind=8), allocatable :: rwork(:)
  !       real (kind=8), allocatable :: workspace_d(:)
  !       integer (kind=4) :: dev_info_d
  !!      integer(kind=cuda_stream_kind) :: stream
  !       integer (kind=4) :: jobz
  !       integer (kind=4) :: uplo=0

  interface
    function hipsolverCreate(handle) bind(c,name = 'hipsolverCreate')
      use iso_c_binding
      use hipsolver_enums
      implicit none
      integer(kind(HIPSOLVER_STATUS_SUCCESS)) :: hipsolverCreate
      type(c_ptr), value :: handle
    end function hipsolverCreate
  end interface

  interface
    function hipsolverDestroy(handle) bind(c,name = 'hipsolverDestroy')
      use iso_c_binding
      use hipsolver_enums
      implicit none
      integer(kind(HIPSOLVER_STATUS_SUCCESS)) :: hipsolverDestroy
      type(c_ptr), value :: handle
    end function hipsolverDestroy
  end interface

  interface
    function hipsolverDgesvd_bufferSize (handle,jobu,jobv,m,n,lwork) &
        bind(c,name = 'hipsolverDgesvd_bufferSize')
      use iso_c_binding
      use hipsolver_enums
      implicit none
      integer(kind(HIPSOLVER_STATUS_SUCCESS)) :: hipsolverDgesvd_bufferSize
      type(c_ptr), value :: handle
      integer(c_signed_char), value :: jobu
      integer(c_signed_char), value :: jobv
      integer(c_int), value :: m
      integer(c_int), value :: n
      integer(c_int), value :: lwork
      !        type(c_ptr), value :: lwork
    end function hipsolverDgesvd_bufferSize
  end interface

  interface
    function hipsolverDgesvd (handle,jobu,jobv,m,n,A,lda,S,U,ldu,V,ldv,work, &
        lwork,rwork,info) bind(c,name = 'hipsolverDgesvd')
      use iso_c_binding
      use hipsolver_enums
      implicit none
      integer(kind(HIPSOLVER_STATUS_SUCCESS)) :: hipsolverDgesvd
      type(c_ptr), value :: handle
      integer(c_signed_char), value :: jobu
      integer(c_signed_char), value :: jobv
      integer(c_int), value :: m
      integer(c_int), value :: n
      type(c_ptr) :: A
      integer(c_int), value :: lda
      type(c_ptr), value :: S
      type(c_ptr), value :: U
      integer(c_int), value :: ldu
      type(c_ptr), value :: V
      integer(c_int), value :: ldv
      type(c_ptr), value :: work
      integer(c_int), value :: lwork
      type(c_ptr), value :: rwork
      type(c_ptr), value :: info
    end function hipsolverDgesvd
  end interface

  interface
    function hipsolverDsyevd_bufferSize (handle,jobz,uplo,n,A,lda,W,lwork) &
        bind(c,name = 'hipsolverDsyevd_bufferSize')
      use iso_c_binding
      use hipsolver_enums
      implicit none
      integer(kind(HIPSOLVER_STATUS_SUCCESS)) :: hipsolverDsyevd_bufferSize
      type(c_ptr), value :: handle
      integer(kind(HIPSOLVER_EIG_MODE_NOVECTOR)), value :: jobz
      integer(kind(HIPSOLVER_FILL_MODE_LOWER)), value :: uplo
      integer(c_int), value :: n
      type(c_ptr), value :: A
      integer(c_int), value :: lda
      type(c_ptr), value :: W
      integer(c_int), value :: lwork
    end function hipsolverDsyevd_bufferSize
  end interface

  interface
    function hipsolverDsyevd(handle,jobz,uplo,n,A,lda,W,work,lwork,info) &
        bind(c,name = 'hipsolverDsyevd')
      use iso_c_binding
      use hipsolver_enums
      implicit none
      integer(kind(HIPSOLVER_STATUS_SUCCESS)) :: hipsolverDsyevd
      type(c_ptr), value :: handle
      integer(kind(HIPSOLVER_EIG_MODE_NOVECTOR)), value :: jobz
      integer(kind(HIPSOLVER_FILL_MODE_LOWER)), value :: uplo
      integer(c_int), value :: n
      type(c_ptr), value :: A
      integer(c_int), value :: lda
      type(c_ptr), value :: W
      type(c_ptr), value :: work
      integer(c_int), value :: lwork
      type(c_ptr), value :: info
    end function hipsolverDsyevd
  end interface

  interface
    function hipsolverDsyevj_bufferSize(handle,jobz,uplo,n,A,lda,W,lwork,params) &
        bind(c,name = 'hipsolverDsyevj_bufferSize')
      use iso_c_binding
      use hipsolver_enums
      implicit none
      integer(kind(HIPSOLVER_STATUS_SUCCESS)) :: hipsolverDsyevj_bufferSize
      type(c_ptr), value :: handle
      integer(kind(HIPSOLVER_EIG_MODE_NOVECTOR)), value :: jobz
      integer(kind(HIPSOLVER_FILL_MODE_LOWER)), value :: uplo
      integer(c_int), value :: n
      type(c_ptr), value :: A
      integer(c_int), value :: lda
      type(c_ptr), value :: W
      !        type(c_ptr), value :: lwork
      integer(c_int), value :: lwork
      type(c_ptr), value :: params
    end function hipsolverDsyevj_bufferSize
  end interface

  interface
    function hipsolverDsyevj(handle,jobz,uplo,n,A,lda,W,work,lwork,info,params) &
        bind(c,name = 'hipsolverDsyevj')
      use iso_c_binding
      use hipsolver_enums
      implicit none
      integer(kind(HIPSOLVER_STATUS_SUCCESS)) :: hipsolverDsyevj
      type(c_ptr), value :: handle
      integer(kind(HIPSOLVER_EIG_MODE_NOVECTOR)), value :: jobz
      integer(kind(HIPSOLVER_FILL_MODE_LOWER)), value :: uplo
      integer(c_int), value :: n
      type(c_ptr), value :: A
      integer(c_int), value :: lda
      type(c_ptr), value :: W
      type(c_ptr), value :: work
      integer(c_int), value :: lwork
      type(c_ptr), value :: info
      type(c_ptr), value :: params
    end function hipsolverDsyevj
  end interface

end module amd_hipsolver
#endif

#ifdef ROCSOLVER
module rocvars
  use iso_c_binding
  type(c_ptr)            :: rocsolver_handle
  integer (kind=4)       :: rocsolver_status
  integer (kind=4)       :: lwork1,lwork2,ndim,mdim
  real (kind=8)          :: tol,residual
  integer (kind=4)       :: max_sweeps,exec_sweeps
  integer (kind=4), parameter :: econ=0
  real (kind=8), allocatable, target :: workspace_d(:)
  integer (kind=4) :: dev_info_d
  !     integer(kind=cuda_stream_kind) :: stream
  integer (kind=4) :: esort,evect
  integer (kind=4) :: uplo=0
  integer (kind=4) :: rocinfo,workmode
  type(c_ptr)      :: workptr
end module rocvars

module amd_rocsolver
  use iso_c_binding
end module amd_rocsolver
#endif

#ifdef MKL
module mkl_solver

  integer (kind=8)    :: lwork1m,lwork2m,lworki,ndimm,mdimm
  real (kind=8),allocatable :: rwork(:)
  real (kind=8),allocatable :: workspace_d(:)
  integer (kind=8), allocatable :: workspace_i(:)

end module mkl_solver
#endif

module gnome_solvers
  enum,bind(c)
    enumerator SOLVER_EISPACK
    enumerator SOLVER_LAPACK
    enumerator SOLVER_MKL
    enumerator SOLVER_CUSOLVER
    enumerator SOLVER_CUSOLVERJ
    enumerator SOLVER_HIPSOLVER
    enumerator SOLVER_HIPSOLVERJ
    enumerator SOLVER_ROCSOLVER
    enumerator SOLVER_ROCSOLVERJ
    enumerator SOLVER_MAGMA
    enumerator SOLVER_SLATE
    enumerator SOLVER_CRAYLIBSCI
  end enum
end module gnome_solvers
