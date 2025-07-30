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

  integer (kind=4) :: me,np,mstr,rocinfo
  integer (kind=4) :: numdev,mydev
  integer (kind=4) :: group_batch,group_heads,group_world
  integer (kind=4), allocatable :: ranks_list(:),ranks_heads(:)

  integer :: nnodes,nrsets,nranks,ncycls,nrnsets,ngpus

  integer :: numacc,numnon,numgrp,numgrn,numidle
  integer :: ngr,meg,npg,comm_group,comm_first,mgr,nalive,nabort
  integer :: mynode,mygroup,myhead,ime,mpibuf
  integer :: iamhead,iamacc,iamactive,nonidle
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
  
  character (len=12) :: machine

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

  integer :: nmol,mstates,nbase,nspin,nbasenoct

  character (len=255) :: root,filinp,filout,filpro,fildbg,filday
  character (len=255) :: filsys,filciv,filvec,filint,filtst,fildet
  character (len=255) :: filone,filtwo,fildat,fillog,filtim,filtmp
  character (len=255) :: filcpr,filarx,filrnk,filcml,filxrx
  integer :: lfninp,lfnout,lfnsys,lfnciv,lfnvec,lfnint,lfndet
  integer :: lfnpro,lfndbg,lfnone,lfntwo,lfndat,lfntim,lfnabt
  integer :: lfnday,lfntst,lfnlog,lfncpr,lfnarx,lfnrnk,lfncml
  integer :: lfnwrn,lfnxrx,lfntmp

  character (len=255) :: user,host,date,time,cwd,command,git_commit
  character (len=255) :: lmodcomp,lmodcompv,lmodmpi,lmodmpiv
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
  real (kind=8), allocatable :: hev(:),hevnoct(:)
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

  integer :: icalc,ipr,ipro,ipvec,idbg,itim,itmp,ires,iint
  integer :: itest,ifault,sv_solver,ev_solver,idevel,idist,labmax
  integer :: ntask,ntaska,nbatch,nbatcha
  integer :: ndbg,mdbg,load,loada
  integer :: iaslvr,jaslvr,inslvr,jnslvr
  integer :: iswsvj,iswevj
  integer :: naccel,nacc0,nacc1,nidle,inpcib,intfil,ncols
  integer :: iday,idipole,itp4,nummps,numgpu,ixpert
  integer :: ins2
  integer :: ibase0,jbase0,idet0,jdet0
  integer :: ndeti,ndetj,nacti,nactj,inacti,inactj
  real (kind=8) :: tau_MO,tau_CI,tau_CI_off,tau_SIN,thresh,thresh_SIN
  real (kind=8) :: tolsvj,tolevj
  integer :: numtasks
  real (kind=8) :: tmax

  integer (kind=4) :: itreq
  integer (kind=8) :: irbuf(4)
  logical :: odbg,odupl,oterm,otreq,otimeout
  logical :: corres
  logical :: lsvcpu,levcpu,lsvtrns

  character (len=1) :: prec
  character (len=128) :: mebfroot,combas
  character (len=16), allocatable :: fragname(:)
  character (len=16), allocatable :: fragstate(:,:)
  character (len=128), allocatable :: fragfile(:)
  character (len=128), allocatable :: vecfile(:)
  character (len=128), allocatable :: detfile(:)

  character (len=128) :: asvd,nsvd,aevd,nevd

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

  ! ntesta/ntestb count electrons on each fragment (closed + open shells) and must match
  ! n1bas is the half-triangular dimension n_bas(n_bas+1)/2
  ! nstdim and mbasel provide maximum array dimensions based on nelecs and nbas.
  ! ijend stores the total number of determinant pairs for a base pair, 
  ! it controls the task loops during parallel work distribution.
  integer :: ntesta,ntestb
  integer :: n1bas,nstdim,mbasel,ijend
  integer (kind=4) :: nelecs
  
  ! va/vb hold MO coefficients for the two fragments
  ! tb is the intermediate matrix S⋅VB used to build the overlap matrix ta
  ! After the MO overlap matrix ta is built, it is copied to a which serves as input to the SVD solver
  ! The SVD produces singular vectors u and wt (transposed right-hand matrix). The transpose is stored in w
  ! diag and sdiag are derived from specific columns of u and w when singular values fall below a threshold.
  ! cdiag/csdiag store the same quantities for later reference
  ! In gronor_tramat, the diagonal vectors and overlap matrices are transformed with va and vb
  ! w1 and w2 are scratch arrays for accumulating intermediate sums
  real (kind=8), allocatable :: va(:,:),vb(:,:),tb(:,:)   
  real (kind=8), allocatable :: ta(:,:),a(:,:)
  real (kind=8), allocatable :: u(:,:),w(:,:),wt(:,:),ev(:)

  real (kind=8), allocatable :: sdiag(:)
  real (kind=8), allocatable, target :: diag(:)
  real (kind=8), allocatable :: bsdiag(:),bdiag(:)
  real (kind=8), allocatable :: csdiag(:),cdiag(:)
  real (kind=8), allocatable :: w1(:),w2(:,:)

  ! veca/vecb collect correlated orbitals for debugging (gronor_cororb.F90)
  real (kind=8), allocatable :: veca(:),vecb(:)

  ! The matrices taa, aaa, aat, tt, and sm are subsequently combined to compute the two-electron contributions in gronor_gntwo
  real (kind=8), allocatable :: taa(:,:)
  real (kind=8), allocatable :: sm(:,:)
  real (kind=8), allocatable :: aaa(:,:),aat(:,:)
  real (kind=8), allocatable :: tt(:,:)


  integer :: ising

  real (kind=8), allocatable :: result(:,:),resultt(:,:)

!  real (kind=8), allocatable :: work(:)
  integer (kind=8) :: len_work_dbl,len_work2_dbl,len_work_int,info

  real (kind=8) :: buffer(17)
!$omp threadprivate(buffer,e2buff,e2summ)

  integer (kind=8) :: numdet,melen,memax,icur,jcur
  integer (kind=4), allocatable :: melist(:,:)
!$omp threadprivate(icur,jcur,melist)
  integer (kind=8), allocatable :: ndxdet(:,:)

  real (kind=8) :: gbmelist

  real (kind=8), allocatable :: c2sum(:,:)
  integer (kind=8), allocatable :: nbdet(:,:)
  real (kind=8), allocatable :: rwork(:)

!$omp threadprivate(a,ta,tb,va,vb,w1,w2,taa,u,w,wt,ev,diag,bdiag,cdiag,bsdiag,csdiag,sdiag,aaa,tt,aat,sm,rwork,nelecs,len_work_dbl,len_work2_dbl,len_work_int,memax)

end module gnome_data

module gnome_integrals

  integer :: intone,int1,int2,mint2,nt(3),myints(7),intndx,jntndx
  integer :: igfr,igto
  real (kind=8), allocatable :: s(:,:)    ! 基函数间的重叠矩阵
  real (kind=8), allocatable :: t(:)      ! 单电子动能积分
  real (kind=8), allocatable :: v(:)      ! 核吸引积分
  real (kind=8), allocatable :: dqm(:,:)  ! 偶极矩(前几列)和四极矩(其余列)积分

#ifdef SINGLEP
  real (kind=4), allocatable :: g(:)!,gg(:)
#else
  real (kind=8), allocatable :: g(:)!,gg(:) ! 双电子排斥积分
#endif

  integer (kind=2), allocatable :: lab(:,:) !每个超索引对应的一对轨道索引
  integer (kind=8), allocatable :: ndx(:)   !索引辅助数组，用于将轨道对转换为g和单电子积分中的位置
  integer (kind=8), allocatable :: ndxk(:)
  integer (kind=4), allocatable :: ndxtv(:)
  integer, allocatable :: ig(:)
  integer :: mbuf,mclab
  integer (kind=4) :: int_comm
  integer (kind=4), allocatable :: new_comm(:)
  integer (kind=4), allocatable :: new_comm1(:)
  integer (kind=4), allocatable :: new_comm2(:,:)
  integer (kind=4) :: int_group,new_group,grp_group
  integer, allocatable :: ntarget(:)

  integer (kind=8) :: mlab,icomm1

  logical :: lcomm1

end module gnome_integrals

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

#ifdef MKL
module mkl_solver

!  integer (kind=8)    :: lwork1m,lwork2m,lworki,ndimm,mdimm
!  real (kind=8),allocatable :: rwork(:)
!  real (kind=8),allocatable :: workspace_d(:)
!  integer (kind=8), allocatable :: workspace_i(:)

end module mkl_solver
#endif

module gnome_solvers
  enum,bind(c)
    enumerator SOLVER_EISPACK
    enumerator SOLVER_MKL
    enumerator SOLVER_MKLD
    enumerator SOLVER_MKLJ
  end enum

  integer (kind=8) :: lwork,liwork,ndimm,mdimm,ndim8,lwork8,liwork8
  integer (kind=4) :: mdim,ndim,ndim4,lwork4,liwork4

  real(kind=8), allocatable :: work(:)
  integer(kind=8), allocatable :: iwork(:)
  real (kind=8),allocatable :: workspace_d(:)
  real (kind=8),allocatable :: workspace2_d(:)
  integer (kind=8), allocatable :: workspace_i(:)
  integer (kind=4), allocatable :: workspace_i4(:)
!$omp threadprivate(workspace_d,workspace2_d,workspace_i,workspace_i4)
!  character*1 :: jobz,uplo
  integer (kind=4) :: jobz,uplo
end module gnome_solvers
