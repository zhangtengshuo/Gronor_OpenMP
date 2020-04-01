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
      character (len=255) :: filsym,filciv,filvec,filint
      character (len=255) :: filf74,filf75,filone,filtwo,fildat
      integer :: lfninp,lfnout,lfnsym,lfnciv,lfnvec,lfnint
      integer :: lfnpro,lfndbg,lfnf74,lfnf75,lfnone,lfntwo,lfndat
      integer :: lfnday

      character (len=255) :: user,host,date,time,cwd,command
      integer :: nprocs

      real (kind=8), allocatable :: civm(:,:), vecsm(:,:,:)
      integer, allocatable :: ioccm(:,:,:), nbasm(:), nactm(:)
      integer, allocatable :: inactm(:), idetm(:)
      real (kind=8), allocatable :: civb(:,:), vecsb(:,:,:)
      integer, allocatable :: ioccb(:,:)
      integer, allocatable ::  ncombv(:,:),inactb(:), nactb(:), idetb(:)
      real (kind=8), allocatable :: hbase(:,:), sbase(:,:)
      integer, allocatable :: nsing(:,:,:)
      integer :: iocch(32)

      end module cidef

      module gnome_parameters

      implicit none

      integer :: iocc(2),iord(2),indx(2),icalc,ipr,ipro,ipvec,idbg,itim
      integer :: iload,ibatch,itest,naccel,nacc0,nacc1,inpcib,intfil
      integer :: iday,idipole,itp4,ione,intadd,nummps,numomp,ixpert
      integer :: ins2
      integer :: ibase0,jbase0,idet0,jdet0
      integer :: ndeti,ndetj,nacti,nactj,inacti,inactj
      real (kind=8) :: thresh, thresh_CI, thresh_SIN

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
      real (kind=8), allocatable :: va(:,:),vb(:,:),tb(:,:)
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

      real (kind=8) :: buffer(7)

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

       type(cusolverDnHandle) :: cusolver_handle
       integer (kind=4)       :: cusolver_status
       integer (kind=4)       :: worksize,worksize2,ndim,mdim
       real(8),allocatable :: rwork(:)
       real(8),allocatable :: workspace_d(:)
       integer (kind=4) :: dev_info_d

       end module cuda_cusolver
#endif
