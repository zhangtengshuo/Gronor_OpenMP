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
!! Main driver for the GronOR application
!! @author  Tjerk P. Straatsma, ORNL
!! @author  Coen de Graaf, URV
!! @date    2016
!>
!> @todo     Maybe we could take out some parts and put them in
!>           separate subroutines. A subroutine for reading the
!>           input would be nice (with some extra checks on
!>           input errors)

#include "gronor_config.fh"

subroutine gronor_main()

  !>    Initialization of the calculation
  !!
  !!     The root name is read from the from the command line,
  !!     and is used to generate all other filenames associated with the calculation.
  !!
  !!     <table>
  !!     <caption id="multi_row">Files</caption>
  !!     <tr><th>File name <th> File Description <th> String <th> Unit <th> Source
  !!     <tr><td>\$(root).inp <td> input file <td> filinp <td> lfninp <td>
  !!     <tr><td>\$(root).out <td> output file <td> filout <td> lfnout <td>
  !!     <tr><td>\$(root).sys <td> system file <td> filsys <td> lfnsys <td>
  !!     <tr><td>\$(root)_nnn.civ <td> CI vectors, one for every nnn=1,mstates
  !!                 <td> filciv <td> lfnciv <td> OpenMolcas Auxiliary
  !!     <tr><td>\$(root)_nnn.vec <td> MO vectors, one for every nnn=1,mstates
  !!                 <td> filvec <td> lfnvec
  !!                 <td> OpenMolcas Auxiliary
  !!     <tr><td>\$(root)_nnn.det <td> CI vectors and MO vectors, one for every nnn=1,mstates
  !!                 These files are a replacement of civ and vec files <td> filvec <td> lfnvec
  !!                 <td> OpenMolcas
  !!     </table>
  !!

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
#ifdef _OPENMP
  use omp_lib
#endif
  !      use iso_c_binding, only : c_loc, c_ptr

#ifdef _OPENACC
  use openacc
  !     use cuda_functions
#endif
  
  implicit none

  external :: gronor_abort
  external :: swatch,timer_init,timer_start,timer_stop
  external :: gronor_timings,gronor_finalize_cml
  external :: gronor_print_dipole_moments,gronor_print_results
  external :: gronor_correlation_energy,gronor_multipoles_nuclear
  external :: gronor_prtmat,gronor_print_matrix
  external :: gronor_worker,gronor_memory_usage
  external :: gronor_master,gronor_read_integrals
  external :: gronor_make_basestate,gronor_assign_managers
  external :: gronor_solver_create_handle
  external :: gronor_results_header_cml,gronor_init_cml
  external :: gronor_env_cml,gronor_gnome_molcas_input
  external :: gronor_read_vectors_and_determinants
  external :: gronor_number_integrals
  external :: gronor_input,gronor_prelude_cml
  external :: getcwd,getlog,hostnm

  external :: MPI_Bcast

  integer (kind=4) :: ierror,ierr,ncount
  integer (kind=8) :: iarg,i,j,jp,idum(55),k,l,iact
  integer :: node
  real (kind=8) :: rdum(6)
  character (len=255) :: string,target,compiler
  logical exist,first_pass

  real(kind=8), external :: timer_wall_total

  integer :: getcpucount
  external :: getcpucount

  integer :: igr,numone,numtwo,maxgrp
  integer (kind=4) :: new

  integer :: ibase,jbase,lnxt,lcur,ksr,nsr(4)

  integer                  :: idet,ndet_rev
  integer, allocatable     :: occu_aux(:)
  integer (kind=1), allocatable :: iocc_aux(:,:)
  real(kind=8),allocatable :: civb_aux(:)
  real(kind=8)             :: thres2
  integer, allocatable     :: nbuf(:,:)

  real (kind=8) :: rint,rndx,rlst,rh,rs,rt,rm(9),c2s

  integer :: ndtot,npl,nrg,igb,ibd,ivc,ibas,ib,lc,nc,mc,nnc

  integer :: l2,nidet,njdet,n

  character (len=3) :: ident(3)
  character (len=4) :: onlabel
  character (len=5) :: version
  character (len=64) :: version_type

  integer :: major,minor

#ifdef USE_POSIXF
  integer*4 len4,ierr4
#endif

#ifdef _OPENACC
  type(c_ptr) :: cpfre, cptot
#endif

  call timer_init()

  call timer_start(99)
  call timer_start(1)
  call timer_start(2)

  bias=0.0d0
  deta=0.0d0

  if(me.eq.mstr) then

    user='                  '
    host='                  '
    date='                  '
    time='                  '
    cwd='                   '

#ifdef USE_POSIXF
    call pxfgetlogin(user,len4,ierr4)
#else
#ifdef IBM
    call hostnm_(host)
    call getenv("USER",user)
#else
#ifdef CRAY
    call hostnm(host)
    call getenv("USER",user)
#else
    call hostnm(host)
    call getlog(user)
#endif
#endif
#endif

    call swatch(date,time)
    call getcwd(cwd)

    !     Read a single string argument 'root' from the command line, that
    !     will be used to generate all other files used to read and write
    !     data associated with this calculation

    !     input file: filinp with logical file number lfninp
    !     output file: filout with logical file number lfnout
    !     system file: filsys with logical file number lfnsys
    !     what is on?: filciv with logical file number lfnciv
    !     what is on?: filvec with logical file number lfnvec

    iarg=0
    call getarg(iarg,command)
    iarg=1
    call getarg(iarg,string)

    root=trim(string)
    if(index(string,'_').gt.0) root=string(1:index(string,'_')-1)
    filinp=trim(string)//'.inp'
    lfninp=5
    filout=trim(string)//'.out'
    lfnout=16
    filsys=trim(root)//'.sys'
    lfnsys=7
    lfnciv=9
    lfnvec=10
    filint=trim(root)//'.int'
    lfnint=11
    lfnone=11
    lfntwo=12
    lfndbg=13
    filpro=trim(string)//'.pro'
    lfnpro=14
    fildat=trim(string)//'.dat'
    lfndat=15
    lfndet=17
    filday=trim(string)//'.day'
    lfnday=18
    filtst=trim(string)//'.tst'
    lfntst=19
    fillog=trim(root)//'.log'
    lfnlog=20
    filcpr=trim(string)//'.cpr'
    lfncpr=21
    filarx=trim(string)//'.arx'
    lfnarx=22
    filrnk=trim(string)//'.rnk'
    lfnrnk=23
    filcml=trim(string)//'.cml'
    lfncml=24
    filtim=trim(string)//'.tim'
    lfntim=25
    lfnabt=26
    lfnwrn=27
    filxrx=trim(string)//'.xrx'
    lfnxrx=28

    open(unit=lfnout,file=trim(filout),form='formatted',status='replace',err=996)

    open(unit=lfnpro,file=trim(filpro),form='formatted',status='unknown',err=996)

    open(unit=lfnday,file=trim(filday),form='formatted',status='replace',err=996)

    open(unit=lfntst,file=trim(filtst),form='formatted',status='unknown',err=996)

    open(unit=lfncpr,file=trim(filcpr),form='unformatted',status='unknown',err=996)

    open(unit=lfnarx,file=trim(filarx),form='formatted',status='unknown',err=996)

    open(unit=lfnxrx,file=trim(filxrx),form='formatted',status='unknown',position='append',err=996)

    open(unit=lfnrnk,file=trim(filrnk),form='formatted',status='replace',err=996)

    open(unit=lfncml,file=trim(filcml),form='formatted',status='replace',err=996)

    call gronor_prelude_cml()

    call timer_stop(99)
    call swatch(date,time)
    write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),'  :  Start of job'
702 format(a8,2x,a8,f12.3,a)
    flush(lfnday)
    call timer_start(99)

!>
!!    <table>
!!    <caption id="multi_row">Input derived variables</caption>
!!    <tr><th> Variable <th> Description
!!    <tr><td> npg <td> number of processes per group (default is 1)
!!    <tr><td> nmol <td> number of molecules
!!    <tr><td> mstates <td> number of molecular states
!!    <tr><td> nbase <td> number of many-electron base (?)
!!    </table>
!!
!!    <table>
!!    <caption id="multi_row">Allocated arrays</caption>
!!    <tr><th> Variable <th> Description
!!    <tr><td> ncombv(nmol,nbase) <td> the combination
!!    <tr><td> nbasm(mstates) <td> number of one-electron basis functions per monomer/state
!!    <tr><td> inactm(mstates) <td> number of inactive orbitals per monomer/state
!!    <tr><td> idetm(mstates) <td> number of determinants per monomer/state
!!    </table>
    
    !
    !     Set defaults
    
    npg=0

    !     Default is one rank per worker group

    mgr=1
    if(npg.eq.0) npg=(np-1)/mgr
    if(np.lt.mgr+1) then
      write(*,'(a,i5,a,i5,a)') 'Error: Size specified (',mgr, &
          ') cannot be supported by available ranks (',np,')'
      call gronor_abort(200,"Error in main")
    endif
    
    nmol=0

    !     mstates will be set depending on Molecule input parameters

    mstates=0

    !     Set ntask and nbatch to undefined

    ntask=-1
    nbatch=-1
    ntaska=-1
    nbatcha=-1
    load=1
    loada=1

    !     Default is not to write a test file

    itest=0

    !     Default label length printed

    labmax=12

    !     Distribution of integrals
    !     0 : MPI_Bcast to all ranks (default)
    !     1 : MPI_Bcast to all nodes, followed by MPI_Bcast to ranks within nodes
    !     2 : MPI_Bcast to all nodes, followed by MPI_iSend/MPI_Recv to ranks within nodes

    idist=0
    !        if(nnodes.ge.128.and.nnodes.le.2048.and.np/nnodes.gt.mgr)
    ! &       idist=2

    idevel=0
    ifault=1

    !     Set solvers undefined

    isolver=-1
    jsolver=-1
    iaslvr=-1
    jaslvr=-1
    inslvr=-1
    jnslvr=-1

    naccel=-1
    nidle=0
    ncols=7
    numgpu=0
    nummps=1
    ixpert=0
    iswsvj=15
    iswevj=15
    nabort=64

#ifdef _OPENACC
    naccel=0
#endif

    !     Default print flag set to medium (20)

    ipr=20

    !     Default is not to turn timing on

    itim=0

    ipro=0
    iday=10
    idbg=0
    nspin=0
    ncorr=0
    mpibuf=168435456
    tau_CI=1.0d-5
    tau_SIN=1.0d-12
    tolsvj=1.0d-07
    tolevj=1.0d-07
    lcpr=.false.

    ngpus=0
    if(numdev.ge.1) ngpus=1
    if(machine.eq.'Juwels  '.and.numdev.eq.1) ngpus=4
    if(machine.eq.'Summit  '.and.numdev.eq.1) ngpus=6
    if(machine.eq.'Snellius'.and.numdev.eq.1) ngpus=4
    if(machine.eq.'Crusher '.and.numdev.eq.1) ngpus=8
    if(machine.eq.'Frontier'.and.numdev.eq.1) ngpus=8

    call gronor_input()

    write(filone,'(a,a,a)') trim(mebfroot),trim(combas),".one"
    
    if(mgr.eq.0) then
      if(npg.eq.0) call gronor_abort(201,"Error in main")
      mgr=(np-1)/npg
    endif

    if(nbatch.eq.-1) nbatch=32
    if(nbatcha.eq.-1) nbatcha=0

    if(ntaska.eq.-1) then
      if(nbatcha.gt.0) then
        ntaska=nbatcha
      else
        ntaska=32
      endif
    endif

    if(ntask.eq.-1) then
      if(nbatch.gt.0) then
        ntask=nbatch
      else
        ntask=32
      endif
    endif

    if(nbatcha.gt.ntaska) nbatcha=ntaska
    if(nbatch.gt.ntask) nbatch=ntask

    !     Sanity check on the task size

    if(numdev.gt.0.and.ntaska.eq.0) call gronor_abort(202,"Error in main")
    if(numdev.eq.0.and.ntask.eq.0) call gronor_abort(202,"Error in main")
    if(ntaska.eq.0.and.ntask.eq.0) call gronor_abort(202,"Error in main")

    !     Unless the user is an expert let's check thresholds

    if(ixpert.eq.0) then
      if(tau_CI.lt.1.0d-6) tau_CI=1.0d-6
      if(tau_CI.gt.1.0d-4) tau_CI=1.0d-4
    endif

    major=_GRONOR_VERSION_MAJOR_
    minor=_GRONOR_VERSION_MINOR_

    if(minor.eq.0.or.minor.gt.12) then
      version_type=" under active development"
    else
      version_type=" official release"
    endif
    write(version,'(i2,a,i2.2)') major,'.',minor

    if(ipr.gt.0) write(lfnout,600) major,minor,trim(version_type)
600 format(/, &
        ' GronOR: Non-Orthogonal Configuration Interaction',//,' Version ',i2,'.',i2.2,a,//, &
        ' T. P. Straatsma',t50,'C. de Graaf',t82,'R. Broer',/, &
        t50,'A. Sanchez',t82,'R. K. Kathir',//, &
        ' National Center for Computational Sciences',t50, &
        'Quantum Chemistry Group',t82, &
        'Department of Theoretical Chemistry',/, &
        ' Oak Ridge National Laboratory',t50, &
        'University Rovira i Virgili',t82, &
        'University of Groningen',/,' Oak Ridge, Tennessee, USA',t50, &
        'Tarragona, Spain',t82, &
        'Groningen, the Netherlands',//, &
        ' Based on GNOME written by R. Broer-Braam, ', &
        'J. Th. van Montfort, and B. Vunderink',//, &
        ' Please cite the following reference when publishing ', &
        'results obtained using GronOR:',//, &
        ' T. P. Straatsma, R. Broer, A. Sanchez-Mansilla,', &
        'C. Sousa, and C.de Graaf',/, &
        ' "GronOR: Scalable and Accelerated Nonorthogonal ', &
        'Configuration Interaction for Molecular Wave ', &
        'Functions"',/, &
        ' Journal of Chemical Theory and Computation, 18, 3549-', &
        '3565 (2022) https://doi.org/10.1021/acs.jctc.2c00266',/)
    
    target=' CPU'
    compiler=' '
#ifdef ACC
#ifdef GPUAMD
    write(target,'(a)') " OpenACC AMD"
#else
    write(target,'(a)') " OpenACC NVIDIA"
#endif
#endif
#ifdef OMPTGT
#ifdef GPUAMD
    write(target,'(a)') " OpenMP target AMD"
#else
    write(target,'(a)') " OpenMP target NVIDIA"
#endif
#endif
    if(ipr.ge.0) write(lfnout,645) trim(target)
#if defined(CRAY)
645 format(//,' Compiled with the CRAY CCE Compiler suite',a)
    compiler='Cray CCE'
#elif defined(IBM)
645 format(//,' Compiled with the IBM XL Compiler suite',a)
    compiler='IBM XL'
#elif defined(PGI)
645 format(//,' Compiled with the PGI Compiler suite',a)
    compiler='PGI'
#elif defined(NVIDIA)
645 format(//,' Compiled with the NVIDIA Compiler suite',a)
    compiler='NVIDIA'
#elif defined(NVHPC)
645 format(//,' Compiled with the NVIDIA NVHPC Compiler suite',a)
    compiler='NVHPC'
#elif defined(GNU)
645 format(//,' Compiled with the GNU Compiler suite',a)
    compiler='GNU'
#elif defined(INTEL)
645 format(//,' Compiled with the INTEL Compiler suite',a)
    compiler='Intel'
#elif defined(FLANG)
645 format(//,' Compiled with the FLANG/CLANG Compiler suite',a)
    compiler='FLANG'
#elif defined(AMD)
645 format(//,' Compiled with the AMD Compiler suite',a)
    compiler='AMD'
#else
645 format(//,' Compiled with unknown compiler suite',a)
#endif
    onlabel='    '
    if(machine.ne.'        ') onlabel=' on '

    if(ipr.ge.20) write(lfnout,601) trim(user),getcpucount(), &
        trim(host),onlabel,trim(machine),nnodes, &
        date(1:8),nrsets, &
        time(1:8),nrnsets,numdev,ngpus,nummps, &
        np,np/nnodes,np/nrsets, &
        ncycls,num_threads
601 format(//, &
        ' User',t30,a,t60,'CPU count',t100,i10,/,/, &
        ' Host',t30,a,a,a,t60,'Number of nodes',t100,i10,/, &
        ' Date',t30,a,t60, &
        'Number of resource sets',t100,i10,/, &
        ' Time',t30,a,t60, &
        'Number of resource sets per node',t100,i10,/, &
        t60,'Number of GPUs per resource set',t100,i10,/, &
        t60,'Number of GPUs per node',t100,i10,/, &
        t60,'Number of MPI ranks per GPU (MPS)',t100,i10,/, &
        t60,'Number of MPI ranks',t100,i10,/, &
        t60,'Number of MPI ranks per node',t100,i10,/, &
        t60,'Number of MPI ranks per resource set',t100,i10,/, &
        t60,'Number of rank assignment cycles',t100,i10,/, &
        t60,'Number of OPENMP threads per MPI rank',t100,i10)
    if(ipr.ge.20.and.naccel.gt.0) write(lfnout,608) naccel
608 format(t60,'Number of accelerated ranks per node limit',t100,i10)
    if(ipr.ge.20.and.nidle.gt.0) write(lfnout,648) nidle
648 format(t60,'Number of idle ranks per node',t100,i10)
    if(ipr.ge.20.and.numdev.gt.0) write(lfnout,602) &
        dble(memfre)/1073.74d06,dble(memtot)/1073.74d06
602 format(t60,'Available memory on device',t90,f24.3,' GB',/, &
        t60,'Total memory on device',t90,f24.3,' GB')
    if(ipr.gt.0) write(lfnout,603) trim(command),trim(cwd), &
        trim(filinp),trim(filsys), &
        trim(filout),trim(filone), &
        trim(mebfroot),trim(combas), &
        trim(mebfroot),trim(combas)
603 format(/,' Command argument',t30,a,/, &
        ' Current working directory',t30,a,//, &
        ' Input file is',t25,a,t60, &
        'System file is',t92,a,/, &
        ' Output file is',t25,a,t60, &
        'One electron integral file is',t92,a,//, &
        ' CI vector file(s) are',t25,a,a,'_lbl.det',/ &
        ' MO vector file(s) are',t25,a,a,'_lbl.vec')

    call swatch(date,time)
    write(lfnarx,401) trim(user)
    write(lfnxrx,401) trim(user)
401 format('User ',a)
    write(lfnarx,402) trim(date)
    write(lfnxrx,402) trim(date)
402 format('Date ',a)
    write(lfnarx,403) trim(time)
    write(lfnxrx,403) trim(time)
403 format('Time ',a)
    write(lfnarx,404) trim(host),onlabel,trim(machine)
    write(lfnxrx,404) trim(host),onlabel,trim(machine)
404 format('Host ',a,a,a)
    write(lfnarx,405) trim(root)
    write(lfnxrx,405) trim(root)
405 format('Root ',a)
    write(lfnarx,406) nnodes,np
    write(lfnxrx,406) nnodes,np
406 format('Nodes',i10,i10)

    call gronor_number_integrals(numone,numtwo)

    if(npg*mgr+1.gt.np) then
      call gronor_abort(203,"Requested ranks exceeds available")
    endif

    allocate(nbasm(mstates))
    allocate(nactm(mstates))
    allocate(inactm(mstates))
    allocate(idetm(mstates))
    
    icalc=0
    ins2=0
    ipvec=0
    idipole=0
    itp4=0
    bias=0.0d0
    corres=.false.

    if(ipr.ge.30) then
      write(lfnout,604) nmol,ifault,mstates,itest,nbase, &
          idevel,idist
604   format(/,' Number of fragments',t40,i10, &
          t60,'Fault tolerance wait time',t100,i10,' sec',/ &
          ' Number of fragment wave functions',t40,i10, &
          t60,'Test option',t100,i10,/, &
          ' Number of many-electron base states',t40,i10, &
          t60,'Development option',t100,i10,/, &
          t60,'Distribution option',t100,i10)
    else
      write(lfnout,605) nmol,mstates,nbase
605   format(/,' Number of fragments',t40,i10,/, &
          ' Number of fragment wave functions',t40,i10,/, &
          ' Number of many-electron base states',t40,i10,/)
    endif
    if(ipr.ge.20) then
      write(lfnout,606) mgr,ntaska,ntask, &
          max(1,nbatcha),max(1,nbatch),loada,load
606   format(' Number of ranks per task',t40,i10,//, &
          ' Task size',t40,i10,i6,/, &
          ' Batch size',t40,i10,i6,/, &
          ' Load balancing',t40,i10,i6)
    endif

    call gronor_read_vectors_and_determinants()

    if(ipr.gt.0) write(lfnout,607) tau_MO,tau_CI,tau_SIN
607 format(/,' Threshold common molecular orbital transformation ',t55,1pe10.3, &
        /,' Threshold product MEBF expansion coefficients ',t55,e10.3, &
        /,' Threshold determining singularities ',t55,e10.3)

    if(ipr.gt.0) write(lfnout,609) nspin+1
609 format(/,' Spin multiplicity (2S+1)',t40,i10)

    call gronor_gnome_molcas_input() ! OpenMolcas
    close(unit=lfnsys)

    call gronor_env_cml(version,version_type)
    call gronor_init_cml
    !     this block is filled in gronor_print_results
    call gronor_results_header_cml

  endif

  if(numdev.eq.0) then
    if(inslvr.lt.0) inslvr=iaslvr
    if(jnslvr.lt.0) jnslvr=jaslvr
  endif
  
  
  
  !     Distribute input data to all processes
  
  if(np.gt.1) then

    if(me.eq.mstr) then
      idum(1)=nmol
      idum(2)=mstates
      idum(3)=nbase
      idum(4)=maxci
      idum(5)=maxvec
      idum(6)=maxnact
      idum(7)=nidle
      idum(8)=0
      idum(9)=0
      idum(10)=0
      idum(11)=0
      idum(12)=0
      idum(13)=icalc
      idum(14)=ipr
      idum(15)=ins2
      idum(16)=ipvec
      idum(17)=itp4
      idum(18)=numfiles
      idum(19)=0
      if(corres) idum(19)=1
      idum(20)=npg
      idum(21)=nnucl
      idum(22)=nbas
      idum(23)=mclab
      idum(24)=ntask
      idum(25)=ntaska
      idum(26)=mgr
      idum(27)=itest
      idum(28)=naccel
      idum(29)=inpcib
      idum(30)=ipro
      idum(31)=int2
      idum(32)=mpibuf
      idum(33)=idbg
      idum(34)=nummps
      idum(35)=labmax
      idum(36)=itim
      idum(37)=iaslvr
      idum(38)=jaslvr
      idum(39)=inslvr
      idum(40)=jnslvr
      idum(41)=ifault
      idum(42)=idevel
      idum(43)=iswsvj
      idum(44)=iswevj
      idum(45)=nbatch
      idum(46)=nbatcha
      idum(47)=nspin
      idum(48)=managers
      idum(49)=mbuf
      idum(50)=idist
      idum(51)=load
      idum(52)=loada
      idum(53)=ngpus
      idum(54)=num_threads
      idum(55)=nabort

      rdum(1)=tau_CI
      rdum(2)=tau_CI
      rdum(3)=tau_SIN
      rdum(4)=tolsvj
      rdum(5)=tolevj
      rdum(6)=potnuc
    endif

    ncount=55
    call MPI_Bcast(idum,ncount,MPI_INTEGER8,mstr,MPI_COMM_WORLD,ierr)
    ncount=6
    call MPI_Bcast(rdum,ncount,MPI_REAL8,mstr,MPI_COMM_WORLD,ierr)
    ncount=255
    call MPI_Bcast(root,ncount,MPI_CHAR,mstr,MPI_COMM_WORLD,ierr)
    ncount=255
    call MPI_Bcast(mebfroot,ncount,MPI_CHAR,mstr,MPI_COMM_WORLD,ierr)
    ncount=255
    call MPI_Bcast(combas,ncount,MPI_CHAR,mstr,MPI_COMM_WORLD,ierr)

    if(me.ne.mstr) then
      nmol=idum(1)
      mstates=idum(2)
      nbase=idum(3)
      maxci=idum(4)
      maxvec=idum(5)
      maxnact=idum(6)
      nidle=idum(7)
      icalc=idum(13)
      ipr=idum(14)
      ins2=idum(15)
      ipvec=idum(16)
      itp4=idum(17)
      numfiles=idum(18)
      corres=idum(19).eq.1
      npg=idum(20)
      nnucl=idum(21)
      nbas=idum(22)
      mclab=idum(23)
      ntask=idum(24)
      ntaska=idum(25)
      mgr=idum(26)
      itest=idum(27)
      naccel=idum(28)
      inpcib=idum(29)
      ipro=idum(30)
      int2=idum(31)
      mpibuf=idum(32)
      idbg=idum(33)
      nummps=idum(34)
      labmax=idum(35)
      itim=idum(36)
      iaslvr=idum(37)
      jaslvr=idum(38)
      inslvr=idum(39)
      jnslvr=idum(40)
      ifault=idum(41)
      idevel=idum(42)
      iswsvj=idum(43)
      iswevj=idum(44)
      nbatch=idum(45)
      nbatcha=idum(46)
      nspin=idum(47)
      managers=idum(48)
      mbuf=idum(49)
      idist=idum(50)
      load=idum(51)
      loada=idum(52)
      ngpus=idum(53)
      num_threads=idum(54)
      nabort=idum(55)

    endif

    if(managers.gt.0) call gronor_assign_managers()

    int1=(nbas*(nbas+1))/2

    if(me.eq.mstr) then
      call timer_stop(99)
      call swatch(date,time)
      write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),'  :  Input broadcasted'
      flush(lfnday)
      call timer_start(99)
    endif
    
    tau_CI=rdum(1)
    tau_CI=rdum(2)
    tau_SIN=rdum(3)
    tolsvj=rdum(4)
    tolevj=rdum(5)
    potnuc=rdum(6)

    if(numdev.gt.0) nbatch=nbatcha

  endif

  lfnone=11
  lfntwo=12
  lfndbg=13
  lfnabt=26
  lfnwrn=27
  if(idbg.gt.0) then
    write(fildbg,1300) me
1300 format('GronOR_',i5.5,'.dbg ')
    open(unit=lfndbg,file=trim(fildbg),form='formatted',status='unknown',err=996)
  endif

  nacc0=naccel
  nacc1=naccel

  if(me.ne.mstr) then
    nmol=idum(1)
    mstates=idum(2)
    nbase=idum(3)
    maxci=idum(4)
    maxvec=idum(5)
    maxnact=idum(6)

    allocate(ncombv(nmol,nbase))
    allocate(nbasm(mstates))
    allocate(nactm(mstates))
    allocate(inactm(mstates))
    allocate(idetm(mstates))
    allocate(spinm(mstates))
    allocate(civm(maxci,mstates))
    allocate(occm_string(maxci,mstates))
    allocate(vecsm(maxvec,maxvec,mstates))
    allocate(ioccm(maxnact,maxci,mstates))
    if(nmol.ge.3) allocate(inter_couplings(nmol-1,nbase))
    icalc=idum(13)
    ipr=idum(14)
    ins2=idum(15)
    ipvec=idum(16)
    itp4=idum(17)
    corres=idum(19).eq.1
  endif

  allocate(nbuf(mstates,5))
  do i=1,mstates
    nbuf(i,1)=nbasm(i)
    nbuf(i,2)=nactm(i)
    nbuf(i,3)=inactm(i)
    nbuf(i,4)=idetm(i)
    nbuf(i,5)=spinm(i)
  enddo
  if(np.gt.0) then
    ncount=maxci*mstates
    call MPI_Bcast(civm,ncount,MPI_REAL8,mstr,MPI_COMM_WORLD,ierr)
    ncount=maxvec*maxvec*mstates
    call MPI_Bcast(vecsm,ncount,MPI_REAL8,mstr,MPI_COMM_WORLD,ierr)
    ncount=maxnact*maxci*mstates
    call MPI_Bcast(ioccm,ncount,MPI_INTEGER8,mstr,MPI_COMM_WORLD,ierr)
    ncount=nmol*nbase
    call MPI_Bcast(ncombv,ncount,MPI_INTEGER8,mstr,MPI_COMM_WORLD,ierr)
    ncount=5*mstates
    call MPI_Bcast(nbuf,ncount,MPI_INTEGER8,mstr,MPI_COMM_WORLD,ierr)
    if(nmol.ge.3)then
      ncount=(nmol-1)*nbase
      call MPI_Bcast(inter_couplings,ncount,MPI_INTEGER8,mstr,MPI_COMM_WORLD,ierr)
    endif
  endif
  do i=1,mstates
    nbasm(i)=nbuf(i,1)
    nactm(i)=nbuf(i,2)
    inactm(i)=nbuf(i,3)
    idetm(i)=nbuf(i,4)
    spinm(i)=nbuf(i,5)
  enddo
  deallocate(nbuf)

  nbasis=0
  mactb=0
  do i=1,nmol
    nbasis=nbasis+nbasm(ncombv(i,1))
    mactb=mactb+nactm(ncombv(i,1))
  enddo

  call timer_stop(2)
  call timer_start(3)

  allocate(vecsb(nbasis,nbasis,nbase))
  allocate(inactb(nbase))
  allocate(nactb(nbase))
  allocate(idetb(nbase))
  inactb = 0
  nactb = 0

  ngr=(np-1)/mgr

  !     On the JFZ Juwels Booster the slurm implementation sets the
  !     device on each ranks based on best affinity. As a result
  !     numdev from acc_get_num_devices returns 1 as each rank
  !     is already associated with the clostest gpu. To get the
  !     correct map2 entries, the ngpus is set to 4 and nummps
  !     multiplied by ngpus.

  numgpu=nummps*ngpus

  allocate(ranks_heads(ngr+1))

  do i=1,np
    if(map2(i,1).gt.0) then
      map2(i,5)=map2(i,1)
      map2(i,2)=1
    else
      map2(i,5)=-map2(i,2)
    endif
  enddo

  map2(mstr+1,5)=0
  map2(mstr+1,2)=1

  !     Limit the number of accelerated ranks per node to naccel

  if(naccel.gt.0) then
    node=-1
    do i=1,np
      node=max(node,map2(i,6))
    enddo
    do i=1,node+1
      k=0
      do j=1,np
        if(map2(j,6).eq.i-1.and.map2(j,5).gt.0) then
          k=k+1
          if(k.gt.naccel) map2(j,5)=-map2(j,2)
        endif
      enddo
    enddo
  endif

  if(nidle.gt.0) then
    node=-1
    do i=1,np
      node=max(node,map2(i,6))
    enddo
    do i=1,node+1
      k=0
      do j=np,1,-1
        if(map2(j,6).eq.i-1.and.map2(j,5).ne.0) then
          k=k+1
          if(k.le.nidle) map2(j,5)=0
        endif
      enddo
    enddo
  endif

  role=idle
  if(map2(me+1,5).ne.0) role=worker
  if(me.eq.mstr) role=master

  numacc=0
  numnon=0
  do i=1,np
    if(map2(i,5).lt.0) map2(i,2)=num_threads
    if(map2(i,5).ge.0) map2(i,2)=1
    if(map2(i,5).gt.0) numacc=numacc+1
    if(map2(i,5).lt.0) numnon=numnon+1
  enddo
  numacc=numacc/mgr
  numnon=numnon/mgr
  maxgrp=numacc+numnon

#ifdef _OPENMP
  call omp_set_num_threads(int(map2(me+1,2),kind=4))
#endif
  
  allocate(thisgroup(mgr+1))
  allocate(allgroups(maxgrp+1,mgr+1))
  allocate(allheads(maxgrp+1))

  ! Array allgroups defines the ranks for each group
  
  ! For each group with index igrp and size igr=1,mgr
  !   allgroups(igrp,1)     : >0: accelerated; <0: CPU-only
  !   allgroups(igrp,igr+1) : rank id of member of group igrp
  
  numgrp=0
  igr=0
  do i=1,np
    if(map2(i,5).gt.0) then
!    if(map2(i,5).gt.0.and.map2(i,8).eq.worker) then
      if(igr.eq.0) numgrp=numgrp+1
      igr=igr+1
      allgroups(numgrp,1)=map2(i,5)
      allgroups(numgrp,igr+1)=i-1
      if(igr.eq.mgr) igr=0
    endif
  enddo
  if(igr.ne.0.and.numgrp.gt.0) numgrp=numgrp-1
  
  igr=0
  do i=1,np
    if(map2(i,5).lt.0) then
!    if(map2(i,5).lt.0.and.map2(i,8).eq.worker) then
      if(igr.eq.0) numgrp=numgrp+1
      igr=igr+1
      allgroups(numgrp,1)=map2(i,5)
      allgroups(numgrp,igr+1)=i-1
      if(igr.eq.mgr) igr=0
    endif
  enddo
  if(igr.ne.0.and.numgrp.gt.0) numgrp=numgrp-1
  
  do i=1,np
    map2(i,3)=0
  enddo
  do i=1,numgrp
    allheads(i)=allgroups(i,2)
    do j=1,mgr
      map2(allgroups(i,j+1)+1,3)=i
    enddo
  enddo
  allheads(numgrp+1)=mstr

  if(numgrp.eq.0) then
    write(*,'(a)') 'Number of groups is zero'
    call gronor_abort(201,"Error in main")
  endif
  do i=1,numgrp
    do j=1,mgr
      if(me.eq.allgroups(i,j+1)) then
        do k=1,mgr+1
          thisgroup(k)=allgroups(i,k)
        enddo
      endif
    enddo
  enddo

  if(me.ne.mstr) then
    iamhead=0
    if(me.eq.thisgroup(2)) iamhead=1
    myhead=thisgroup(2)
    numdev=thisgroup(1)
    if(numdev.lt.0) numdev=0
  else
    myhead=mstr
    numdev=-1
  endif

  mygroup=0
  do i=1,numgrp
    ranks_heads(i)=allgroups(i,2)
    do j=1,mgr
      if(me.eq.allgroups(i,j+1)) mygroup=i
    enddo
  enddo
  ranks_heads(numgrp+1)=mstr


  new=-1

  call MPI_Comm_Group(MPI_COMM_WORLD,group_world,ierr)
  
  call MPI_Group_incl(group_world,int(numgrp+1,kind=4),ranks_heads,group_heads,ierr)
  if(iamhead.eq.1) then
    call MPI_Group_Rank(group_heads,new,ierr)
  endif

  if(me.eq.mstr) then
    nrg=0
    do j=1,numgrp
      do i=2,mgr+1
        nrg=max(nrg,map2(allgroups(j,i)+1,4))
      enddo
    enddo
    k=5
    if(nrg.lt.10000) k=4
    if(nrg.lt.1000) k=3
    if(nrg.lt.100) k=2
    if(nrg.lt.10) k=1
    l=5
    if(np.lt.10000) l=4
    if(np.lt.1000) l=3
    if(np.lt.100) l=2
    if(np.lt.10) l=1
    npl=100/(k+l+3)
    write(lfnrnk,730)
730 format(//,' MPI Rank Assignments')
734 format(' ')
735 format(/,' Accelerated Groups on ',i8,' ranks',/)
733 format(/,' Non-Accelerated Groups on ',i8,' ranks',/)
    k=1
    l=min(k+npl,numacc)
    if(numacc.gt.0) then
      write(lfnrnk,735) numacc
      do while (k.le.l)
        if(mgr.gt.1) write(lfnrnk,734)
        do i=2,mgr+1
          if(nrg.lt.10) then
            if(np.lt.10) then
              write(lfnrnk,6611) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6612) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6613) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6614) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6615) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          elseif(nrg.lt.100) then
            if(np.lt.10) then
              write(lfnrnk,6621) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6622) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6623) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6624) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6625) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          elseif(nrg.lt.1000) then
            if(np.lt.10) then
              write(lfnrnk,6631) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6632) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6633) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6634) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6635) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          elseif(nrg.lt.10000) then
            if(np.lt.10) then
              write(lfnrnk,6641) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6642) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6643) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6644) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6645) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          else
            if(np.lt.10) then
              write(lfnrnk,6651) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6652) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6653) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6654) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6655) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          endif
        enddo
        k=l+1
        l=min(k+npl,numacc)
      enddo
      write(lfnrnk,734)
    endif
    flush(lfnrnk)
    if(numnon.gt.0) then
      write(lfnrnk,733) numnon
      k=l+1
      l=min(k+npl,numgrp)
      do while (k.le.l)
        if(mgr.gt.1) write(lfnrnk,734)
6611    format(5x,20(i1,':',i1,'  '))
6612    format(5x,20(i1,':',i2,'  '))
6613    format(5x,20(i1,':',i3,'  '))
6614    format(5x,20(i1,':',i4,'  '))
6615    format(5x,20(i1,':',i5,'  '))
6621    format(5x,20(i2,':',i1,'  '))
6622    format(5x,20(i2,':',i2,'  '))
6623    format(5x,20(i2,':',i3,'  '))
6624    format(5x,20(i2,':',i4,'  '))
6625    format(5x,20(i2,':',i5,'  '))
6631    format(5x,20(i3,':',i1,'  '))
6632    format(5x,20(i3,':',i2,'  '))
6633    format(5x,20(i3,':',i3,'  '))
6634    format(5x,20(i3,':',i4,'  '))
6635    format(5x,20(i3,':',i5,'  '))
6641    format(5x,20(i4,':',i1,'  '))
6642    format(5x,20(i4,':',i2,'  '))
6643    format(5x,20(i4,':',i3,'  '))
6644    format(5x,20(i4,':',i4,'  '))
6645    format(5x,20(i4,':',i5,'  '))
6651    format(5x,20(i5,':',i1,'  '))
6652    format(5x,20(i5,':',i2,'  '))
6653    format(5x,20(i5,':',i3,'  '))
6654    format(5x,20(i5,':',i4,'  '))
6655    format(5x,20(i5,':',i5,'  '))
        do i=2,mgr+1
          if(nrg.lt.10) then
            if(np.lt.10) then
              write(lfnrnk,6611) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6612) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6613) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6614) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6615) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          elseif(nrg.lt.100) then
            if(np.lt.10) then
              write(lfnrnk,6621) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6622) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6623) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6624) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6625) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          elseif(nrg.lt.1000) then
            if(np.lt.10) then
              write(lfnrnk,6631) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6632) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6633) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6634) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6635) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          elseif(nrg.lt.10000) then
            if(np.lt.10) then
              write(lfnrnk,6641) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6642) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6643) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6644) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6645) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          else
            if(np.lt.10) then
              write(lfnrnk,6651) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.100) then
              write(lfnrnk,6652) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.1000) then
              write(lfnrnk,6653) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            elseif(np.lt.10000) then
              write(lfnrnk,6654) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            else
              write(lfnrnk,6655) (map2(allgroups(j,i)+1,4),allgroups(j,i),j=k,l)
            endif
          endif
        enddo
        k=l+1
        l=min(k+npl,numgrp)
      enddo
    endif
    flush(lfnrnk)
  endif

  j=0
  jp=-1
  do i=1,np
    if(map2(i,4).ne.jp) j=0
    if(map2(i,5).gt.0.and.map2(i,1).gt.0) then
      map2(i,7)=mod(j,map2(i,1))
    else
      map2(i,7)=-1
    endif
    j=j+1
    jp=map2(i,4)
  enddo

  if(me.eq.mstr) then
    if(managers.eq.0) then
      write(lfnrnk,6666)
6666  format(//,' Rank map ND=NumDev DI=DevId Nt=NumThr Ac=Accel',//, &
          '    Rank  ND DI NT   Group    RSet Ac    Node  ', &
          '    Rank  ND DI NT   Group    RSet Ac    Node  ', &
          '    Rank  ND DI NT   Group    RSet Ac    Node',/)
      nc=3
      lc=np/nc
      i=np
      mc=mod(i,nc)
      if(mc.gt.0) lc=lc+1
      do i=1,lc
        if(mc.gt.0.and.i.eq.lc) nc=mc
        write(ident(1),'(a3)') '   '
        write(ident(2),'(a3)') '   '
        write(ident(3),'(a3)') '   '
        nnc=nc
        if(i+(nc-1)*lc.gt.np) nnc=nnc-1
        do j=1,nnc
          if(map2(i+(j-1)*lc,5).gt.0) write(ident(j),'(i3)') map2(i+(j-1)*lc,7)
        enddo
        write(lfnrnk,6667) (i+(j-1)*lc-1,map2(i+(j-1)*lc,1),ident(j), &
            (map2(i+(j-1)*lc,k),k=2,6),j=1,nnc)
6667    format(3(i8,':',i3,a3,i3,2i8,i3,i8,2x))
      enddo
    else
      write(lfnrnk,6668)
6668  format(//,' Rank map ND=NumDev DI=DevId Nt=NumThr Ac=Accel',//, &
          '    Rank  ND DI NT   Group    RSet Ac    Node    AcID Role Manager  ', &
          '    Rank  ND DI NT   Group    RSet Ac    Node    AcID Role Manager  ',/)
      nc=2
      lc=np/nc
      i=np
      mc=mod(i,nc)
      if(mc.gt.0) lc=lc+1
      do i=1,lc
        if(mc.gt.0.and.i.eq.lc) nc=mc
        write(ident(1),'(a3)') '   '
        write(ident(2),'(a3)') '   '
        write(ident(3),'(a3)') '   '
        nnc=nc
        if(i+(nc-1)*lc.gt.np) nnc=nnc-1
        do j=1,nnc
          if(map2(i+(j-1)*lc,5).gt.0) write(ident(j),'(i3)') map2(i+(j-1)*lc,7)
        enddo
        write(lfnrnk,6669) (i+(j-1)*lc-1,map2(i+(j-1)*lc,1),ident(j), &
            (map2(i+(j-1)*lc,k),k=2,9),j=1,nnc)
6669    format(2(i8,':',i3,a3,i3,2i8,i3,2i8,i5,i8,2x))
      enddo
    endif
    flush(lfnrnk)
  endif

  if(me.eq.mstr) then
    flush(lfnout)
    flush(lfnrnk)
  endif

  iamacc=0
  iamactive=0
  if(me.eq.mstr) iamactive=1
  if(me.ne.mstr.and.thisgroup(1).gt.0) then
    do i=1,mgr
      if(me.eq.thisgroup(i+1)) iamacc=1
    enddo
  endif
  if(thisgroup(1).ne.0) then
    do i=1,mgr
      if(me.eq.thisgroup(i+1)) iamactive=1
    enddo
  endif

  numidle=0
  do i=1,np
    if(map2(i,5).lt.0.and.(map2(i,3).eq.0.or.ntask.eq.0)) numidle=numidle+1
    if(map2(i,5).eq.0.and.i.ne.mstr+1) numidle=numidle+1
  enddo

  if(iamacc.eq.1) then

    ! Set ntask to the value for accelerated ranks

    ntask=ntaska

    !     only accelerated ranks need to set device

    !     On the JFZ Juwels Booster the slurm implementation sets the
    !     device on each ranks based on best affinity. As a result
    !     no acc_set_device_num should be executed here.

#ifdef _OPENACC
    if(numdev.gt.1) then
      mydev=map2(me+1,7)
      if(mydev.ge.0) then
#ifdef GPUAMD
        call acc_set_device_num(mydev,ACC_DEVICE_AMD)
#else
        call acc_set_device_num(mydev,ACC_DEVICE_NVIDIA)
#endif
        cpfre=c_loc(memfre)
        cptot=c_loc(memtot)
        !     istat=cudaMemGetInfo(cpfre,cptot)
        memavail=memfre
        if(idbg.gt.0) then
          call swatch(date,time)
          write(lfndbg,1301) date(1:8),time(1:8),' Device set to ',mydev,' of ',numdev
1301      format(a,1x,a,1x,a,i3,a,i3)
          flush(lfndbg)
        endif
      else
        iamacc=0
      endif
    endif
#endif

    call gronor_solver_create_handle()

  endif

#ifdef CUSOLVERJ
  if(iaslvr.lt.0) iaslvr=2
#endif
#ifdef CUSOLVER
  if(iaslvr.lt.0) iaslvr=1
  if(jaslvr.lt.0) jaslvr=1
#endif
#ifdef MKL
  if(inslvr.lt.0) inslvr=1
  if(jnslvr.lt.0) jnslvr=1
#endif
  if(inslvr.lt.0) inslvr=0
  if(jnslvr.lt.0) jnslvr=0
  if(iaslvr.lt.0) iaslvr=0
  if(jaslvr.lt.0) jaslvr=0

  if(iaslvr.lt.0) iaslvr=0
  if(jaslvr.lt.0) jaslvr=0
  if(inslvr.lt.0) inslvr=0
  if(jnslvr.lt.0) jnslvr=0

  if(iamacc.eq.1) then
    isolver=SOLVER_EISPACK
    if(iaslvr.eq.0) isolver=SOLVER_EISPACK
    if(iaslvr.eq.1) isolver=SOLVER_CUSOLVER
    if(iaslvr.eq.2) isolver=SOLVER_CUSOLVERJ
    jsolver=SOLVER_EISPACK
    if(jaslvr.eq.0) jsolver=SOLVER_EISPACK
    if(jaslvr.eq.1) jsolver=SOLVER_CUSOLVER
    if(jaslvr.eq.2) jsolver=SOLVER_CUSOLVERJ
  else
    isolver=SOLVER_EISPACK
    if(inslvr.eq.0) isolver=SOLVER_EISPACK
    if(inslvr.eq.1) isolver=SOLVER_MKL
    jsolver=SOLVER_EISPACK
    if(jnslvr.eq.0) jsolver=SOLVER_EISPACK
    if(jnslvr.eq.1) jsolver=SOLVER_MKL
  endif

  if(me.eq.mstr.and.ipr.ge.20) then
    write(lfnout,610)
610 format(/,' Linear algebra solvers',/)
    if(numacc.gt.0) then
      if(iaslvr.eq.SOLVER_EISPACK) write(lfnout,611)
611   format(' Accelerated ranks use EISPACK SVD on CPU')
      if(iaslvr.eq.SOLVER_CUSOLVER) write(lfnout,612)
612   format(' Accelerated ranks use CUSOLVER QR gesvd')
      if(iaslvr.eq.SOLVER_CUSOLVERJ) write(lfnout,613)
613   format(' Accelerated ranks use CUSOLVER Jacobi gesvdj')

      if(jaslvr.eq.SOLVER_EISPACK) write(lfnout,614)
614   format(' Accelerated ranks use EISPACK TRED2/TQL on CPU')
      if(jaslvr.eq.SOLVER_CUSOLVER) write(lfnout,615)
615   format(' Accelerated ranks use CUSOLVER QR syevd')
      if(jaslvr.eq.SOLVER_CUSOLVER) write(lfnout,616)
616   format(' Accelerated ranks use CUSOLVER Jacobi syevj')

    endif
    if(inslvr.eq.SOLVER_EISPACK) write(lfnout,617)
617 format(' Non-accelerated ranks use EISPACK SVD')
    if(inslvr.eq.SOLVER_MKL) write(lfnout,618)
618 format(' Non-accelerated ranks use MKL QR dgesvd')
    if(jnslvr.eq.SOLVER_EISPACK) write(lfnout,619)
619 format(' Non-accelerated ranks use EISPACK TRED2/TQL on CPU')
    if(jnslvr.eq.SOLVER_MKL) write(lfnout,620)
620 format(' Non-accelerated ranks use MKL QR dsyevd')
  endif

  if(me.eq.mstr) numdev=0

  call timer_stop(3)
  call timer_start(7)

  if(me.eq.mstr) then
    call timer_stop(99)
    call swatch(date,time)
    write(lfnday,1707) date(1:8),time(1:8),timer_wall_total(99), &
        '  :  Start of base state generation'
1707 format(a8,2x,a8,f12.3,a)
    flush(lfnday)
    call timer_start(99)
  endif
  
  maxcoef=0.0
  numact = 0
  allocate(alldets(nbase))
  ! First pass through make_basestate to determine maxcib
  first_pass = .true.
  do i=1,nbase

    if(idbg.ge.50) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8), &
          ' entering make_basestate for ibase ',i
      flush(lfndbg)
    endif

    call gronor_make_basestate(i,first_pass)

    if(idbg.ge.50) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8), &
          ' returned from make_basestate for ibase ',i
      flush(lfndbg)
    endif

    if(me.eq.mstr) then
      call timer_stop(99)
      call swatch(date,time)
      write(lfnday,1706) date(1:8),time(1:8),timer_wall_total(99), &
          '  :  First pass for base state  ',i,' completed '
1706  format(a8,2x,a8,f12.3,a,i4,a)
      flush(lfnday)
      call timer_start(99)
    endif

    if(idbg.ge.50) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8), &
          ' Base state completed for ibase ',i
      flush(lfndbg)
    endif
  enddo

  ! Now that maxcib is known, civb and iocc can be allocated
  allocate(civb(maxcib,nbase))
  allocate(iocc(maxcib,nbase,numact))

  ! Second pass through make_basestate to actually fill civb and iocc
  first_pass = .false.
  do i=1,nbase

    if(idbg.ge.50) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8), &
          ' entering make_basestate for ibase ',i
      flush(lfndbg)
    endif

    call gronor_make_basestate(i,first_pass)

    if(idbg.ge.50) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8), &
          ' returned from make_basestate for ibase ',i
      flush(lfndbg)
    endif

    if(me.eq.mstr) then
      call timer_stop(99)
      call swatch(date,time)
      write(lfnday,1706) date(1:8),time(1:8),timer_wall_total(99), &
          '  :  Second pass for base state ',i,' completed '
      flush(lfnday)
      call timer_start(99)
    endif

    if(idbg.ge.50) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8), &
          ' Base state completed for ibase ',i
      flush(lfndbg)
    endif
  enddo
  
  thres2=tau_CI / maxcoef
  allocate(civb_aux(maxcib))
  allocate(occu_aux(maxcib))
  allocate(iocc_aux(maxcib,numact))
  do ibase=1, nbase
    ndeti=idetb(ibase)
    ndet_rev=0
    civb_aux=0.0
    occu_aux=0
    do idet=1,idetb(ibase)
      if(abs(civb(idet,ibase)).gt.thres2) then
        ndet_rev=ndet_rev + 1
        civb_aux(ndet_rev)=civb(idet,ibase)
        do iact=1,nactb(ibase)
          iocc_aux(ndet_rev,iact)=iocc(idet,ibase,iact)
        enddo
      endif
    enddo
    do idet=1,maxcib
      civb(idet,ibase)=0.0
      do iact=1,numact
        iocc(idet,ibase,iact)=0
      enddo
    enddo
    do idet=1,ndet_rev
      civb(idet,ibase)=civb_aux(idet)
      do iact=1,numact
        iocc(idet,ibase,iact)=iocc_aux(idet,iact)
      enddo
    enddo
    idetb(ibase)=ndet_rev
    if(me.eq.mstr) then
      call timer_stop(99)
      call swatch(date,time)
      write(lfnday,1710) date(1:8),time(1:8),timer_wall_total(99), &
          '  :  Base state ',ibase, &
          ' coefficient list reduced from    ', &
          alldets(ibase),' to ',ndet_rev
1710  format(a8,2x,a8,f12.3,a,i4,a,i16,a,i16)
      flush(lfnday)
      call timer_start(99)
    endif
  enddo
  deallocate(civb_aux,occu_aux)
  deallocate(iocc_aux)

  do i=1,nbase
    if(idbg.gt.50) then
      ndeti=idetb(i)
      nacti=nactb(i)
      write(lfndbg,1601) i,ndeti
1601  format(/,' Basestate ',i4,' wavefunction has ',i8,' determinants:',/)
      do ib=1,ndeti
        write(lfndbg,1602) ib,civb(ib,i),(iocc(ib,i,k),k=1,nacti)
1602    format(i5,f15.5,32i3)
      enddo
      write(lfndbg,1603) i,nactb(i)+inactb(i),nbas
1603  format(/,' Basestate ',i4,' vector ',2i8)
      do ivc=1,nactb(i)+inactb(i)
        write(lfndbg,1604) '(',ivc,')',(vecsb(ibas,ivc,i),ibas=1,nbas)
1604    format(a2,i3,a1,(t9,10f12.8))
      enddo
      flush(lfndbg)
    endif
  enddo

  call timer_stop(7)
  call timer_start(8)

  if(me.eq.mstr) then
    call timer_stop(99)
    call swatch(date,time)
    write(lfnday,1708) date(1:8),time(1:8),timer_wall_total(99),'  :  Start of ME list generation'
1708 format(a8,2x,a8,f12.3,a)
    flush(lfnday)
    call timer_start(99)
  endif

  allocate(c2sum(nbase,nbase),nbdet(nbase,nbase))

  numdet=0
  memax=0
  icur=0
  jcur=0
  ndtot=0
  do ibase=1,nbase
    do jbase=1,ibase
      melen=0
      ndeti=idetb(ibase)
      ndetj=idetb(jbase)
      ibd=0
      c2s=0.0d0
      if(ibase.eq.jbase) then
        ndtot=ndtot+ndeti*(ndeti+1)/2
        do i=1,ndeti
          do j=i,ndeti
            if(dabs(civb(i,ibase)*civb(j,jbase)).lt.tau_CI) exit
            numdet=numdet+1
            melen=melen+1
            ibd=ibd+1
            if(i.eq.j) then
              c2s=c2s+(civb(i,ibase)*civb(j,jbase))**2
            else
              c2s=c2s+2.0*(civb(i,ibase)*civb(j,jbase))**2
            endif
          enddo
        enddo
      else
        ndtot=ndtot+ndeti*ndetj
        do i=1,ndeti
          do j=1,ndetj
            if(dabs(civb(i,ibase)*civb(j,jbase)).lt.tau_CI) exit
            numdet=numdet+1
            melen=melen+1
            ibd=ibd+1
            c2s=c2s+(civb(i,ibase)*civb(j,jbase))**2
          enddo
        enddo
      endif
      nbdet(ibase,jbase)=ibd
      nbdet(jbase,ibase)=ibd
      c2sum(ibase,jbase)=c2s
      c2sum(jbase,ibase)=c2s
      memax=max(memax,melen)
    enddo
  enddo

  if(me.eq.mstr) then
    call timer_stop(99)
    call swatch(date,time)
    write(lfnday,1703) date(1:8),time(1:8),timer_wall_total(99), &
        '  :  ME list dimension reduced from   ',ndtot,' to ',numdet
1703 format(a8,2x,a8,f12.3,a,16x,i16,a,i16)
    flush(lfnday)
    call timer_start(99)
  endif

  allocate(ndxdet(nbase,nbase))

  numdet=0
  do ibase=1,nbase
    do jbase=1,ibase
      ndeti=idetb(ibase)
      ndetj=idetb(jbase)
      if(ibase.eq.jbase) then
        do i=1,ndeti
          do j=i,ndeti
            if(dabs(civb(i,ibase)*civb(j,jbase)).lt.tau_CI) exit
            numdet=numdet+1
            melen=melen+1
          enddo
        enddo
      else
        do i=1,ndeti
          do j=1,ndetj
            if(dabs(civb(i,ibase)*civb(j,jbase)).lt.tau_CI) exit
            numdet=numdet+1
            melen=melen+1
          enddo
        enddo
      endif
      ndxdet(ibase,jbase)=numdet
    enddo
  enddo

  if(me.ne.mstr) then
    lnxt=0
    lcur=0
    do ibase=1,nbase
      do jbase=1,ibase
        lcur=ndxdet(ibase,jbase)
        ndxdet(ibase,jbase)=lnxt
        lnxt=lcur
      enddo
    enddo
  endif

  if(me.eq.mstr) then
    call timer_stop(99)
    call swatch(date,time)
    write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),'  :  Base states generated'
    flush(lfnday)
    call timer_start(99)
  endif

  allocate(hbase(nbase,nbase),sbase(nbase,nbase),tbase(nbase,nbase))
  allocate(dqbase(nbase,nbase,9))      
  allocate(hev(nbase))
  allocate(nsing(nbase,nbase,5))
  allocate(bpdone(nbase,nbase))

  if(me.eq.mstr) then
    do i=1,nbase
      do j=1,nbase
        bpdone(i,j)=.false.
        hbase(i,j)=0.0d0
        sbase(i,j)=0.0d0
        tbase(i,j)=0.0d0            
        do ksr=1,9
          dqbase(i,j,ksr)=0.0d0
        enddo
        do ksr=1,4
          nsing(i,j,ksr)=-1
          nsing(j,i,ksr)=-1
        enddo
      enddo
    enddo
    rewind(lfncpr)
    do k=1,nbase*nbase
      read(lfncpr,end=299) i,j,rh,rs,rt,nsr,(rm(ksr),ksr=1,9)
      if(i.gt.0.and.j.gt.0) then
        hbase(i,j)=rh
        sbase(i,j)=rs
        tbase(i,j)=rt
        hbase(j,i)=rh
        sbase(j,i)=rs
        tbase(j,i)=rt
        do ksr=1,4
          nsing(i,j,ksr)=nsr(ksr)
          nsing(j,i,ksr)=nsr(ksr)
        enddo
        do ksr=1,9
          dqbase(i,j,ksr)=rm(ksr)
          dqbase(j,i,ksr)=rm(ksr)
        enddo
        bpdone(i,j)=.true.
        bpdone(j,i)=.true.
      endif
    enddo
299 continue
    backspace(lfncpr)
  endif

  !     Define the process groups

  if(npg.le.0) npg=1
  if(npg.gt.np) npg=np-numidle-1
  ngr=np/npg

  if(numidle.gt.0) then
    npg=(np-numidle-1)/mgr
    ngr=(np-numidle-1)/npg
  endif

  if(me.eq.mstr.and.ipr.ge.0) then
    write(lfnout,621) np,numgrp,np-numgrp*mgr-1,mgr,numacc,numnon
621 format(/,' Rank Distribution',//, &
        ' Number of ranks',t40,i10,t60, &
        'Number of rank groups',t100,i10,/, &
        ' Number of idle ranks',t40,i10,t60, &
        'Number of ranks per group',t100,i10,/, &
        ' Number of accelerated ranks',t40,i10,/, &
        ' Number of non-accelerated ranks',t40,i10)
  endif
  ! to be removed
  flush(lfnout)

  call timer_stop(8)
  call timer_start(9)

  if(me.eq.mstr) then
    call timer_stop(99)
    call swatch(date,time)
    write(lfnday,1702) date(1:8),time(1:8),timer_wall_total(99),'  :  Start reading integrals'
1702 format(a8,2x,a8,f12.3,a)
    flush(lfnday)
    call timer_start(99)
  endif

  !     Read the integrals from lfnint

  call gronor_read_integrals()

  if(me.eq.mstr) then
    call timer_stop(99)
    call swatch(date,time)
    write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),'  :  Reading of integrals completed'
    flush(lfnday)
    call timer_start(99)
  endif

  call timer_stop(9)
  
  if(me.eq.mstr.and.ipr.ge.0) then
#ifdef SINGLEP
    rint=dble(4*int2)*1.073741824d-9
#else
    rint=dble(8*int2)*1.073741824d-9
#endif
    rndx=(dble(mlab*4)+dble(mlab*4))*1.073741824d-09
    rlst=dble(2*8*numdet)*1.073741824d-9
    igb=1
    if(max(rint,rndx,rlst).lt.1.0d0) then
      igb=0
      rint=rint*1024.0
      rndx=rndx*1024.0
      rlst=rlst*1024.0
    endif
    if(ipr.ge.3) then
      if(igb.eq.1) then
        write(lfnout,622) nbas,rndx
        write(lfnout,623) int1
        write(lfnout,624) int2,rint
        write(lfnout,641) numdet,rlst
622     format(/,' Number of basisfunctions',t50,i16,t80, &
            ' Size of index arrays',t115,f8.3,' GB')
623     format(' Number of one-electron integrals',t50,i16,t123,' GB')
624     format(' Number of two-electron integrals',t50,i16,t80, &
            ' Size of two-electron integrals',t115,f8.3,' GB')
641     format(' Number of determinant pairs',t50,i16,t80, &
            ' Size of determinant pair list',t115,f8.3,' GB')
      else
        write(lfnout,625) nbas,rndx
        write(lfnout,626) int1
        write(lfnout,627) int2,rint
        write(lfnout,642) numdet,rlst
625     format(/,' Number of basisfunctions',t50,i16,t80, &
            ' Size of index arrays',t115,f8.3,' MB')
626     format(' Number of one-electron integrals',t50,i16,t123,' MB')
627     format(' Number of two-electron integrals',t50,i16,t80, &
            ' Size of two-electron integrals',t115,f8.3,' MB')
642     format(' Number of determinant pairs',t50,i16,t80, &
            ' Size of determinant pair list',t115,f8.3,' MB')
      endif
    else
      write(lfnout,628) nbas
      write(lfnout,629) int1
      write(lfnout,630) int2
628   format(/,' Number of basisfunctions',t50,i16)
629   format(' Number of one-electron integrals',t50,i16)
630   format(' Number of two-electron integrals',t50,i16)
    endif
  endif

  if(me.eq.mstr.and.ipr.gt.0) then
    write(lfnout,631) nbase,nbase*(nbase+1)/2
631 format(' Number of base states',t50,i16,/, &
        ' Number of unique Hamiltonian matrix elements',t50,i16)
  endif

  if(me.eq.mstr.and.ipr.ge.30) then
#ifdef SINGLEP
    write(lfnout,639)
639 format(/,' Integrals are used in single precision')
#else
    write(lfnout,640)
640 format(/,' Integrals are used in double precision')
#endif
  endif
  if(me.eq.mstr) flush(lfnout)

  if(me.eq.mstr) then
    allocate(numrecs(np))
    do i=1,np
      numrecs(i)=0
    enddo
  endif

  if(me.ne.mstr) call timer_start(98)
  call timer_start(4)

  if(idbg.gt.0) then
    call swatch(date,time)
#ifdef GPUAMD
    write(lfndbg,'(a,1x,a,1x,a,11i5)') date(1:8),time(1:8), &
        ' AMD    ',numdev,mydev,iamacc,nummps,numgpu,(map2(me+1,i),i=1,5)
#else
    write(lfndbg,'(a,1x,a,1x,a,11i5)') date(1:8),time(1:8), &
        ' NVIDIA ',numdev,mydev,iamacc,nummps,numgpu,(map2(me+1,i),i=1,5)
#endif
    flush(lfndbg)
  endif

  if(me.eq.mstr) then
    if(idbg.gt.0) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),' Calling GronOR_master'
      flush(lfndbg)
    endif
    call gronor_memory_usage()
    call gronor_master()
  else
    if(idbg.gt.0) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,1x,a,2i5)') date(1:8),time(1:8), &
 &         ' iamactive, iamacc=',iamactive,iamacc
      flush(lfndbg)
    endif
    if(iamactive.eq.1) then

      l2=0
      mnact=0
      mvec=0
      do ibase=1,nbase
        do jbase=1,ibase
          nidet=idetb(ibase)
          njdet=idetb(jbase)
          if(ibase.eq.jbase) then
            l2=max(l2,nidet*(nidet+1)/2)
          else
            l2=max(l2,nidet*njdet)
          endif
        enddo
        mnact=max(mnact,nactb(ibase))
        mvec=max(mvec,nactb(ibase)+inactb(ibase))
      enddo

      nelecs=0
      nveca=0
      n=0
      do ibase=1,nbase
        nveca=max(nveca,inactb(ibase)+nactb(ibase))
        n=2*inactb(ibase)
        do iact=1,nactb(ibase)
          n=n+iabs(int(iocc(1,ibase,iact),kind=kind(n)))
        enddo
        nelecs=max(nelecs,n)
      enddo

      nvecb=nveca
      nstdim=max(1,nelecs*nelecs,nbas*(nbas+1)/2)
      mbasel=max(nelecs,nbas)

      allocate(ioccup(mnact,2))
      allocate(vec(mvec,mbasel,2))
      allocate(vtemp(mvec,mbasel,2))
      allocate(itemp(21),ioccn(20,2))
      allocate(st(nstdim))

      allocate(va(nveca,mbasel))
      allocate(vb(nvecb,mbasel))
      allocate(veca(mbasel))
      allocate(vecb(mbasel))
      allocate(ta(mbasel,max(mbasel,nveca)))
      allocate(taa(mbasel,max(mbasel,nveca)))
      allocate(aaa(mbasel,max(mbasel,nveca)))
      allocate(tb(mbasel,nvecb))
      allocate(s12d(mbasel,mbasel))
      allocate(tt(mbasel,max(mbasel,nveca)))
      allocate(aat(mbasel,max(mbasel,nveca)))
      allocate(sm(mbasel,max(mbasel,nveca)))

      allocate(a(nelecs,nelecs))
      allocate(u(nelecs,nelecs))
      allocate(w(nelecs,nelecs))
      allocate(ev(nelecs))
      allocate(sdiag(max(nelecs,nbas,mbasel)))
      allocate(diag(max(nelecs,nbas,mbasel)))
      allocate(bsdiag(max(nelecs,nbas,mbasel)))
      allocate(bdiag(max(nelecs,nbas,mbasel)))
      allocate(csdiag(max(nelecs,nbas,mbasel)))
      allocate(cdiag(max(nelecs,nbas,mbasel)))

      if(nbatch.lt.0) then
        allocate(prefac(ntask))
        allocate(diagl(ntask,nbas))
        allocate(bsdiagl(ntask,nbas))
        allocate(bdiagl(ntask,nbas))
        allocate(csdiagl(ntask,nbas))
        allocate(sml(ntask,nbas,nbas))
        allocate(tal(ntask,nbas,nbas))
        allocate(ttl(ntask,nbas,nbas))
        allocate(aatl(ntask,nbas,nbas))
        allocate(aaal(ntask,nbas,nbas))
        allocate(tatl(ntask,nbas,nbas))
        allocate(prefac0(1))
        allocate(sm0(1,1,1))
        allocate(ta0(1,1,1))
        allocate(tt0(1,1,1))
        allocate(aat0(1,1,1))
        allocate(aaa0(1,1,1))
        allocate(prefac1(1))
        allocate(diag1(1,1))
        allocate(bsdiag1(1,1))
        allocate(bdiag1(1,1))
        allocate(csdiag1(1,1))
        allocate(sm1(1,1,1))
        allocate(ta1(1,1,1))
        allocate(tt1(1,1,1))
        allocate(aat1(1,1,1))
        allocate(aaa1(1,1,1))
      elseif(nbatch.eq.0) then
        allocate(diagl(1,1))
        allocate(bsdiagl(1,1))
        allocate(bdiagl(1,1))
        allocate(csdiagl(1,1))
        allocate(sml(1,1,1))
        allocate(tal(1,1,1))
        allocate(ttl(1,1,1))
        allocate(aatl(1,1,1))
        allocate(aaal(1,1,1))
        allocate(tatl(1,1,1))
        allocate(prefac(1))
        allocate(prefac0(1))
        allocate(sm0(1,1,1))
        allocate(ta0(1,1,1))
        allocate(tt0(1,1,1))
        allocate(aat0(1,1,1))
        allocate(aaa0(1,1,1))
        allocate(prefac1(1))
        allocate(diag1(1,1))
        allocate(bsdiag1(1,1))
        allocate(bdiag1(1,1))
        allocate(csdiag1(1,1))
        allocate(sm1(1,1,1))
        allocate(ta1(1,1,1))
        allocate(tt1(1,1,1))
        allocate(aat1(1,1,1))
        allocate(aaa1(1,1,1))
      else
        allocate(prefac0(nbatch))
        allocate(prefac1(nbatch))
        if(iamacc.eq.0) then
          allocate(sm0(nbatch,nbas,nbas))
          allocate(ta0(nbatch,nbas,nbas))
          allocate(tt0(nbatch,nbas,nbas))
          allocate(aat0(nbatch,nbas,nbas))
          allocate(aaa0(nbatch,nbas,nbas))
          allocate(diag1(nbatch,nbas))
          allocate(bsdiag1(nbatch,nbas))
          allocate(bdiag1(nbatch,nbas))
          allocate(csdiag1(nbatch,nbas))
          allocate(sm1(nbatch,nbas,nbas))
          allocate(ta1(nbatch,nbas,nbas))
          allocate(tt1(nbatch,nbas,nbas))
          allocate(aat1(nbatch,nbas,nbas))
          allocate(aaa1(nbatch,nbas,nbas))
        else
          allocate(sm0(nbas,nbas,nbatch))
          allocate(ta0(nbas,nbas,nbatch))
          allocate(tt0(nbas,nbas,nbatch))
          allocate(aat0(nbas,nbas,nbatch))
          allocate(aaa0(nbas,nbas,nbatch))
          allocate(diag1(nbas,nbatch))
          allocate(bsdiag1(nbas,nbatch))
          allocate(bdiag1(nbas,nbatch))
          allocate(csdiag1(nbas,nbatch))
          allocate(sm1(nbas,nbas,nbatch))
          allocate(ta1(nbas,nbas,nbatch))
          allocate(tt1(nbas,nbas,nbatch))
          allocate(aat1(nbas,nbas,nbatch))
          allocate(aaa1(nbas,nbas,nbatch))
        endif
        allocate(diagl(1,1))
        allocate(bsdiagl(1,1))
        allocate(bdiagl(1,1))
        allocate(csdiagl(1,1))
        allocate(sml(1,1,1))
        allocate(tal(1,1,1))
        allocate(ttl(1,1,1))
        allocate(aatl(1,1,1))
        allocate(aaal(1,1,1))
        allocate(tatl(1,1,1))
        allocate(prefac(1))
      endif

      allocate(w1(max(nelecs,nbas,mbasel)))
      allocate(w2(max(nelecs,nbas,mbasel),max(nelecs,nbas,mbasel)))

      allocate(temp(nelecs,nelecs),rwork(nelecs))

      if(iamacc.eq.1) then
        if(idbg.gt.0) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a,i12)') date(1:8),time(1:8),' mint2= ',mint2
          flush(lfndbg)
        endif
#ifdef ACC
!$acc data copyin(g,lab,ndx,t,v,dqm,ndxtv,s) &
!$acc& create(a,ta,tb,w1,w2,taa,u,w,ev,temp,rwork) &
!$acc& create(diag,bdiag,cdiag,bsdiag,csdiag,sdiag,aaa,tt,aat,sm) &
!$acc& create(diagl,bdiagl,bsdiagl,csdiagl,sml,aaal,ttl,aatl,tatl,tal) &
!$acc& create(sm0,aaa0,tt0,aat0,ta0,ta1) &
!$acc& create(diag1,bdiag1,bsdiag1,csdiag1,sm1,aaa1,tt1,aat1)
#endif
#ifdef OMPTGT
!$omp target data map(to:g,lab,ndx,t,v,dqm,ndxtv,s) &
!$omp& map(alloc:a,ta,tb,w1,w2,taa,u,w,ev,temp,rwork) &
!$omp& map(alloc:diag,bdiag,cdiag,bsdiag,csdiag,sdiag,aaa,tt,aat,sm) &
!$omp& map(alloc:diagl,bdiagl,bsdiagl,csdiagl,sml,aaal,ttl,aatl,tatl) &
!$omp& map(alloc:tal,aaa0,tt0,aat0,ta0,ta1) &
!$omp& map(alloc:diag1,bdiag1,bsdiag1,csdiag1,sm1,aaa1,tt1,aat1)
#endif
        if(idbg.gt.0) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),' Calling GronOR_worker'
          flush(lfndbg)
        endif
        call gronor_memory_usage()
        if(managers.eq.0) then
          if(role.eq.worker) call gronor_worker()
          if(role.eq.idle) call gronor_idle()
!          call gronor_worker()
        else
          if(role.eq.worker) call gronor_worker()
          if(role.eq.manager) call gronor_manager()
          if(role.eq.idle) call gronor_idle()
        endif
#ifdef ACC
!$acc end data
#endif
#ifdef OMPTGT
!$omp end target data
#endif
      elseif(ntask.ne.0) then
        lwork=10
        allocate(work(lwork))
        call gronor_memory_usage()
        if(managers.eq.0) then
          if(role.eq.worker) call gronor_worker()
          if(role.eq.idle) call gronor_idle()
!          call gronor_worker()
        else
          if(role.eq.worker) call gronor_worker()
          if(role.eq.manager) call gronor_manager()
          if(role.eq.idle) call gronor_idle()
        endif
      endif

      if(nbatch.lt.0) then
        deallocate(tal,tatl,aaal,aatl,ttl,sml)
        deallocate(csdiagl,bdiagl,bsdiagl,diagl)
        deallocate(prefac)
      elseif(nbatch.eq.0) then
        deallocate(tal,tatl,aaal,aatl,ttl,sml)
        deallocate(csdiagl,bdiagl,bsdiagl,diagl)
        deallocate(prefac)
      else
        deallocate(ta0,aaa0,aat0,tt0,sm0)
        deallocate(prefac0)
        deallocate(ta1,aaa1,aat1,tt1,sm1)
        deallocate(csdiag1,bdiag1,bsdiag1,diag1)
        deallocate(prefac1)
      endif
      deallocate(cdiag,csdiag,bdiag,bsdiag,diag,sdiag,ev,w,u,a)
      deallocate(aat,tt,sm)
      deallocate(s12d,tb,aaa,taa,ta,vecb,veca,vb,va)
      deallocate(s,st)
      deallocate(temp)
    endif
  endif

  call timer_stop(4)
  call timer_stop(98)
  
  
  if(me.eq.mstr.and.ipr.ge.30) then
    key='     '
    header='Hamiltonian Matrix (unnormalized)'
    call gronor_print_matrix(lfnout,0,0,0,header,key, &
        mebfLabels,mebfLabel,.false.,1.0d0,hbase,nbase,7,6,lablen.le.labmax)
    header='Overlap Matrix (unnormalized)'
    call gronor_print_matrix(lfnout,0,0,0,header,key, &
        mebfLabels,mebfLabel,.false.,1.0d0,sbase,nbase,7,6,lablen.le.labmax)

    
    
    
    write(lfnout,632)
632 format(//,' Hamiltonian Matrix unnormalized')
    call gronor_prtmat(lfnout,hbase,nbase)
    write(lfnout,635)
635 format(//,' Overlap Matrix unnormalized')
    call gronor_prtmat(lfnout,sbase,nbase)
  endif
  
  if(me.eq.mstr) then
    allocate(nociwf(nbase,nbase))
    do i=1,nbase
      do j=1,i-1
        hbase(i,j)=(1.0d0/(dsqrt(sbase(i,i))*dsqrt(sbase(j,j))))*hbase(i,j)
        do k=1,9
          dqbase(i,j,k)=(1.0d0/(dsqrt(sbase(i,i))*dsqrt(sbase(j,j))))*dqbase(i,j,k)
        enddo
        sbase(i,j)=(1.0d0/(dsqrt(sbase(i,i))*dsqrt(sbase(j,j))))*sbase(i,j)
        hbase(j,i)=hbase(i,j)
        do k=1,9
          dqbase(j,i,k)=dqbase(i,j,k)
        enddo
        sbase(j,i)=sbase(i,j)
      enddo
    enddo
    do i=1,nbase
      hbase(i,i)=(1.0d0/sbase(i,i))*hbase(i,i)
      do j=1,9
        dqbase(i,i,j)=(1.0d0/sbase(i,i))*dqbase(i,i,j)
      enddo
      sbase(i,i)=(1.0d0/sbase(i,i))*sbase(i,i)
    enddo
    call gronor_multipoles_nuclear()
    call gronor_print_results(hbase)
    if(ncorr.eq.1) then
      allocate( hcorr(nbase,nbase) )
      hcorr=0.0
      call gronor_correlation_energy( &
          nwt,mstates,ecorr,nbase,nmol,ncombv,hbase,sbase,hcorr)
      if(ipr.le.40) write(lfnout,636)
636   format(//,' After applying a shift to the diagonal of H')
      call gronor_print_results(hcorr)
      deallocate( hcorr )
    endif
    call gronor_print_dipole_moments()
    
    inquire(file=trim(fillog),exist=EXIST)
    open(unit=lfnlog,file=trim(fillog),form='formatted',status='unknown',position='append',err=993)
    call timer_stop(99)
    call swatch(date,time)
    if(.not.exist) then
      write(lfnlog,800)
800   format( &
          80x,'      T     ',30x,'       ',/, &
          80x,'      h     ',30x,'      S',/, &
          80x,'      r     ',30x,'F  B  t',/, &
          80x,'      e    D',30x,'r  a  a',/, &
          80x,'M  G  a  M i',30x,'g  s  t',/ &
          80x,'P  P  d  g s',30x,'m  e  e',/, &
          '  Date     Time      Setup        Main       ', &
          'Total  Nodes  Ranks    Acc nonAcc  S  U  s  r t ', &
          'Solvers   Task      Batch    s  s  s  tau_MO   ', &
          'tau_CI  User Jobname Command Host Compiler',/)
    endif
    write(lfnlog,801) date(1:8),time(1:8), &
        timer_wall_total(99)-timer_wall_total(98), &
        timer_wall_total(98),timer_wall_total(99), &
        nnodes,np,numacc,numnon,nummps,numgpu,num_threads,mgr, &
        idist,iaslvr,jaslvr,inslvr,jnslvr, &
        ntaska,ntask,max(1,nbatcha),max(1,nbatch), &
        nmol,nbase,mstates,tau_MO,tau_CI, &
        trim(user),trim(string), &
        trim(command),trim(host),trim(compiler),trim(target)
801 format(a8,1x,a8,f9.3,2f12.3,4i7,4i3,5i2,4i5,3i3, &
        1pe9.2,e9.2,1x,a,1x,a,1x,a,1x,a,1x,a,1x,a)
    write(lfnlog,803) (hbase(i,i),i=1,nbase)
    if(nbase.gt.1) then
      write(lfnlog,803) &
          (27.2114d0*(hbase(i,i)-hbase(1,1)),i=1,nbase)
      !        endif
      write(lfnlog,803) (hev(i),i=1,nbase)
      !     Close the cml file after closing the remaining open tags
      !        if(nbase.gt.1) then
      write(lfnlog,803) (27.2114d0*(hev(i)-hev(1)),i=1,nbase)
    endif
    call gronor_finalize_cml
    close(unit=lfncml,status='keep')
    !
803 format(5x,10f20.10)
    write(lfnlog,804)
804 format(' ')
    close(unit=lfnlog,status='keep')
    call timer_start(99)
    deallocate(nociwf)
  endif

  if(role.ne.idle) then
    deallocate(t,v,dqm,ndxtv)
    deallocate(lab,ndx,ig,g)
  endif
  deallocate(idetb,nactb,inactb,vecsb,civb,iocc)
  deallocate(ioccm,vecsm,civm,occm_string)
  deallocate(idetm,inactm,nactm,nbasm,ncombv)
  if(nmol.ge.3)deallocate(inter_couplings)

  deallocate(hbase,sbase,tbase,dqbase,hev,nsing,bpdone)

  deallocate(ndxdet)
  
  deallocate(c2sum,nbdet)

  call timer_stop(1)

  if(me.eq.mstr) then
    call timer_stop(99)
    call swatch(date,time)
    write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99), &
        '  :  End of Hamiltonian calculation'
    flush(lfnday)
    call timer_start(99)

    open(unit=lfntim,file=trim(filtim),form='formatted',status='unknown',err=992)

  endif

  call gronor_timings(lfnout,lfnday,lfntim)

  deallocate(map2)
  
  if(me.eq.mstr) then
    flush(lfnout)
    flush(lfnday)
    flush(lfntim)
  endif

  if(me.eq.mstr) then
    write(lfnrnk,4444)
4444 format(//,' Number of tasks completed per rank',/)
    nc=7
    lc=np/nc
    i=np
    j=nc
    mc=mod(i,j)
    if(mc.gt.0) lc=lc+1
    do i=1,lc
      nnc=nc
      do j=1,lc
        if(i+(nnc-1)*lc.gt.np) nnc=nnc-1
      enddo
      write(lfnrnk,4445) (i+(j-1)*lc-1,numrecs(i+(j-1)*lc),j=1,nnc)
4445  format(10(i6,': ',i8,3x))
    enddo
    flush(lfnrnk)
  endif

  if(me.eq.mstr) deallocate(numrecs)

  !      close(unit=lfndbg,status='keep')
  if(me.eq.mstr) then
    if(ipro.ge.3) then
      close(unit=lfnpro,status='keep')
    else
      close(unit=lfnpro,status='delete')
    endif
    close(unit=lfnday,status='keep')
    if(itest.eq.0) then
      close(unit=lfntst,status='delete')
    else
      close(unit=lfntst,status='keep')
    endif
    if(lcpr) then
      close(unit=lfncpr,status='keep')
    else
      close(unit=lfncpr,status='delete')
    endif
    close(unit=lfnarx,status='keep')
    close(unit=lfnxrx,status='keep')
    close(unit=lfncml,status='keep')
    if(itim.eq.0) then
      close(unit=lfntim,status='delete')
    else
      close(unit=lfntim,status='keep')
    endif

    flush(lfnout)

    close(unit=lfnrnk,status='keep')
  endif

  if(me.eq.mstr) then
    if(np.gt.nabort.and.(nalive+numidle+1.ne.np.or.otimeout)) then
      ierror=0
      ierr=0
      if(idbg.gt.0) then
        write(lfndbg,'(a)') "Issuing mpi_abort"
        flush(lfndbg)
      endif
      call swatch(date,time)
      if(ipr.ge.0) write(lfnout,637) trim(date),trim(time)
637   format(/,' Completion/abort of run ',2a10,/)
      flush(lfnout)
      close(unit=lfnout,status='keep')
      call mpi_abort(MPI_COMM_WORLD,ierror,ierr)
    else
      call swatch(date,time)
      if(ipr.ge.0) write(lfnout,638) trim(date),trim(time)
638   format(/,' Completion/finalize of run ',2a10,/)
      flush(lfnout)
      close(unit=lfnout,status='keep')
    endif
  endif
  
  if(idbg.gt.0) then
    write(lfndbg,'(a)') "Issuing mpi_finalize"
    flush(lfndbg)
  endif

  call mpi_finalize(ierr)

  if(idbg.gt.0) then
    call swatch(date,time)
    write(lfndbg,'(a,1x,a,a,a)') date(1:8),time(1:8),' Closing ',trim(fildbg)
    flush(lfndbg)
    close(unit=lfndbg,status='keep')
  endif

  return
992 write(lfnout,983) trim(filtim)
  call gronor_abort(204,trim(filtim))
993 write(lfnout,984) trim(fillog)
  call gronor_abort(205,trim(fillog))
994 write(lfnout,984) trim(fildbg)
  call gronor_abort(206,trim(fildbg))
995 write(lfnout,985) trim(filsys)
  call gronor_abort(207,trim(filsys))
996 write(lfnout,986) trim(filout)
  call gronor_abort(208,trim(filout))
997 write(lfnout,987) trim(filvec)
  call gronor_abort(209,trim(filvec))
998 write(lfnout,988) trim(filciv)
  call gronor_abort(210,trim(filciv))
983 format('Unable to open timings file ',a)
984 format('Unable to open debug file ',a)
985 format('Unable to open system file ',a)
986 format('Unable to open output file ',a)
987 format('Unable to open vects file ',a)
988 format('Unable to open civec file ',a)
end subroutine gronor_main

#ifdef IBM
integer function getcpucount()
  getcpucount=0
  return
end function getcpucount
#endif
