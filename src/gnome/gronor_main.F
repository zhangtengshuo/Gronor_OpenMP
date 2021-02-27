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

#include "tgt_guard.fh"
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
!!     <tr><td>\$(root).sym <td> symmetry of MOs for different sets,
!!                 coordinates, basis set <td> filsym <td> lfnsym <td> OpenMolcas Auxiliary
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

#ifdef CUSOLVER
      use cusolverDn
      use cuda_cusolver
#endif
#ifdef _OPENACC
      use openacc
c     use cuda_functions
#endif
      
      implicit none

      integer (kind=4) :: ierror,ierr
      integer :: iarg,i,j,idum(51),k,l
      integer,allocatable :: jspinmol(:)
      integer :: node
      real (kind=8) :: rdum(5)
      character (len=255) :: string

      real(kind=8), external :: timer_wall_total
      
      integer :: getcpucount
      external :: getcpucount

      integer :: igr,numone,numtwo,maxgrp,nk,ii,il
      integer (kind=4) :: new

      integer :: ibase,jbase,lnxt,lcur,ksr,nsr(4)

      integer                  :: idet,length,ndet_rev
      integer, allocatable     :: sort_occ(:),occu_aux(:)
      real(kind=8),allocatable :: sort_coef(:),civb_aux(:)
      real(kind=8)             :: thres2,max_coef

      real (kind=8) :: rint,rlab,rndx,rh,rs,rt

      integer :: ndtot,npl,nrg,igb

      character (len=64) :: version_type

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

      if(me.eq.master) then

        user='                  '
        host='                  '
        date='                  '
        time='                  '
        cwd='                   '

#ifdef USE_POSIXF
        call pxfgetlogin(user,len4,ierr4)
#else
        call getlog(user)
        call hostnm(host)
        call swatch(date,time)
        call getcwd(cwd)
#endif

!     Read a single string argument 'root' from the command line, that
!     will be used to generate all other files used to read and write
!     data associated with this calculation

!     input file: filinp with logical file number lfninp
!     output file: filout with logical file number lfnout
!     symmetry file: filsym with logical file number lfnsym
!     what is on?: filciv with logical file number lfnciv
!     what is on?: filvec with logical file number lfnvec

        iarg=0
        call getarg(iarg,command)
        iarg=1
        call getarg(iarg,string)

        root=string(1:inp_strlen(string))
        if(index(string,'_').gt.0) root=string(1:index(string,'_')-1)
        filinp=string(1:inp_strlen(string))//'.inp'
        lfninp=5
        filout=string(1:inp_strlen(string))//'.out'
        lfnout=16
        filsym=root(1:inp_strlen(root))//'.sym'
        lfnsym=7
        lfnciv=9
        lfnvec=10
        filint=root(1:inp_strlen(root))//'.int'
        filone=root(1:inp_strlen(root))//'.one'
        lfnint=11
        lfnone=11
        lfntwo=12
        lfndbg=13
        filpro=string(1:inp_strlen(string))//'.pro'
        lfnpro=14
        fildat=string(1:inp_strlen(string))//'.dat'
        lfndat=15
        lfndet=17
        filday=string(1:inp_strlen(string))//'.day'
        lfnday=18
        filtst=string(1:inp_strlen(string))//'.tst'
        lfntst=19
        fillog=root(1:inp_strlen(root))//'.log'
        lfnlog=20
        filcpr=string(1:inp_strlen(string))//'.cpr'
        lfncpr=21
        filarx=string(1:inp_strlen(string))//'.arx'
        lfnarx=22
        filrnk=string(1:inp_strlen(string))//'.rnk'
        lfnrnk=23
        filcml=string(1:inp_strlen(string))//'.cml'
        lfncml=24

        open(unit=lfnout,file=filout(1:inp_strlen(filout)),
     &       form='formatted',status='unknown',err=996)

        open(unit=lfnpro,file=filpro(1:inp_strlen(filpro)),
     &       form='formatted',status='unknown',err=996)

        open(unit=lfnday,file=filday(1:inp_strlen(filday)),
     &       form='formatted',status='unknown',err=996)
        
        open(unit=lfntst,file=filtst(1:inp_strlen(filtst)),
     &       form='formatted',status='unknown',err=996)
        
        open(unit=lfncpr,file=filcpr(1:inp_strlen(filcpr)),
     &       form='unformatted',status='unknown',err=996)
        
        open(unit=lfnarx,file=filarx(1:inp_strlen(filarx)),
     &       form='formatted',status='unknown',err=996)
        
        open(unit=lfnrnk,file=filrnk(1:inp_strlen(filrnk)),
     &       form='formatted',status='unknown',err=996)
        
        call timer_stop(99)
        call swatch(date,time)
        write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  Start of job'
 702    format(a8,2x,a8,f12.3,a)
        flush(lfnday)
        call timer_start(99)

!----------------------
      
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
          write(*,'(a,i5,a,i5,a)') 'Error: Size specified (',mgr,
     &         ') cannot be supported by available ranks (',np,')'
          call errquit(me,2,"Error in main")
        endif

        nmol=0

!     mstates will be set depending on Molecule input parameters

        mstates=0

!     Default batch size is 32
        
        ntask=32

!     Set ntask and nbatch to undefined

        ntask=-1
        nbatch=-1

!     Default is not to write a test file
        
        itest=0

        idevel=0
        ifault=1

!     Default solvers are those included in the code
        
        isolver=0
        jsolver=0
        iaslvr=0
        jaslvr=0
        inslvr=0
        jnslvr=0

        naccel=-1
        ncols=7
        nummps=1
        numomp=1
        ixpert=0
        iswsvj=15
        iswevj=15

!     labels   0 : read labels from the integral file
!     1 : use Molcas integral ordering
        
        labels=1
        
#ifdef _OPENACC
        naccel=0
#endif

!     Default print flag set to medium (20)
        
        ipr=20

!     Default is not to turn timing on

        itim=0

        ipro=-1
        iday=10
        idbg=0
        nspin=0
        mpibuf=168435456
        thresh=1.0d-10          ! What is the difference between thresh
        thresh_CI=1.0d-5        !   and tresh_CI? Simplify?
        thresh_SIN=1.0d-12
        tolsvj=1.0d-07
        tolevj=1.0d-07
        lcpr=.false.

        call gronor_input()
        
        if(mgr.eq.0) then
          if(npg.eq.0) call errquit(me,2,"Error in main")
          mgr=(np-1)/npg
        endif

        if(ntask.eq.-1) then
          if(nbatch.gt.0) then
            ntask=nbatch
          else
            ntask=32
          endif
          nbatch=0
        else
          if(ntask.le.0) ntask=32
        endif
        
        if(nbatch.gt.ntask) nbatch=ntask

!     Unless the user is an expert let's check thresholds
        
        if(ixpert.eq.0) then
          thresh=min(thresh,1.0d-05)
          thresh_CI=min(thresh_CI,1.0d-05)
        endif

!------------------------------

        
        if(_GRONOR_VERSION_MINOR_.eq.0.or.
     &       _GRONOR_VERSION_MINOR_.gt.12) then
          version_type=" under active development"
        else
          version_type=" official release"
        endif
        
        if(ipr.gt.0) write(lfnout,600) _GRONOR_VERSION_MAJOR_,
     &       _GRONOR_VERSION_MINOR_,trim(version_type)
 600    format(/,
     &       ' GronOR: Non-Orthogonal Configuration Interaction',//,
     &       ' Version ',i2,'.',i2.2,a,//,
     &       ' T. P. Straatsma',t50,'C. de Graaf',t82,'R. Broer',/,
     &       t50,'A. Sanchez',t82,'R. K. Kathir',//,
     &       ' National Center for Computational Sciences',t50,
     &       'Quantum Chemistry Group',t82,
     &       'Department of Theoretical Chemistry',/,
     &       ' Oak Ridge National Laboratory',t50,
     &       'University Rovira i Virgili',t82,
     &       'University of Groningen',/,
     &       ' Oak Ridge, Tennessee, USA',t50,
     &       'Tarragona, Spain',t82,
     &       'Groningen, the Netherlands',//,
     &       ' Based on GNOME written by R. Broer-Braam, ',
     &       'J. Th. van Montfort, and B. Vunderink',//,
     &       ' Please cite the following reference when publishing ',
     &       'results obtained using GronOR:',//,
     &       ' T. P. Straatsma, R. Broer, S. Faraji,',
     &       'R. W. A. Havenith, L. E. Aguilar Suarez, R. K. Kathir, ',
     &       'M. Wibowo, and C.de Graaf',/,
     &       ' "GronOR: Massively parallel and GPU-accelerated ',
     &       'non-orthogonal configuration interaction for large ',
     &       'molecular systems"',/,
     &       ' Journal of Chemical Physics, 152 (2020), ',
     &       'https://doi.org/10.1063/1.5141358',/)

        if(ipr.ge.20) write(lfnout,601) trim(user),getcpucount(),
     &       trim(host),nnodes,
     &       date(1:8),nrsets,
     &       time(1:8),nrnsets,numdev,
     &       np,np/nnodes,np/nrsets,
     &       ncycls,num_threads
 601    format(//,
     &       ' User',t30,a,t60,'CPU count',t100,i10,/,/,
     &       ' Host',t30,a,t60,'Number of nodes',t100,i10,/,
     &       ' Date',t30,a,t60,
     &       'Number of resource sets',t100,i10,/,
     &       ' Time',t30,a,t60,
     &       'Number of resource sets per node',t100,i10,/,
     &       t60,'Number of GPUs per resource set',t100,i10,/,
     &       t60,'Number of MPI ranks',t100,i10,/,
     &       t60,'Number of MPI ranks per node',t100,i10,/,
     &       t60,'Number of MPI ranks per resource set',t100,i10,/,
     &       t60,'Number of rank assignment cycles',t100,i10,/,
     &       t60,'Number of OPENMP threads',t100,i10)
        if(ipr.ge.20.and.numdev.gt.0) write(lfnout,602) memfre,memtot
 602    format(t60,'Available memory on device',t90,i20,/,
     &       t60,'Total memory on device',t90,i20)
        if(ipr.gt.0) write(lfnout,603) trim(command),trim(cwd),
     &       filinp(1:inp_strlen(filinp)),filsym(1:inp_strlen(filsym)),
     &       filout(1:inp_strlen(filout)),filone(1:inp_strlen(filone)),
     &       root(1:inp_strlen(root)),
     &       root(1:inp_strlen(root))
 603    format(/,' Command argument',t30,a,/,
     &       ' Current working directory',t30,a,//,
     &       ' Input file is',t25,a,t60,
     &       'Symmetry file is',t92,a,/,
     &       ' Output file is',t25,a,t60,
     &       'One electron integral file is',t92,a,//,
     &       ' CI vector file(s) are',t25,a,'_nnn.det',/
     &       ' MO vector file(s) are',t25,a,'_nnn.vec')

        write(lfnarx,401) trim(user),trim(host),root(1:inp_strlen(root))
 401    format('*** ',a,' on ',a,' : ',a)
        write(lfnarx,402) nnodes,np
 402    format('Nodes',i10,i10)


!----------------------------
  
!----------------
        
        call gronor_number_integrals(numone,numtwo)
        
        if(npg*mgr+1.gt.np) then
          call errquit(me,2,"Requested ranks exceeds available")
        endif
        
        allocate(nbasm(mstates))
        allocate(nactm(mstates))
        allocate(inactm(mstates))
        allocate(idetm(mstates))
        
        iocc(1)=1
        iocc(2)=2
        iord(1)=0
        iord(2)=0
        indx(1)=0
        indx(2)=0
        icalc=0
        ins2=0
        ipvec=0
        idipole=0
        itp4=0
        ione=0
        bias=0.0d0
        corres=.false.

        if(ipr.ge.30) then
          write(lfnout,604) nmol,ifault,mstates,itest,nbase,idevel
 604      format(/,' Number of fragments',t40,i10,
     &         t60,'Fault tolerance wait time',t100,i10,' sec',/
     &         ' Number of fragment wave functions',t40,i10,
     &         t60,'Test option',t100,i10,/,
     &         ' Number of many-electron base states',t40,i10,
     &         t60,'Development option',t100,i10)
        else
          write(lfnout,605) nmol,mstates,nbase
 605      format(/,' Number of fragments',t40,i10,/,
     &       ' Number of fragment wave functions',t40,i10,/,
     &       ' Number of many-electron base states',t40,i10,/)
        endif
        if(ipr.ge.20) then
          write(lfnout,606) ntask,mgr,nbatch,numomp
 606      format(' Task size',t40,i10,
     &         t60,'Number of ranks per task',t100,i10,/,
     &         ' Batch size',t40,i10,
     &         t60,'Number of OpenMP threads per rank',t100,i10)
          if(numdev.gt.0) write(lfnout,638) nummps
 638      format(t60,'Over-subscription of accelerators (MPS)',t100,i10)
        endif

        if(ipr.gt.0) write(lfnout,607) thresh_CI
 607    format(/,' Threshold for expension coefficients',t40,1pe12.3)
        if(ipr.ge.30) write(lfnout,608) thresh_SIN
 608    format(' Threshold for singularity',t40,1pe12.3,/)

        if(ipr.gt.0) write(lfnout,609) nspin+1
 609    format(/,' Spin multiplicity (2S+1)',t40,i10)
        
        call gronor_read_vectors_and_determinants()
        
        call gronor_gnome_molcas_input() ! OpenMolcas
        close(unit=lfnsym)

      endif

      if(numdev.eq.0) then
        if(inslvr.lt.0) inslvr=iaslvr
        if(jnslvr.lt.0) jnslvr=jaslvr
      endif
      
!     Distribute input data to all processes

      if(np.gt.1) then
        
        if(me.eq.master) then
          idum(1)=nmol
          idum(2)=mstates
          idum(3)=nbase
          idum(4)=maxci
          idum(5)=maxvec
          idum(6)=maxnact
          idum(7)=iocc(1)
          idum(8)=iocc(2)
          idum(9)=iord(1)
          idum(10)=iord(2)
          idum(11)=indx(1)
          idum(12)=indx(2)
          idum(13)=icalc
          idum(14)=ipr
          idum(15)=ins2
          idum(16)=ipvec
          idum(17)=itp4
          idum(18)=ione
          idum(19)=0
          if(corres) idum(19)=1
          idum(20)=npg
          idum(21)=nnucl
          idum(22)=nbas
          idum(23)=0
          idum(24)=ntask
          idum(25)=mgr
          idum(26)=itest
          idum(27)=naccel
          idum(28)=inpcib
          idum(29)=ipro
          idum(30)=0
          idum(31)=mpibuf
          idum(32)=idbg
          idum(33)=nummps
          idum(34)=numomp
          idum(35)=itim
          idum(36)=iaslvr
          idum(37)=jaslvr
          idum(38)=inslvr
          idum(39)=jnslvr
          idum(40)=ifault
          idum(41)=idevel
          idum(42)=iswsvj
          idum(43)=iswevj
          idum(44)=nbatch
          idum(45)=nspin
          idum(46)=labels
          
          rdum(1)=thresh
          rdum(2)=thresh_CI
          rdum(3)=thresh_SIN
          rdum(4)=tolsvj
          rdum(5)=tolevj
        endif
        
        call MPI_Bcast(idum,46,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
        call MPI_Bcast(rdum,5,MPI_REAL8,master,MPI_COMM_WORLD,ierr)
        
        nmol   =idum(1)
        mstates=idum(2)
        nbase  =idum(3)
        maxci  =idum(4)
        maxvec =idum(5)
        maxnact=idum(6)
        iocc(1)=idum(7)
        iocc(2)=idum(8)
        iord(1)=idum(9)
        iord(2)=idum(10)
        indx(1)=idum(11)
        indx(2)=idum(12)
        icalc  =idum(13)
        ipr    =idum(14)
        ins2   =idum(15)
        ipvec  =idum(16)
        itp4   =idum(17)
        ione   =idum(18)
        corres =idum(19).eq.1
        npg    =idum(20)
        nnucl  =idum(21)
        nbas   =idum(22)

        ntask =idum(24)
        mgr    =idum(25)
        itest  =idum(26)
        naccel =idum(27)
        inpcib =idum(28)
        ipro   =idum(29)

        mpibuf =idum(31)
        idbg   =idum(32)
        nummps =idum(33)
        numomp =idum(34)
        itim   =idum(35)
        iaslvr =idum(36)
        jaslvr =idum(37)
        inslvr =idum(38)
        jnslvr =idum(39)
        ifault =idum(40)
        idevel =idum(41)
        iswsvj =idum(42)
        iswevj =idum(43)
        nbatch =idum(44)
        nspin  =idum(45)
        labels =idum(46)
        
        if(me.eq.master) then
          call timer_stop(99)
          call swatch(date,time)
          write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &         '  :  Input broadcasted'
          flush(lfnday)
          call timer_start(99)
        endif
        
        thresh=rdum(1)
        thresh_CI=rdum(2)
        thresh_SIN=rdum(3)
        tolsvj=rdum(4)
        tolevj=rdum(5)
        
      endif

#ifdef _OPENMP
      if(numomp.gt.0) num_threads=numomp
      if(num_threads.gt.0) call omp_set_num_threads(num_threads)
#endif

      lfnone=11
      lfntwo=12
      lfndbg=13
      if(idbg.gt.0) then
        write(fildbg,1300) me
 1300   format('GronOR_',i5.5,'.dbg ')
        open(unit=lfndbg,file=fildbg(1:inp_strlen(fildbg)),
     &       form='formatted',status='unknown',err=996)
      endif
      
      nacc0=naccel
      nacc1=naccel
      
      if(me.ne.master) then
        nmol   =idum(1)
        mstates=idum(2)
        nbase  =idum(3)
        maxci  =idum(4)
        maxvec =idum(5)
        maxnact=idum(6)
        
        allocate(ncombv(nmol,nbase))
        allocate(nbasm(mstates))
        allocate(nactm(mstates))
        allocate(inactm(mstates))
        allocate(idetm(mstates))
        allocate(civm(maxci,mstates))
        allocate(vecsm(maxvec,maxvec,mstates))
        allocate(ioccm(maxnact,maxci,mstates))
        
        iocc(1)=idum(7)
        iocc(2)=idum(8)
        iord(1)=idum(9)
        iord(2)=idum(10)
        indx(1)=idum(11)
        indx(2)=idum(12)
        icalc  =idum(13)
        ipr    =idum(14)
        ins2   =idum(15)
        ipvec  =idum(16)
        itp4   =idum(17)
        ione   =idum(18)
        corres =idum(19).eq.1
      endif
      if(np.gt.0) then
        call MPI_Bcast(civm,maxci*mstates,MPI_REAL8,master,
     &       MPI_COMM_WORLD,ierr)
        call MPI_Bcast(vecsm,maxvec*maxvec*mstates,MPI_REAL8,master,
     &       MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ioccm,maxnact*maxci*mstates,MPI_INTEGER8,master,
     &       MPI_COMM_WORLD,ierr)
        call MPI_Bcast(ncombv,nmol*nbase,MPI_INTEGER8,master,
     &       MPI_COMM_WORLD,ierr)
        call MPI_Bcast(nbasm,mstates,MPI_INTEGER8,master,
     &       MPI_COMM_WORLD,ierr)
        call MPI_Bcast(nactm,mstates,MPI_INTEGER8,master,
     &       MPI_COMM_WORLD,ierr)
        call MPI_Bcast(inactm,mstates,MPI_INTEGER8,master,
     &       MPI_COMM_WORLD,ierr)
        call MPI_Bcast(idetm,mstates,MPI_INTEGER8,master,
     &       MPI_COMM_WORLD,ierr)
      endif
      
      nbasis=0
      mactb=0
      do i=1,nmol
        nbasis=nbasis+nbasm(ncombv(i,1))
        mactb=mactb+nactm(ncombv(i,1))
      enddo
      allocate(jspinmol(mstates))
      do i=1,mstates
        jspinmol(i)=0
        do j=1,nactm(i)
          if(abs((ioccm(j,1,i))).eq.1) jspinmol(i)=jspinmol(i)+
     &         ioccm(j,1,i)
        enddo
      enddo

      call timer_stop(2)
      call timer_start(3)

      allocate(vecsb(nbasis,nbasis,nbase))
      allocate(inactb(nbase))
      allocate(nactb(nbase))
      allocate(idetb(nbase))      

      call gronor_determine_nci()
      call gronor_determine_maxcib()

      allocate(ioccb(maxcib,nbase))
      allocate(civb(maxcib,nbase))

      ngr=(np-1)/mgr

      allocate(ranks_heads(ngr+1))

      master=np-1

      if(ncycls.gt.1) then
        j=0
        do i=1,np
          if(map2(i,5).gt.0) j=j+1
          if(i-1.eq.master) map2(i,5)=0
        enddo
        do i=1,min(np,nummps*j)
          map2(i,5)=map2(i,1)
          if(i-1.eq.master) map2(i,5)=0
        enddo
      else
        node=map2(1,4)
        j=nummps*map2(1,1)
        do i=1,np
          map2(i,2)=numomp
          if(map2(i,4).eq.node) then
            if(j.gt.0) then
              map2(i,5)=map2(i,1)
              if(naccel.lt.0) map2(i,5)=-map2(i,2)
              j=j-1
            else
              map2(i,5)=-map2(i,2)
            endif
          else
            node=map2(i,4)
            j=nummps*map2(i,1)
            if(j.gt.0) then
              map2(i,5)=map2(i,1)
              if(naccel.lt.0) map2(i,5)=-map2(i,2)
              j=j-1
            else
              map2(i,5)=-map2(i,2)
            endif
          endif
          if(i-1.eq.master) map2(i,5)=0
        enddo
      endif
      
      if(me.eq.master) then
        numacc=0
        numnon=0
        do i=1,np
          if(i-1.eq.master) map2(i,5)=0
          if(map2(i,5).gt.0) numacc=numacc+1
          if(map2(i,5).lt.0) numnon=numnon+1
        enddo
        numacc=numacc/mgr
        numnon=numnon/mgr
        maxgrp=numacc+numnon
      endif
      call MPI_Bcast(maxgrp,1,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      allocate(thisgroup(mgr+1))
      allocate(allgroups(maxgrp+1,mgr+1))
      allocate(allheads(maxgrp+1))
      if(me.eq.master) then
        numgrp=0
        igr=0
        do i=1,np
          if(map2(i,5).gt.0) then
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
            if(igr.eq.0) numgrp=numgrp+1
            igr=igr+1
            allgroups(numgrp,1)=map2(i,5)
            allgroups(numgrp,igr+1)=i-1
            if(igr.eq.mgr) igr=0
          endif
        enddo

        if(igr.ne.0.and.numgrp.gt.0) numgrp=numgrp-1
        do i=1,numgrp
          allheads(i)=allgroups(i,2)
          do j=1,mgr
            map2(allgroups(i,j+1)+1,3)=i
          enddo
        enddo
        allheads(numgrp+1)=master
      endif
      call MPI_Bcast(numgrp,1,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(map2,5*np,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(allgroups,int((maxgrp+1)*(mgr+1),kind=4),
     &     MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(allheads,int(maxgrp,kind=4),MPI_INTEGER8,master,
     &     MPI_COMM_WORLD,ierr)
      if(numgrp.eq.0) then
        write(*,'(a)') 'Number of groups is zero'
        call errquit(me,2,"Error in main")
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
      
      if(me.ne.master) then
        iamhead=0
        if(me.eq.thisgroup(2)) iamhead=1
        myhead=thisgroup(2)
        numdev=thisgroup(1)
        if(numdev.lt.0) numdev=0
      else
        myhead=master
        numdev=0
      endif
      
      mygroup=0
      do i=1,numgrp
        ranks_heads(i)=allgroups(i,2)
        do j=1,mgr
          if(me.eq.allgroups(i,j+1)) mygroup=i
        enddo
      enddo
      ranks_heads(numgrp+1)=master
      
      new=-1
      
      call MPI_Comm_Group(MPI_COMM_WORLD,group_world,ierr)
      
      call MPI_Group_incl(group_world,int(numgrp+1,kind=4),ranks_heads,
     &     group_heads,ierr)
      if(iamhead.eq.1) then
        call MPI_Group_Rank(group_heads,new,ierr)
      endif
      
      if(me.eq.master) then
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
 730    format(//,' MPI Rank Assignments')
 734    format(' ')
 735    format(/,' Accelerated Groups',/)
 733    format(/,' Non-Accelerated Groups',/)
        k=1
        l=min(k+npl,numacc)
        if(numacc.gt.0) then
          write(lfnrnk,735)
          do while (k.le.l)
            if(mgr.gt.1) write(lfnrnk,734)
            do i=2,mgr+1                
              if(nrg.lt.10) then
                if(np.lt.10) then
                  write(lfnrnk,6611)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6612)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6613)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6614)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6615)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              elseif(nrg.lt.100) then
                if(np.lt.10) then
                  write(lfnrnk,6621)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6622)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6623)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6624)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6625)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              elseif(nrg.lt.1000) then
                if(np.lt.10) then
                  write(lfnrnk,6631)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6632)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6633)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6634)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6635)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              elseif(nrg.lt.10000) then
                if(np.lt.10) then
                  write(lfnrnk,6641)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6642)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6643)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6644)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6645)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              else
                if(np.lt.10) then
                  write(lfnrnk,6651)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6652)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6653)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6654)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6655)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              endif
            enddo
            k=l+1
            l=min(k+npl,numacc)
          enddo
          write(lfnrnk,734)
        endif
        if(numnon.gt.0) then
          write(lfnrnk,733)
          k=l+1
          l=min(k+npl,numgrp)
          do while (k.le.l)
            if(mgr.gt.1) write(lfnrnk,734)
 6611       format(5x,20(i1,':',i1,'  '))
 6612       format(5x,20(i1,':',i2,'  '))
 6613       format(5x,20(i1,':',i3,'  '))
 6614       format(5x,20(i1,':',i4,'  '))
 6615       format(5x,20(i1,':',i5,'  '))
 6621       format(5x,20(i2,':',i1,'  '))
 6622       format(5x,20(i2,':',i2,'  '))
 6623       format(5x,20(i2,':',i3,'  '))
 6624       format(5x,20(i2,':',i4,'  '))
 6625       format(5x,20(i2,':',i5,'  '))
 6631       format(5x,20(i3,':',i1,'  '))
 6632       format(5x,20(i3,':',i2,'  '))
 6633       format(5x,20(i3,':',i3,'  '))
 6634       format(5x,20(i3,':',i4,'  '))
 6635       format(5x,20(i3,':',i5,'  '))
 6641       format(5x,20(i4,':',i1,'  '))
 6642       format(5x,20(i4,':',i2,'  '))
 6643       format(5x,20(i4,':',i3,'  '))
 6644       format(5x,20(i4,':',i4,'  '))
 6645       format(5x,20(i4,':',i5,'  '))
 6651       format(5x,20(i5,':',i1,'  '))
 6652       format(5x,20(i5,':',i2,'  '))
 6653       format(5x,20(i5,':',i3,'  '))
 6654       format(5x,20(i5,':',i4,'  '))
 6655       format(5x,20(i5,':',i5,'  '))
            do i=2,mgr+1
              if(nrg.lt.10) then
                if(np.lt.10) then
                  write(lfnrnk,6611)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6612)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6613)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6614)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6615)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              elseif(nrg.lt.100) then
                if(np.lt.10) then
                  write(lfnrnk,6621)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6622)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6623)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6624)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6625)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              elseif(nrg.lt.1000) then
                if(np.lt.10) then
                  write(lfnrnk,6631)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6632)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6633)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6634)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6635)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              elseif(nrg.lt.10000) then
                if(np.lt.10) then
                  write(lfnrnk,6641)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6642)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6643)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6644)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6645)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              else
                if(np.lt.10) then
                  write(lfnrnk,6651)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.100) then
                  write(lfnrnk,6652)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.1000) then
                  write(lfnrnk,6653)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                elseif(np.lt.10000) then
                  write(lfnrnk,6654)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                else
                  write(lfnrnk,6655)
     &                 ((map2(allgroups(j,i)+1,4),
     &                 allgroups(j,i)),j=k,l)
                endif
              endif
            enddo
            k=l+1
            l=min(k+npl,numgrp)
          enddo
        endif
      endif

      if(me.eq.master) then
        flush(lfnout)
        flush(lfnrnk)
      endif
      
      iamacc=0
      iamactive=0
      if(me.eq.master) iamactive=1
      if(me.ne.master.and.thisgroup(1).gt.0) then
        do i=1,mgr
          if(me.eq.thisgroup(i+1)) iamacc=1
        enddo
      endif
      if(thisgroup(1).ne.0) then
        do i=1,mgr
          if(me.eq.thisgroup(i+1)) iamactive=1
        enddo
      endif
      
      if(iamacc.eq.1) then

!     only accelerated ranks need to set device

#ifdef _OPENACC
        if(numdev.gt.1) then
          mydev=mod(me,numdev)
          call acc_set_device_num(mydev,ACC_DEVICE_NVIDIA)
          cpfre=c_loc(memfre)
          cptot=c_loc(memtot)
c     istat=cudaMemGetInfo(cpfre,cptot)
          memavail=memfre
        endif
#endif
        
!     only accelerated ranks need to define cusolver handles
        
#ifdef CUSOLVER
        cusolver_status = cusolverDnCreate(cusolver_handle)
        if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &       print *, me,'cusolver_handle creation failed'

!     if(numdev.gt.1) then
!     cpfre=c_loc(memfre)
!     cptot=c_loc(memtot)
!     istat=cudaMemGetInfo(cpfre,cptot)
!     endif

#endif

        
      else
        
!     only non-accelerated ranks will use OpenMP threads

#ifdef MKL

!     if solvers were not specified use MKL solver for svd only
        
        isolver=inslvr
        jsolver=jnslvr
        if(isolver.lt.0) isolver=1
        if(jsolver.lt.0) jsolver=0
        if(inslvr.lt.0) inslvr=1
        if(jnslvr.lt.0) jnslvr=0

#endif
        
      endif

!     if solvers were not specified use the default ones
      
      if(isolver.lt.0) isolver=0
      if(jsolver.lt.0) jsolver=0
      if(iaslvr.lt.0) iaslvr=0
      if(jaslvr.lt.0) jaslvr=0
      if(inslvr.lt.0) inslvr=0
      if(jnslvr.lt.0) jnslvr=0
      
      
      if(me.eq.master.and.ipr.ge.20) then
        write(lfnout,610)
 610    format(/,' Linear algebra solvers',/)
        if(numdev.gt.0) then
          if(iaslvr.eq.0) write(lfnout,611)
 611      format(' Accelerated ranks use GronOR SVD on CPU')
#ifdef CUSOLVER
          if(iaslvr.eq.1) write(lfnout,612)
 612      format(' Accelerated ranks use CUSOLVER QR gesvd')
#endif
#ifdef CUSOLVERJ
          if(iaslvr.eq.1) write(lfnout,612)
 612      format(' Accelerated ranks use CUSOLVER QR gesvd')
          if(iaslvr.eq.2) write(lfnout,613)
 613      format(' Accelerated ranks use CUSOLVER Jacobi gesvdj')
#endif
          if(jaslvr.eq.0) write(lfnout,614)
 614      format(' Accelerated ranks use GronOR TRED2/TQL on CPU')
#ifdef CUSOLVER
          if(jaslvr.eq.1) write(lfnout,615)
 615      format(' Accelerated ranks use CUSOLVER QR syevd')
#endif
#ifdef CUSOLVERJ
          if(jaslvr.eq.1) write(lfnout,615)
 615      format(' Accelerated ranks use CUSOLVER QR syevd')
          if(jaslvr.eq.2) write(lfnout,616)
 616      format(' Accelerated ranks use CUSOLVER Jacobi syevj')
#endif
        endif
        if(inslvr.eq.0) write(lfnout,617)
 617    format(' Non-accelerated ranks use GronOR SVD')
#ifdef MKL
        if(inslvr.eq.1) write(lfnout,618)
 618    format(' Non-accelerated ranks use MKL QR dgesvd')
#endif
        if(jnslvr.eq.0) write(lfnout,619)
 619    format(' Non-accelerated ranks use GronOR TRED2/TQL on CPU')
#ifdef CUSOLVER
        if(jnslvr.eq.1) write(lfnout,620)
 620    format(' Non-accelerated ranks use MKL QR dsyevd')
#endif
      endif

      call timer_stop(3)
      call timer_start(7)
      
      if(me.eq.master) then
        call timer_stop(99)
        call swatch(date,time)
        write(lfnday,1707) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  Start of base state generation'
 1707   format(a8,2x,a8,f12.3,a)
        flush(lfnday)
        call timer_start(99)
      endif
      
      max_coef = 0.0
      
      do i=1,nbase

        call gronor_make_basestate(i)
        
        allocate(sort_coef(maxcib))
        allocate(sort_occ(maxcib))
        sort_coef = 0.0
        sort_occ = 0
        do idet = 1, idetb(i)
          sort_coef(idet) = civb(idet,i)
          sort_occ(idet)  = ioccb(idet,i)
        enddo
        length = idetb(i)
        call quicksort(sort_coef,sort_occ,length)
        do idet = 1, idetb(i)
          civb(idet,i)  = sort_coef(idet)
          ioccb(idet,i) = sort_occ(idet)
        enddo
        deallocate(sort_coef,sort_occ)
        if(max_coef.lt.abs(civb(1,i))) then
          max_coef = abs(civb(1,i))
        endif
        
        if(me.eq.master) then
          call timer_stop(99)
          call swatch(date,time)
          write(lfnday,1706) date(1:8),time(1:8),timer_wall_total(99),
     &         '  :  Base state ',i,' completed '
 1706     format(a8,2x,a8,f12.3,a,i4,a)
          flush(lfnday)
          call timer_start(99)
        endif
        
      enddo

      thres2 = thresh_CI / max_coef
      allocate(civb_aux(maxcib))
      allocate(occu_aux(maxcib))
      do ibase = 1, nbase
        ndeti=idetb(ibase)
        ndet_rev = 0
        civb_aux = 0.0
        occu_aux = 0
        do idet = 1, idetb(ibase)
          if(abs(civb(idet,ibase)).gt.thres2) then
            ndet_rev = ndet_rev + 1
            civb_aux(ndet_rev) = civb(idet,ibase)
            occu_aux(ndet_rev) = ioccb(idet,ibase)
          endif
        enddo
        civb(:,ibase) = 0.0
        ioccb(:,ibase) = 0
        do idet = 1, ndet_rev
          civb(idet,ibase) = civb_aux(idet)
          ioccb(idet,ibase) = occu_aux(idet)
        enddo
        idetb(ibase) = ndet_rev
        if(me.eq.master) then
          call timer_stop(99)
          call swatch(date,time)
          write(lfnday,1710) date(1:8),time(1:8),timer_wall_total(99),
     &         '  :  Base state ',ibase,
     &         ' coefficient list reduced from    ',
     &         ndeti,' to ',ndet_rev
 1710     format(a8,2x,a8,f12.3,a,i4,a,i16,a,i16)
          flush(lfnday)
          call timer_start(99)
        endif
      enddo
      deallocate(civb_aux,occu_aux)

      call timer_stop(7)
      call timer_start(8)
      
      if(me.eq.master) then
        call timer_stop(99)
        call swatch(date,time)
        write(lfnday,1708) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  Start of ME list generation'
 1708   format(a8,2x,a8,f12.3,a)
        flush(lfnday)
        call timer_start(99)
      endif
      
      numdet=0
      ndtot=0
      do ibase=1,nbase
        do jbase=1,ibase
          ndeti=idetb(ibase)
          ndetj=idetb(jbase)
          if(ibase.eq.jbase) then
            ndtot=ndtot+ndeti*(ndeti+1)/2
            do i=1,ndeti
              do j=i,ndeti
                if(dabs(civb(i,ibase)*civb(j,jbase)).lt.thresh_CI) exit
                numdet=numdet+1
              enddo
            enddo
          else
            ndtot=ndtot+ndeti*ndetj
            do i=1,ndeti
              do j=1,ndetj
                if(dabs(civb(i,ibase)*civb(j,jbase)).lt.thresh_CI) exit
                numdet=numdet+1
              enddo
            enddo
          endif
        enddo
      enddo

      if(me.eq.master) then
        call timer_stop(99)
        call swatch(date,time)
        write(lfnday,1703) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  ME list dimension reduced from   ',
     &       ndtot,' to ',numdet
 1703   format(a8,2x,a8,f12.3,a,16x,i16,a,i16)
        flush(lfnday)
        call timer_start(99)
      endif

      allocate(ndxdet(nbase,nbase))
      allocate(lstdet(numdet,2))

      numdet=0
      do ibase=1,nbase
        do jbase=1,ibase
          ndeti=idetb(ibase)
          ndetj=idetb(jbase)
          if(ibase.eq.jbase) then
            do i=1,ndeti
              do j=i,ndeti
                if(dabs(civb(i,ibase)*civb(j,jbase)).lt.thresh_CI) exit
                numdet=numdet+1
                lstdet(numdet,1)=i
                lstdet(numdet,2)=j
              enddo
            enddo
          else
            do i=1,ndeti
              do j=1,ndetj
                if(dabs(civb(i,ibase)*civb(j,jbase)).lt.thresh_CI) exit
                numdet=numdet+1
                lstdet(numdet,1)=i
                lstdet(numdet,2)=j
              enddo
            enddo
          endif
          ndxdet(ibase,jbase)=numdet
        enddo
      enddo

!     Following Bcasts should not be necessary, but added for safety

!     call MPI_Bcast(ndxdet,nbase*nbase,MPI_INTEGER8,master,
!     &     MPI_COMM_WORLD,ierr)
!     call MPI_Bcast(lstdet,numdet*2,MPI_INTEGER8,master,
!     &     MPI_COMM_WORLD,ierr)

      if(me.ne.master) then
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
      
      if(me.eq.master) then
        call timer_stop(99)
        call swatch(date,time)
        write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  Base states generated'
        flush(lfnday)
        call timer_start(99)
      endif
      
      allocate(hbase(nbase,nbase),sbase(nbase,nbase),tbase(nbase,nbase))
      allocate(nsing(nbase,nbase,5))
      allocate(bpdone(nbase,nbase))

      if(me.eq.master) then
        do i=1,nbase
          do j=1,nbase
            bpdone(i,j)=.false.
            hbase(i,j)=0.0d0
            sbase(i,j)=0.0d0
            tbase(i,j)=0.0d0
            do ksr=1,4
              nsing(i,j,ksr)=-1
              nsing(j,i,ksr)=-1
            enddo
          enddo
        enddo
        rewind(lfncpr)
        do k=1,nbase*nbase
          read(lfncpr,end=299) i,j,rh,rs,rt,nsr
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
            bpdone(i,j)=.true.
            bpdone(j,i)=.true.
          endif
        enddo
 299    continue
      endif
      
!     Define the process groups
      
      if(npg.le.0) npg=1
      if(npg.gt.np) npg=np
      ngr=np/npg
      
      if(me.eq.master.and.ipr.ge.0) then
        write(lfnout,621) np,numgrp,np-numgrp*mgr-1,mgr
 621    format(/,' Process Distribution',//,
     &       ' Number of processes',t40,i10,t60,
     &       'Number of process groups',t100,i10,/,
     &       ' Number of idle processes',t40,i10,t60,
     &       'Number of processes per group',t100,i10)
      endif
      

      call timer_stop(8)
      call timer_start(9)

!     Read the integrals from lfnint
      
      call gronor_read_integrals()
      
      if(me.eq.master) then
        call timer_stop(99)
        call swatch(date,time)
        write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  Integrals read and distributed'
        flush(lfnday)
        call timer_start(99)
      endif
      
      call timer_stop(9)
      
      if(me.eq.master.and.ipr.ge.0) then
#ifdef SINGLEP
        rint=dble(4*int2)*1.0d-9
#else
        rint=dble(8*int2)*1.0d-9
#endif
        if(labels.eq.0) then
          rlab=dble(8*int2)*1.0d-9
          rndx=0.0d0
        else
          rlab=0.0d0
          rndx=(dble(mlab*4)+dble(mlab*4))*1.0d-09
        endif
        igb=1
        if(max(rint,rlab,rndx).lt.1.0d0) then
          igb=0
          rint=rint*1024.0
          rlab=rlab*1024.0
          rndx=rndx*1024.0
        endif
        if(ipr.ge.20) then
          if(labels.eq.0) then
            if(igb.eq.1) then
              write(lfnout,622) nbas,rndx
              write(lfnout,623) int1,rlab
              write(lfnout,624) int2,rint
 622          format(/,' Number of basisfunctions',t50,i16,t80,
     &             ' Size of index arrays',t115,f8.3,' GB')
 623          format(' Number of one-electron integrals',t50,i16,t80,
     &             ' Size of integral labels',t115,f8.3,' GB')
 624          format(' Number of two-electron integrals',t50,i16,t80,
     &             ' Size of two-electron integrals',t115,f8.3,' GB')
            else
              write(lfnout,625) nbas,rndx
              write(lfnout,626) int1,rlab
              write(lfnout,627) int2,rint
 625          format(/,' Number of basisfunctions',t50,i16,t80,
     &             ' Size of index arrays',t115,f8.3,' MB')
 626          format(' Number of one-electron integrals',t50,i16,t80,
     &             ' Size of integral labels',t115,f8.3,' MB')
 627          format(' Number of two-electron integrals',t50,i16,t80,
     &             ' Size of two-electron integrals',t115,f8.3,' MB')
            endif
          else
            write(lfnout,628) nbas
            write(lfnout,629) int1
 628        format(/,' Number of basisfunctions',t50,i16)
 629        format(' Number of one-electron integrals',t50,i16)
            if(igb.eq.1) then
              write(lfnout,624) int2,rint
            else
              write(lfnout,627) int2,rint
            endif
          endif
        else
          write(lfnout,628) nbas
          write(lfnout,629) int1
          write(lfnout,630) int2
 630      format(' Number of two-electron integrals',t50,i16)
        endif
      endif
      
      if(me.eq.master.and.ipr.gt.0) then
        write(lfnout,631) nbase,nbase*(nbase+1)/2
 631    format(' Number of base states',t50,i16,/,
     &       ' Number of unique Hamiltonian matrix elements',t50,i16)
      endif
      
      if(me.eq.master.and.ipr.ge.1) then
#ifdef SINGLEP
        write(lfnout,639)
 639    format(/,' Integrals are used in single precision')
#else
        write(lfnout,640)
 640    format(/,' Integrals are used in double precision')
#endif
      endif 
      if(me.eq.master) flush(lfnout)
      
      if(me.eq.master) then
        allocate(numrecs(np))
        do i=1,np
          numrecs(i)=0
        enddo
      endif

      call timer_start(98)
      
      if(me.eq.master) then
        call timer_stop(99)
        call swatch(date,time)
        write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  Start Hamiltonian calculation'
        flush(lfnday)
        call timer_start(99)
      endif
      
      if(idbg.gt.0) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,10i5)') date(1:8),time(1:8),
     &       ' NVIDIA ',numdev,mydev,iamacc,nummps,
     &       (map2(me+1,i),i=1,5)
        flush(lfndbg)
      endif

      if(me.eq.master) then
        if(idbg.gt.0) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &         ' Calling GronOR_master'
          flush(lfndbg)
        endif
        call gronor_master()
      else
        if(idbg.gt.0) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a,2i5)') date(1:8),time(1:8),
     &         ' iamactive, iamacc = ',iamactive,iamacc
          flush(lfndbg)
        endif
        if(iamactive.eq.1) then
          if(iamacc.eq.1) then
            if(idbg.gt.0) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,1x,a,i12)') date(1:8),time(1:8),
     &             ' mint2= ',mint2
              flush(lfndbg)
            endif            
!_ACCTGT_($acc data copyin(g,lab,ndx,t,v,ndxtv))
!_OMPTGT_($omp target data map(to:g,lab,ndx,t,v,ndxtv))
            if(idbg.gt.0) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &             ' Calling GronOR_worker'
              flush(lfndbg)
            endif
            call gronor_worker()
!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)
          else
            lwork=10
            allocate(work(lwork))
            call gronor_worker()
          endif
        endif
      endif

      call timer_stop(98)
      
      if(idbg.gt.0) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,a,a)') date(1:8),time(1:8),
     &       ' Closing ',trim(fildbg)
        flush(lfndbg)
        close(unit=lfndbg,status='keep')
      endif

      if(me.eq.master.and.ipr.ge.20) then
        write(lfnout,632)
 632    format(//,' Hamiltonian Matrix unnormalized',/)

        nk=nbase/ncols
        if(mod(nbase,ncols).ne.0) nk=nk+1
        do k=1,nk
          ii=(k-1)*ncols+1
          il=min(nbase,k*ncols)
          write(lfnout,633) (i,i=ii,il)
          do j=1,nbase
            write(lfnout,634) j,(hbase(i,j),i=ii,il)
          enddo
          write(lfnout,*) ' '
        enddo

        nk=nbase/ncols
        if(mod(nbase,ncols).ne.0) nk=nk+1
        do k=1,nk
          ii=(k-1)*ncols+1
          il=min(nbase,k*ncols)
          write(lfnout,635)
 635      format(//,' Overlap Matrix unnormalized',/)
          write(lfnout,633) (i,i=ii,il)
 633      format(6x,7(6x,i8,6x))
          do j=1,nbase
            write(lfnout,634) j,(sbase(i,j),i=ii,il)
 634        format(i5,1x,10f20.10)
          enddo
          write(lfnout,*) ' '
        enddo
      endif

      if(me.eq.master) then
        do i=1,nbase
          do j=1,i-1
            hbase(i,j)=(1.0d0/(dsqrt(sbase(i,i))*
     &           dsqrt(sbase(j,j))))*hbase(i,j)
            sbase(i,j)=(1.0d0/(dsqrt(sbase(i,i))*
     &           dsqrt(sbase(j,j))))*sbase(i,j)
            hbase(j,i)=hbase(i,j)
            sbase(j,i)=sbase(i,j)
          enddo
        enddo
        do i=1,nbase
          hbase(i,i)=(1.0d0/sbase(i,i))*hbase(i,i)
          sbase(i,i)=(1.0d0/sbase(i,i))*sbase(i,i)
        enddo
        call gronor_print_results(hbase,sbase,nbase)
        if(ncorr.eq.1) then
          allocate( hcorr(nbase,nbase) )
          hcorr = 0.0
          call gronor_correlation_energy(
     &         nwt,mstates,ecorr,nbase,nmol,ncombv,hbase,sbase,hcorr)
          if(ipr.le.40) write(lfnout,636)
 636      format(//,'After applying a shift to the diagonal of H')
          call gronor_print_results(hcorr,sbase,nbase)
          deallocate( hcorr )
        endif
        open(unit=lfnlog,file=fillog(1:inp_strlen(fillog)),
     &       form='formatted',status='unknown',access='append',err=993)
        call timer_stop(99)
        call swatch(date,time)
        write(lfnlog,801) date(1:8),time(1:8),
     &       timer_wall_total(99)-timer_wall_total(98),
     &       timer_wall_total(98),timer_wall_total(99),
     &       nnodes,np,numacc,numnon,nummps,numomp,mgr,
     &       iaslvr,jaslvr,inslvr,jnslvr,
     &       ntask,nbatch,nmol,nbase,mstates,thresh,
     &       user(1:inp_strlen(user)),string(1:inp_strlen(string))
 801    format(a8,1x,a8,3f12.3,4i7,3i3,4i2,i5,i2,3i5,e12.3,1x,a,1x,a)
        write(lfnlog,803) (hbase(i,i),i=1,nbase)
        if(nbase.gt.1) then
          write(lfnlog,803) (hbase(i,i)-hbase(1,1),i=1,nbase)
        endif
 803    format(5x,8f20.10)
        write(lfnlog,804)
 804    format(' ')
        close(unit=lfnlog,status='keep')
        call timer_start(99)
      endif

      deallocate(t,v,ndxtv)
      deallocate(lab,ndx,ig,g)
      deallocate(jspinmol)
      deallocate(idetb,nactb,inactb,ioccb,vecsb,civb)
      deallocate(ioccm,vecsm,civm)
      deallocate(idetm,inactm,nactm,nbasm,ncombv)
      deallocate(map2)

      deallocate(hbase,sbase,tbase,nsing,bpdone)

      deallocate(lstdet,ndxdet)

      call timer_stop(1)

      if(me.eq.master) then
        call timer_stop(99)
        call swatch(date,time)
        write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  End of Hamiltonian calculation'
        flush(lfnday)
        call timer_start(99)
      endif

      if(me.eq.master) call gronor_cml()

      call gronor_timings(lfnout,lfnday)

      if(me.eq.master) deallocate(numrecs)

      close(unit=lfndbg,status='keep')
      if(me.eq.master) then
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
        
        call swatch(date,time)
        write(*,637) trim(date),trim(time)
        if(ipr.ge.0) write(lfnout,637) trim(date),trim(time)
 637    format(/,' Completion of run ',2a10,/)
        
        close(unit=lfnrnk,status='keep')
        close(unit=lfnout,status='keep')
      endif

      if(me.eq.master.and.nalive+1.ne.np) then
        call mpi_abort(MPI_COMM_WORLD,ierror,ierr)
      else
        call mpi_finalize(ierr)
      endif
      
      return
 993  write(lfnout,984) fillog(1:inp_strlen(fillog))
      call errquit(me,2,"Error in main")
 994  write(lfnout,984) filsym(1:inp_strlen(fildbg))
      call errquit(me,2,"Error in main")
 995  write(lfnout,985) filsym(1:inp_strlen(filsym))
      call errquit(me,2,"Error in main")
 996  write(lfnout,986) filout(1:inp_strlen(filout))
      call errquit(me,2,"Error in main")
 997  write(lfnout,987) filvec(1:inp_strlen(filvec))
      call errquit(me,2,"Error in main")
 998  write(lfnout,988) filciv(1:inp_strlen(filciv))
      call errquit(me,2,"Error in main")
 984  format('Unable to open debug file ',a)
 985  format('Unable to open symmetry file ',a)
 986  format('Unable to open output file ',a)
 987  format('Unable to open vects file ',a)
 988  format('Unable to open civec file ',a)
      end
      
