
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
c      use cuda_functions
#endif
      
      implicit none

      integer (kind=4) :: ierror,ierr,istat,icnt
      integer :: iarg,i,j,idum(51),k,l,len,nn,ifac,imul,scounter
      integer,allocatable :: iocc_tmp(:)
      integer,allocatable :: jspinmol(:)
      integer :: maxspinmul
      integer :: mdum, node
      real (kind=8) :: sdum,ctime,wtime,rmat(2),rmatt(2),rdum(5)
      character (len=255) :: string, filnam, card

      real(kind=8), external :: timer_wall_total
      
      integer :: l1
      logical :: osame

      integer :: getcpucount
      external :: getcpucount

      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl
      integer :: nmax,kbloc

      integer :: igr, numone, numtwo
      integer (kind=4) :: new,newme,newgrp,newrnk
      integer :: maxcib1
      real (kind=8), allocatable :: tc(:,:), hbt(:)

      character (len=4) :: token

      integer :: maxgrp
      integer :: nk,ii,il,ik

      integer :: ibase, jbase, ib, jb, l2, lnxt, lcur, ksr, ln, lsum
      integer :: l3, l4, ltot, ndsum, ioff, nsr(4)
      integer (kind=8), allocatable :: numd(:), ndxd(:), lstdet0(:,:)

      integer jdet, kdet, kdetbatch, ldet, mdet, ndet

      integer                  :: idet,length,ndet_rev
      integer, allocatable     :: sort_occ(:),occu_aux(:)
      real(kind=8),allocatable :: sort_coef(:),civb_aux(:)
      real(kind=8)             :: thres2,max_coef

      real (kind=8) :: rint,rlab,rndx,rh,rs,rt

      integer :: li, ndtot
      real(kind=8) :: p, q

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
       lfnf74=74
       filf74=string(1:inp_strlen(string))//'.f74'
       lfnf75=75
       filf75=string(1:inp_strlen(string))//'.f75'

       open(unit=lfnout,file=filout(1:inp_strlen(filout)),
     &  form='formatted',status='unknown',err=996)

       open(unit=lfnpro,file=filpro(1:inp_strlen(filpro)),
     &  form='formatted',status='unknown',err=996)

       open(unit=lfnday,file=filday(1:inp_strlen(filday)),
     &      form='formatted',status='unknown',err=996)
       
       open(unit=lfntst,file=filtst(1:inp_strlen(filtst)),
     &  form='formatted',status='unknown',err=996)
       
       open(unit=lfncpr,file=filcpr(1:inp_strlen(filcpr)),
     &  form='unformatted',status='unknown',err=996)
       
       open(unit=lfnarx,file=filarx(1:inp_strlen(filarx)),
     &  form='formatted',status='unknown',err=996)

       call timer_stop(99)
       call swatch(date,time)
       write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  Start of job'
  702  format(a8,2x,a8,f12.3,a)
       flush(lfnday)
       call timer_start(99)
           
       write(lfnout,600)
  600  format(/,
     &  ' GronOR: Non-Orthogonal Configuration Interaction program',//,
     &      4x,'T. P. Straatsma',t53,
     &      'C. de Graaf',t85,
     &      'R. W. A. Havenith and R. K. Kathir',/,
     &      4x,'National Center for Computational Sciences',t53,
     &      'Quantum Chemistry Group',t85,
     &      'Department of Theoretical Chemistry',/,
     &      4x,'Oak Ridge National Laboratory',t53,
     &      'University Rovira i Virgili',t85,
     &      'University of Groningen',/,
     &      4x,'Oak Ridge, Tennessee, USA',t53,
     &      'Tarragona, Spain',t85,
     &      'Groningen, the Netherlands',//,
     &      4x,'Based on GNOME written by R. Broer-Braam, ',
     &      'J. Th. van Montfort, and B. Vunderink',//,
     &      4x,'Please cite the following reference when publishing ',
     &      'results obtained using GronOR:',//,
     &      4x,'T. P. Straatsma, R. Broer, S. Faraji,',
     &      'R. W. A. Havenith, L. E. Aguilar Suarez, R. K. Kathir, ',
     &      'M. Wibowo, and C.de Graaf',/,
     &      4x,'"GronOR: Massively parallel and GPU-accelerated ',
     &      'non-orthogonal configuration interaction for large ',
     &      'molecular systems"',/,
     &      4x,'Journal of Chemical Physics, 152 (2020), ',
     &      'https://doi.org/10.1063/1.5141358',/)

       write(lfnout,601) trim(user),getcpucount(),
     &      trim(host),nnodes,
     &      date(1:8),nrsets,
     &      time(1:8),nrnsets,numdev,
     &      np,np/nnodes,np/nrsets,
     &      ncycls,num_threads,trim(command),trim(cwd),
     &      filinp(1:inp_strlen(filinp)),filsym(1:inp_strlen(filsym)),
     &      filout(1:inp_strlen(filout)),filone(1:inp_strlen(filone)),
     &      root(1:inp_strlen(root)),memfre,
     &      root(1:inp_strlen(root)),memtot

 601   format(//,
     &      ' User',t30,a,t60,'CPU count',t100,i10,/,/,
     &      ' Host',t30,a,t60,'Number of nodes',t100,i10,/,
     &      ' Date',t30,a,t60,
     &      'Number of resource sets',t100,i10,/,
     &      ' Time',t30,a,t60,
     &      'Number of resource sets per node',t100,i10,/,
     &      t60,'Number of GPU accelerators per resource set',t100,i10,/,
     &      t60,'Number of MPI ranks',t100,i10,/,
     &      t60,'Number of MPI ranks per node',t100,i10,/,
     &      t60,'Number of MPI ranks per resource set',t100,i10,/,
     &      t60,'Number of rank assignment cycles',t100,i10,/,
     &      t60,'Number of OPENMP threads (OMP_NUM_THREADS)',t100,i10,//,
     &      ' Command argument',t30,a,/,
     &      ' Current working directory',t30,a,//,
     &      ' Input file is',t25,a,t60,
     &      'Symmetry file is',t92,a,/,
     &      ' Output file is',t25,a,t60,
     &      'One electron integral file is',t92,a,//,
     &      ' CI vector file(s) are',t25,a,'_nnn.civ',t60,
     &      'Available memory on device',t90,i20,/,
     &      ' MO vector file(s) are',t25,a,'_nnn.vec',t60,
     &      'Total memory on device',t90,i20,//)

       write(lfnarx,401) trim(user),trim(host),root(1:inp_strlen(root))
 401   format('*** ',a,' on ',a,' : ',a)
       write(lfnarx,402) nnodes,np
 402   format('Nodes',i10,i10)

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
!Set defaults

       npg=0

!     Default is one rank per worker group

       mgr=1
       if(npg.eq.0) npg=(np-1)/mgr
       if(np.lt.mgr+1) then
         write(*,'(a,i5,a,i5,a)') 'Error: Size specified (',mgr,
     &        ') cannot be supported by available ranks (',np,')'
         call errquit()
       endif

       nmol=0

!     mstates will be set depending on Molecule input parameters

       mstates=0
       iload=1

!     Default batch size is 32
       
       nbatch=32

!     Desult is to process integrals for each task entry

       lbatch=0

!     Default is to write a test file
       
       itest=1

!     Development option idevel 0: original makebasestate; 1: most recent makebasestate
!     This may be reset depending on the file check

       idevel=1
       ifault=1

!     Default solvers are those included in the code
       
       isolver=0
       jsolver=0
       iaslvr=0
       jaslvr=0
       inslvr=0
       jnslvr=0

       intadd=0
       naccel=-1
       ncols=7
       nummps=1
       numomp=1
       ixpert=0
       iswsvj=15
       iswevj=15

!     labels   0 : read labels from the integral file
!              1 : use Molcas integral ordering
       
       labels=1
       
#ifdef _OPENACC
       naccel=0
#endif
       nqcp=3

!     Default print flag set to 2
       
       ipr=2

!     Default is to turn timing on

       itim=1

       ipro=-1
       iday=10
       idbg=0
       nspin=0
       mpibuf=168435456
       thresh=1.0d-10
       thresh_CI=1.0d-5
       thresh_SIN=1.0d-12
       tolsvj=1.0d-07
       tolevj=1.0d-07
       lcpr=.false.

       call gronor_input()
      
       if(mgr.eq.0) then
         if(npg.eq.0) call errquit()
         mgr=(np-1)/npg
       endif
      
!     Unless the user is an expert let's check thresholds
      
      if(ixpert.eq.0) then
        thresh=min(thresh,1.0d-05)
        thresh_CI=min(thresh_CI,1.0d-05)
      endif

!     Check if vec/civ or det files are available

      call gronor_check_vectors_and_determinants()

!     For now, we report on the algorithm for generating the MEBFs
!     Furthermore, it need to ensured that reading the root_xxx.det
!     files activates the new makebasestate. Reading civ files must
!     be accompanied by going through the original makebasestate
      
      if(idevel .eq. 0) then
        write(lfnout,'(/,A,/)')
     &     ' Original makebasestate is used to contruct the MEBFs'
        vecdet = .false.
      else
        write(lfnout,'(/,2A,/)')
     &       ' New makebasestate with general spin coupling is',
     &       ' used to contruct the MEBFs'
        vecdet = .true.
      endif
      
      call gronor_number_integrals(numone,numtwo)
      
!     Coen: Warning, numone and numtwo are read from a dummy entry. The real
!           number of two electron integrals is on the next record of the 1-el
!           integral file. Now, they always take the values 8421 and -1. 
!           When npg = 0, something very strange will come out (I guess).
       
      if(npg.eq.0) then
        npg=numtwo/intsrnk
        if(mod(numtwo,intsrnk).gt.0) npg=npg+1
      endif
      
      if(mgr.eq.0) then
        if(iload.eq.0) mgr=np/npg
        if(iload.eq.1) mgr=(np-1)/npg
      else
        if(npg*mgr+1.gt.np) then
          write(lfnout,983)
 983      format(' Requested number of ranks exceeds available')
          call errquit()
        endif
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
      
      if(ifault.eq.0) then
        write(lfnout,627) nummps
 627    format(' No fault tolerance',t60,'
     &       Over-subscription of accelerators (MPS)',t100,i10)
      else
        write(lfnout,625) ifault,nummps
 625    format(' Fault tolerance wait time',t40,i10,' sec',t60,
     &       'Over-subscription of accelerators (MPS)',t100,i10)
      endif
      write(lfnout,626) itest,mgr,idevel,nmol,nbatch,mstates,
     &     lbatch,nbase
 626  format(' Test option',t40,i10,t60,
     &     'Integrals divided over tasks number of ranks',t100,i10,/,
     &     ' Development option',t40,i10,t60,
     &     'Number of molecules (nmol)',t100,i10,/,
     &     ' Batch size',t40,i10,t60,
     &     'Number of molecular states (mstates)',t100,i10,/,
     &     ' Batched option',t40,i10,t60,
     &     'Number of many-electron states (nbase)',t100,i10,/)

      write(lfnout,623) thresh_CI, thresh_SIN
 623  format(' Threshold for CI coefficients',t40,1pe12.3,/,
     &     ' Threshold for ISING',t40,1pe12.3,/)

      write(lfnout,636) nspin+1
 636  format(/,' Spin state (2S+1)',t40,i10,/)
      
      if(idevel .ne. 0) then
        
!     new, experimental code
        call gronor_read_vectors_and_determinants()
      
      else
        
!     standard code begins here, no extra indent
!     to keep everything as it is now

      maxci=0
      maxvec=0
      maxnact=0
      
      if(ipr.ge.1) write(lfnout,604)
 604  format(/,' Dimensions',//,
     &     '  istate   ndetm   nactm  ninatm   nbasm   files',/)
      
      do i=1,mstates
        write(filciv,200) root(1:inp_strlen(root)),i
 200    format(a,'_',i3.3,'.civ')
        len=inp_strlen(filciv)
        open(unit=lfnciv,file=filciv(1:len),form='formatted',
     &       status='old',err=998)
        read(lfnciv,*) idetm(i),nactm(i),inactm(i)
        close(unit=lfnciv)
        write(filvec,201) root(1:inp_strlen(root)),i
 201    format(a,'_',i3.3,'.vec')
        len=inp_strlen(filvec)
        open(unit=lfnvec,file=filvec(1:len),form='formatted',
     &       status='old',err=998)
        read(lfnvec,*) nbasm(i)
        close(unit=lfnvec)
        maxci=max(maxci,idetm(i))
        maxvec=max(maxvec,nbasm(i))
        maxnact=max(maxnact,nactm(i))
        if(ipr.ge.1)
     &       write(lfnout,605) i,idetm(i),nactm(i),inactm(i),nbasm(i),
     &       filciv(1:inp_strlen(filciv)),filvec(1:inp_strlen(filvec))
 605    format(5i8,3x,a,', ',a)
      enddo
      
       if(ipr.ge.1) write(lfnout,606) maxci,maxnact,maxvec
 606   format(8x,'  ------  ------',8x,'  ------',/,
     &  '    max:',i8,i8,8x,i8)

       allocate(civm(maxci,mstates))
       allocate(vecsm(maxvec,maxvec,mstates))
       allocate(ioccm(maxnact,maxci,mstates))
       
       write(lfnout,613) (i,i=1,nbase)
 613   format(/,' Molecular states included in this calculation',//,
     &      ' State        : ',20i4)
       write(lfnout,624)
 624   format(/)
       do i=1,nmol
         write(lfnout,617) i,(ncombv(i,j),j=1,nbase)
 617     format(' Molecule',i4,' : ',20i4)
       enddo
       
       if(ipr.ge.3) write(lfnout,607)
 607   format(//,'  istate   idetm           ci coefficient')
       
!     Coen Oct17_19: changed idum for iocc_tmp, idum is used later for
!     something completely different.
       
       allocate( iocc_tmp(maxnact) )
       do i=1,mstates
         iocc_tmp = 0
         write(filciv,200) root(1:inp_strlen(root)),i
         len=inp_strlen(filciv)
         open(unit=lfnciv,file=filciv(1:len),form='formatted',
     &        status='old',err=998)
         read(lfnciv,*)
         do j=1,idetm(i)
           read(lfnciv,1002) nn,civm(j,i),(iocc_tmp(k),k=1,nactm(i))
 1002      format(i7,f25.14,6x,20(i2,x))
           if(ipr.ge.3) then
             if(j.eq.1) then
               write(lfnout,608) i,nn,civm(j,i),
     &              (iocc_tmp(k),k=1,nactm(i))
 608           format(/,2i8,f25.14,6x,20(i2.2,1x))
             else
               write(lfnout,609) nn,civm(j,i),
     &              (iocc_tmp(k),k=1,nactm(i))
 609           format(8x,i8,f25.14,6x,20(i2.2,1x))
             endif
           endif
           do k=1,nactm(i)
             if(iocc_tmp(k).eq.11) then
               ioccm(k,j,i)=2
             elseif(iocc_tmp(k).eq.1) then
               ioccm(k,j,i)=-1
             elseif(iocc_tmp(k).eq.10) then
               ioccm(k,j,i)=1
             elseif(iocc_tmp(k).eq.0) then
               ioccm(k,j,i)=0
             else
               print *,(iocc_tmp(l),l=1,nactm(i))
               stop 'weird occupation'
             endif
           enddo
         enddo
         close(unit=lfnciv)
         
         write(filvec,201) root(1:inp_strlen(root)),i
         len=inp_strlen(filvec)
         open(unit=lfnvec,file=filvec(1:len),form='formatted',
     &        status='old',err=998)
         read(lfnvec,*)
         if(ipr.ge.4) then
           write(lfnout,2002) filvec(1:len),nbasm(i)
 2002      format(/,' Vector File ',a,//,
     &          ' Number of basis functions ',i8,/)
         endif
         do j=1,nbasm(i)
           read(lfnvec,1003) (vecsm(k,j,i),k=1,nbasm(i))
 1003      format(4F18.14)
           if(ipr.ge.4) then
             write(lfnout,2003) j,(vecsm(k,j,i),k=1,nbasm(i))
 2003        format(i5,5f25.14,/,(5x,5f25.14))
           endif
         enddo
         close(unit=lfnvec)
         flush(lfnout)
       enddo
       deallocate (iocc_tmp)
       write(filsym,202) root(1:inp_strlen(root))
 202   format(a,'.sym')
       len=inp_strlen(filsym)
       open(unit=lfnsym,file=filsym(1:len),form='formatted',
     &      status='old',err=995)
       
       if(ipr.ge.3) write(lfnout,610) filsym(1:inp_strlen(filsym))
 610   format(/,' System information read from ',a,/)

       
      endif
      
!     new and old code come together again      

!
!     Coen: remove SYMOL from the list of QC programs
!     if (nqcp.eq.1) then
!     call gronor_gnome_symmetry_input(lfnsym,lfnout,ipr)
!     elseif (nqcp.eq.2) then
       
       if ( nqcp .eq. 2 ) then
         call gronor_gnome_symmetry_input_gamess(lfnsym,lfnout,ipr)
       elseif (nqcp.eq.3) then
         call gronor_gnome_molcas_input(lfnsym,lfnout,lfnarx,ipr) ! OpenMolcas
       endif
       close(unit=lfnsym)
       
      else                      ! belonging to "if(inp_read()) then" (I guess)
        lfnout=6
        write(filout,650) me
 650    format('gronor_',i5.5,'.dbg')
        if(ipr.ge.10) then
          open(unit=lfnout,file=filout(1:inp_strlen(filout)),
     &         form='formatted',status='unknown',err=996)
        endif
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
          idum(23)=iload
          idum(24)=nbatch
          idum(25)=mgr
          idum(26)=itest
          idum(27)=naccel
          idum(28)=inpcib
          idum(29)=ipro
          idum(30)=intadd
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
          idum(44)=lbatch
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
        iload  =idum(23)
        nbatch =idum(24)
        mgr    =idum(25)
        itest  =idum(26)
        naccel =idum(27)
        inpcib =idum(28)
        ipro   =idum(29)
        intadd =idum(30)
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
        lbatch =idum(44)
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
        call MPI_Bcast(nqcp,1,MPI_INTEGER8,master,
     &       MPI_COMM_WORLD,ierr)
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
     &           ioccm(j,1,i)
         enddo
       enddo


       if(idevel .eq. 0) then
!     again no indent to keep the code as it is
         
      maxspinmul=0
      do i=1,nmol
         if(maxspinmul.lt.jspinmol(ncombv(i,1))) then
            maxspinmul=jspinmol(ncombv(i,1))
         endif
      enddo
      maxspinmul=maxspinmul+1
      maxcib=0
      do i=1,nbase
        maxcib1=1
        do j=1,nmol
          if(mod(jspinmol(ncombv(j,i))+1,2).eq.0) then
            scounter=((jspinmol(ncombv(j,i))+1)/2)-1
          else
            scounter=jspinmol(ncombv(j,i))/2
          endif
          maxcib1=maxcib1*idetm(ncombv(j,i))
          imul=jspinmol(ncombv(j,i))+1
          do k=1,scounter
            maxcib1=maxcib1*((imul)*(imul))/nactm(ncombv(i,1))
            imul=imul-2
          enddo
        enddo
        maxcib=max(maxcib,maxcib1)
      enddo
      maxcib=2*maxspinmul*maxcib
      if (me.eq.master) then
        write(lfnout,629) maxcib
 629    format(/,' Recommended value of maxcib is ',i10)
        flush(lfnout)
      endif
      if(inpcib.ne.0) maxcib=inpcib
      
      allocate(civb(maxcib,nbase))
      allocate(vecsb(nbasis,nbasis,nbase))
      allocate(ioccb(maxcib,nbase))
      allocate(inactb(nbase))
      allocate(nactb(nbase))
      allocate(idetb(nbase))
      

      else
!     alternative code starts here

        allocate(vecsb(nbasis,nbasis,nbase))
        allocate(inactb(nbase))
        allocate(nactb(nbase))
        allocate(idetb(nbase))      

        call gronor_determine_nci()
        call gronor_determine_maxcib()

        allocate(ioccb(maxcib,nbase))
        allocate(civb(maxcib,nbase))
        if (me.eq.master) then
          write(lfnout,637) maxcib
 637      format(/,' Calculated value of maxcib is ',i10)
          flush(lfnout)
          call timer_stop(99)
          call swatch(date,time)
          write(lfnday,1704) date(1:8),time(1:8),timer_wall_total(99),
     &         '  :  Maxcib calculated value is           ',maxcib
 1704     format(a8,2x,a8,f12.3,a,i16)
          flush(lfnday)
          call timer_start(99)
        endif

      endif
!     old and new come together again
      
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
        call errquit()
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
        write(lfnout,630)
 630    format(//,' MPI Rank Assignments')
 635    format(/,' Accelerated Groups')
        if(ipr.lt.2) then
        k=1
        l=min(k+14,numacc)
        if(numacc.gt.0) then
          write(lfnout,635)
          do while (k.le.l)
            write(lfnout,634)
 634        format(/)
            write(lfnout,631) (j,j=k,l)
 631        format(10x,15i8)
            write(lfnout,632) (' ----',j=k,l)
 632        format(10x,15a8)
            do i=2,mgr+1
              write(lfnout,631) (allgroups(j,i),j=k,l)
            enddo
            k=l+1
            l=min(k+14,numacc)
          enddo
          write(lfnout,634)
        endif
        if(numnon.gt.0) then
          write(lfnout,633)
 633      format(/,' Non-Accelerated Groups')
          k=l+1
          l=min(k+14,numgrp)
          do while (k.le.l)
            write(lfnout,634)
            write(lfnout,631) (j,j=k,l)
            write(lfnout,632) (' ----',j=k,l)
            do i=2,mgr+1
              write(lfnout,631) (allgroups(j,i),j=k,l)
            enddo
            k=l+1
            l=min(k+14,numgrp)
          enddo
        endif
       else
        k=1
        l=min(k+7,numacc)
        if(numacc.gt.0) then
          write(lfnout,635)
          do while (k.le.l)
            write(lfnout,634)
            write(lfnout,1631) (j,j=k,l)
 1631       format(10x,8i16)
 1633       format(10x,8(i14,'/',i1))
 2633       format(10x,8(i13,'/',i2))
 3633       format(10x,8(i12,'/',i3))
 4633       format(10x,8(i11,'/',i4))
 5633       format(10x,8(i8,'/',i7))
            write(lfnout,1632) (' ------------',j=k,l)
 1632        format(10x,8a16)
            do i=2,mgr+1
              if(np.lt.10) then
                write(lfnout,1633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              elseif(np.lt.100) then
                write(lfnout,2633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              elseif(np.lt.1000) then
                write(lfnout,3633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              elseif(np.lt.10000) then
                write(lfnout,4633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              else
                write(lfnout,5633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              endif
            enddo
            k=l+1
            l=min(k+7,numacc)
          enddo
          write(lfnout,634)
        endif
        if(numnon.gt.0) then
          write(lfnout,633)
          k=l+1
          l=min(k+7,numgrp)
          do while (k.le.l)
            write(lfnout,634)
            write(lfnout,1631) (j,j=k,l)
            write(lfnout,1632) (' ------------',j=k,l)
            do i=2,mgr+1
              if(np.lt.10) then
                write(lfnout,1633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              elseif(np.lt.100) then
                write(lfnout,2633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              elseif(np.lt.1000) then
                write(lfnout,3633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              elseif(np.lt.10000) then
                write(lfnout,4633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              else
                write(lfnout,5633) ((allgroups(j,i),
     &               map2(allgroups(j,i)+1,4)),j=k,l)
              endif
            enddo
            k=l+1
            l=min(k+7,numgrp)
          enddo
        endif
       endif
      endif

      if(me.eq.master) flush(lfnout)
      
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
        if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)
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
      
      
      if(me.eq.master) then
        write(lfnout,622) iaslvr,jaslvr,inslvr,jnslvr,numomp
 622    format(//,' Solver options accelerated',t40,2i5,/,
     &       ' Solver options non-accelerated',t40,2i5,t60,
     &       'Number of OpenMP threads per MPI (OMP)',t100,i10)
        write(lfnout,611) nbase,nmol,nbasis
 611    format(//,' Dimensions',//,
     &       ' Number of states (nbase)',t40,i10,t60,
     &       'Number of molecules (nmol)', t100,i10,/,
     &       ' Number of basis functions (nbasis)',t40,i10)
        flush(lfnout)
      endif
      
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

        if(idevel.eq.0) then
          call gronor_old_makebasestate
     &         (lfnout,civm,vecsm,ioccm,nbasm,nactm,inactm,
     &         idetm,maxci,maxvec,maxnact,mstates,civb(1,i),
     &         vecsb(1,1,i),ioccb(1,i),inactb(i),nactb(i),
     &         idetb(i),maxcib,nbasis,mactb,
     &         nmol,ncombv(1,i),nspin,thresh_CI)
        else
          if(ipr.gt.2.and.me.eq.master) then
            write(lfnout,*) ' Calling the new make_basestate'
            flush(lfnout)
          endif
          call gronor_make_basestate(i)
!        call gronor_makebasestate_remove
!     &       (lfnout,civm,vecsm,ioccm,nbasm,nactm,inactm,
!     &       idetm,maxci,maxvec,maxnact,mstates,civb(1,i),vecsb(1,1,i),
!     &       ioccb(1,i),inactb(i),nactb(i),idetb(i),maxcib,nbasis,mactb,
!     &         nmol,ncombv(1,i),nspin,thresh_CI)
        
        endif
        
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
        if ( max_coef .lt. abs(civb(1,i)) ) then
          max_coef = abs(civb(1,i))
        endif
        
        if(me.eq.master.and.ipr.ge.3) then
          write(lfnout,612) i
 612      format(/,' Base state',i4,' completed')
          flush(lfnout)
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
          if ( abs(civb(idet,ibase)) .gt. thres2 ) then
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

! Following Bcasts should not be necessary, but added for safety

!      call MPI_Bcast(ndxdet,nbase*nbase,MPI_INTEGER8,master,
!     &     MPI_COMM_WORLD,ierr)
!      call MPI_Bcast(lstdet,numdet*2,MPI_INTEGER8,master,
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
        write(lfnout,618) np,numgrp,np-numgrp*mgr-1,mgr
 618    format(/,' Process Distribution',//,
     &       ' Number of processes',t40,i10,t60,
     &       'Number of process groups',t100,i10,/,
     &       ' Number of idle processes',t40,i10,t60,
     &       'Number of processes per group',t100,i10)
      endif
      
      call timer_stop(2)
      call timer_start(3)
      
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
        write(lfnout,614) nbas,rndx
        write(lfnout,615) int1,rlab
        write(lfnout,616) int2,rint
 614    format(' Number of basisfunctions',t50,i16,
     &       t80,' Size of index arrays',t115,f8.3,' GB')
 615    format(' Number of one-electron integrals',t50,i16,
     &       t80,' Size of integral labels',t115,f8.3,' GB')
 616    format(' Number of two-electron integrals',t50,i16,
     &       t80,' Size of two-electron integrals',t115,f8.3,' GB')
      endif
      
      call timer_stop(3)
      call timer_start(4)
      
      if(me.eq.master.and.ipr.ge.5) then
        write(lfnout,628) nbase,nbase*(nbase+1)/2
 628    format(' Number of base states',t50,i16,/,                       &
     &       ' Number of base state elements',t50,i16)
      endif
      
      if(me.eq.master) flush(lfnout)
      
      if(me.eq.master) then
        allocate(numrecs(np))
        do i=1,np
          numrecs(i)=0
        enddo
      endif
      
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
     &         ' NVIDIA ',numdev,mydev,iamacc,nummps,
     &           (map2(me+1,i),i=1,5)
          flush(lfndbg)
        endif

      call timer_start(98)
      if(me.eq.master.and.iload.eq.1) then
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

      if(me.eq.master) then
       write(lfnout,674)
 674   format(//,' Hamiltonian Matrix unnormalised',/)
       open(unit=lfnf74,file=filf74,form='unformatted',status='unknown')
       write(74)nbase
       write(74)((hbase(j,i),i=1,j),j=1,nbase)
       write(74)((sbase(j,i),i=1,j),j=1,nbase)
       close(unit=lfnf74,status='keep')

       nk=nbase/ncols
       if(mod(nbase,ncols).ne.0) nk=nk+1
       do k=1,nk
         ii=(k-1)*ncols+1
         il=min(nbase,k*ncols)
         write(lfnout,671) (i,i=ii,il)
         do j=1,nbase
           write(lfnout,672) j,(hbase(i,j),i=ii,il)
         enddo
         write(lfnout,*) ' '
       enddo

       nk=nbase/ncols
       if(mod(nbase,ncols).ne.0) nk=nk+1
       do k=1,nk
         ii=(k-1)*ncols+1
         il=min(nbase,k*ncols)
         write(lfnout,675)
 675     format(//,' Overlap Matrix unnormalised',/)
         write(lfnout,671) (i,i=ii,il)
 671     format(6x,7(6x,i8,6x))
         do j=1,nbase
           write(lfnout,672) j,(sbase(i,j),i=ii,il)
 672       format(i5,1x,10f20.10)
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
        if ( ncorr .eq. 1 ) then
          allocate( hcorr(nbase,nbase) )
          hcorr = 0.0
          call gronor_correlation_energy(
     &         nwt,mstates,ecorr,nbase,nmol,ncombv,hbase,sbase,hcorr)
          write(lfnout,*)
          write(lfnout,*)
          write(lfnout,*) 'After applying a shift to the diagonal of H'
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
     &       nbatch,lbatch,nmol,nbase,mstates,thresh,
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

      call timer_stop(4)
      call timer_stop(1)

      if(me.eq.master) then
       call timer_stop(99)
       call swatch(date,time)
       write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  End of Hamiltonian calculation'
       flush(lfnday)
       call timer_start(99)
      endif

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
       write(*,602) trim(date),trim(time)
       write(lfnout,602) trim(date),trim(time)
 602   format(/,' Completion of run ',2a10,/)
       
       close(unit=lfnout,status='keep')
      else
!       if(ipr.lt.5) then
!        close(unit=lfnout,status='delete')
!       else
!        close(unit=lfnout,status='keep')
!       endif
!       close(unit=lfnpro,status='delete')
      endif

      if(me.eq.master.and.nalive+1.ne.np) then
        call mpi_abort(MPI_COMM_WORLD,ierror,ierr)
      else
        call mpi_finalize(ierr)
      endif
      
      return
 993  write(lfnout,984) fillog(1:inp_strlen(fillog))
      call errquit()
 994  write(lfnout,984) filsym(1:inp_strlen(fildbg))
      call errquit()
 995  write(lfnout,985) filsym(1:inp_strlen(filsym))
      call errquit()
 996  write(lfnout,986) filout(1:inp_strlen(filout))
      call errquit()
 997  write(lfnout,987) filvec(1:inp_strlen(filvec))
      call errquit()
 998  write(lfnout,988) filciv(1:inp_strlen(filciv))
      call errquit()
 984  format('Unable to open debug file ',a)
 985  format('Unable to open symmetry file ',a)
 986  format('Unable to open output file ',a)
 987  format('Unable to open vects file ',a)
 988  format('Unable to open civec file ',a)
      end
      
