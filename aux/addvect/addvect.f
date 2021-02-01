*************************************************************************
*                                                                       *
*     program addvect                                                   *
*                                                                       *
*     Goal of the program is to generate one formatted                  *
*     vector file out of two seperate fragments                         *
*                                                                       *
*                                                                       *
*-----------------------------------------------------------------------*
*                                                                       *
*     written by:                                                       *
*     Coen de Graaf                                                     *
*     University of Groningen, The Netherlands, Jan. 1996               *
*                                                                       *
*                                                                       *
* Jan 2018:   Adapted for new vector format, only the baab ordering     *
* March 2020: Removed the old vector format; generalized the way in     *
*             which the different orderings are treated; simplified     *
*             the input; added the possibility to split a file in 2     *
*             fragments                                                 *
*                                                                       *
* INPUT DESCRIPTION (free format)                                       *
*                                                                       *
* CARD 1       title for vector file of fragment (A+B)                  *
* CARD 2       filename of vectors fragment A                           *
* CARD 3       filename of vectors fragment B                           *
* CARD 4       order of the vectors of the two fragments                *
*              baab --> occupied fragmnet B, occupied fragment A,       *
*                 virtuals fragment A, virtuals fragment B              *
*              abba, abab, aabb, etc.                                   *
* CARD 5       'add' or 'split' -  superimpose A and B or split into    *
*              A and B. Coefficients are read from INPORB and vecfile1  *
*              and vecfile2 are used to determine the number of basis   *
*              functions and occupied orbitals of the fragments in the  *
*              latter case.                                             *
*                                                                       *
* no special requirements for compilation, just plain gfortran will     *
* do the job (or any other standard compiler), no need for linking.     *
*                                                                       *
*************************************************************************
      module fileInfo
      implicit none
      character(len=80)     :: title
      character(len=25)     :: vecfile1,vecfile2
      end module

      module genInfo
      implicit none
      integer,allocatable    :: nBas(:)
      integer,allocatable    :: nOcc(:)
      integer                :: nIrrep,mxBas,totOrb
      end module

      module decisions
      implicit none
      logical                :: add,symmetry_is_lowered
      end module

      module coefficients
      implicit none
      real(kind=8),allocatable  :: coeff(:,:,:)
      end module
      
* -- End of modules ----

      program main
      use decisions
      implicit none

      character (len=4)    :: order

      call read_input(order)
      call read_vectorfiles
      if ( add ) then
        call superimpose(order)
        write(*,'(A50,A4,A9)') 
     &    'Orbital files of A and B have been merged using a  ',
     &     order,' ordering'
      else
        call split(order)
                write(*,'(A39,A4,A9)')
     &    'Orbital file has been split assuming a ',
     &     order,' ordering'
      end if

      end program main



      subroutine read_input(order)
      use fileInfo
      use decisions
      implicit none
      character(len=4)      :: order        
      character(len=5)      :: split_or_add

      read(5,'(A80)') title
      read(5,'(A25)') vecfile1
      read(5,'(A25)') vecfile2
      read(5,'(A4)') order
      read(5,'(A5)') split_or_add
      call capitalize(order)
      call capitalize(split_or_add)

      if ( title(1:1) .ne. '*') title = '*'//title
      if ( title(2:2) .ne. ' ') title = '* '//title(2:80)

      add = .true.
      if ( split_or_add(1:5) .eq. 'SPLIT' ) add = .false.

      end subroutine read_input



      subroutine read_vectorfiles
      use fileInfo
      use decisions
      implicit none

      character (len=25)  :: filename

      filename = 'INPORB'
      open (10,file='INPORB')
      filename = vecfile1
      open (11,file=filename,status='old',err=401)
      filename = vecfile2
      open (12,file=filename,status='old',err=401)
      
      call read_info
      if ( add ) then
        call read_coeff_add
      else
        call read_coeff_split 
      end if
      call read_label

      goto 402
 401  write(*,*) 'Something went wrong trying to open ',filename
      stop

 402  continue
      end subroutine read_vectorfiles
      


      subroutine read_info
      use genInfo
      use decisions
      implicit none
      character(len=80)     :: dummy
      character(len=132)    :: line
      integer               :: idummy,i,nSym,nBasA,nBasB

      symmetry_is_lowered = .false.

      call find_mark('#INFO ',11)
      read(11,'(A132)') dummy
      read(11,'(2I8)') idummy,nIrrep
      allocate( nBas(2*nIrrep) )
      read(11,'(8I8)') (nBas(i),i=1,nIrrep)
      call find_mark('#INFO ',12)
      read(12,'(A132)') dummy
      read(12,'(A132)') dummy
      read(12,'(8I8)') (nBas(i),i=nIrrep+1,2*nIrrep)
      mxBas = maxval(nBas)
      totOrb = sum(nBas)

      if ( .not. add ) then
        call find_mark('#INFO ',10)
        read(10,'(A132)') dummy
        read(10,'(2I8)') idummy,nSym
        if ( nIrrep .ne. nSym ) then
          if ( nSym .eq. 1 ) then
            symmetry_is_lowered = .true.
          else
            write(*,*) 'Going from ',nIrrep,' to ',nSym,'irreps'
     &                ,' has not been implemented'
            stop
          end if 
        end if
      end if

      if ( .not. symmetry_is_lowered ) nSym = nIrrep

      if ( nSym .gt. 1 ) then
        write(*,*) 'Number of basis functions per irrep'
      else
        write(*,*) 'Number of basis functions'
      end if
      if ( symmetry_is_lowered ) then
        nBasA = 0
        nBasB = 0
        do i = 1, nIrrep
          nbasA = nBasA + nBas(i)
          nBasB = nBasB + nBas(i+nIrrep)
        end do
        write(*,'(A,I5)')'fragment A: ',nBasA
        write(*,'(A,I5)')'         B: ',nbasB
      else
        write(*,'(A,8I5)')'fragment A: ',(nBas(i),i=1,nSym)
        write(*,'(A,8I5)')'         B: ',(nBas(i),i=nSym+1,2*nSym)
      end if
      write(*,*)

      end subroutine read_info


      subroutine read_label
      use genInfo
      use decisions
      implicit none

      character, allocatable :: label(:)
      character(len=80)      :: dummy
      integer                :: i,j,countOcc,nOccA,nOccB
      
      allocate( label(mxBas) ) 
      allocate( nOcc(2*nIrrep) )
      nOcc = 0
      call find_mark('#INDEX',11)
      do i = 1, nIrrep
        if ( nBas(i) .gt. 0 ) then
          label = ' '
          read(11,'(A80)') dummy
          read(11,'(2x,10A)')(label(j),j=1,nBas(i))
        end if
        nOcc(i) = countOcc(label,mxBas)
      end do
      call find_mark('#INDEX',12)
      do i = 1, nIrrep
        if ( nBas(i+nIrrep) .gt. 0 ) then
          label = ' '
          read(12,'(A80)') dummy
          read(12,'(2x,10A)')(label(j),j=1,nBas(i+nIrrep))
        end if
        nOcc(i+nIrrep) = countOcc(label,mxBas)
      end do

      write(*,*) 'Occupied orbitals'
      if ( symmetry_is_lowered ) then
        nOccA = 0
        nOccB = 0
        do i = 1, nIrrep
          nOccA = nOccA + nOcc(i)
          nOccB = nOccB + nOcc(i+nIrrep)
        end do
        write(*,'(A,8I5)') 'fragment A: ',nOccA
        write(*,'(A,8I5)') '         B: ',nOccB
      else
        write(*,'(A,8I5)') 'fragment A: ',(nOcc(i),i=1,nIrrep)
        write(*,'(A,8I5)') '         B: ',(nOcc(i),i=nIrrep+1,2*nIrrep)
      end if
      write(*,*)

      end subroutine read_label


      function countOcc(label,n) result(nOcc)
      implicit none

      character, intent(in)   :: label(n)
      integer                 :: nOcc,j,n

      nOcc = 0
      do j = 1,n
        if (label(j) .eq. 'f') nOcc = nOcc + 1     ! frozen in RASSCF
        if (label(j) .eq. 'i') nOcc = nOcc + 1     ! inactive
        if (label(j) .eq. '1') nOcc = nOcc + 1     ! active (ras1)
        if (label(j) .eq. '2') nOcc = nOcc + 1     ! active (ras2)
        if (label(j) .eq. '3') nOcc = nOcc + 1     ! active (ras3)
      end do

      end function countOcc


      subroutine read_coeff_add
      use genInfo
      use coefficients
      implicit none

      integer            :: i,j,k
      character(len=80)  :: dummy

      allocate( coeff(nIrrep,totOrb,totOrb) )
      coeff = 0.0
      call find_mark('#ORB  ',11)
      do i = 1 , nIrrep
        if ( nBas(i) .gt. 0 ) then
          do j = 1, nBas(i)
            read(11,'(A80)') dummy
            read(11,'(5E22.14)') (coeff(i,j,k),k = 1,nBas(i))
          end do
        end if
      end do
      
      call find_mark('#ORB  ',12)
      do i = 1 , nIrrep
        if ( nBas(i + nIrrep) .gt. 0 ) then
          do j = 1, nBas(i + nIrrep)
            read(12,'(A80)') dummy
            read(12,'(5E22.14)') (coeff(i,j + nBas(i),k + nBas(i)),
     &                                   k = 1,nBas(i + nIrrep))
          end do
        end if
      end do

      end subroutine read_coeff_add


      subroutine read_coeff_split
      use genInfo
      use decisions
      use coefficients
      implicit none

      integer            :: i,j,k,n,nSym
      character(len=80)  :: dummy

      allocate( coeff(nIrrep,totOrb,totOrb) )
      coeff = 0.0
      call find_mark('#ORB  ',10)
      if ( symmetry_is_lowered ) then
        nSym = 1
      else
        nSym = nIrrep
      endif
      do i = 1 , nSym
        if ( symmetry_is_lowered ) then
          n = totOrb
        else
          n = nBas(i) + nBas(i+nIrrep)
        end if
        if ( n .gt. 0 ) then
          do j = 1, n
            read(10,'(A80)') dummy
            read(10,'(5E22.14)') (coeff(i,j,k),k = 1,n)
          end do
        end if
      end do

      end subroutine read_coeff_split



      subroutine find_mark(mark,lUnit)
      implicit none
      integer            :: lUnit,l
      character(len=6)   :: mark
      character(len=132) :: line

      l = len(trim(mark))
      rewind(lUnit)
 101  read(lUnit,'(A132)') line
*      write(*,*) line
      if ( line(1:l) .ne. mark(1:l) ) goto 101
      
      end subroutine find_mark



      subroutine superimpose(order)
      use fileInfo
      use genInfo
      use coefficients
      implicit none

      character (len=4)        :: order
      integer                  :: nBasTot,i,j,k
      logical                  :: occA_done,occB_done
      real(kind=8),allocatable :: occNu(:)
      character,allocatable    :: dummy_label(:)

      open (9,file='AddOrb',status='unknown')
      write(9,'(A11)') '#INPORB 2.2'
      write(9,'(A5)')  '#INFO'
      write(9,'(A80)') title
      write(9,'(3I8)') 0,nIrrep,0
      write(9,'(8I8)') (nBas(i) + nBas(i+nIrrep),i=1,nIrrep)
      write(9,'(8I8)') (nBas(i) + nBas(i+nIrrep),i=1,nIrrep)
      write(9,'(A4)') '#ORB'


      do i = 1, nIrrep
        occA_done = .false.
        occB_done = .false.
        nBasTot = nBas(i) + nBas(i + nIrrep)

        if ( order(1:1) .eq. 'A' ) then
          occA_done = .true.
          if (nOcc(i) .gt. 0 ) then
            do j = 1, nOcc(i)
              write(9,601) '* Occupied A: orbital ',j,' of irrep ',i
              write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
            end do
          end if
        else
          occB_done = .true.
          if ( nOcc(i + nIrrep) .gt. 0 ) then
            do j = 1, nOcc(i + nIrrep)
              write(9,601) '* Occupied B: orbital ',j,' of irrep ',i
              write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
            end do
          end if
        end if

        if ( order(2:2) .eq. 'A' ) then
          if ( occA_done ) then
            if ( nBas(i) - nOcc(i) .gt. 0 ) then
              do j = nOcc(i) + 1, nBas(i)
                write(9,601) '* Virtual A: orbital ',j,' or irrep ',i
                write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
              end do
            end if
          else
            occA_done = .true.
            if ( nOcc(i) .gt. 0 ) then
              do j = 1, nOcc(i)
                write(9,601) '* Occupied A: orbital ',j,' or irrep ',i
                write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
              end do
            end if
          end if
        else
          if ( occB_done ) then
            if ( nBas(i + nIrrep) - nOcc(i + nIrrep) .gt. 0 ) then
              do j = nOcc(i+nIrrep) + 1, nBas(i+nIrrep) 
                write(9,601) '* Virtual B: orbital ',j,' or irrep ',i
                write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
              end do
            end if
          else
            occB_done = .true.
            if ( nOcc(i + nIrrep ) .gt. 0 ) then
              do j = 1, nOcc(i + nIrrep)
                write(9,601) '* Occupied B: orbital ',j,' or irrep ',i
                write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
              end do
            end if
          end if
        end if

        if ( order(3:3) .eq. 'A' ) then
          if ( occA_done ) then
            if ( nBas(i) - nOcc(i) .gt. 0 ) then
              do j = nOcc(i) + 1, nBas(i)
                write(9,601) '* Virtual A: orbital ',j,' or irrep ',i
                write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
              end do
            end if
          else
            occA_done = .true.
            if ( nOcc(i) .gt. 0 ) then
              do j = 1, nOcc(i)
                write(9,601) '* Occupied A: orbital ',j,' or irrep ',i
                write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
              end do
            end if
          end if
        else
          if ( occB_done ) then
            if ( nBas(i + nIrrep) - nOcc(i + nIrrep) .gt. 0 ) then
              do j = nOcc(i+nIrrep) + 1, nBas(i+nIrrep)
                write(9,601) '* Virtual B: orbital ',j,' or irrep ',i
                write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
              end do
            end if
          else
            occB_done = .true.
            if ( nOcc(i + nIrrep ) .gt. 0 ) then
              do j = 1, nOcc(i + nIrrep)
                write(9,601) '* Occupied B: orbital ',j,' or irrep ',i
                write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
              end do
            end if
          end if
        end if

        if ( order(4:4) .eq. 'A' ) then
          if ( nBas(i) - nOcc(i) .gt. 0 ) then
            do j = nOcc(i) + 1, nBas(i)
              write(9,601) '* Virtual A: orbital ',j,' or irrep ',i
              write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
            end do      
          end if
        else
          if ( nBas(i + nIrrep) - nOcc(i + nIrrep) .gt. 0 ) then
            do j = nOcc(i+nIrrep) + 1, nBas(i+nIrrep)
              write(9,601) '* Virtual B: orbital ',j,' or irrep ',i
              write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
            end do
          end if
        endif
      end do
 601  format(2(A,I4))

* These are just dummy arrays, not reflecting real occupation numbers, nor the space
* to which the orbital belongs
      allocate( occNu(totOrb) )
      write(9,'(A4)')'#OCC'
      do i = 1, nIrrep
        nBasTot = nBas(i) + nBas(i + nIrrep)
        if ( nBasTot .gt. 0 ) then
          occNu = 0.0
          do j = 1, nOcc(i) + nOcc(i + nIrrep)
            occNu(j) = 2.0
          end do 
          write(9,'(5E22.14)') (occNu(j),j=1,nBasTot)
        end if
      end do
      allocate( dummy_label(totOrb) )
      write(9,'(A6)')'#INDEX'
      do i = 1, nIrrep
        write(9,'(A12)')'* 1234567890'
        nBasTot = nBas(i) + nBas(i + nIrrep)
        if ( nBasTot .gt. 0 ) then
          dummy_label = 's'
          do j = 1, nOcc(i) + nOcc(i + nIrrep)
            dummy_label(j) = 'i'
          end do
          do j = 1, nBasTot, 10
            if ( j + 9 .le. nBasTot) then
              write(9,602)mod(int(j/10),10),(dummy_label(k),k=j,j+9)
            else
              write(9,602)mod(int(j/10),10),(dummy_label(k),k=j,nBasTot)
            end if
          end do
        end if
      end do
 602  format(I1,x,10A1)

      end subroutine superimpose

      subroutine capitalize(string)
      implicit none
      integer      :: i
      character(*) string

      do i = 1, len(string)
        if (ichar(string(i:i)).gt.96) then
          string(i:i) = char(ichar(string(i:i))-32)
        endif
      end do
      return
      end subroutine capitalize


      subroutine split(order)
      use fileInfo
      use genInfo
      use decisions
      use coefficients
      implicit none

      character (len=4)              :: order
      integer                        :: nBasTot,i,j,k,jj,start,kk,nSym
      logical                        :: occA_done,occB_done
      real(kind=8),allocatable       :: occNu(:)
      character,allocatable          :: dummy_label(:)
      character(len=80),dimension(2) :: vectit,filename


      filename(1) = 'INPORB_A'
      filename(2) = 'INPORB_B'
      vectit(1) = '* Orbitals fragment A from splitting A--B'
      vectit(2) = '* Orbitals fragment B from splitting A--B'

      if ( symmetry_is_lowered ) then
        nSym = 1
        j = 0
        k = 0
        jj = 0
        kk = 0
        do i = 1, nIrrep
          j  = j  + nBas(i)
          k  = k  + nBas(i+nIrrep)
          jj = jj + nOcc(i)
          kk = kk + nOcc(i+nIrrep)
        end do
        nBas(1) = j
        nBas(2) = k
        nOcc(1) = jj
        nOcc(2) = kk
      else
        nSym = nIrrep
      end if
      do j = 1,2
        jj = j + 12
        start = (j - 1)*nSym
        open(jj,file=filename(j),status='unknown')
        write(jj,'(A11)') '#INPORB 2.2'
        write(jj,'(A5)')  '#INFO'
        write(jj,'(A80)') vectit(j)
        write(jj,'(3I8)') 0,nSym,0
        write(jj,'(8I8)') (nBas(i),i=start+1,start+nSym)
        write(jj,'(8I8)') (nBas(i),i=start+1,start+nSym)
        write(jj,'(A4)') '#ORB'
      end do
      
      do i = 1, nSym
        occA_done = .false.
        occB_done = .false.
        start = 0
        nBasTot = nBas(i) + nBas(i+nSym)

        if ( order(1:1) .eq. 'A' ) then
          occA_done = .true.
          if ( nBas(i) .gt. 0 ) then
            do j = 1, nOcc(i)
              write(13,601)'* Occupied A: orbital ',j,' of irrep ',i
              write(13,'(5E22.14)')(coeff(i,j,k),k=1,nBas(i))
            end do
          end if
          start = nOcc(i)
        else
          occB_done = .true.
          if ( nBas(i+nSym) .gt. 0 ) then
            do j = 1, nOcc(i+nSym)
              write(14,601)'* Occupied B: orbital ',j,' of irrep ',i
              write(14,'(5E22.14)')(coeff(i,j,k),k=nBas(i)+1,nBasTot)
            end do
          end if
          start = nOcc(i+nSym)
        end if

        if ( order(2:2) .eq. 'A' ) then
          if ( occA_done ) then
            if ( nBas(i) .gt. 0 ) then
              do j = nOcc(i) + 1, nBas(i)
                write(13,601)'* Virtual A: orbital ',j,' of irrep ',i
                write(13,'(5E22.14)')(coeff(i,j,k),k=1,nBas(i))
              end do
            end if
            start = nBas(i)
          else
            occA_done = .true.
            if ( nBas(i) .gt. 0 ) then
              do j = 1, nOcc(i)
                write(13,601)'* Occupied A: orbital ',j,' of irrep ',i
                write(13,'(5E22.14)')(coeff(i,j+start,k),k=1,nBas(i))
              end do
            end if
            start = start + nOcc(i) 
          end if
        else
          if ( occB_done ) then
            if ( nBas(i+nSym) .gt. 0 ) then
              do j = nOcc(i+nSym) + 1, nBas(i+nSym)
                write(14,601)'* Virtual B: orbital ',j,' of irrep ',i
                write(14,'(5E22.14)')(coeff(i,j,k),k=nBas(i)+1,nBasTot)
              end do
            end if
            start = nBas(i+nSym)
          else
            occB_done = .true.
            if ( nBas(i+nSym) .gt. 0 ) then
              do j = 1, nOcc(i+nSym)
                write(14,601)'* Occupied B: orbital ',j,' of irrep ',i
                write(14,'(5E22.14)')(coeff(i,j+start,k),
     &                                         k=nBas(i)+1,nBasTot)
              end do
            end if
            start = start + nOcc(i+nSym)
          end if
        end if

        if ( order(3:3) .eq. 'A' ) then
          if ( occA_done ) then
            if ( nBas(i) .gt. 0 ) then
              do j = 1, nBas(i) - nOcc(i)
                write(13,601)'* Virtual A: orbital ',j+nOcc(i),
     &                                           ' of irrep ',i
                write(13,'(5E22.14)')(coeff(i,j+start,k),k=1,nBas(i))
              end do
            end if
            start = start + nBas(i) - nOcc(i)
          else
            if ( nBas(i) .gt. 0 ) then
              do j = 1, nOcc(i)
                write(13,601)'* Occupied A: orbital ',j,' of irrep ',i
                write(13,'(5E22.14)')(coeff(i,j+start,k),k=1,nBas(i))
              end do
            end if
            start = start + nOcc(i)
          end if
        else
          if ( occB_done ) then
            if ( nBas(i+nSym) .gt. 0 ) then
              do j = 1, nBas(i+nSym) - nOcc(i+nSym)
                write(14,601)'* Virtual B: orbital ',j+nOcc(i+nSym),
     &                                                 ' of irrep ',i
                write(14,'(5E22.14)')(coeff(i,j+start,k),
     &                                       k=nBas(i)+1,nBasTot)
              end do
            end if
            start = start + nBas(i+nSym) - nOcc(i+nSym)
          else
            if ( nBas(i+nSym) .gt. 0 ) then
              do j = 1, nOcc(i+nSym)
                write(14,601)'* Occupied B: orbital ',j,' of irrep ',i
                write(14,'(5E22.14)')(coeff(i,j+start,k),
     &                                       k=nBas(i)+1,nBasTot)
              end do
            end if
            start = start + nOcc(i+nSym)
          end if
        end if

        if ( order(4:4) .eq. 'A' ) then
          if ( nBas(i) .gt. 0 ) then
            do j = 1, nBas(i) - nOcc(i)
              write(13,601)'* Virtual A: orbital ',j+nOcc(i),
     &                                                   ' of irrep ',i
              write(13,'(5E22.14)')(coeff(i,j+start,k),k=1,nBas(i))
            end do
          end if
        else
          if ( nBas(i+nSym) .gt. 0 ) then
            do j = 1, nBas(i+nSym) - nOcc(i+nSym)
              write(14,601)'* Virtual B: orbital ',j+nOcc(i+nSym),
     &                                                    ' of irrep ',i
              write(14,'(5E22.14)')(coeff(i,j+start,k),
     &                                     k=nBas(i)+1,nBasTot)
            end do
          end if
        end if
      end do
 601  format(2(A,I4))

* These are just dummy arrays, not reflecting real occupation numbers, nor the space
* to which the orbital belongs
      allocate ( occNu(maxval(nBas)) )
      allocate ( dummy_label(maxval(nBas)) )
      do j = 1,2
        jj = j + 12
        start = (j - 1)*nSym
        write(jj,'(A4)')'#OCC'
        write(jj,'(A)')'* OCCUPATION NUMBERS'
        do i = 1, nSym
          nBasTot = nBas(i + start)
          if ( nBasTot .gt. 0 ) then
            occNu = 0.0
            do k = 1, nOcc(i+start) 
              occNu(k) = 2.0
            end do
            write(jj,'(5E22.14)') (occNu(k),k=1,nBasTot)
          end if
        end do
        write(jj,'(A6)')'#INDEX'
        do i = 1, nSym   
          write(jj,'(A12)')'* 1234567890'
          nBasTot = nBas(i + start)
          if ( nBasTot .gt. 0 ) then
            dummy_label = 's'
            do k = 1, nOcc(i + start)
              dummy_label(k) = 'i'
            end do
            do k = 1, nBasTot, 10
              if ( k + 9 .le. nBasTot) then
                write(jj,602)mod(int(k/10),10),
     &                             (dummy_label(kk),kk=k,k+9)
              else
                write(jj,602)mod(int(k/10),10),
     &                             (dummy_label(kk),kk=k,nBasTot)
              end if
            end do
          end if
        end do
      end do
 602  format(I1,x,10A1)



      end subroutine split 
