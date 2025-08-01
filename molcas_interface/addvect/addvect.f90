!*************************************************************************
!*                                                                       *
!*     program addvect                                                   *
!*                                                                       *
!*     Goal of the program is to generate one formatted                  *
!*     vector file out of two seperate fragments                         *
!*                                                                       *
!*                                                                       *
!*-----------------------------------------------------------------------*
!*                                                                       *
!*     written by:                                                       *
!*     Coen de Graaf                                                     *
!*     University of Groningen, The Netherlands, Jan. 1996               *
!*                                                                       *
!*                                                                       *
!* Jan 2018:   Adapted for new vector format, only the baab ordering     *
!* March 2020: Removed the old vector format; generalized the way in     *
!*             which the different orderings are treated; simplified     *
!*             the input; added the possibility to split a file in 2     *
!*             fragments                                                 *
!* April 2022: Added the possibility to divide the vectors in three      *
!*             subspaces: inactive - active - virtual                    *
!*                                                                       *
!* INPUT DESCRIPTION (free format)                                       *
!*                                                                       *
!* CARD 1       title for vector file of fragment (A+B)                  *
!* CARD 2       filename of vectors fragment A                           *
!* CARD 3       filename of vectors fragment B                           *
!* CARD 4       order of the vectors of the two fragments                *
!*              baab --> occupied fragmnet B, occupied fragment A,       *
!*                 virtuals fragment A, virtuals fragment B              *
!*              abba, abab, aabb, etc.                                   *
!*                or                                                     *
!*              abbaab --> inactive A, inactive B, active B, active B,   *
!*                         virtual A, virtual B                          *
!*              aabbab, ababab, etc.                                     *
!* CARD 5       'add' or 'split' -  superimpose A and B or split into    *
!*              A and B. Coefficients are read from INPORB and vecfile1  *
!*              and vecfile2 are used to determine the number of basis   *
!*              functions and occupied orbitals of the fragments in the  *
!*              latter case.                                             *
!*                                                                       *
!* no special requirements for compilation, just plain gfortran will     *
!* do the job (or any other standard compiler), no need for linking.     *
!*                                                                       *
!*************************************************************************
module addvect_fileInfo
  implicit none
  character(len=80)     :: title
  character(len=25)     :: vecfile1,vecfile2
end module

module addvect_genInfo
  implicit none
  integer,allocatable    :: nBas(:)
  integer,allocatable    :: nOcc(:),nAct(:)
  integer                :: nIrrep,mxBas,totOrb
end module

module addvect_decisions
  implicit none
  logical                :: add,symmetry_is_lowered
end module

module addvect_coefficients
  implicit none
  real(kind=8),allocatable  :: coeff(:,:,:)
end module
      
!* -- End of modules ----

program main
use addvect_decisions
implicit none

external :: addvect_read_input,addvect_read_vectorfiles
external :: addvect_superimpose,addvect_split
external :: addvect_superimpose2

character (len=6)    :: order

call addvect_read_input(order)
call addvect_read_vectorfiles
if (len(trim(order)).eq.4) then
  if ( add ) then
    call addvect_superimpose(order)
    write(*,'(A50,A4,A9)')'Orbital files of A and B have been merged using a  ',order,' ordering'
  else
    call addvect_split(order)
    write(*,'(A39,A4,A9)')'Orbital file has been split assuming a ',order,' ordering'
  end if
  else  
  if ( add ) then
    call addvect_superimpose2(order)
    write(*,'(A50,A6,A9)')'Orbital files of A and B have been merged using a  ',order,' ordering'
  else
    write(*,*) 'not yet implemented'
  endif
endif

end program main



subroutine addvect_read_input(order)
use addvect_fileInfo
use addvect_decisions
implicit none

external :: addvect_capitalize

character(len=6)      :: order        
character(len=5)      :: split_or_add

read(5,'(A80)') title
read(5,'(A25)') vecfile1
read(5,'(A25)') vecfile2
read(5,*) order
read(5,'(A5)') split_or_add
call addvect_capitalize(order)
call addvect_capitalize(split_or_add)

if ( title(1:1) .ne. '*') title = '*'//title(1:79)
if ( title(2:2) .ne. ' ') title = '* '//title(2:80)

add = .true.
if ( split_or_add(1:5) .eq. 'SPLIT' ) add = .false.

if ( len(trim(order)) .eq. 4 ) then
  write(*,*) 'ordering occupied - virtuals'
elseif ( len(trim(order)) .eq. 6 ) then
  write(*,*) 'ordering inactive - active - virtuals'
else
  write(*,*) trim(order),' not implemented'
endif 

end subroutine addvect_read_input



subroutine addvect_read_vectorfiles
use addvect_fileInfo
use addvect_decisions
implicit none

external :: addvect_read_info,addvect_read_coeff_add
external :: addvect_read_coeff_split,addvect_read_label

character (len=25)  :: filename

filename = 'INPORB'
open (10,file='INPORB')
filename = vecfile1
open (11,file=filename,status='old',err=401)
filename = vecfile2
open (12,file=filename,status='old',err=401)

call addvect_read_info
if ( add ) then
  call addvect_read_coeff_add
else
  call addvect_read_coeff_split 
end if
call addvect_read_label

goto 402
401  write(*,*) 'Something went wrong trying to open ',filename
stop

402  continue
end subroutine addvect_read_vectorfiles
      


subroutine addvect_read_info
use addvect_genInfo
use addvect_decisions
implicit none

external :: addvect_find_mark

character(len=80)     :: dummy
integer               :: idummy,i,nSym,nBasA,nBasB

symmetry_is_lowered = .false.

call addvect_find_mark('#INFO ',11)
read(11,'(A132)') dummy
read(11,'(2I8)') idummy,nIrrep
allocate( nBas(2*nIrrep) )
read(11,'(8I8)') (nBas(i),i=1,nIrrep)
call addvect_find_mark('#INFO ',12)
read(12,'(A132)') dummy
read(12,'(A132)') dummy
read(12,'(8I8)') (nBas(i),i=nIrrep+1,2*nIrrep)
mxBas = maxval(nBas)
totOrb = sum(nBas)

if ( .not. add ) then
  call addvect_find_mark('#INFO ',10)
  read(10,'(A132)') dummy
  read(10,'(2I8)') idummy,nSym
  if ( nIrrep .ne. nSym ) then
    if ( nSym .eq. 1 ) then
      symmetry_is_lowered = .true.
    else
      write(*,*) 'Going from ',nIrrep,' to ',nSym,'irreps has not been implemented'
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

end subroutine addvect_read_info


subroutine addvect_read_label
use addvect_genInfo
use addvect_decisions
implicit none

external :: addvect_find_mark

character, allocatable :: label(:)
character(len=80)      :: dummy
integer                :: i,j,countOcc,nOccA,nOccB
integer                :: nActA,nActB,countAct

allocate( label(mxBas) ) 
allocate( nOcc(2*nIrrep) )
allocate( nAct(2*nIrrep) )
nOcc = 0
call addvect_find_mark('#INDEX',11)
do i = 1, nIrrep
  if ( nBas(i) .gt. 0 ) then
    label = ' '
    read(11,'(A80)') dummy
    read(11,'(2x,10A)')(label(j),j=1,nBas(i))
  end if
  nOcc(i) = countOcc(label,mxBas)
  nAct(i) = countAct(label,mxBas)
end do
call addvect_find_mark('#INDEX',12)
do i = 1, nIrrep
  if ( nBas(i+nIrrep) .gt. 0 ) then
    label = ' '
    read(12,'(A80)') dummy
    read(12,'(2x,10A)')(label(j),j=1,nBas(i+nIrrep))
  end if
  nOcc(i+nIrrep) = countOcc(label,mxBas)
  nAct(i+nIrrep) = countAct(label,mxBas)
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

write(*,*) 'Active orbitals'
if ( symmetry_is_lowered ) then
  nActA = 0
  nActB = 0
  do i = 1, nIrrep
    nActA = nActA + nAct(i)
    nActB = nActB + nAct(i+nIrrep)
  end do
  write(*,'(A,8I5)') 'fragment A: ',nActA
  write(*,'(A,8I5)') '         B: ',nActB
else
  write(*,'(A,8I5)') 'fragment A: ',(nAct(i),i=1,nIrrep)
  write(*,'(A,8I5)') '         B: ',(nAct(i),i=nIrrep+1,2*nIrrep)
end if
write(*,*)



end subroutine addvect_read_label


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



function countAct(label,n) result(nAct)
implicit none

character, intent(in)   :: label(n)
integer                 :: nAct,j,n

nAct = 0
do j = 1,n
  if (label(j) .eq. '1') nAct = nAct + 1     ! active (ras1)
  if (label(j) .eq. '2') nAct = nAct + 1     ! active (ras2)
  if (label(j) .eq. '3') nAct = nAct + 1     ! active (ras3)
end do

end function countAct



subroutine addvect_read_coeff_add
use addvect_genInfo
use addvect_coefficients
implicit none

external :: addvect_find_mark

integer            :: i,j,k
character(len=80)  :: dummy

allocate( coeff(nIrrep,totOrb,totOrb) )
coeff = 0.0
call addvect_find_mark('#ORB  ',11)
do i = 1 , nIrrep
  if ( nBas(i) .gt. 0 ) then
    do j = 1, nBas(i)
      read(11,'(A80)') dummy
      read(11,'(5E22.14)') (coeff(i,j,k),k = 1,nBas(i))
    end do
  end if
end do

call addvect_find_mark('#ORB  ',12)
do i = 1 , nIrrep
  if ( nBas(i + nIrrep) .gt. 0 ) then
    do j = 1, nBas(i + nIrrep)
      read(12,'(A80)') dummy
      read(12,'(5E22.14)') (coeff(i,j + nBas(i),k + nBas(i)),k = 1,nBas(i + nIrrep))
    end do
  end if
end do

end subroutine addvect_read_coeff_add


subroutine addvect_read_coeff_split
use addvect_genInfo
use addvect_decisions
use addvect_coefficients
implicit none

external :: addvect_find_mark

integer            :: i,j,k,n,nSym
character(len=80)  :: dummy

allocate( coeff(nIrrep,totOrb,totOrb) )
coeff = 0.0
call addvect_find_mark('#ORB  ',10)
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

end subroutine addvect_read_coeff_split



subroutine addvect_find_mark(mark,lUnit)
implicit none
integer            :: lUnit,l
character(len=6)   :: mark
character(len=132) :: line

l = len(trim(mark))
rewind(lUnit)
101  read(lUnit,'(A132)') line
if ( line(1:l) .ne. mark(1:l) ) goto 101
      
end subroutine addvect_find_mark



subroutine addvect_superimpose(order)
use addvect_fileInfo
use addvect_genInfo
use addvect_coefficients
implicit none

character (len=6)        :: order
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
if (order .eq. 'ABAB') then
  write(9,'(a)')'* Occupied and virtual orbitals of the fragments'
  write(9,'(8I8)')(nOcc(i),i=1,2*nIrrep)
  write(9,'(8I8)')(nBas(i) - nOcc(i),i=1,2*nIrrep)
end if
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

!* These are just dummy arrays, not reflecting real occupation numbers, nor the space
!* to which the orbital belongs
allocate( occNu(totOrb) )
write(9,'(A4)')'#OCC'
write(9,'(A20)')'* OCCUPATION NUMBERS'
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

end subroutine addvect_superimpose


subroutine addvect_superimpose2(order)
use addvect_fileInfo
use addvect_genInfo
use addvect_coefficients
implicit none

character (len=6)        :: order
integer                  :: nBasTot,i,j,k,jj
logical                  :: inactA_done,inactB_done
logical                  :: actA_done,actB_done
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
  inactA_done = .false.
  inactB_done = .false.
  actA_done = .false.
  actB_done = .false.
  nBasTot = nBas(i) + nBas(i + nIrrep)

  if ( order(1:1) .eq. 'A' ) then
    inactA_done = .true.
    if (nOcc(i) - nAct(i) .gt. 0 ) then
      do j = 1, nOcc(i) - nAct(i)
        write(9,601) '* Inactive A: orbital ',j,' of irrep ',i
        write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
      end do
    end if
  else
    inactB_done = .true.
    if ( nOcc(i + nIrrep) - nAct(i + nIrrep).gt. 0 ) then
      do j = 1, nOcc(i + nIrrep)-nAct(i + nIrrep)
        write(9,601) '* Inactive B: orbital ',j,' of irrep ',i
        write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
      end do
    end if
  end if

  if ( order(2:2) .eq. 'A' ) then
    if ( inactA_done ) then
      actA_done = .true.
      if ( nAct(i) .gt. 0 ) then
        do j = 1, nAct(i)
          jj = nOcc(i) - nAct(i) + j
          write(9,601) '* Active A: orbital ',jj,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,jj,k),k=1,nBasTot)
        end do
      end if
    else
      inactA_done = .true.
      if ( nOcc(i) - nAct(i) .gt. 0 ) then
        do j = 1, nOcc(i) - nAct(i)
          write(9,601) '* Inactive A: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
        end do
      end if
    end if
  else
    if ( inactB_done ) then
      actB_done = .true.
      if ( nAct(i + nIrrep) .gt. 0 ) then
        do j = 1, nAct(i + nIrrep) 
          jj = nOcc(i + nIrrep) - nAct(i + nIrrep) + j
          write(9,601) '* Active B: orbital ',jj,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,jj+nBas(i),k),k=1,nBasTot)
        end do
      end if
    else
      inactB_done = .true.
      if ( nOcc(i + nIrrep ) - nAct(i + nIrrep) .gt. 0 ) then
        do j = 1, nOcc(i + nIrrep) - nAct(i + nIrrep)
          write(9,601) '* Inactive B: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
        end do
      end if
    end if
  end if

  if ( order(3:3) .eq. 'A' ) then
    if ( inactA_done .and. actA_done) then
      if ( nBas(i) - nOcc(i) .gt. 0 ) then
        do j = nOcc(i) + 1, nBas(i)
          write(9,601) '* Virtual A: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
        end do
      end if
    elseif ( inactA_done ) then
      actA_done = .true.
      if ( nAct(i) .gt. 0 ) then
        do j = 1, nAct(i) 
          jj = nOcc(i) - nAct(i) + j
          write(9,601) '* Active A: orbital ',jj,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,jj,k),k=1,nBasTot)
        end do
      end if
    else
      inactA_done = .true.
      if ( nOcc(i) - nAct(i) .gt. 0 ) then
        do j = 1, nOcc(i) - nAct(i)
          write(9,601) '* Inactive A: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
        end do
      end if
    end if
  else
    if ( inactB_done .and. actB_done ) then
      if ( nBas(i + nIrrep) - nOcc(i + nIrrep) .gt. 0 ) then
        do j = nOcc(i+nIrrep) + 1, nBas(i+nIrrep)
          write(9,601) '* Virtual B: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
        end do
      end if
    elseif ( inactB_done ) then
      actB_done = .true.
      if ( nAct(i + nIrrep) .gt. 0 ) then
        do j = 1, nAct(i + nIrrep)
          jj = nOcc(i + nIrrep) - nAct(i+nIrrep) + j
          write(9,601) '* Active B: orbital ',jj,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,jj+nBas(i),k),k=1,nBasTot)
        end do
      end if
    else
      inactB_done = .true.
      if ( nOcc(i + nIrrep) - nAct(i + nIrrep) .gt. 0 ) then
        do j = 1, nOcc(i + nIrrep) - nAct(i + nIrrep)
          write(9,601) '* Inactive B: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
        end do
      end if
    end if
  end if

  if ( order(4:4) .eq. 'A' ) then
    if ( inactA_done .and. actA_done) then
      if ( nBas(i) - nOcc(i) .gt. 0 ) then
        do j = nOcc(i) + 1, nBas(i)
          write(9,601) '* Virtual A: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
        end do
      end if
    elseif ( inactA_done ) then
      actA_done = .true.
      if ( nAct(i) .gt. 0 ) then
        do j = 1, nAct(i)
          jj = nOcc(i) - nAct(i) + j
          write(9,601) '* Active A: orbital ',jj,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,jj,k),k=1,nBasTot)
        end do
      end if
    else
      inactA_done = .true.
      if ( nOcc(i) - nAct(i) .gt. 0 ) then
        do j = 1, nOcc(i) - nAct(i)
          write(9,601) '* Inactive A: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
        end do
      end if
    end if
  else
    if ( inactB_done .and. actB_done ) then
      if ( nBas(i + nIrrep) - nOcc(i + nIrrep) .gt. 0 ) then
        do j = nOcc(i+nIrrep) + 1, nBas(i+nIrrep)
          write(9,601) '* Virtual B: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
        end do
      end if
    elseif ( inactB_done ) then
      actB_done = .true.
      if ( nAct(i + nIrrep) .gt. 0 ) then
        do j = 1, nAct(i + nIrrep)
          jj = nOcc(i + nIrrep) - nAct(i+nIrrep) + j
          write(9,601) '* Active B: orbital ',jj,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,jj+nBas(i),k),k=1,nBasTot)
        end do
      end if
    else
      inactB_done = .true.
      if ( nOcc(i + nIrrep) - nAct(i + nIrrep) .gt. 0 ) then
        do j = 1, nOcc(i + nIrrep) - nAct(i + nIrrep)
          write(9,601) '* Inactive B: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
        end do
      end if
    end if
  end if


  if ( order(5:5) .eq. 'A' ) then
    if ( actA_done ) then
      if ( nBas(i) - nOcc(i) .gt. 0 ) then
        do j = nOcc(i) + 1, nBas(i)
          write(9,601) '* Virtual A: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
        end do      
      end if
    else
      actA_done = .true.
      if ( nAct(i) .gt. 0 ) then
        do j = 1, nAct(i)
          jj = nOcc(i) - nAct(i) + j
          write(9,601) '* Active A: orbital ',jj,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,jj,k),k=1,nBasTot)
        end do
      end if
    end if 
  else
    if ( actB_done ) then
      if ( nBas(i + nIrrep ) - nOcc(i + nIrrep ) .gt. 0 ) then
        do j = nOcc(i + nIrrep ) + 1, nBas(i + nIrrep )
          write(9,601) '* Virtual B: orbital ',j,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
        end do
      end if
    else
      actB_done = .true.
      if ( nAct(i + nIrrep ) .gt. 0 ) then
        do j = 1, nAct(i + nIrrep )
          jj = nOcc(i + nIrrep ) - nAct(i + nIrrep ) + j
          write(9,601) '* Active A: orbital ',jj,' or irrep ',i
          write(9,'(5E22.14)') (coeff(i,jj+nBas(i),k),k=1,nBasTot)
        end do
      end if
    end if
  endif

  if ( order(6:6) .eq. 'A' ) then
    if ( nBas(i) - nOcc(i) .gt. 0 ) then
      do j = nOcc(i) + 1, nBas(i)
        write(9,601) '* Virtual A: orbital ',j,' or irrep ',i
        write(9,'(5E22.14)') (coeff(i,j,k),k=1,nBasTot)
      end do
    end if
  else
    if ( nBas(i + nIrrep ) - nOcc(i + nIrrep ) .gt. 0 ) then
      do j = nOcc(i + nIrrep ) + 1, nBas(i + nIrrep )
        write(9,601) '* Virtual B: orbital ',j,' or irrep ',i
        write(9,'(5E22.14)') (coeff(i,j+nBas(i),k),k=1,nBasTot)
      end do
    end if
  endif
end do
601  format(2(A,I4))

!* These are just dummy arrays, not reflecting real occupation numbers, nor the space
!* to which the orbital belongs
allocate( occNu(totOrb) )
write(9,'(A4)')'#OCC'
write(9,'(A20)')'* OCCUPATION NUMBERS'
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

end subroutine addvect_superimpose2


subroutine addvect_capitalize(string)
implicit none
integer      :: i
character(*) string

do i = 1, len(string)
  if (ichar(string(i:i)).gt.96) then
    string(i:i) = char(ichar(string(i:i))-32)
  endif
end do
return
end subroutine addvect_capitalize


subroutine addvect_split(order)
use addvect_fileInfo
use addvect_genInfo
use addvect_decisions
use addvect_coefficients
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
          write(14,'(5E22.14)')(coeff(i,j+start,k),k=nBas(i)+1,nBasTot)
        end do
      end if
      start = start + nOcc(i+nSym)
    end if
  end if

  if ( order(3:3) .eq. 'A' ) then
    if ( occA_done ) then
      if ( nBas(i) .gt. 0 ) then
        do j = 1, nBas(i) - nOcc(i)
          write(13,601)'* Virtual A: orbital ',j+nOcc(i),' of irrep ',i
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
          write(14,601)'* Virtual B: orbital ',j+nOcc(i+nSym),' of irrep ',i
          write(14,'(5E22.14)')(coeff(i,j+start,k),k=nBas(i)+1,nBasTot)
        end do
      end if
      start = start + nBas(i+nSym) - nOcc(i+nSym)
    else
      if ( nBas(i+nSym) .gt. 0 ) then
        do j = 1, nOcc(i+nSym)
          write(14,601)'* Occupied B: orbital ',j,' of irrep ',i
          write(14,'(5E22.14)')(coeff(i,j+start,k),k=nBas(i)+1,nBasTot)
        end do
      end if
      start = start + nOcc(i+nSym)
    end if
  end if

  if ( order(4:4) .eq. 'A' ) then
    if ( nBas(i) .gt. 0 ) then
      do j = 1, nBas(i) - nOcc(i)
        write(13,601)'* Virtual A: orbital ',j+nOcc(i),' of irrep ',i
        write(13,'(5E22.14)')(coeff(i,j+start,k),k=1,nBas(i))
      end do
    end if
  else
    if ( nBas(i+nSym) .gt. 0 ) then
      do j = 1, nBas(i+nSym) - nOcc(i+nSym)
        write(14,601)'* Virtual B: orbital ',j+nOcc(i+nSym),' of irrep ',i
        write(14,'(5E22.14)')(coeff(i,j+start,k),k=nBas(i)+1,nBasTot)
      end do
    end if
  end if
end do
601  format(2(A,I4))

!* These are just dummy arrays, not reflecting real occupation numbers, nor the space
!* to which the orbital belongs
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
          write(jj,602)mod(int(k/10),10),(dummy_label(kk),kk=k,k+9)
        else
          write(jj,602)mod(int(k/10),10),(dummy_label(kk),kk=k,nBasTot)
        end if
      end do
    end if
  end do
end do
602  format(I1,x,10A1)

end subroutine addvect_split 
