! ======================================================================
module orthoFMO_geninfo
implicit none
integer                       :: occA,occB,virtA,virtB
integer                       :: nBas
real(kind=8),allocatable      :: vec(:,:),occNu(:)
character(len=1),allocatable  :: orbLabel(:)
end module orthoFMO_geninfo
! ======================================================================


program lowdin_ortho
use orthoFMO_geninfo
implicit none

external :: orthoFMO_read_info,orthoFMO_getAtomicOverlap,orthoFMO_lowdin
external :: orthoFMO_dump_vectors

integer                   :: i,j,k,l,luone
real(kind=8)              :: overlap(4,4),ortho(4,4)
real(kind=8),allocatable  :: FMO(:,:),orthoFMO(:,:),sAO(:,:)
character(len=12)         :: filename

call orthoFMO_read_info

allocate(sAO(nBas,nBas))
allocate(FMO(4,nBas))
allocate(orthoFMO(4,nBas))

filename = 'ONEINT'
luone = 35
call orthoFMO_getAtomicOverlap(filename,luone,nBas,sAO)

FMO = 0.0d0
orthoFMO = 0.0d0

j = occA + occB - 2
do i = 1, 4
  FMO(i,:) = vec(j+i,:)
end do

overlap = 0.0d0
do i = 1, 4
  do j = 1, i
    do k = 1, nBas
      do l = 1, nBas
        overlap(i,j) = overlap(i,j) + FMO(i,k) * FMO(j,l) * sAO(k,l)
      end do
    end do
    overlap(j,i) = overlap(i,j)
    write(*,'(2i4,f14.8)') i,j,overlap(i,j)
  end do
end do
write(*,*)

call orthoFMO_lowdin(4,overlap,ortho)
orthoFMO = matmul(transpose(ortho),FMO)

overlap = 0.0d0
do i = 1, 4
  do j = 1, i
    do k = 1, nBas
      do l = 1, nBas
        overlap(i,j) = overlap(i,j) + orthoFMO(i,k) * orthoFMO(j,l) * sAO(k,l)
      end do
    end do
    overlap(j,i) = overlap(i,j)
    write(*,'(2i4,f14.8)') i,j,overlap(i,j)
  end do
end do

j = occA + occB - 2
do i = 1, 4
  vec(i+j,:) = orthoFMO(i,:)
end do

call orthoFMO_dump_vectors

end program lowdin_ortho

! ======================================================================

subroutine orthoFMO_read_info
use orthoFMO_geninfo
implicit none

external :: orthoFMO_find_mark

character (len=132)           :: line
integer                       :: dummy,nSym,i,j,k

open (9,file='AddOrb',status='old')

call orthoFMO_find_mark('#INFO ',9)
read(9,'(A132)') line
read(9,*) dummy, nSym
if ( nSym .ne. 1 ) then
  write(*,*)'Program not (yet) prepared for symmetry'
  write(*,*)'come back later'
  stop
end if
nBas = 0
read(9,*) nBas
read(9,*) dummy
read(9,*) line
read(9,'(2I8)') occA,occB
read(9,'(2I8)') virtA,virtB

allocate(vec(nBas,nBas))
allocate(occNu(nBas))
allocate(orbLabel(nBas))

call orthoFMO_find_mark('#ORB  ',9)
k = 0
do i = 1, occA-1
  k = k + 1
  read(9,'(A132)') line
  read(9,'(5E22.14)') (vec(i,j),j=1,nBas)
end do
k = k + 1
write(*,'(a,2i5)')'homo A',k,occA+occB-1
read(9,'(A132)') line
read(9,'(5E22.14)') (vec(occA+occB-1,j),j=1,nBas)
do i = 1, occB-1
  k = k + 1
  read(9,'(A132)') line
  read(9,'(5E22.14)') (vec(i+occA-1,j),j=1,nBas)
end do
k = k + 1
write(*,'(a,2i5)')'homo B',k,occA+occB
read(9,'(A132)') line
read(9,'(5E22.14)') (vec(occA+occB,j),j=1,nBas)
k = k + 1
write(*,'(a,2i5)')'lumo A',k,occA+occB+1
read(9,'(A132)') line
read(9,'(5E22.14)') (vec(occA+occB+1,j),j=1,nBas)
do i = 2, virtA
  k = k + 1
  read(9,'(A132)') line
  read(9,'(5E22.14)') (vec(i+occA+occB+1,j),j=1,nBas)
end do
k = k + 1
write(*,'(a,2i5)')'lumo B',k,occA+occB+2
read(9,'(A132)') line
read(9,'(5E22.14)') (vec(occA+occB+2,j),j=1,nBas)
do i = 1, virtB-1
  k = k + 1
  read(9,'(A132)') line
  read(9,'(5E22.14)') (vec(i+occA+occB+virtA+1,j),j=1,nBas)
end do
write(*,*)

call orthoFMO_find_mark('#OCC  ',9)
read(9,'(A132)') line
read(9,'(5E22.14)') (occNu(j),j=1,nBas)

call orthoFMO_find_mark('#INDEX',9)
read(9,'(A132)') line
read(9,'(2x,10A)')(orbLabel(j),j=1,nBas)

end subroutine orthoFMO_read_info

! ======================================================================

subroutine orthoFMO_lowdin(nbase,sbase,slow)
implicit none

external :: dsyev

integer, intent(in)        :: nbase
integer                    :: info,j

real(kind=8), intent(in)   :: sbase(nbase,nbase)
real(kind=8), intent(out)  :: slow (nbase,nbase)
real(kind=8), allocatable  :: work(:)

real(kind=8)               :: u(nbase,nbase)
real(kind=8)               :: ut(nbase,nbase)
real(kind=8)               :: ev(nbase)

! diagonalize S matrix: 
u=sbase
allocate(work(4*nbase))
work = 0.0
info = 0
call dsyev('V','L',nbase,u,nbase,ev,work,4*nbase,info)
if ( info .ne. 0 ) then
  write(*,*)'Something went wrong in dsyev in lowdin '
  write(*,*) 'info = ',info
  stop
end if
deallocate(work)

! S(diag)=U'SU
ut=transpose(u)
slow=matmul(ut,matmul(sbase,u))

! S(diag)^-1/2
do j=1,nbase
  slow(j,j)=1.0d0/dsqrt(slow(j,j))
enddo

! S^-1/2 = U S(diag)^-1/2 U'
slow=matmul(u,matmul(slow,ut))

end subroutine orthoFMO_lowdin

! ======================================================================
 
subroutine orthoFMO_find_mark(mark,lUnit)
implicit none
integer            :: lUnit,l
character(len=6)   :: mark
character(len=132) :: line

l = len(trim(mark))
rewind(lUnit)
101 read(lUnit,'(A132)') line
if ( line(1:l) .ne. mark(1:l) ) goto 101

end subroutine orthoFMO_find_mark

! ======================================================================

subroutine orthoFMO_getAtomicOverlap(filename,luOne,n,sAO)
implicit none

external :: namerun,opnone,rdone,clsone

integer,intent(in)                :: n,luOne
integer                           :: iCounter,iComponent
integer                           :: iRC,iOpt,iSymLbl,j,k

real (kind=8)                     :: s( n * (n + 1 ) / 2 + 4)
real (kind=8),intent(out)         :: sAO(n,n)

character (len=12),intent(in)     :: filename

call NameRun('RUNFILE')
s = 0.0
sAO = 0.0
iRc=-1
iOpt=0
Call OpnOne(iRC,iOpt,filename,LuOne)
if (iRC.ne.0) write(6,*)'Something went wrong opening ',filename
iRC =  0
iOpt = 2
iComponent = 1
iSymLbl = 1
Call RdOne(iRC,iOpt,'Mltpl  0',iComponent,s,iSymLbl)
iCounter = 1
do j = 1, n
  do k = 1, j
    sAO(j,k) = s(iCounter)
    sAO(k,j) = s(iCounter)
    iCounter = iCounter + 1
  end do
end do
iOpt = 0
Call ClsOne(iRc,iOpt)
return
end subroutine orthoFMO_getAtomicOverlap

! ======================================================================

subroutine orthoFMO_dump_vectors
use orthoFMO_geninfo
implicit none

integer    :: i,j,k


open(12,file='FMOORB')
write(12,'(A11)')'#INPORB 2.2'
write(12,'(A5)') '#INFO'
write(12,'(A36)')'* Lowdin orthogonalized frontier MOs'
write(12,'(A24)')'       0       1       0'
write(12,'(I8)') nBas
write(12,'(I8)') nBas
write(12,'(A4)')'#ORB'
do i = 1, nBas
  write(12,'(A14,I4)')'* ORBITAL    1',i
  write(12,'(5E22.14)')(vec(i,j),j=1,nBas)
end do
write(12,'(A4)')'#OCC'
write(12,'(A20)')'* OCCUPATION NUMBERS'
write(12,'(5E22.14)')(occNu(i),i=1,nBas)
write(12,'(A6)')'#INDEX'
write(12,'(A12)')'* 1234567890'
do j = 1, nBas, 10
  if ( j + 9 .le. nBas) then
    write(12,102)mod(int(j/10),10),(orbLabel(k),k=j,j+9)
  else
    write(12,102)mod(int(j/10),10),(orbLabel(k),k=j,nBas)
  end if
end do
102 format(I1,x,10A1)
close(12)
return
end subroutine orthoFMO_dump_vectors

! ======================================================================

