!================================================================================!

module dipro_data
implicit none
integer                    :: nSym,maxBas
integer,allocatable        :: nBas(:)
real(kind=8),allocatable   :: phiAB(:,:,:),epsilonAB(:,:)
real(kind=8),allocatable   :: sAO(:,:,:)
end module dipro_data

!================================================================================!

program dipro
use dipro_data
implicit none

external :: NameRun,Get_iScalar,Get_iArray
external :: dipro_getAtomicOverlap,dipro_read_FMO,dipro_read_vec

integer                    :: i,j,k,l
integer                    :: symHOMO_A,symHOMO_B
integer                    :: symLUMO_A,symLUMO_B
real(kind=8)               :: epsHOMO_A,epsHOMO_B
real(kind=8)               :: epsLUMO_A,epsLUMO_B
real(kind=8)               :: J_AB,t_AB,S_AB,gammaA,gammaB
real(kind=8),allocatable   :: HOMO_A(:),LUMO_A(:)
real(kind=8),allocatable   :: HOMO_B(:),LUMO_B(:)

write(*,'(a)') 'Electron and hole transport based on:'
write(*,'(a)') 'B. Baumeier, J. Kirkpatrick, D. Andrienko'
write(*,'(a)') 'PCCP 2010, 12, 11103. DOI: 10.1039/c002337j'
write(*,*)

call NameRun('RUNFILE')
call Get_iScalar('nSym',nSym)
allocate(nBas(nSym))
call Get_iArray('nBas',nBas,nSym) 
maxBas = maxval(nBas)
allocate(sAO(nSym,maxBas,maxBas))
allocate(phiAB(nSym,maxBas,maxBas))
allocate(epsilonAB(nsym,maxBas))
allocate(HOMO_A(maxBas))
allocate(HOMO_B(maxBas))
allocate(LUMO_A(maxBas))
allocate(LUMO_B(maxBas))

call dipro_getAtomicOverlap
write(*,*) ' Fragment A'
call dipro_read_FMO('INPORB_A',HOMO_A,LUMO_A,epsHOMO_A,epsLUMO_A,symHOMO_A,symLUMO_A)
write(*,*)
write(*,*) ' Fragment B'
call dipro_read_FMO('INPORB_B',HOMO_B,LUMO_B,epsHOMO_B,epsLUMO_B,symHOMO_B,symLUMO_B)
write(*,*)
call dipro_read_vec('INPORB_AB')
! First the hole transport (HOMO on A and B)
J_AB = 0.0d0
S_AB = 0.0d0
if (symHOMO_A .eq. symHOMO_B) then
  i = symHOMO_A
  do j = 1, nBas(i)
    gammaA = 0.0d0
    gammaB = 0.0d0
    do k = 1, nBas(i)
      do l = 1, nBas(i)
        gammaA = gammaA + HOMO_A(k) * phiAB(i,j,l) * sAO(i,k,l)
        gammaB = gammaB + HOMO_B(k) * phiAB(i,j,l) * sAO(i,k,l)
      end do
      S_AB = S_AB + HOMO_A(j) * HOMO_B(k) * sAO(i,j,k)
    end do
    J_AB = J_AB + epsilonAB(i,j) * gammaA * gammaB
  end do
end if
t_AB = (J_AB - 0.5*(epsHOMO_A + epsHOMO_B) * S_AB) / (1 - S_AB**2)
write(*,*) '==> Hole transport'
write(*,'(a,t15,F12.4,a)') 'J_AB = ',J_AB*27211.386,' meV'
write(*,'(a,t15,F12.4)') 'S_AB = ',S_AB
write(*,'(a,t15,2F12.4,a)') 'e_A, e_B = ',epsHOMO_A,epsHOMO_B,' Eh'
write(*,'(a,t15,F12.4,a)') 't_AB = ',t_AB*27211.386,' meV'
write(*,*)
! Next the electron transport (LUMO on A and B)
J_AB = 0.0d0
S_AB = 0.0d0
if (symLUMO_A .eq. symLUMO_B) then
  i = symLUMO_A
  do j = 1, nBas(i)
    gammaA = 0.0d0
    gammaB = 0.0d0
    do k = 1, nBas(i)
      do l = 1, nBas(i)
        gammaA = gammaA + LUMO_A(k) * phiAB(i,j,l) * sAO(i,k,l)
        gammaB = gammaB + LUMO_B(k) * phiAB(i,j,l) * sAO(i,k,l)
      end do
      S_AB = S_AB + LUMO_A(j) * LUMO_B(k) * sAO(i,j,k)
    end do
    J_AB = J_AB + epsilonAB(i,j) * gammaA * gammaB
  end do
end if
t_AB = (J_AB - 0.5*(epsLUMO_A + epsLUMO_B) * S_AB) / (1 - S_AB**2)
write(*,*) '==> Electron transport'
write(*,'(a,t15,F12.4,a)') 'J_AB = ',J_AB*27211.386,' meV'
write(*,'(a,t15,F12.4)') 'S_AB = ',S_AB
write(*,'(a,t15,2F12.4,a)') 'e_A, e_B = ',epsLUMO_A,epsLUMO_B,' Eh'
write(*,'(a,t15,F12.4,a)') 't_AB = ',t_AB*27211.386,' meV'

end program dipro

! ==========================================================================

subroutine dipro_getAtomicOverlap
use dipro_data
implicit none

external :: OpnOne,RdOne,ClsOne

integer                           :: luOne,ltriangle,nBasTot
integer                           :: iCounter,iComponent
integer                           :: iRC,iOpt,iSymLbl,i,j,k

real (kind=8),allocatable         :: s(:)

character(len=6)                  :: filename

luOne = 20
nBasTot = 0
do i = 1, nSym
  nBasTot = nBasTot + nBas(i)
end do
ltriangle = (nBasTot * (nBasTot+1))/2   
allocate(s(ltriangle+4))
filename = 'ONEINT'
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
do i = 1, nSym
  do j = 1, nBas(i)
    do k = 1, j
      sAO(i,j,k) = s(iCounter)
      sAO(i,k,j) = s(iCounter)
      iCounter = iCounter + 1
    end do
  end do
end do
iOpt = 0
Call ClsOne(iRc,iOpt)
deallocate(s)

end subroutine dipro_getAtomicOverlap

! ==========================================================================

subroutine dipro_read_FMO(filename,HOMO,LUMO,epsHOMO,epsLUMO,symHOMO,symLUMO)
use dipro_data
implicit none

character(len=8),intent(in)     :: filename
real(kind=8),intent(out)        :: HOMO(maxbas),LUMO(maxBas)
real(kind=8),intent(out)        :: epsHOMO,epsLUMO
integer,intent(out)             :: symHOMO,symLUMO

integer                         :: luvec,i,j,k
integer,allocatable             :: nOcc(:)
real(kind=8),allocatable        :: phi(:,:,:)
real(kind=8),allocatable        :: eps(:,:)
character(len=1),allocatable    :: orblabel(:,:)
character(len=6)                :: mark
character(len=255)              :: line

allocate(nOcc(nSym))
allocate(phi(nSym,maxBas,maxBas))
allocate(eps(nSym,maxBas))
allocate(orblabel(nSym,maxBas))
nOcc     = 0
phi      = 0.0d0
eps      = 0.0d0
orbLabel = ' '
epsHOMO  = -1000.0
epsLUMO  =  1000.0
HOMO     = 0.0d0
LUMO     = 0.0d0

luvec = 12
open(luvec,file=filename,status='old')

mark = '#INDEX'
46 read(luvec,'(A132)') line
if (line(1:6).ne.mark) goto 46
do i = 1, nSym
  read(luvec,'(A132)') line
  if (nBas(i) .ne. 0) then
    read(luvec,'(2x,10A)')(orbLabel(i,j),j=1,nBas(i))
  end if
  do j = 1, nBas(i)
    if (orbLabel(i,j) .eq. 'f') nOcc(i) = nOcc(i) + 1     ! frozen in RASSCF
    if (orbLabel(i,j) .eq. 'i') nOcc(i) = nOcc(i) + 1     ! inactive
    if (orbLabel(i,j) .eq. '1') nOcc(i) = nOcc(i) + 1     ! active (ras1)
    if (orbLabel(i,j) .eq. '2') nOcc(i) = nOcc(i) + 1     ! active (ras2)
    if (orbLabel(i,j) .eq. '3') nOcc(i) = nOcc(i) + 1     ! active (ras3)
  end do
end do
rewind(luvec)

mark = '#ONE'
47 read(luvec,'(A132)') line
if (line(1:4).ne.mark) goto 47
read(luvec,'(A132)') line
do i = 1, nSym
  if (nBas(i) .ne. 0) then
    read(luvec,'(10E12.4)') (eps(i,j),j=1,nBas(i))
  end if
  if (eps(i,nOcc(i)) .gt. epsHOMO) then
    epsHOMO = eps(i,nOcc(i))
    symHOMO = i
  end if
  if (eps(i,nOcc(i)+1) .lt. epsLUMO) then
    epsLUMO = eps(i,nOcc(i)+1)
    symLUMO = i
  end if
end do
!epsLUMO = eps(1,nOcc(1)+2)
write(*,'(2(a,i3),a,f8.4)') 'HOMO: irrep',symHOMO,' number ',nOcc(symHOMO),' epsilon: ',epsHOMO
!write(*,'(2(a,i3),a,f8.4)') 'LUMO: irrep',symLUMO,' number ',nOcc(symLUMO)+2,' epsilon: ',epsLUMO
write(*,'(2(a,i3),a,f8.4)') 'LUMO: irrep',symLUMO,' number ',nOcc(symLUMO)+1,' epsilon: ',epsLUMO
rewind(luvec)

mark = '#ORB'
48 read(luvec,'(A132)') line
if (line(1:4).ne.mark) goto 48
do i = 1, nSym
  if (nBas(i) .ne. 0) then
    do j = 1, nBas(i)
      read(luvec,'(A132)') line
      read(luvec,'(5E22.14)') (phi(i,j,k),k=1,nBas(i))
    end do
  end if
end do
do k = 1, nBas(symHOMO)
  HOMO(k) = phi(symHOMO,nOcc(symHOMO),k)
end do
do k = 1, nBas(symLUMO)
  LUMO(k) = phi(symLUMO,nOcc(symLUMO)+1,k)
!  LUMO(k) = phi(symLUMO,nOcc(symLUMO)+2,k)
end do

close(luvec)
deallocate(nOcc,phi,eps,orblabel)

end subroutine dipro_read_FMO


! ==========================================================================

subroutine dipro_read_vec(filename)
use dipro_data
implicit none

integer                      :: luvec,i,j,k
character(len=4)             :: mark
character(len=9),intent(in)  :: filename
character(len=255)           :: line

luvec = 12
open(luvec,file=filename,status='old')

mark = '#ORB'
47 read(luvec,'(A132)') line
if (line(1:4).ne.mark) goto 47
do i = 1, nSym
  if (nBas(i) .ne. 0) then
    do j = 1, nBas(i)
      read(luvec,'(A132)') line
      read(luvec,'(5E22.14)') (phiAB(i,j,k),k=1,nBas(i))
    end do
  end if
end do
rewind(luvec)

mark = '#ONE'
48 read(luvec,'(A132)') line
if (line(1:4).ne.mark) goto 48
read(luvec,'(A132)') line
do i = 1, nSym
  if (nBas(i) .ne. 0) then
    read(luvec,*) (epsilonAB(i,j),j=1,nBas(i))
  end if
end do

close(luvec)

end subroutine dipro_read_vec

! ==========================================================================
