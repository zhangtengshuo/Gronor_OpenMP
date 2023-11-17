program michl_SF_model
implicit none

external :: DANAME,WR_MOTRA_Info,dDAFILE,DANAME_MF,iDAFILE,dsyev

integer, parameter        :: nToc = 64
integer, parameter        :: nTraToc = 106
integer, parameter        :: mxSym = 1, maxBfn = 10000
integer, parameter        :: labelSize = 2 * maxBfn * 4
integer, parameter        :: nTraBuf = 9600
integer                   :: iToc(nToc)
integer                   :: iTraToc(nTraToc)
integer                   :: luOne,luTra,iAd30,iAd50
integer                   :: i,j,k,l,ltuvx,n,nt,nOrbtt
integer                   :: nOrb,nBas,nFrozen,nDMO,nSym
integer                   :: nTimes,lastInts,iRc,lwork
integer                   :: projMax1,projMax2,projMax3

real(kind=8), allocatable :: fock(:),kinInt(:)
real(kind=8), allocatable :: Heff(:,:),coef(:,:)
real(kind=8), allocatable :: work(:),eigenvalues(:),projection(:)
real(kind=8)              :: potnuc,Eone,Etwo
real(kind=8)              :: traBuf(nTraBuf)
real(kind=8)              :: H(5,5),coupling1,coupling2,S(3,3)

character(len=1)          :: basLabel(labelSize)

logical                   :: huge_listing

huge_listing = .false.
H = 0.0d0

! open the file with the one-electron integrals (TRAONE)
luOne = 30
call DANAME(luOne,'TRAONE')

! read info from TRAONE
iAd30 = 0
call WR_MOTRA_Info(luOne,2,iAd30,iToc, nToc, potNuc,nSym, nBas, nOrb,     &
                           nFrozen, nDMO, mxSym, basLabel, labelSize)
!write(*,*)'Info from TRAONE'
!write(*,'(A,I4)')'nBas    :',nBas
!write(*,'(A,I4)')'nFrozen :',nFrozen
!write(*,'(A,I4)')'nOrb    :',nOrb
!write(*,'(A,I4)')'nDeleted:',nDMO
!write(*,*)
nOrbtt = (nOrb * (nOrb + 1)) / 2 

if (nOrb .ne. 4) then
  write(*,*) ' * * * * * * *   W A R N I N G   * * * * * * *'
  write(*,*) '     program only works correctly for four'
  write(*,*) '         frontier molecular orbitals'
  write(*,*) 
endif

! read one-electron integrals
allocate ( fock(nOrbtt) )
allocate ( kinInt(nOrbtt) )
iAd30 = iToc(2)
call dDAFILE(luOne,2,fock,nOrbtt,iAd30)
iAd30 = iToc(3)
call dDAFILE(luOne,2,kinInt,nOrbtt,iAd30)

if (huge_listing) then
  write(*,*) 'fock and Kinetic'
  do i = 1, nOrbtt
    write(*,'(i5,3F14.8)') i,fock(i),kinInt(i)
  end do
endif

! open the file with the two-electron integrals (TRAINT)
luTra = 50
call DANAME_MF(lutra,'TRAINT')

! read table of contents of TRAINT
iAd50 = 0
call iDAFILE(luTra,2,iTraToc,nTraToc,iad50)

! calculate number of two-electron integrals, ltuvx
ltuvx = 0
nt = (nOrb * (nOrb + 1)) / 2
do i = 1, nOrb
  do j = 1, i
    do k = 1, nt
      ltuvx = ltuvx + 1
    end do
    nt = nt - 1
  end do
end do

nTimes = int(ltuvx/nTraBuf)
lastInts = ltuvx - ntimes*nTraBuf

do i = 1, nTimes
  call dDAFILE(luTra,2,traBuf,nTraBuf,iad50)
  if (huge_listing) then
    do n = 1, nTraBuf
      write(*,'(I6,F22.15)') (i-1)*nTraBuf+n, traBuf(n)
    end do
  end if
end do
call dDAFILE(luTra,2,traBuf,lastInts,iad50)
if (huge_listing) then
  do n = 1, lastInts
    write(*,'(I8,F22.15)') nTimes*nTraBuf+n, traBuf(n)
  end do
endif

Eone = 2*fock(1)+2*fock(3)
Etwo = traBuf(1)+traBuf(20)+4*traBuf(3)-2*traBuf(11)
write(*,*) Eone,Etwo,potNuc
write(*,'(a,f20.12,a)') 'E(S0) = ',potNuc+Eone+Etwo,' Eh'
write(*,*)

Eone = 2*fock(1)+fock(3)+fock(10)
Etwo = traBuf(1)+2*traBuf(3)+2*traBuf(10)+traBuf(27)-traBuf(11)-traBuf(46)+traBuf(50)
H(1,1) = potNuc+Eone+Etwo

Eone = fock(1)+fock(6)+2*fock(3)
Etwo = traBuf(6)+2*traBuf(3)+2*traBuf(23)+traBuf(20)-traBuf(11)-traBuf(35)+traBuf(28)
H(2,2) = potNuc+Eone+Etwo

Eone = fock(1)+fock(3)+fock(6)+fock(10)
Etwo = traBuf(3)+traBuf(6)+traBuf(10)+traBuf(23)+traBuf(27)+traBuf(45)-traBuf(28)-traBuf(50)
H(3,3) = potNuc+Eone+Etwo

Eone = fock(1)+2*fock(3)+fock(10)
Etwo = 2*traBuf(3)+traBuf(10)+traBuf(20)+2*traBuf(27)-traBuf(11)-traBuf(50)
H(4,4) = potNuc+Eone+Etwo

Eone = 2*fock(1)+fock(3)+fock(6)
Etwo = traBuf(1)+2*traBuf(6)+2*traBuf(3)+traBuf(23)-traBuf(28)-traBuf(11)
H(5,5) = potNuc+Eone+Etwo

H(1,2) = 2*traBuf(32)-traBuf(18)
H(2,1) = H(1,2)

H(1,3) = sqrt(3.0d0/2.0d0)*(traBuf(14)-traBuf(48))
H(2,3) = sqrt(3.0d0/2.0d0)*(traBuf(16)-traBuf(39))
H(3,1) = H(1,3)
H(3,2) = H(2,3)

H(1,4) = -fock(2)-traBuf(2)-traBuf(12)-traBuf(19)+2*traBuf(47)
H(4,1) = H(1,4)
H(2,4) = fock(9)+traBuf(9)+2*traBuf(26)-traBuf(38)+traBuf(31)
H(4,2) = H(2,4)
H(1,5) = fock(9)+traBuf(38)+2*traBuf(9)-traBuf(31)+trabuf(26)
H(5,1) = H(1,5)
H(2,5) = -fock(2)+2*traBuf(29)-traBuf(2)-traBuf(15)-traBuf(12)
H(5,2) = H(2,5)

H(3,4) = -sqrt(3.0d0/2.0d0)*(fock(5)+traBuf(5)-traBuf(13)+traBuf(40)+traBuf(22))
H(4,3) = H(3,4)
H(3,5) = -sqrt(3.0d0/2.0d0)*(fock(7)+traBuf(7)+traBuf(42)+traBuf(24)-traBuf(17))
H(5,3) = H(3,5)

H(4,5) = 2*traBuf(37)-traBuf(18)
H(5,4) = H(4,5)

write(*,*) 'Raw Hamiltonian'
write(*,*) '---------------'
write(*,'(14x,a)')'S0S1            S1S0            T1T1            D+D-            D-D+'
write(*,'(a,5F16.8)') 'S0S1 ',H(1,:)
write(*,'(a,5F16.8)') 'S1S0 ',H(2,:)
write(*,'(a,5F16.8)') 'T1T1 ',H(3,:)
write(*,'(a,5F16.8)') 'D+D- ',H(4,:)
write(*,'(a,5F16.8)') 'D-D+ ',H(5,:)
write(*,*)
write(*,'(a,f8.2,a)') '< S0S1 | H | T1T1 > :',H(1,3)*27211.4,' meV'
write(*,'(a,f8.2,a)') '< S1S0 | H | T1T1 > :',H(2,3)*27211.4,' meV'
write(*,*)
write(*,*)
write(*,*)
write(*,*)
!coupling of S0S1 +/- S1S0 with T1T1
write(*,*) 'S0S1 +/- S1S0 coupling with 1TT'
write(*,*) '-------------------------------'
write(*,*)
lwork = 8
allocate (work(lwork))
allocate (eigenvalues(2))
allocate (coef(2,2))
coef(1,1) = H(1,1)
coef(2,2) = H(2,2)
coef(1,2) = H(1,2)
coef(2,1) = H(2,1)
work = 0.0d0
eigenvalues = 0.0d0
iRc = 1
call dsyev('V','L',2,coef,2,eigenvalues,work,lwork,iRc)
write(*,'(a,2(f11.6,a))') 'Phi1 =',coef(1,1),' |S0S1> +',coef(2,1),' |S1S0>'
write(*,'(a,2(f11.6,a))') 'Phi2 =',coef(1,2),' |S0S1> +',coef(2,2),' |S1S0>'
write(*,'(a,f14.6)') 'E(Phi1) = ',eigenvalues(1)
write(*,'(a,f14.6)') 'E(Phi2) = ',eigenvalues(2)
write(*,'(a,f16.8,a)') '<Phi1 | H | T1T1 > :',(coef(1,1)*H(1,3) + coef(2,1)*H(2,3))*27211.4,' meV'
write(*,'(a,f16.8,a)') '<Phi2 | H | T1T1 > :',(coef(1,2)*H(1,3) + coef(2,2)*H(2,3))*27211.4,' meV'
deallocate(coef,eigenvalues,work)
write(*,*)
write(*,*)
write(*,*)
write(*,*)
!CT enhanced
write(*,*) 'Total coupling (CT enhanced)'
write(*,*) '----------------------------'
write(*,*)
allocate(coef(5,3))
coef=0.0d0
lwork = 16
allocate (work(lwork))
allocate (eigenvalues(4))
allocate (Heff(4,4))
Heff(1,1) = H(1,1)
Heff(2,2) = H(2,2)
Heff(3,3) = H(4,4)
Heff(4,4) = H(5,5)

Heff(1,2) = H(1,2)
Heff(1,3) = H(1,4)
Heff(1,4) = H(1,5)
Heff(2,3) = H(2,4)
Heff(2,4) = H(2,5)
Heff(3,4) = H(4,5)

Heff(2,1) = H(2,1)
Heff(3,1) = H(4,1)
Heff(4,1) = H(5,1)
Heff(3,2) = H(4,2)
Heff(4,2) = H(5,2)
Heff(4,3) = H(5,4)


write(*,*)
write(*,*) 'Block1 Hamiltonian'
write(*,*) '---------------'
write(*,'(14x,a)')'S0S1            S1S0            D+D-            D+D-'
write(*,'(a,5F16.8)') 'S0S1 ',Heff(1,:)
write(*,'(a,5F16.8)') 'S1S0 ',Heff(2,:)
write(*,'(a,5F16.8)') 'D+D- ',Heff(3,:)
write(*,'(a,5F16.8)') 'D-D+ ',Heff(4,:)
write(*,*)

work = 0.0d0
eigenvalues = 0.0d0
iRc = 1
call dsyev('V','L',4,Heff,4,eigenvalues,work,lwork,iRc)
do i = 1, 4
  write(*,'(i4,f14.6,4f11.6)')i,eigenvalues(i),Heff(:,i)
enddo
allocate(projection(4))
projection = 0.0d0
do i = 1, 4
  do j = 1, 2
    projection(i) = projection(i) + Heff(j,i)**2
  end do
end do
projMax1 = maxloc(projection,1)
projection(projMax1) = 0.0d0
projMax2 = maxloc(projection,1)
write(*,'(11x,a)')'|S0S1>     |S1S0>     |D+D->     |D-D+>'
write(*,'(a,4f11.6)') 'Phi1: ',Heff(:,projMax1)
write(*,'(a,4f11.6)') 'Phi2: ',Heff(:,projMax2)
write(*,'(a,f14.6)') 'E(Phi1) = ',eigenvalues(ProjMax1)
write(*,'(a,f14.6)') 'E(Phi2) = ',eigenvalues(projMax2)
coef(1,1) = Heff(1,projMax1)
coef(2,1) = Heff(2,projMax1)
coef(4,1) = Heff(3,projMax1)
coef(5,1) = Heff(4,projMax1)
coef(1,2) = Heff(1,projMax2)
coef(2,2) = Heff(2,projMax2)
coef(4,2) = Heff(3,projMax2)
coef(5,2) = Heff(4,projMax2)
deallocate(eigenvalues,projection,work,Heff)
lwork = 12
allocate(eigenvalues(3))
allocate(work(lwork))
allocate(Heff(3,3))
Heff(1,1) = H(3,3)
Heff(2,2) = H(4,4)
Heff(3,3) = H(5,5)
Heff(1,2) = H(3,4)
Heff(1,3) = H(3,5)
Heff(2,3) = H(4,5)
Heff(2,1) = H(4,3)
Heff(3,1) = H(5,3)
Heff(3,2) = H(5,4)

write(*,*)
write(*,*) 'Block2 Hamiltonian'
write(*,*) '---------------'
write(*,'(14x,a)')'T1T1            D+D-            D+D-'
write(*,'(a,5F16.8)') 'T1T1 ',Heff(1,:)
write(*,'(a,5F16.8)') 'D+D- ',Heff(2,:)
write(*,'(a,5F16.8)') 'D-D+ ',Heff(3,:)
write(*,*)


work = 0.0d0
eigenvalues = 0.0d0
iRc = 1
call dsyev('V','L',3,Heff,3,eigenvalues,work,lwork,iRc)
do i = 1, 3
  write(*,'(i4,f14.6,3f11.6)')i,eigenvalues(i),Heff(:,i)
enddo
allocate(projection(3))
projection = 0.0d0
do i = 1, 3
  projection(i) = Heff(1,i)**2
end do
projMax3 = maxloc(projection,1)
write(*,*)
write(*,'(11x,a)')'|T1T1>     |D+D->     |D-D+>'
write(*,'(a,3f11.6)') 'Phi3: ',Heff(:,projMax3)
write(*,'(a,f14.6)') 'E(Phi3) = ',eigenvalues(ProjMax3)
coef(3,3) = Heff(1,ProjMax3)
coef(4,3) = Heff(2,ProjMax3)
coef(5,3) = Heff(3,ProjMax3)
S=0.0d0
Heff=0.0d0
do i = 1, 3
  do j = 1, i
    do k = 1, 5
      do l = 1, 5
        Heff(i,j) = Heff(i,j) + coef(k,i)*coef(l,j)*H(k,l)
      end do
      S(i,j) = S(i,j) + coef(k,i)*coef(k,j)
    end do
    Heff(j,i) = Heff(i,j)
    S(j,i) = S(i,j)
  end do
end do
write(*,*)
write(*,'(a)')'Hamiltonian and overlap in the basis of the new MEBFs'
do i = 1,3
  write(*,'(3f14.6,8x,3f14.6)')Heff(i,:),S(i,:)
end do
write(*,*)
coupling1 = (Heff(3,1)-0.5*(Heff(1,1)+Heff(3,3))*S(3,1))/(1-S(3,1)**2)
coupling2 = (Heff(3,2)-0.5*(Heff(2,2)+Heff(3,3))*S(3,2))/(1-S(3,2)**2)
write(*,'(a,f16.8,a)') '<Phi1 | H | Phi3 > :',coupling1*27211.4,' meV'
write(*,'(a,f16.8,a)') '<Phi2 | H | Phi3 > :',coupling2*27211.4,' meV'
write(*,*)
end program michl_SF_model

