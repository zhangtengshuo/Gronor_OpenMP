! *****************************************************************
! ***   locactive                                               ***
! ***   Program to localize inactive orbitals                   ***
! ***     by projection of a model vector                       ***
! ***       E. Bordas,  R. Caballol, Coen de Graaf              ***
! ***    *Departament de Quimica  Fisica i Inorganica           ***
! ***     Universitat Rovira i Virgili, Tarragona (Spain)       ***
! ***     caballol@quimica.urv.es   coen@urv.net                ***
! ***                                                           ***
! ***                                                           ***
! *** - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ***
! ***                                                           ***
! ***    May 2021: rewritten (coen.degraaf@urv.cat)             ***
! ***    February 2023: Added symmetry                          ***
! *****************************************************************

module locactive_general_info
implicit none
integer                               :: nSym,nBasMax,nBasTot,nProjMax
integer, parameter                    :: lblL = 6                ! LenIN  in openMolcas
integer, parameter                    :: lblL8 = lblL + 8        ! LenIN8 in openMolcas
integer, allocatable                  :: nBas(:),nFroz(:),nInac(:)
integer, allocatable                  :: nRas1(:),nRas2(:),nRas3(:)
integer, allocatable                  :: nDble(:)
real(kind=8),allocatable              :: molOrb(:,:,:),occNu(:,:)
real(kind=8),allocatable              :: sAO(:,:,:),orthoProjMO(:,:,:)
character (len = 1 ) , allocatable    :: orbLabel(:,:)
character (len = lblL8), allocatable  :: basFnLbl(:,:)
end module locactive_general_info

module locactive_input_data
implicit none
logical                    :: debug,ymult,auto
integer                    :: printlevel
integer,allocatable        :: multAct(:,:),nProj(:)
real(kind=8),allocatable   :: modelVector(:,:,:)
end module locactive_input_data

! ===============================================================================

program locactive
implicit none

external  :: locative_read_info,locactive_read_input,locactive_getAtomicOverlap
external  :: locactive_getBasisLabels,locactive_construct_modelVectors
external  :: loactive_projection,loactive_dumpVectors

write(*,*) 'Projection of the active orbitals'
write(*,*) 'based on the method published in'
write(*,*) 'Bordas et al., Chem. Phys. 309, 259–269 (2005).'
write(*,*)

call locactive_read_info
call locactive_read_input
call locactive_getAtomicOverlap
call locactive_getBasisLabels

call locactive_construct_modelVectors
call locactive_projection
call locactive_dump_Vectors

end program locactive

! ===============================================================================


subroutine locactive_read_info
use locactive_general_info
use locactive_input_data
implicit none

character (len = 132)         :: line
character (len = 6)           :: mark
integer                       :: dummy,j,k,iSym

open (9,file='INPORB',status='old')

mark = '#INFO'
44 read(9,'(A132)') line
if (line(1:5).ne.mark) goto 44
read(9,'(A132)') line
read(9,*) dummy, nSym
allocate(nBas(nSym))
read(9,*) (nBas(iSym),iSym=1,nSym)
nBasMax = maxval(nBas)
nBasTot = sum(nBas)

allocate(nProj(nSym))
allocate(molOrb(nSym,nBasMax,nBasMax))
allocate(occNu(nSym,nBasMax))
allocate(orbLabel(nSym,nBasMax))
nProj = 0
molorb = 0.0d0
occNu = 0.0d0
orbLabel =''

mark = '#ORB'
45 read(9,'(A132)') line
if (line(1:4).ne.mark) goto 45
do iSym = 1, nSym
  if (nBas(iSym) .ne. 0) then
    do j = 1, nBas(iSym)
      read(9,'(A132)') line
      read(9,'(5E22.14)') (molOrb(iSym,j,k),k=1,nBas(iSym))
    end do
  end if
end do

mark = '#OCC'
46 read(9,'(A132)') line
if (line(1:4).ne.mark) goto 46
read(9,'(A132)') line
do iSym = 1, nSym
  if (nBas(iSym) .ne. 0 ) then
    read(9,'(5E22.14)') (occNu(iSym,j),j=1,nBas(iSym))
  end if
end do

mark = '#INDEX'
47 read(9,'(A132)') line
if (line(1:6).ne.mark) goto 47
allocate (nFroz(nSym))
allocate (nInac(nSym))
allocate (nRas1(nSym))
allocate (nRas2(nSym))
allocate (nRas3(nSym))
allocate (nDble(nSym))
nFroz = 0
nInac = 0
nRas1 = 0
nRas2 = 0
nRas3 = 0
nDble = 0
do iSym = 1, nSym
  if (nBas(iSym) .ne. 0) then
    read(9,'(A132)') line
    read(9,'(2x,10A)')(orbLabel(iSym,j),j=1,nBas(iSym))
    do j = 1, nBas(iSym)
      if (orbLabel(iSym,j) .eq. 'f') nFroz(iSym) = nFroz(iSym) + 1     ! frozen
      if (orbLabel(iSym,j) .eq. 'i') nInac(iSym) = nInac(iSym) + 1     ! inactive
      if (orbLabel(iSym,j) .eq. '1') nRas1(iSym) = nRas1(iSym) + 1     ! active (ras1)
      if (orbLabel(iSym,j) .eq. '2') nRas2(iSym) = nRas2(iSym) + 1     ! active (ras2)
      if (orbLabel(iSym,j) .eq. '3') nRas3(iSym) = nRas3(iSym) + 1     ! active (ras3)
    end do
    nDble(iSym) = nFroz(iSym) + nInac(iSym)
  end if
end do
close(9)
return
end subroutine locactive_read_info

! ===============================================================================

subroutine locactive_read_input
use locactive_input_data
use locactive_general_info
implicit none

external  :: locactive_capitalize,locactive_locate

integer, parameter                   :: nKeys = 5
integer                              :: iSym,iProj,iBas,jj,iKey
character (len=4)                    :: key
character (len=4), dimension(nKeys)  :: keyword      
character (len=132)                  :: line

logical                              :: all_ok = .true.
logical, dimension(nKeys)            :: hit = .false.

data keyword /'PROJ','MULT','PRIN','DEBU','NOMU'/


debug = .false.
printLevel = 5
ymult = .true.
auto = .true.

do while (all_ok)
  read(5,*,iostat=jj) line
  line = adjustl(line)
  key = line(1:4)
  call locactive_capitalize(key)
  do iKey = 1, nKeys
    if ( key .eq. keyword(iKey) ) hit(iKey) = .true.
  end do
  if (  jj .lt. 0 ) all_ok = .false.
end do

do iKey = 1, nKeys
  if ( hit(iKey) ) then
    select case(iKey)
      case(1)
        call locactive_locate('PROJ')
        read(5,*) (nProj(iSym),iSym=1,nSym)
        if (sum(nProj) .gt. sum(nRas2)) then
          write(*,*) 'Number of model vectors cannot be larger',          &
                     'than the number of active orbitals'
          stop
        end if
        nProjMax = maxval(nProj)
        allocate( modelVector(iSym,nProjMax,nBasMax) )
        modelVector = 0.0d0
        do iSym = 1, nSym
          do iProj = 1, nProj(iSym)
            read(5,*)(modelVector(iSym,iProj,iBas),iBas=1,nBas(iSym))
          end do
        end do
      case(2)
        call locactive_locate('MULT')
        ymult = .true.
        auto = .false.
        allocate( multAct(nSym,nProjMax) )
        multAct = 0
        do iSym = 1, nSym
          if (nProj(iSym) .ne. 0) then
            read(5,*)(multAct(iSym,iProj),iProj=1,nProj(iSym))
          end if
        end do
      case(3)
        call locactive_locate('PRIN')
        read(5,*) printLevel
      case(4)
        debug = .true.
      case(5)
        ymult = .false.
        auto = .false.
        if (hit(2)) then
          write(*,*) 'MULTiply and NOMUltiply are mutually exclusive'
          stop
        end if
    end select
  end if
end do
return
end subroutine locactive_read_input

! ===============================================================================

subroutine locactive_construct_modelVectors
use locactive_general_info
use locactive_input_data
implicit none

integer                       :: iSym,iProj,i,j,k,start,finish,iAct
real(kind = 8), allocatable   :: mVec(:,:,:),norm(:,:),testVec(:,:,:)
real(kind = 8), allocatable   :: testnorm(:,:)
real(kind = 8)                :: maxnorm

if ( ymult ) then
  allocate( mVec(nSym,nProjMax,nBasMax) )
  if ( auto ) then
    allocate(multAct(nSym,nProjMax))
    multAct = 0
    do iSym = 1, nSym
      allocate( testnorm(nProj(iSym),nRas2(iSym)) )
      allocate( testVec(nProj(iSym),nRas2(iSym),nBas(iSym)) )
      testnorm = 0
      do iProj = 1, nProj(iSym)
        do iAct = 1, nRas2(iSym)
          do j = 1, nBas(iSym)
            testvec(iProj,iAct,j) = molOrb(iSym,nDble(iSym)+nRas1(iSym)+iAct,j) &
                                          * modelVector(iSym,iProj,j)
          end do
          do j = 1, nBas(iSym)
            do k = 1, nBas(iSym)
              testnorm(iProj,iAct) = testnorm(iProj,iAct) + testVec(iProj,iAct,j) & 
                                      * testVec(iProj,iAct,k) * sAO(iSym,j,k)
            end do
          end do
        end do
      end do
      do iProj = 1, nProj(iSym)
        maxnorm = 0
        do iAct = 1, nRas2(iSym)
          if ( testnorm(iProj,iAct).gt.maxnorm ) then
            maxnorm = testnorm(iProj,iAct)
            multAct(iSym,iProj) = nDble(iSym)+nRas1(iSym)+iAct
          end if
        end do
        testnorm(:,multAct(iSym,iProj)-(nDble(iSym)+nRas1(iSym))) = 0.0
      end do
      deallocate(testnorm,testVec)
    end do
    if ( printLevel .ge. 10 ) then
      write(*,'(a)')'model vectors are multiplied with the following active orbitals'
      do iSym = 1, nSym
        if (nSym .gt. 1) write(*,'(a,i3)') 'Symmetry ',iSym
        write(*,'(a,20i4)') 'model vectors  : ',(iProj,iProj=1,nProj(iSym))
        write(*,'(a,20i4)') 'active orbitals: ',(multAct(iSym,iProj),iProj=1,nProj(iSym))
      end do
      write(*,*)
    endif
  endif
  mVec = 0.0d0
  do iSym = 1, nSym
    do iProj = 1, nProj(iSym)
      k = multAct(iSym,iProj)
      do j = 1, nBas(iSym)
        mVec(iSym,iProj,j) = modelVector(iSym,iProj,j) * molOrb(iSym,k,j)
      end do
    end do
  end do
  if ( printLevel .ge. 15 ) then
    write(*,*) 'Model vectors after multiplication'
    write(*,*)
    write(*,*) '                     Model vector    Active Orbital    Product'
    do iSym = 1, nSym
      do iProj = 1, nProj(iSym)
        write(*,'(a,2i3)') 'model vector ',iSym,iProj
        k = multAct(iSym,iProj)
        do j = 1, nBas(iSym)
          write(*,'(A,3F18.8)') basFnLbl(iSym,j),modelVector(iSym,iProj,j),     &
                                  molOrb(iSym,k,j),mVec(iSym,iProj,j)
        end do
        write(*,*)
      end do
    end do
  end if
  modelVector = mVec
  deallocate(mVec)
else
  if ( printlevel .ge. 15 ) then
    write(*,*)
    write(*,*) '   Model vector(s)'
    do iSym = 1, nSym
      do i = 1, int(nProj(iSym)/10)+1
        start = (i-1) * 10 + 1
        finish = min(nProj(iSym),i*10)
        do k = 1, nBas(iSym)
          write(*,'(A,10F18.8)') basFnLbl(iSym,k),                              &
                             (modelVector(iSym,j,k),j=start,finish)
        end do
        write(*,*)
      end do
    end do
  end if
endif
! Check the norm of the model vectors
allocate( norm(nSym,nProjMax) )
do iSym = 1, nSym
  do iProj = 1, nProj(iSym)
    norm(iSym,iProj) = 0.0
    do j = 1 , nBas(iSym)
      do k = 1, nBas(iSym)
        norm(iSym,iProj) = norm(iSym,iProj) + modelVector(iSym,iProj,j) *  &
             modelVector(iSym,iProj,k) * sAO(iSym,j,k)
      end do
    end do
    if (abs(norm(iSym,iProj)) .lt. 1e-3) then
      write(*,*) 'Warning: norm of model vector ',iProj,                      &
          'of symmetry', iSym, 'is smaller than 1e-3.'
      write(*,*) ' Carefully check the input.'
    end if
  end do
  write(*,'(a,i3)')' Norm of the model vectors of symmetry',iSym
  do iProj=1,nProj(iSym)
    write(*,'(I4,F18.8)') iProj,norm(iSym,iProj)
  end do
end do
deallocate(norm)
return
end subroutine locactive_construct_modelVectors

! ===============================================================================

subroutine locactive_projection
use locactive_general_info
use locactive_input_data
implicit none

external :: locactive_lowdin

integer                   :: i,j,k,l,iProj,start,finish,iSym
real(kind=8)              :: norm,invsqnorm,maxoverlap
real(kind=8),allocatable  :: aux(:,:),projector(:,:),projMO(:,:,:)
real(kind=8),allocatable  :: overlap(:,:),ortho(:,:)

allocate(aux(nBasMax,nBasMax) )
allocate(projector(nBasMax,nBasMax))
allocate(orthoProjMO(nSym,nProjMax,nBasMax))
orthoProjMO = 0.0
do iSym = 1, nSym
  aux = 0.0
  projector = 0.0
  write(*,*)
  if (printLevel .gt. 5 .or. debug) then
    write(*,'(a,i3)') 'Projection in symmetry : ', iSym
  endif
  do i = nDble(iSym) + nRas1(iSym) + 1, nDble(iSym) + nRas1(iSym) + nRas2(iSym)
    do j = 1, nBas(iSym)
      do k = 1, nBas(iSym)
        aux(j,k) = aux(j,k) + molorb(iSym,i,j) * molorb(iSym,i,k)
      end do
    end do
  end do
  if ( debug ) then
    write(*,*) 'aux'
    do i =1 , nBas(iSym)
      write(*,'(5F18.8)') (aux(i,j),j=1,nBas(iSym))
    end do
  end if
  do j = 1, nBas(iSym)
    do k = 1, nBas(iSym)
      do l = 1, nBas(iSym)
        projector(j,k) = projector(j,k) + aux(j,l) * sAO(iSym,k,l)
      end do
    end do
  end do 
  if ( debug ) then
    write(*,*) 'projector'
    do i =1 , nBas(iSym)
      write(*,'(5F18.8)') (projector(i,j),j=1,nBas(iSym))
    end do
  end if
  allocate(projMO(nSym,nProj(iSym),nBasMax))
  projMO = 0.0
  do iProj = 1, nProj(iSym)
    do j = 1, nBas(iSym)
      do k = 1, nBas(iSym)
        projMO(iSym,iProj,j) = projMO(iSym,iProj,j) + projector(j,k) *           &
                                     modelVector(iSym,iProj,k)
      end do
    end do
  end do 
  if ( debug ) then
    write(*,*) 'projected'
    do i =1 , nProj(iSym)
      write(*,'(5F18.8)') (projMO(iSym,i,j),j=1,nBas(iSym))
    end do
  end if
  do iProj = 1, nProj(iSym)
    norm = 0.0
    do j = 1, nBas(iSym)
      do k = 1, nBas(iSym)
        norm = norm + projMO(iSym,iProj,j) * projMO(iSym,iProj,k) * sAO(iSym,j,k)
      end do
    end do
    invsqnorm = 1 / sqrt(norm)
    do j = 1, nBas(iSym)
      projMO(iSym,iProj,j) = projMO(iSym,iProj,j) * invsqnorm
    end do
  end do
  if ( debug ) then
    write(*,*) 'projected and normalized'
    do i =1 , nProj(iSym)
      write(*,'(5F18.8)') (projMO(iSym,i,j),j=1,nBas(iSym))
    end do
  end if
  allocate( overlap(nProj(iSym),nProj(iSym)) )
  overlap = 0.0
  do iProj = 1, nProj(iSym)
    do j = 1, iProj
      do k = 1, nBas(iSym)
        do l = 1, nBas(iSym)
          overlap(iProj,j) = overlap(iProj,j) + projMO(iSym,iProj,k) *                &
             projMO(iSym,j,l) * sAO(iSym,k,l)
        end do
      end do
      overlap(j,iProj) = overlap(iProj,j)
    end do
  end do
  if ( printLevel .ge. 5 .and. nProj(iSym) .gt. 1 ) then
    write(*,*)
    write(*,*) 'Overlap matrix of the projected model vectors after normalization'
    do iProj = 1, nProj(iSym)
      write(*,'(10F15.8)') (overlap(iProj,j),j=1,nProj(iSym))
    end do
    maxoverlap = 0.0
    do iProj = 1, nProj(iSym)
      do j = 1, iProj-1
        if (abs(overlap(iProj,j)) .gt. abs(maxoverlap)) maxoverlap = overlap(iProj,j)
      end do
    end do
    write(*,'(A,F14.6)') 'Largest off-diagonal element: ',maxoverlap 
  end if
  allocate( ortho(nProj(iSym),nProj(iSym)) )
  ortho = 0.0
  call locactive_lowdin (nProj(iSym),overlap,ortho)
  if ( debug ) then
    write(*,'(a,2i5)') 'Size of ortho : ',size(ortho,1),size(ortho,2)
    write(*,'(a,3i5)') 'Size of projMO: ',size(projMO,1),size(projMO,2),size(projMO,3)
  endif
  orthoProjMO(iSym,:,:) = matmul( transpose(ortho),projMO(iSym,:,:) )
  if ( printLevel .ge. 25 .and. nProj(iSym) .gt. 1 ) then
    overlap = 0.0
    do i = 1, nProj(iSym)
      do j = 1, i
        do k = 1, nBas(iSym)
          do l = 1, nBas(iSym)
            overlap(i,j) = overlap(i,j) + orthoProjMO(iSym,i,k) *                  &
               orthoProjMO(iSym,j,l) *sAO(iSym,k,l)
          end do
        end do
        overlap(j,i) = overlap(i,j)
      end do
    end do
    write(*,*) 'overlap after lowdin'
    do i = 1, nProj(iSym)
      write(*,'(6F18.8)') (overlap(i,j),j=1,nProj(iSym))
    end do
  end if
  if ( printLevel .gt. 10 ) then
    write(*,*) 
    write(*,*) 'Orthonormalized projections'
    do i = 1, int(nProj(iSym)/10)+1
      start = (i-1) * 10 + 1
      finish = min(nProj(iSym),i*10)
      do k = 1, nBas(iSym)
        write(*,'(A,10F18.8)') basFnLbl(iSym,k),(orthoProjMO(iSym,j,k),j=start,finish)
      end do
      write(*,*)
    end do
  end if
  deallocate(projMO,overlap,ortho)
end do
end subroutine locactive_projection

! ===============================================================================

subroutine locactive_dump_Vectors
use locactive_general_info
use locactive_input_data
implicit none

integer  :: i,j,k,iProj,iSym

open(12,file='PROJORB')
write(12,'(A11)')'#INPORB 2.2'
write(12,'(A5)') '#INFO'
write(12,'(A27)')'* Projected active orbitals'
write(12,'(3I8)')0,nSym,0
write(12,'(8I8)') (nBas(iSym),iSym = 1, nSym)
write(12,'(8I8)') (nBas(iSym),iSym = 1, nSym)
write(12,'(A4)')'#ORB'
do iSym = 1, nSym
  do i = 1, nDble(iSym)
    write(12,'(A10,2I4)')'* ORBITAL ',iSym,i
    write(12,'(5E22.14)')(molOrb(iSym,i,j),j=1,nBas(iSym))
  end do
  do iProj = 1, nProj(iSym)
    write(12,'(A10,2I4)')'* ORBITAL ',iSym,nDble(iSym) + iProj
    write(12,'(5E22.14)')(orthoProjMO(iSym,iProj,j),j=1,nBas(iSym))
  end do
  do i = nDble(iSym) + nProj(iSym) + 1, nBas(iSym)
    write(12,'(A10,2I4)')'* ORBITAL ',iSym,i
    write(12,'(5E22.14)')(molOrb(iSym,i,j),j=1,nBas(iSym))
  end do
end do
write(12,'(A4)')'#OCC'
write(12,'(A20)')'* OCCUPATION NUMBERS'
do iSym = 1, nSym
  write(12,'(5E22.14)')(occNu(iSym,i),i=1,nBas(iSym))
end do
write(12,'(A6)')'#INDEX'
do iSym = 1, nSym
  if (nBas(iSym) .ne. 0) then
    write(12,'(A12)')'* 1234567890'
    do j = 1, nBas(iSym), 10
      if ( j + 9 .le. nBas(iSym)) then
        write(12,102)mod(int(j/10),10),(orbLabel(iSym,k),k=j,j+9)
      else
        write(12,102)mod(int(j/10),10),(orbLabel(iSym,k),k=j,nBas(iSym))
      end if
    end do
  end if
end do
102 format(I1,x,10A1)       
close(12)
return
end subroutine locactive_dump_Vectors

! ===============================================================================

subroutine locactive_getAtomicOverlap
use locactive_general_info
implicit none

external :: NameRun,OpnOne,RdOne,ClsOne

integer                           :: iCounter,iComponent,iSym
integer                           :: iRC,iOpt,iSymLbl,j,k,luone

real (kind=8)                     :: s(nBasTot*(nBasTot+1)/2 + 4)
character (len=12)                :: filename

call NameRun('RUNFILE')
s = 0.0
allocate( sAO(nSym,nBasMax,nBasMax) )
sAO = 0.0
iRc=-1
iOpt=0
filename = 'ONEINT'
luOne = 77
Call OpnOne(iRC,iOpt,filename,LuOne)
if (iRC.ne.0) write(6,*)'Something went wrong opening ',filename
iRC =  0
iOpt = 2
iComponent = 1
iSymLbl = 1
Call RdOne(iRC,iOpt,'Mltpl  0',iComponent,s,iSymLbl)
iCounter = 1
do iSym = 1, nSym
  do j = 1, nBas(iSym)
    do k = 1, j
      sAO(iSym,j,k) = s(iCounter)
      sAO(iSym,k,j) = s(iCounter)
      iCounter = iCounter + 1
    end do
  end do
end do
iOpt = 0
Call ClsOne(iRc,iOpt)
return
end subroutine locactive_getAtomicOverlap

! ===============================================================================

subroutine locactive_getBasisLabels
use locactive_general_info
implicit none

external :: get_cArray

integer              :: iCounter,iSym,j
character(len=lblL8) :: label(nBasTot)

allocate( basFnLbl(nSym,nBasMax) )   
basFnLbl = ''

Call Get_cArray('Unique Basis Name',label,nBasTot*LblL8)

iCounter = 1
do iSym = 1, nSym
  do j = 1, nBas(iSym)
    basfnLbl(iSym,j) = label(iCounter)
    iCounter = iCounter + 1
  end do
end do

return
end subroutine locactive_getBasisLabels

! ===============================================================================

subroutine locactive_capitalize(string)
implicit none
integer      :: i
character(*) string

do i = 1, len(string)
  if (ichar(string(i:i)).gt.96) then
    string(i:i) = char(ichar(string(i:i))-32)
  endif
end do
return
end subroutine locactive_capitalize

! ===============================================================================

subroutine locactive_locate(string)
implicit none

external :: locactive_capitalize

character(4)   ::  string,string2
character(132) ::  line
rewind(5)
40 read(5,*) line
line = adjustl(line)
string2=line(1:4)
call locactive_capitalize(string2)
if (string2.ne.string) goto 40
return
end subroutine locactive_locate

! ===============================================================================

subroutine locactive_lowdin(nbase,sbase,slow)
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

end subroutine locactive_lowdin

! ===============================================================================
