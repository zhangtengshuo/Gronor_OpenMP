!****************************************************************************!
!                                                                            !
!     program for overlapping fragments                                      !
!                                                                            !
! Corresponding orbitals of fragment A and B are determined and written      !
!   in SUPERORB. For corresponding orbitals with a singular value larger     !
!   than 'thrs' (default = 0.2), only the (A+B)/2 orbital is written. In     !
!   this way, the overlap between A and B is eliminated.                     !
!   The program can also eliminate basis functions for atoms that are        !
!   present in the fragments but do not appear in the supermolecule.         !
!   Typically, hydrogens that were used to saturate bond that were broken    !
!   in the definition of the fragments from the supermolecule.               !
!                                                                            !
! INPUT                                                                      !
!   - INPORB1     :   Orbitals of fragment A                                 !
!   - INPORB2     :   Orbitals of fragment B                                 !
!      overlapping atoms should be the last ones in fragment A and the       !
!      first ones in fragment B                                              !
!   - RUNFILE     :   RunFile from the supermolecule                         !
!   - ONEINT      :   One-electron integrals of the supermolecule            !
!   - GUESSORB    :   Any vector file of the supermolecule (e.g. GssOrb)     !
!      Only used to write the virtual orbitals in SUPERORB                   !
!   - vecdetA     :   CI vector of fragment A                                ! 
!   - vecdetB     :   CI vector of fragment B                                !
!                                                                            !
!  OUTPUT                                                                    !
!   - SUPERORB       : Vector file with the superimposed fragments           !
!   - CORRORB.A      : Corresponding orbitals of fragment A                  !
!   - CORRORB.B      : Corresponding orbitals of fragment B                  !
!   - SUPERORB.ORTHO : Vector file with the orthogonalized superimposed      !
!                      fragments                                             !
!   - detAB          : CI vector of the supermolecule                        !
!                                                                            !
!      The atoms that are to be eliminated from the supermolecule            !
!      should be labeled with the letter 'Q'                                 !
!                                                                            !
! To be linked to molcas libraries for accessing ONEINT and RUNFILE, the     !
!    program also uses the 'dgesvd' routine from lapack                      !
!                                                                            !
! Input file has three keywords: 'SVDThreshold', 'CITHreshold' and 'SPIN'    ! 
!      SVDThreshold  :   Threshold for considering an orbital to be part     !
!                        of the overlapping fragment                         !  
!      CITHreshold   :   Threshold for writing a determinant on the det file !
!      SPIN          :   Target spin of the complete molecule                !
!         (input is case-insensitive and only the first four characters      !
!          of the keywords are being read)                                   !
!                                                                            !
!    Coen de Graaf, Sept. 2019, beta version                                 !
!        the program seems to work, but more testing is highly desirable     ! 
!                                                                            !
!    January 2021, added the possibility to use RASSCF fragment wave         !
!        functions                                                           ! 
!                                                                            !
!   February 2021, orthogonalized superimposed orbitals are also written     !
!                                                                            !
!   May 2021: implementation of the superposition of the CI vectors of       !
!             the two fragments, Aitor Mansilla (URV)                        !
!                                                                            !
!   February 2022: -fixed incorrect number of inactive orbitals on det file  !  
!                  -removed commnand line input, replaced with keyword-      !
!                   based input from a file                                  !
!                  -introduced threshold for writing a det on the det file   !
!                  -duplicate dets are eliminated                            !
!                                                                            !
!   December 2024: -Introduced the possiblity to 'freeze' the core orbitals  !
!                  -Normalization of the corresponding orbitals, important   !
!                   when dummy atoms are removed                             !
!                  -Deactivated orthogonalized superorb                      !
!                                                                            !
!                                                                            !
!****************************************************************************!

module ovlfdets_data
implicit none

integer                  :: target_spin,nci,maxcib,maxnact,nSym
integer,dimension(2)     :: inactm,nactm,idetm,spinm
integer,allocatable      :: ioccm(:,:,:),nFrozen1(:),nFrozen2(:)
real(kind=8),allocatable :: civm(:,:)
real(kind=8)             :: ciThreshold,svdThreshold
logical                  :: normalize

end module ovlfdets_data

program overlapping_fragments
use ovlfdets_data
implicit none

external :: ovlf_readin,NameRun,Get_iArray,Get_cArray,OpnOne,RdOne,ClsOne
external :: dgesvd,ovlf_dets

integer                       :: iSym
integer                       :: lDim,offset1,offset2
integer                       :: i,j,k,kk,l,m1,m2,m3,jOrb,first,last
integer                       :: LuOne,iRC,iOpt,iComponent,iSymLbl
integer                       :: nBasTot,iCounter
integer                       :: lTriangle
integer                       :: dummy,start,finish
integer                       :: maxInact,nBasFinal,maxFrozen
integer, allocatable          :: nBas1(:),nBas2(:),nBas(:)
integer, allocatable          :: nInact1(:),nInact2(:)
integer, allocatable          :: nAct1(:),nAct2(:)
integer, allocatable          :: nRas1_1(:),nRas1_2(:),nRas2_1(:),nRas2_2(:),nRas3_1(:),nRas3_2(:)
integer, allocatable          :: n_elim(:),to_be_averaged(:),to_be_averaged_frozen(:)

real (kind = 8), allocatable  :: c1(:,:,:),c2(:,:,:),c3(:,:,:),c4(:,:,:)
real (kind = 8), allocatable  :: c_1(:,:),c_2(:,:)
real (kind = 8), allocatable  :: c1c(:,:),c2c(:,:)
real (kind = 8), allocatable  :: occNu1(:,:),occNu2(:,:)
real (kind = 8), allocatable  :: occNu3(:,:)
real (kind = 8), allocatable  :: s(:),sAO(:,:,:),sAO_mod(:,:,:),sMO(:,:)
real (kind = 8), allocatable  :: U(:,:),VT(:,:),V(:,:)
real (kind = 8), allocatable  :: aux(:,:),work(:)
real (kind = 8), allocatable  :: sigma(:)
real (kind = 8)               :: norm

character (len = 132)         :: line,title
character (len = 6)           :: mark
character (len = 14) , allocatable    :: basLabel(:),basLabel_mod(:,:)
character (len = 1 ) , allocatable    :: orbLabel(:),orblabel12(:,:)

logical , allocatable         :: eliminate(:,:)
logical                       :: debug,new

integer,allocatable       :: nOcc(:)
integer                   :: nB,nAct,nInact,nFrozen

debug = .false.

!...Open the RUNFILE and retrieve the number of irreps
call NameRun('RUNFILE')
call Get_iScalar('nSym',nSym)
allocate(nFrozen1(nSym))
allocate(nFrozen2(nSym))

!...Read the input file
call ovlf_readin
write(*,'(A,F8.3)') 'Running ovlf with threshold        : ',svdThreshold
write(*,'(A,I3)')   'for a state with spin multiplicity : ',target_spin
write(*,*)
target_spin = target_spin - 1

!...Open the files where the corresponding orbitals will be stored
open (12,file='SUPERORB',status='unknown')
open (13,file='CORRORB.A',status='unknown')
open (14,file='CORRORB.B',status='unknown')

!...Retrieving info from  INPORB1
open (9,file='INPORB1',status='old')
read(9,'(A132)') line
write(12,'(A132)') line
write(13,'(A132)') line
read(9,'(A132)') line
write(12,'(A132)') line
write(13,'(A132)') line
read(9,'(A132)') title
write(6,*) "Title vector file 1 : ",title
write(12,'(A32)')'* Superimposed fragment orbitals'
write(13,'(A35)')'* Corresponding orbitals fragment A'
read(9,*) dummy, nSym, dummy
write(12,'(3I8)')0,nSym,0
write(13,'(3I8)')0,nSym,0
allocate ( nBas1(nSym) )
nBas1 = 0
read(9,*)(nBas1(iSym),iSym=1,nSym)
write(6,'(A,8I5)')"Number of basis functions : ",(nBas1(iSym),iSym=1,nSym)
write(13,'(8I8)')(nBas1(iSym),iSym = 1, nSym)
write(13,'(8I8)')(nBas1(iSym),iSym = 1, nSym)
write(13,'(A4)')'#ORB'
lDim = maxval(nBas1)
allocate( c1(nSym,lDim,lDim) )
allocate( occNu1(nSym,lDim) )
c1 = 0.0
occNu1 = 0.0
mark = '#ORB'
44 read(9,'(A132)') line
if (line(1:4).ne.mark) goto 44
do iSym = 1, nSym
  if (nBas1(iSym).ne.0) then
    do j = 1, nBas1(iSym)
      read(9,'(A132)') line
      read(9,'(5E22.14)') (c1(iSym,j,k),k=1,nBas1(iSym))
    end do
  endif
  if ( debug ) then
    write(6,*) 'MO coefficients, irrep ',iSym
    do j = 1, nBas1(isym)
      write(6,'(20F10.5)')(c1(iSym,j,k),k=1,nBas1(iSym))
    end do
  end if
end do
rewind(9)
mark = '#OCC'
45 read(9,'(A132)') line
if (line(1:4).ne.mark) goto 45
read(9,'(A132)') line
do iSym = 1, nSym
  if (nBas1(iSym).ne.0) then
    read(9,'(5E22.14)') (occNu1(iSym,j),j=1,nBas1(iSym))
  endif
end do
rewind(9)
mark = '#INDEX'
allocate( orbLabel(lDim) )
allocate(  nInact1(nSym) )
allocate(  nRas1_1(nSym) )
allocate(  nRas2_1(nSym) )
allocate(  nRas3_1(nSym) )
allocate(    nAct1(nSym) )
nInact1 = 0
nRas1_1 = 0
nRas2_1 = 0
nRas3_1 = 0
46 read(9,'(A132)') line
if (line(1:6).ne.mark) goto 46
read(9,'(A132)') line
do iSym = 1, nSym
  orbLabel = ' '
  if ( nBas1(iSym) .ne. 0 ) then
    read(9,'(2x,10A)')(orbLabel(j),j=1,nBas1(iSym))
  end if
  do j = 1, nBas1(iSym)
    if (orbLabel(j) .eq. 'i') nInact1(iSym) = nInact1(iSym) + 1
    if (orbLabel(j) .eq. '1') nRas1_1(iSym) = nRas1_1(iSym) + 1
    if (orbLabel(j) .eq. '2') nRas2_1(iSym) = nRas2_1(iSym) + 1
    if (orbLabel(j) .eq. '3') nRas3_1(iSym) = nRas3_1(iSym) + 1
  end do
  nInact1(iSym) = nInact1(iSym) - nFrozen1(iSym)
  nAct1(iSym) = nRas1_1(iSym) + nRas2_1(iSym) + nRas3_1(iSym)
end do
write(6,'(A,8I5)') '          frozen orbitals : ',(nFrozen1(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '        inactive orbitals : ',(nInact1(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '          active orbitals : ',(nAct1(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '          RAS1   orbitals : ',(nRas1_1(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '          RAS2   orbitals : ',(nRas2_1(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '          RAS3   orbitals : ',(nRas3_1(iSym),iSym=1,nSym)
write(6,*)
deallocate( orbLabel )
close(9)

!...Retrieving info from  INPORB2
open (10,file='INPORB2',status='old')
read(10,'(A132)') line
write(14,'(A132)') line
read(10,'(A132)') line
write(14,'(A132)') line
read(10,'(A132)') title
write(6,*) "Title vector file 2 : ",title
write(14,'(A35)')'* Corresponding orbitals fragment B'
read(10,*) dummy,nSym,dummy
write(14,'(3I8)')0,nSym,0
allocate ( nBas2(nSym) )
nBas2 = 0
read(10,*) (nBas2(iSym),iSym=1,nSym)
write(6,'(A,8I5)')"Number of basis functions : ",(nBas2(iSym),iSym=1,nSym)
write(14,'(8I8)')(nBas2(iSym),iSym = 1, nSym)
write(14,'(8I8)')(nBas2(iSym),iSym = 1, nSym)
write(14,'(A4)')'#ORB'
lDim = maxval(nBas2)
allocate( orbLabel(lDim) )
allocate( c2(nSym,lDim,lDim) )
allocate( occNu2(nSym,lDim) )
c2 = 0.0
occNu2 = 0.0
mark = '#ORB'
54 read(10,'(A132)') line
if (line(1:4).ne.mark) goto 54
do iSym = 1, nSym
  if (nBas2(iSym).ne.0) then
    do j = 1, nBas2(iSym)
      read(10,'(A132)') line
      read(10,'(5E22.14)') (c2(iSym,j,k),k=1,nBas2(iSym))
    end do
  endif
  if ( debug ) then
    write(6,*) 'MO coefficients, irrep ',iSym
    do j = 1, nBas2(isym)
      write(6,'(20F10.5)')(c2(iSym,j,k),k=1,nBas2(iSym))
    end do
  end if
end do
rewind(10)
mark = '#OCC'
55 read(10,'(A132)') line
if (line(1:4).ne.mark) goto 55
read(10,'(A132)') line
do iSym = 1, nSym
  if (nBas2(iSym).ne.0) then
    read(10,'(5E22.14)') (occNu2(iSym,j),j=1,nBas2(iSym))
  endif
end do
rewind(10)
mark = '#INDEX'
56 read(10,'(A132)') line
if (line(1:6).ne.mark) goto 56
read(10,'(A132)') line
allocate(  nInact2(nSym) )
allocate(  nRas1_2(nSym) )
allocate(  nRas2_2(nSym) )
allocate(  nRas3_2(nSym) )
allocate(    nAct2(nSym) )
nInact2 = 0
nRas1_2 = 0
nRas2_2 = 0
nRas3_2 = 0
do iSym = 1, nSym
  orbLabel = ' '
  if ( nBas1(iSym) .ne. 0 ) then
    read(10,'(2x,10A)')(orbLabel(j),j=1,nBas2(iSym))
  end if
  do j = 1, nBas2(iSym)
    if (orbLabel(j) .eq. 'i') nInact2(iSym) = nInact2(iSym) + 1
    if (orbLabel(j) .eq. '1') nRas1_2(iSym) = nRas1_2(iSym) + 1
    if (orbLabel(j) .eq. '2') nRas2_2(iSym) = nRas2_2(iSym) + 1
    if (orbLabel(j) .eq. '3') nRas3_2(iSym) = nRas3_2(iSym) + 1
  end do
  nAct2(iSym) = nRas1_2(iSym) + nRas2_2(iSym) + nRas3_2(iSym)
  nInact2(iSym) = nInact2(iSym) - nFrozen2(iSym)
end do
write(6,'(A,8I5)') '          frozen orbitals : ',(nFrozen2(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '        inactive orbitals : ',(nInact2(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '          active orbitals : ',(nAct2(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '          RAS1   orbitals : ',(nRas1_2(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '          RAS2   orbitals : ',(nRas2_2(iSym),iSym=1,nSym)
write(6,'(A,8I5)') '          RAS3   orbitals : ',(nRas3_2(iSym),iSym=1,nSym)
write(6,*)
close(10)

!...Read some more info from the RunFile (number of basis functions and basis labels)
allocate ( nBas(nSym) )
allocate ( nOcc(nSym) )
nBas = 0.0
nOcc = 0.0
!Call NameRun('RUNFILE')
Call Get_iArray('nBas',nBas,nSym)
nBasTot = sum(nBas)
lTriangle = ( nBasTot * ( nBasTot + 1 ) ) / 2 + 4
lDim = maxval(nBas)
allocate ( basLabel(nBasTot) )
allocate ( basLabel_mod(nSym,lDim) )
allocate ( eliminate(nSym,lDim) )
allocate ( n_elim(nSym),to_be_averaged(nSym),to_be_averaged_frozen(nSym) )
to_be_averaged = 0
to_be_averaged_frozen = 0
n_elim = 0
Call Get_cArray('Unique Basis Names',basLabel,14*nBasTot)
iCounter = 0
do iSym = 1, nSym
  l = 0
  do j = 1, nBas(iSym)
    iCounter = iCounter + 1
    if (basLabel(iCounter)(1:1) .eq. 'Q') then
      eliminate(iSym,j) = .True.
      n_elim(iSym) = n_elim(iSym) + 1 
    else
      eliminate(iSym,j) = .False.
      l = l + 1
      basLabel_mod(iSym,l) = basLabel(iCounter)
    end if
  end do
end do
write(12,'(8I8)')(nBas(iSym)-n_elim(nSym),iSym = 1, nSym)
write(12,'(8I8)')(nBas(iSym)-n_elim(nSym),iSym = 1, nSym)
write(12,'(A4)')'#ORB'

!...Open the OneInt file and read the overlap matrix of the ao-basis
allocate ( s(lTriangle) )
s = 0.0
allocate ( sAO(nSym,lDim,lDim) )
sAO = 0.0
LuOne = 77
iRc=-1
iOpt=0
Call OpnOne(iRC,iOpt,'ONEINT',LuOne)
if (iRC.ne.0) write(6,*) 'Something went wrong opening ONEINT'
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
  if ( debug ) then
    write(6,*) 'AO overlap'
    do j = 1,nBas(iSym)
      write(6,'(20F10.5)')(sAO(iSym,j,k),k=1,nBas(iSym))
    end do
    write(6,*)
  end if
end do
deallocate( s )
iOpt = 0
Call ClsOne(iRc,iOpt)
!...Construct the AO overlap matrix of the supermolecule after removing the dummy atoms
allocate( sAO_mod(nSym,lDim,lDim) )
do iSym = 1, nSym
  l = 0
  do j = 1, nBas(iSym)
    if ( .not. eliminate(iSym,j) ) then
      l = l + 1
      i = 0
      do k = 1, nBas(iSym)
        if ( .not. eliminate(iSym,k) ) then
          i = i + 1
          sAO_mod(iSym,l,i) = sAO(iSym,j,k)
        endif
      end do
    endif
  end do
end do
!...Open GUESSORB of the supermolecule for some reasonable virtual orbitals
open(11,file='GUESSORB',status='old')
allocate( c3(nSym,lDim,lDim) )
allocate( c4(nSym,lDim,lDim) )
allocate( occNu3(nSym,lDim) )
c3 = 0.0
c4 = 0.0
occNu3 = 0.0
mark = '#ORB'
64 read(11,'(A132)') line
if (line(1:4).ne.mark) goto 64
do iSym = 1, nSym
  if (nBas(iSym).ne.0) then
    do j = 1, nBas(iSym)
      read(11,'(A132)') line
      read(11,'(5E22.14)') (c3(iSym,j,k),k=1,nBas(iSym))
    end do
  endif
end do
close(11)
!...End of data collection

!...Replace the first vectors by the active orbitals of the
!     fragments to save them for the final vector file SUPERORB
do iSym = 1, nSym
  jOrb = 0
  offset1 = nInact1(iSym)+nFrozen1(iSym)
  offset2 = nInact2(iSym)+nFrozen2(iSym)
  do j = 1, nRas1_1(iSym)
    jOrb = jOrb + 1
    c3(iSym,jOrb,:) = 0.0d0
    do k = 1, nBas1(iSym)
      c3(iSym,jOrb,k) = c1(iSym,j+offset1,k)
    end do
  end do
  do j = 1, nRas1_2(iSym)
    jOrb = jOrb + 1
    c3(iSym,jOrb,:) = 0.0
    do k = 1, nBas2(iSym)
      c3(iSym,jOrb,k+nBas(iSym)-nBas2(iSym)) = c2(iSym,j+offset2,k)
    end do
  end do
  offset1 = offset1 + nRas1_1(iSym)
  offset2 = offset2 + nRas1_2(iSym)
  do j = 1, nRas2_1(iSym)
    jOrb = jOrb + 1
    c3(iSym,jOrb,:) = 0.0
    do k = 1, nBas1(iSym)
      c3(iSym,jOrb,k) = c1(iSym,j+offset1,k)
    end do
  end do
  do j = 1, nRas2_2(iSym)
    jOrb = jOrb + 1
    c3(iSym,jOrb,:) = 0.0
    do k = 1, nBas2(iSym)
      c3(iSym,jOrb,k+nBas(iSym)-nBas2(iSym)) = c2(iSym,j+offset2,k)
    end do
  end do
  offset1 = offset1 + nRas2_1(iSym)
  offset2 = offset2 + nRas2_2(iSym)
  do j = 1, nRas3_1(iSym)
    jOrb = jOrb + 1
    c3(iSym,jOrb,:) = 0.0
    do k = 1, nBas1(iSym)
      c3(iSym,jOrb,k) = c1(iSym,j+offset1,k)
    end do
  end do
  do j = 1, nRas3_2(iSym)
    jOrb = jOrb + 1
    c3(iSym,jOrb,:) = 0.0
    do k = 1, nBas2(iSym)
      c3(iSym,jOrb,k+nBas(iSym)-nBas2(iSym)) = c2(iSym,j+offset2,k)
    end do
  end do
end do
!...Corresponding orbital transformation in the frozen orbitals
do iSym = 1, nSym
  m1 = nFrozen1(iSym) + nFrozen2(iSym)
  if (m1 .eq. 0) cycle       !...No frozen in this irrep, try the next one
  m2 = nBas(iSym)
  m3 = nBas(iSym) - nBas2(iSym)
  if (debug) write (6,*) 'm1 m2 m3 : ',m1,m2,m3
  allocate ( aux(m1,m2) )
  allocate ( sMO(m1,m1) )
  allocate (   U(m1,m1) )
  allocate (  VT(m1,m1) )
  allocate (   V(m1,m1) )
  allocate (  sigma(m1) )
  allocate ( work(5*m1) )
  allocate ( c_1(m1,m2) )
  allocate ( c_2(m1,m2) )
  allocate ( c1c(m1,m2) )
  allocate ( c2c(m1,m2) )
  c_1 = 0.0
  c_2 = 0.0
  c1c = 0.0
  c2c = 0.0
!...Adding zeros to c1 and c2 to fit the dimension of the "super molecule"
  do j = 1, nFrozen1(iSym)
    do k = 1, nBas1(iSym)
      c_1(j,k) = c1(iSym,j,k)
    end do
  end do
  do j = 1, nFrozen2(iSym)
    do k = 1, nBas2(iSym)
      c_2(j+nFrozen1(iSym),k + m3) = c2(iSym,j,k)
    end do
  end do
  if ( debug ) then
    write(6,*)'c_1 '
    do j = 1, m1
      write(*,*) j
      write(6,'(20F10.5)') (c_1(j,k),k=1,m2)
    end do
    write(6,*)'c_2 '
    do j = 1, m1
      write(*,*) j
      write(6,'(20F10.5)') (c_2(j,k),k=1,m2)
    end do
  end if
  aux = 0.0
  sMO = 0.0
  do j = 1, m1
    do k = 1, m2
      do l = 1, m2
        aux(j,k) = aux(j,k) + c_2(j,l) * sAO(iSym,l,k)
      end do
    end do
  end do
  do j = 1, m1
    do k = 1, m1
      do l = 1, m2
        sMO(j,k) = sMO(j,k) + aux(j,l) * c_1(k,l)
      end do
    end do
  end do
  if (debug) then
    write(*,*) 'frozen MO overlap'
    do j = 1, m1
      write(6,'(20F10.5)')(sMO(j,k),k=1,m1)
    end do
    write(6,*)
  end if
! Singular value decomposition of sMO
  call dgesvd('A','A',m1,m1,sMO,m1,sigma,U,m1,VT,m1,work,5*m1,iRC)
  V = Transpose(VT)
  if ( debug ) then
    write(6,*)'Singular Value Decomposition'
    write(6,*) 'left sigular matrix, U'
    do j = 1, m1
      write(6,'(20F10.5)') (U(j,k),k = 1, m1)
    end do
    write(6,*) 'right sigular matrix, V'
    do j = 1, m1
      write(6,'(20F10.5)') (V(j,k),k = 1, m1)
    end do
  end if
  write(6,*) 'Singular values, sigma'
  write(6,'(10F12.6)') (sigma(j),j = 1, m1)
  write(6,*)
  do j = 1, m1
    do k = 1, m2
      do l = 1, m1
        c1c(j,k) = c1c(j,k) + c_1(l,k) * V(l,j)
        c2c(j,k) = c2c(j,k) + c_2(l,k) * U(l,j)
      end do
    end do
  end do
  if ( debug ) then
    write(6,*)
    write(6,*) 'Corresponding frozen orbitals'
    write(6,*) 'transformed c1'
    do j = 1, m1
      write(6,'(20F10.5)') (c1c(j,k),k = 1, m2)
    end do
    write(6,*) 'transformed c2'
    do j = 1, m1
      write(*,'(20F10.5)') (c2c(j,k),k = 1, m2)
    end do
  end if
!...Dump the corresponding orbitals on CORRORB.A and CORRORB.B
  do j = 1, nFrozen1(iSym)
    write(13,'(A,2I5)') '* ORBITAL',iSym,j
    write(13,'(5E22.14)')(c1(iSym,j,k),k = 1, nBas1(iSym))
    occNu1(iSym,j) = 2.0
  end do
  do j = 1, nFrozen2(iSym)
    write(14,'(A,2I5)') '* ORBITAL',iSym,j
    write(14,'(5E22.14)')(c2(iSym,j,k),k = 1, nBas2(iSym))
    occNu2(iSym,j) = 2.0
  end do
!...Remove the basis functions belonging to the dummy atoms
  l = 0
  do k = 1, nBas(iSym)
    if (.not.eliminate(iSym,k)) then
      l = l + 1
      do j = 1, nFrozen1(iSym)
        c1c(j,l) = c1c(j,k)
      end do
      do j = 1, nFrozen2(iSym)
        c2c(j,l) = c2c(j,k)
      end do
      do j = 1, nFrozen1(iSym)+nFrozen2(iSym)
        c_1(j,l) = c_1(j,k)
        c_2(j,l) = c_2(j,k)
      end do
    endif
  end do
!...Construct the frozen orbitals that belong to the overlapping fragment
  nBasFinal = nBas(isym) - n_elim(iSym)
  maxFrozen = max(nFrozen1(iSym),nFrozen2(iSym))
  jOrb = 1
  to_be_averaged_frozen(iSym) = 0
  do j = 1, maxFrozen
    if ( sigma(jOrb) .gt. svdThreshold ) then
      to_be_averaged_frozen(iSym) = to_be_averaged_frozen(iSym) + 1
      do k = 1, nBasFinal
        c4(iSym,jOrb,k) = (c1c(j,k) + c2c(j,k)) / sqrt(2+2*sigma(jOrb))
      end do
      jOrb = jOrb + 1
    else
      if ( j .le. nFrozen1(iSym) ) then
        c4(iSym,jOrb,:) = c1c(j,:)
        jOrb = jOrb + 1
      endif
      if ( j .le. nFrozen2(iSym) ) then
        c4(iSym,jOrb,:) = c2c(j,:)
        jOrb = jOrb + 1
      endif 
    endif
  enddo
  sMO = 0.0d0
  do j = 1, to_be_averaged_frozen(iSym)
    do k = 1, nBasFinal
      do l = 1, nBasFinal
        aux(j,k) = aux(j,k) + c4(iSym,j,l) * sAO_mod(iSym,l,k)
      enddo
    enddo
  enddo
  do j = 1, to_be_averaged_frozen(iSym)
    do k = 1, nFrozen1(iSym)
      do l = 1, nBasFinal
        sMO(j,k) = sMO(j,k) + aux(j,l) * c_1(k,l)
      end do
    end do
    do k = nFrozen1(iSym) + 1, nFrozen1(iSym) + nFrozen2(iSym)
      do l = 1, nBasFinal
        sMO(j,k) = sMO(j,k) + aux(j,l) * c_2(k,l)
      enddo
    enddo
  enddo
!...Select the other frozen orbitals from the non-transformed ones by overlap
!   with the corresponding frozen orbitals of the overlapping fragment. This
!   might require an improved selection criterion instead of the numeric that we
!   have now. For example, exclude the ones with the largest overlap
  jOrb = to_be_averaged_frozen(iSym)
  do j = 1, nFrozen1(iSym)
    new = .true.
    do k = 1, to_be_averaged_frozen(iSym)
      if (abs(sMO(k,j)) .gt. 0.3) new = .false.
    end do
    if (new) then
      jOrb = jOrb +1 
      c4(iSym,jOrb,:) = c_1(j,:)
    endif
  end do
  do j = nFrozen1(iSym) + 1, nFrozen1(iSym) + nFrozen2(iSym)
    new = .true.
    do k = 1, to_be_averaged_frozen(iSym)
      if (abs(sMO(k,j)) .gt. 0.3) new = .false.
    end do
    if (new) then
      jOrb = jOrb +1 
      c4(iSym,jOrb,:) = c_2(j,:)
    endif
  end do
  if (debug) then
    do i = 1, nFrozen1(iSym) + nFrozen2(iSym) - to_be_averaged_frozen(iSym)
      write(*,*) i
      do k = 1, nBasFinal
        if (abs(c4(iSym,i,k)) .gt. 5e-2) write(*,'(a,F15.6)')basLabel_mod(iSym,k),c4(iSym,i,k)
      end do
    end do
  endif
    
! deallocate for next symmetry
  deallocate ( aux )
  deallocate ( sMO )
  deallocate (   U )
  deallocate (  VT )
  deallocate (   V )
  deallocate ( sigma )
  deallocate ( work )
  deallocate ( c_1 )
  deallocate ( c_2 )
  deallocate ( c1c )
  deallocate ( c2c )
end do
 
!...Next, the inactive orbitals 
!...Construct overlap matrix of the MOs, sMO = (c2)^T x sAO x c1 
do iSym = 1, nSym
  m1 = nInact1(iSym)+nInact2(iSym)
  m2 = nBas(iSym)
  m3 = nBas(iSym)-nBas2(iSym)
  if (debug) write (6,*) 'm1 m2 m3 : ',m1,m2,m3
  allocate ( aux(m1,m2) )
  allocate ( sMO(m1,m1) )
  allocate (   U(m1,m1) )
  allocate (  VT(m1,m1) )
  allocate (   V(m1,m1) )
  allocate (  sigma(m1) )
  allocate ( work(5*m1) )
  allocate ( c_1(m1,m2) )
  allocate ( c_2(m1,m2) )
  allocate ( c1c(m1,m2) )
  allocate ( c2c(m1,m2) )
  c_1 = 0.0
  c_2 = 0.0
  c1c = 0.0
  c2c = 0.0
! Adding zeros to c1 and c2 to fit the dimension of the "super molecule"
  do j = 1, nInact1(iSym)
    do k = 1, nBas1(iSym)
      c_1(j,k) = c1(iSym,j+nFrozen1(iSym),k)
    end do
  end do
  do j = 1, nInact2(iSym)
    do k = 1, nBas2(iSym)
      c_2(j+nInact1(iSym),k + m3) = c2(iSym,j+nFrozen2(iSym),k)
    end do
  end do
  if ( debug ) then
    write(6,*)'c_1 '
    do j = 1, m1
      write(6,'(20F10.5)') (c_1(j,k),k=1,m2)
    end do
    write(6,*)'c_2 '
    do j = 1, m1
      write(6,'(20F10.5)') (c_2(j,k),k=1,m2)
    end do
  end if
  aux = 0.0
  sMO = 0.0
  do j = 1, m1
    do k = 1, m2
      do l = 1, m2
        aux(j,k) = aux(j,k) + c_2(j,l) * sAO(iSym,l,k) 
      end do
    end do
  end do
  do j = 1, m1
    do k = 1, m1
      do l = 1, m2
        sMO(j,k) = sMO(j,k) + aux(j,l) * c_1(k,l)
      end do
    end do
  end do
  if (debug) then
    write(*,*) 'MO overlap'
    do j = 1, m1
      write(6,'(20F10.5)')(sMO(j,k),k=1,m1)
    end do
    write(6,*)
  end if
! Singular value decomposition of sMO
  call dgesvd('A','A',m1,m1,sMO,m1,sigma,U,m1,VT,m1,work,5*m1,iRC)
  V = Transpose(VT)
  if ( debug ) then
    write(6,*)'Singular Value Decomposition'
    write(6,*) 'left sigular matrix, U'
    do j = 1, m1
      write(6,'(20F10.5)') (U(j,k),k = 1, m1)
    end do
    write(6,*) 'right sigular matrix, V'
    do j = 1, m1
      write(6,'(20F10.5)') (V(j,k),k = 1, m1)
    end do
  end if
  write(6,*) 'Singular values, sigma'
  write(6,'(10F12.6)') (sigma(j),j = 1, m1)
  write(6,*)
  do j = 1, m1
    do k = 1, m2
      do l = 1, m1
        c1c(j,k) = c1c(j,k) + c_1(l,k) * V(l,j)
        c2c(j,k) = c2c(j,k) + c_2(l,k) * U(l,j)
      end do
    end do
  end do
  if ( debug ) then
    write(6,*)
    write(6,*) 'Corresponding orbitals'
    write(6,*) 'transformed c1'
    do j = 1, m1
      write(6,'(20F10.5)') (c1c(j,k),k = 1, m2)
    end do
    write(6,*) 'transformed c2'
    do j = 1, m1
      write(*,'(20F10.5)') (c2c(j,k),k = 1, m2)
    end do
  end if
! Dump the corresponding orbitals on CORRORB.A and CORRORB.B
  do j = 1, nInact1(iSym)
    write(13,'(A,2I5)') '* ORBITAL',iSym,j+nFrozen1(iSym)
    write(13,'(5E22.14)')(c1c(j,k),k = 1, nBas1(iSym))
    occNu1(iSym,j+nFrozen1(iSym)) = sigma(j)
  end do
  do j = nFrozen1(iSym)+nInact1(iSym)+1,nBas1(iSym)
    write(13,'(A,2I5)') '* ORBITAL',iSym,j
    write(13,'(5E22.14)')(c1(iSym,j,k),k = 1, nBas1(iSym))
  end do
  do j = 1, nInact2(iSym)
    write(14,'(A,2I5)') '* ORBITAL',iSym,j+nFrozen2(iSym)
    write(14,'(5E22.14)')(c2c(j,k),k = m3+1, m2)
    occNu2(iSym,j+nFrozen2(iSym)) = sigma(j)
  end do
  do j = nFrozen2(iSym)+nInact2(iSym)+1, nBas2(iSym)
    write(14,'(A,2I5)') '* ORBITAL',iSym,j
    write(14,'(5E22.14)')(c2(iSym,j,k),k = 1, nBas2(iSym))
  end do
! Check that we have <phi_A,i | phi_B,j > =  sigma_ij (only non-zero for i = j)
  if ( debug ) then
    aux = 0.0
    do j = 1, m1 
      do k = 1, m2 
        do l = 1, m2 
          aux(j,k) = aux(j,k) + c2c(j,l) * sAO(iSym,l,k)
        end do
      end do
    end do
    sMO = 0.0
    do j = 1, m1 
      do k = 1, m1 
        do l = 1, m2
          sMO(j,k) = sMO(j,k) + aux(j,l) * c1c(k,l)
        end do
      end do
    end do
    write(6,*) 'Corresponding orbitals overlap'
    do j = 1, m1 
      write(6,'(20F10.5)')(sMO(j,k),k=1,m1)
    end do
    write(*,*)
  end if
! Before writing the corresponding orbitals to SUPERORB (unit 12), we first
! eliminate (if needed) the basis function(s) associated to atoms not present
! in the final supermolecule (typically saturating hydrogens)

  if ( debug ) then
    write(*,*) 'basis functions to be eliminated :'
    do k = 1, nBas(iSym)
      write(*,*) k,eliminate(iSym,k)
    end do
  endif

  l = 0
  if ( sum(n_elim) .gt. 0 ) then
    write(*,'(I4,A)') sum(n_elim),' basis functions have been eliminated'
  end if
  if ( normalize ) write(*,*) 'Corresponding orbitals of the supermolecule are normalized'
  do k = 1, nBas(iSym)
    if ( .not.eliminate(iSym,k) ) then
      l = l + 1
      do j = 1, nInact1(iSym)
        c1c(j,l) = c1c(j,k)
      end do
      do j = 1, nInact2(iSym)
        c2c(j,l) = c2c(j,k)
      end do
      do j = 1, nBas(iSym)
        c3(iSym,j,l) = c3(iSym,j,k)
      end do
    end if
  end do

  if ( debug ) then
    write(*,*) 'c1c and c2c after eliminating some basis functions'
    do k = 1, nBas(iSym)-n_elim(iSym)
      write(*,'(I4,2x,A,30F15.8)')k,basLabel_mod(iSym,k),(c1c(j,k),c2c(j,k),j=1,min(nInact1(iSym),nInact2(iSym)))
    end do
    write(*,*) 'norm of c1c'
    do i = 1, nInact1(iSym)
      norm = 0.0d0
      do j = 1, nBas(iSym)-n_elim(iSym)
        do k = 1, nBas(iSym)-n_elim(iSym)
          norm = norm + c1c(i,j) * c1c(i,k) * sAO_mod(iSym,j,k)
        end do
      end do
      write(*,*) i,norm
    end do
    write(*,*) 'norm of c2c'
    do i = 1, nInact2(iSym)
      norm = 0.0d0
      do j = 1, nBas(iSym)-n_elim(iSym)
        do k = 1, nBas(iSym)-n_elim(iSym)
          norm = norm + c2c(i,j) * c2c(i,k) * sAO_mod(iSym,j,k)
        end do
      end do
      write(*,*) i,norm
    end do
    write(*,*) 'How does c3 look like?'
    do k = 1, nBas(iSym)-n_elim(iSym)
      write(*,'(I4,2x,A,30F15.8)')k,basLabel_mod(iSym,k),(c3(iSym,j,k),j=1,nAct1(iSym)+nAct2(iSym))
    end do
  endif

! dump corresponding orbitals in orbital file of the supermolecule 
! taking the average for those orbitals with a sigma larger than the threshold
  nBasFinal = nBas(iSym) - n_elim(iSym)
  nFrozen = nFrozen1(iSym) + nFrozen2(iSym) - to_be_averaged_frozen(iSym)
  nAct = nAct1(iSym) + nAct2(iSym)
  to_be_averaged(iSym) = 0
  do j = 1, m1
    if ( sigma(j) .ge. svdThreshold ) to_be_averaged(iSym) = to_be_averaged(iSym) + 1
  end  do
  nOcc(isym) = m1 - to_be_averaged(iSym) + nAct1(iSym) + nAct2(iSym)
  if (nFrozen .ne. 0) then
    write (6,'(I4,A,I4,A)') nFrozen,' frozen orbitals and', &
                   m1 - to_be_averaged(iSym),' inactive orbitals will be dumped on SUPERORB'
  else
    write (6,'(I4,A)') m1 - to_be_averaged(iSym),' inactive orbitals will be dumped on SUPERORB'
  end if
  write(6,*) 'Please, check!!  If not correct, change --svdThreshold-- and run again'
  write(6,*)
  maxInact = max(nInact1(iSym),nInact2(iSym))
  jOrb = 1
  do j = 1, maxInact
    if ( sigma(jOrb) .gt. svdThreshold ) then
      do k = 1, nBasFinal
        c4(iSym,jOrb+nFrozen,k) = (c1c(j,k) + c2c(j,k)) / sqrt(2+2*sigma(jOrb))
      end do
      jOrb = jOrb + 1
    else
      if (j .le. nInact1(iSym)) then
        c4(iSym,jOrb+nFrozen,:) = c1c(j,:)
        jOrb = jOrb + 1
      endif
      if (j .le. nInact2(iSym)) then
        c4(iSym,jOrb+nFrozen,:) = c2c(j,:)
        jOrb = jOrb + 1
      endif
    endif
  end do
  do j = 1, nAct
    c4(iSym,jOrb+nFrozen,:) = c3(iSym,j,:)
    jOrb = jOrb + 1
  end do
  do j = jOrb+nFrozen, nBasFinal
    c4(iSym,j,:) = c3(iSym,j,:)
  end do
  if ( normalize ) then
    do j = 1, nBasFinal
      norm = 0.0d0
      do k = 1, nBasFinal
        do l = 1, nBasFinal
          norm = norm + c4(iSym,j,k) * c4(iSym,j,l) * sAO_mod(iSym,k,l)
        end do
      end do
      if (debug) write(*,'(i4,F15.8)') j,norm
      norm = 1.0d0/sqrt(norm)
      do k = 1, nBasFinal
        c4(iSym,j,k) = c4(iSym,j,k) * norm
      end do
      if ( debug ) then
        norm = 0.0d0
        do k = 1, nBasFinal
          do l = 1, nBasFinal
            norm = norm + c4(iSym,j,k) * c4(iSym,j,l) * sAO_mod(iSym,k,l)
          end do
        end do
        write(*,'(i4,f14.7)') j,norm
      endif
    end do
  endif
  do j = 1, nBasFinal
    write(12,'(A,2I5)')'* ORBITAL',iSym,j
    write(12,'(5E22.14)')(c4(iSym,j,k),k=1,nBasFinal)
  end do
! deallocate for next symmetry
  deallocate ( aux )
  deallocate ( sMO )
  deallocate (   U )
  deallocate (  VT )
  deallocate (   V )
  deallocate ( sigma )
  deallocate ( work )
  deallocate ( c_1 )
  deallocate ( c_2 )
  deallocate ( c1c )
  deallocate ( c2c )
end do
! construct the occupation numbers of the supermolecule
do iSym = 1, nSym
  start = 1
  finish = nFrozen1(iSym)+nFrozen2(iSym)+nInact1(iSym)+nInact2(iSym)-to_be_averaged_frozen(iSym)-to_be_averaged(iSym)
  do i = start, finish
    occNu3(iSym,i) = 2.0
  end do
  start = finish + 1
  finish = start + nAct1(iSym) - 1
  iCounter = 0
  do i = start, finish
    iCounter = iCounter + 1
    occNu3(iSym,i) = occNu1(iSym,nFrozen1(iSym)+nInact1(iSym)+iCounter)
  end do
  start = finish + 1
  finish = start + nAct2(iSym) - 1
  iCounter = 0
  do i = start,finish
    iCounter = iCounter + 1
    occNu3(iSym,i) = occNu2(iSym,nFrozen2(iSym)+nInact2(iSym)+iCounter)
  end do
  do i = finish + 1, nBas(iSym) - n_elim(iSym)
    occNu3(iSym,i) = 0.0
  end do
end do
! writing the occupation numbers
write(6,'(A,8I4)')'Number of basis functions of the supermolecule: ',   &
                   (nBas(iSym)-n_elim(iSym),iSym = 1, nSym)   
write(12,'(A4)')'#OCC'
write(12,'(A20)') '* OCCUPATION NUMBERS'
write(13,'(A4)')'#OCC'
write(13,'(A20)') '* OCCUPATION NUMBERS'
write(14,'(A4)')'#OCC'
write(14,'(A20)') '* OCCUPATION NUMBERS'
do iSym = 1, nSym
  write(12,'(5E22.14)')(occNu3(iSym,j),j=1,nBas(iSym)-n_elim(iSym))
  write(13,'(5E22.14)')(occNu1(iSym,j),j=1,nBas1(iSym))
  write(14,'(5E22.14)')(occNu2(iSym,j),j=1,nBas2(iSym))
end do
!...Construct the INDEX for the supermolecule
allocate(orblabel12(nSym,maxval(nBas)))
orblabel12 = ' '
do iSym = 1, nSym
  first = 1
  last = nFrozen1(iSym)+nFrozen2(iSym)
  do j = first, last
    orblabel12(iSym,j) = 'f'
  end do
  first = last + 1
  last = last + nInact1(iSym)+nInact2(iSym) - to_be_averaged(iSym)
  do j = first, last
    orblabel12(iSym,j) = 'i'
  end do
  first = last + 1
  last = last + nAct1(iSym)+nAct2(iSym)
  do j = first, last
    orblabel12(iSym,j) = '2'
  end do
  first = last + 1
  do j = first, nBas(iSym) - n_elim(iSym)
    orblabel12(iSym,j) = 's'
  end do
end do
write(12,'(A6)')'#INDEX'
do iSym = 1, nSym
  nB = nBas(iSym) - n_elim(iSym)
  write(12,'(A12)')'* 1234567890'
  if (nB .gt. 0 .and. nOcc(iSym) .gt. 0) then
    do k = 1, nB, 10
      if ( k + 9 .le. nB) then
        write(12,602)mod(int(k/10),10),(orblabel12(iSym,kk),kk=k,k+9)
      else
        write(12,602)mod(int(k/10),10),(orblabel12(iSym,kk),kk=k,nB)
      end if
    end do
  end if
end do
602 format(I1,x,10A1)
nInact = sum(nFrozen1) + sum (nFrozen2) - sum(to_be_averaged_frozen)    &
       + sum(nInact1) + sum(nInact2) - sum(to_be_averaged) 
if (ciThreshold.ne.0.0) call ovlf_dets(nInact)
 
! clean up, deallocate all the rest
close(12)
close(13)
close(14)
deallocate (nBas)
deallocate (nBas1)
deallocate (nBas2)
deallocate (nOcc)
deallocate (c1)
deallocate (c2)
deallocate (c3)
deallocate (occNu1)
deallocate (occNu2)
deallocate (occNu3)
deallocate (nInact1)
deallocate (nInact2)
deallocate (nAct1)
deallocate (nAct2)
deallocate (sAO)
deallocate (n_elim,to_be_averaged)
deallocate (eliminate)
deallocate (basLabel)
deallocate (orbLabel,orblabel12)
end program overlapping_fragments

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

subroutine ovlf_lowdin(nbase,sbase,slow)
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


end subroutine ovlf_lowdin

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

subroutine ovlf_dets(inact_ovlf)
use ovlfdets_data

implicit none

external :: ovlf_read_dets,ovlf_determine_nci,ovlf_determine_maxcib,ovlf_sminop
external :: ovlf_clebsch_gordan

integer                  :: lfnovlf
integer                  :: i,j,k,kk
integer                  :: det1,det2,iact,jact,ms1,ms2
integer                  :: ndets,ndet_large
integer                  :: nact_ovlf,inact_ovlf
integer, allocatable     :: iocc_tmp(:,:)
integer, allocatable     :: occ1_ms(:,:,:),occ2_ms(:,:,:)
integer, allocatable     :: ndet1_ms(:),ndet2_ms(:)
integer, allocatable     :: iocc_ovlf(:,:)
real(kind=8)             :: prod,dnorm,ovlf_clebsch_gordan
real(kind=8),allocatable :: cicoef_tmp(:)
real(kind=8),allocatable :: ci1_ms(:,:),ci2_ms(:,:)
real(kind=8),allocatable :: civ_ovlf(:)
character (len=40)       :: determinant


call ovlf_read_dets
call ovlf_determine_nci
call ovlf_determine_maxcib

allocate(ndet1_ms(spinm(1)+1))
allocate(ndet2_ms(spinm(2)+1))
allocate(ci1_ms(nci,spinm(1)+1))
allocate(ci2_ms(nci,spinm(2)+1))
allocate(occ1_ms(maxnact,nci,spinm(1)+1))
allocate(occ2_ms(maxnact,nci,spinm(2)+1))

ndet1_ms = 0
ndet2_ms = 0
ci1_ms = 0.0
ci2_ms = 0.0
occ1_ms = 0
occ2_ms = 0


write(*,*)
write(*,'(A,I8)')'Number of determinants frag 1 : ',idetm(1)
write(*,'(A,I8)')'                       frag 2 : ',idetm(2)

! generate M_S components of fragment 1
allocate(iocc_tmp(nactm(1),idetm(1)))
allocate(cicoef_tmp(idetm(1)))
do j = 1, idetm(1)
  cicoef_tmp(j) = civm(j,1)
  do k = 1, nactm(1)
    iocc_tmp(k,j) = ioccm(k,j,1)
  end do
end do
call ovlf_sminop(cicoef_tmp,ci1_ms,iocc_tmp,occ1_ms,idetm(1),     &
             ndet1_ms,spinm(1),nactm(1),nci,maxnact)
deallocate(cicoef_tmp)
deallocate(iocc_tmp)
! generate M_S components of fragment 2
allocate(iocc_tmp(nactm(2),idetm(2)))
allocate(cicoef_tmp(idetm(2)))
do j = 1, idetm(2)
  cicoef_tmp(j) = civm(j,2)
  do k = 1, nactm(2)
    iocc_tmp(k,j) = ioccm(k,j,2)
  end do
end do
call ovlf_sminop(cicoef_tmp,ci2_ms,iocc_tmp,occ2_ms,idetm(2),     &
             ndet2_ms,spinm(2),nactm(2),nci,maxnact)
deallocate(cicoef_tmp)
deallocate(iocc_tmp)

nact_ovlf = nactm(1)+nactm(2)
allocate(civ_ovlf(maxcib))
allocate(iocc_ovlf(nact_ovlf,maxcib))
iocc_ovlf=0
k = 0
do i = 1, spinm(1) + 1
  ms1  = spinm(1) - 2*(i-1)
  do j = 1, spinm(2) + 1
    ms2 = spinm(2) - 2*(j-1)
    if ( ms1 + ms2 .eq. target_spin ) then
      do det1 = 1, ndet1_ms(i)
        do det2 = 1, ndet2_ms(j)
          dnorm = ovlf_clebsch_gordan(spinm(1)/2.d0,spinm(2)/2.d0,        &
               ms1/2.d0,ms2/2.d0,target_spin/2.d0,target_spin/2.d0)
          prod = ci1_ms(det1,i)*ci2_ms(det2,j)
          k = k + 1
          civ_ovlf(k) = prod * dnorm
          kk = 0
          do iact = 1, nactm(1)
            kk = kk + 1
            iocc_ovlf(kk,k) = occ1_ms(iact,det1,i)
          end do
          do jact = 1, nactm(2)
            kk = kk + 1
            iocc_ovlf(kk,k) = occ2_ms(jact,det2,j)
          end do
        end do
      end do
    end if
  end do
end do

ndets = k
dnorm = 0.0
do i = 1, ndets
  dnorm = dnorm + civ_ovlf(i)**2
end do
dnorm = 1/dsqrt(dnorm)
do i = 1, ndets
  civ_ovlf(i) = civ_ovlf(i) * dnorm
end do

!     Write to a new detfile

write(*,'(A,I12)')'Number of determinants supermolecule : ',ndets

ndet_large = 0
lfnovlf = 23
open(lfnovlf,file='detAB')
write(lfnovlf,'(i4)') inact_ovlf
do i=1,ndets
  if ( abs(civ_ovlf(i)) .gt. ciThreshold ) then
    ndet_large = ndet_large + 1
    determinant=''
    do j=1,nact_ovlf
      if(iocc_ovlf(j,i).eq.2)then
        determinant=trim(determinant)//'2'
      else if(iocc_ovlf(j,i).eq.1)then
        determinant=trim(determinant)//'a'
      else if(iocc_ovlf(j,i).eq.-1)then
        determinant=trim(determinant)//'b'
      else if(iocc_ovlf(j,i).eq.0)then
        determinant=trim(determinant)//'0'
      endif
    enddo

    write(lfnovlf,'(e15.8,6x,A)')civ_ovlf(i),determinant
  end if
enddo
write(*,'(I8,A,E10.3)')ndet_large, ' with |c1*c2| larger than ',ciThreshold
close(lfnovlf)
deallocate(civm)
deallocate(ioccm)

end subroutine ovlf_dets

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ! 

subroutine ovlf_read_dets
use ovlfdets_data
implicit none

external :: ovlf_quicksort

integer                          :: i,j,k,dummy,lfndetA,lfndetB
integer                          :: istat,maxci
integer,dimension(2)             :: idetm_raw
real(kind=8)                     :: sdum
real(kind=8),allocatable         :: civm_raw(:,:)
character (len=80)               :: dumstr
character (len=255), allocatable :: iocc_raw(:,:),iocc(:,:)

lfndetA = 38
lfndetB = 39
open(lfndetA,file='vecdetA',err=998)
open(lfndetB,file='vecdetB',err=999)
read(lfndetA,*) dummy
read(lfndetB,*) dummy
idetm_raw = 0
idetm = 0
do
  read(lfndetA,*,iostat=istat)
  if ( istat .ne. 0 ) then
    rewind(lfndetA)
    exit
  else
    idetm_raw(1) = idetm_raw(1) + 1
  end if
end do
do
  read(lfndetB,*,iostat=istat)
  if ( istat .ne. 0 ) then
    rewind(lfndetB)
    exit
  else
    idetm_raw(2) = idetm_raw(2) + 1
  end if
end do
read(lfndetA,*) inactm(1)
read(lfndetA,*) sdum,dumstr
nactm(1)=len_trim(dumstr)
read(lfndetB,*) inactm(2)
read(lfndetB,*) sdum,dumstr
nactm(2)=len_trim(dumstr)
maxci = 0
maxnact = 0
do i=1,2
  maxci=max(maxci,idetm_raw(i))
  maxnact=max(maxnact,nactm(i))
enddo
allocate(civm(maxci,2))
allocate(iocc(maxci,2))
allocate(civm_raw(maxci,2))
allocate(iocc_raw(maxci,2))
allocate(ioccm(maxnact,maxci,2))
rewind(lfndetA)
rewind(lfndetB)

read(lfndetA,*) dummy
do j=1,idetm_raw(1)
  read(lfndetA,*) civm_raw(j,1),iocc_raw(j,1)
  iocc_raw(j,1)=adjustl(iocc_raw(j,1))
enddo
if (idetm_raw(1).gt.1) then
  call ovlf_quicksort(civm_raw(:,1),iocc_raw(:,1),idetm_raw(1))
end if
idetm(1) = 1
iocc(1,1) = iocc_raw(1,1)
civm(1,1) = civm_raw(1,1)
do i = 2, idetm_raw(1)
  if ( iocc_raw(i,1) .ne. iocc_raw(i-1,1) ) then
    idetm(1) = idetm(1) + 1
    iocc(idetm(1),1) = iocc_raw(i,1)
    civm(idetm(1),1) = civm_raw(i,1)
  else
    civm(idetm(1),1) = civm(idetm(1),1) + civm_raw(i,1)
  end if
end do
close(lfndetA)

read(lfndetB,*) dummy
do j=1,idetm_raw(2)
  read(lfndetB,*) civm_raw(j,2),iocc_raw(j,2)
  iocc_raw(j,2)=adjustl(iocc_raw(j,2))
enddo
if (idetm_raw(2).gt.1) then
  call ovlf_quicksort(civm_raw(:,2),iocc_raw(:,2),idetm_raw(2))
end if
idetm(2) = 1
iocc(1,2) = iocc_raw(1,2)
civm(1,2) = civm_raw(1,2)
do i = 2, idetm_raw(2)
  if ( iocc_raw(i,2) .ne. iocc_raw(i-1,2) ) then
    idetm(2) = idetm(2) + 1
    iocc(idetm(2),2) = iocc_raw(i,2)
    civm(idetm(2),2) = civm_raw(i,2)
  else
    civm(idetm(2),2) = civm(idetm(2),2) + civm_raw(i,2)
  end if
end do
close(lfndetB)

do i=1,2
  do j=1,idetm(i)
    dumstr=trim(adjustl(iocc(j,i)))
    do k=1,nactm(i)
      if(dumstr(k:k).eq.'2') then
        ioccm(k,j,i)=2
      elseif(dumstr(k:k).eq.'a') then
        ioccm(k,j,i)=1
      elseif(dumstr(k:k).eq.'b') then
        ioccm(k,j,i)=-1
      elseif(dumstr(k:k).eq.'0') then
        ioccm(k,j,i)=0
      else
        write(*,601) dumstr
601 format(' Weird occupation in ',a)
        write(*,602) j,i
602 format(' Determinant ',i4,' of state ',i4)
        stop 'Inconsistent occupation'
      endif
    enddo
  enddo
enddo

deallocate(civm_raw)
deallocate(iocc_raw)
deallocate(iocc)

return
998  write(*,*) 'vecdetA'
stop 'Error reading determinants'
999  write(*,*) 'vecdetB'
stop 'Error reading determinants'

end subroutine ovlf_read_dets 

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

subroutine ovlf_determine_nci
use ovlfdets_data

implicit none

external :: ovlf_pascal

integer     :: ovlf_pascal,ms
integer     :: state,det,act,alpha

integer,allocatable  :: spin(:)
integer,allocatable  :: ndets(:)


allocate( spin(2) )
do state = 1, 2
  spin(state) = 0                          ! spin of the monomer function
  do act = 1, nactm(state)
    if ( abs(ioccm(act,1,state)) .eq. 1 ) then
      spin(state) = spin(state) + ioccm(act,1,state)
    end if
  end do
end do
allocate( ndets(2) )
nci = 0
do state = 1, 2
  ndets(state) = 0
  if ( spin(state) .eq. 0 ) then
    ndets(state) = idetm(state)
  else
    do det = 1, idetm(state)
      alpha = 0
      do act = 1, nactm(state)
        if ( ioccm(act,det,state) .eq. 1 ) alpha = alpha + 1
      end do
      do ms = 1, spin(state) + 1               ! from ms_max to -ms_max
        ndets(state) = ndets(state) + ovlf_pascal(alpha+1,ms)
      end do
    end do
  end if
  nci = max(nci,ndets(state))
end do

deallocate(spin)
deallocate(ndets)
return
end subroutine ovlf_determine_nci

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = !

subroutine ovlf_determine_maxcib
use ovlfdets_data
implicit none

external :: ovlf_sminop

integer              :: iact,i,j,k
integer              :: det1,det2,cib
integer              :: ms1,ms2,offsetm,mspin_states
integer              :: minspin,maxspin
integer,allocatable  :: mspin(:),ndet_mspin(:,:),occ_mspin(:,:,:)
integer,allocatable  :: iocc_temp(:,:)

real (kind=8), allocatable  :: ci_mspin(:,:)
real (kind=8), allocatable  :: cicoef_temp(:)


spinm=0
do i=1,2
  do iact=1,nactm(i)
    if(abs(ioccm(iact,1,i)).eq.1) then
      spinm(i)=spinm(i)+ioccm(iact,1,i)
    endif
  enddo
enddo
maxcib=0

minspin=abs(spinm(1)-spinm(2) )
maxspin=    spinm(1)+spinm(2)
if (( maxspin.lt.target_spin ).or.( minspin.gt.target_spin))then
  write(*,602) spinm(1),spinm(2),target_spin
602 format(' Spins of fragments incompatible with total spin',/              &
           ' Spin(1) = ',i4,/,' Spin(2) = ',i4,/,' cannot be coupled')
  stop 'Incompatible spin state'
end if

cib=0
if((target_spin.eq.0) .and. (spinm(1).eq.0) .and. (spinm(2).eq.0)) then
  do det1=1,idetm(1)
    do det2=1,idetm(2)
      cib=cib+1
    enddo
  enddo
else
  mspin_states=spinm(1)+1+spinm(2)+1
  allocate(mspin(mspin_states)) ! ms values
  allocate(ndet_mspin(2,mspin_states))
  allocate(ci_mspin(nci,mspin_states))
  allocate(occ_mspin(maxnact,nci,mspin_states)) ! not really needed here, only for compatibility with ovlf_sminop call
  ndet_mspin=0
  offsetm=0
  do i=1,2
    allocate(iocc_temp(nactm(i),idetm(i)))
    allocate(cicoef_temp(idetm(i)))
    iocc_temp=0
    cicoef_temp=0
    do j=1,idetm(i)
      cicoef_temp(j)=civm(j,i)
      do k=1,nactm(i)
        iocc_temp(k,j)=ioccm(k,j,i)
      end do
    end do
    call ovlf_sminop(cicoef_temp,ci_mspin(1,offsetm+1),iocc_temp,occ_mspin(1,1,offsetm+1),     &
                idetm(i),ndet_mspin(i,:),spinm(i),nactm(i),nci,maxnact)
    deallocate(cicoef_temp)
    deallocate(iocc_temp)
    mspin(offsetm+1)=spinm(i)
    do j=2,spinm(i)+1
      mspin(offsetm+j)=mspin(offsetm+j-1)-2
    enddo
    offsetm=offsetm+spinm(i)+1
  enddo

  do j=1,spinm(1)+1
    ms1 =spinm(1)-2*(j-1)
    do k=1,spinm(2)+1
      ms2=spinm(2)-2*(k-1)
      if(ms1+ms2.eq.target_spin) then
        do det1=1,ndet_mspin(1,j)
          do det2=1,ndet_mspin(2,k)
            cib=cib+1
          enddo
        enddo
      endif
    enddo
  enddo
  deallocate(mspin,ndet_mspin,ci_mspin,occ_mspin)
endif
if(cib.eq.0) then
  stop 'Incompatible spin state'
endif

maxcib=max(cib,maxcib)

return
end subroutine ovlf_determine_maxcib

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ! 

subroutine ovlf_sminop(cicoef_in,cicoef_out,occ_in,occ_out,        &
                       ndetin,ndetout,spin,nact,nci,maxnact)

implicit none

external ovlf_check_duplicate

integer          :: nact,spin,ndetin,nci,maxnact
integer          :: ndetout(spin+1),occ_in(nact,ndetin)
integer          :: occ_out(maxnact,nci,spin+1)
integer          :: i,j,k,l,m,idup
integer,allocatable :: checklist(:,:)
real (kind=8)    :: cicoef_in(ndetin),dnorm
real (kind=8)    :: cicoef_out(nci,spin+1)
logical          :: duplicated

ndetout=0
ndetout(1)=ndetin
do i=1,ndetout(1)
   cicoef_out(i,1)=cicoef_in(i)
   do j=1,nact
      occ_out(j,i,1)=occ_in(j,i)
   enddo
enddo
! generate ms states with the s minus operator
do i=2,spin+1
   do j=1,ndetout(i-1)
      do k=1,nact
         if(occ_out(k,j,i-1).eq.1)then
            allocate(checklist(nact,ndetout(i)))
            do l=1,ndetout(i)
               do m=1,nact
                  checklist(m,l)=occ_out(m,l,i)
               enddo
            enddo
            call ovlf_check_duplicate(checklist,occ_out(:,j,i-1),k,nact,ndetout(i),duplicated,idup)
            deallocate(checklist)
            if(.not.duplicated)then
               ndetout(i)=ndetout(i)+1
               do l=1,nact
                  occ_out(l,ndetout(i),i)=occ_out(l,j,i-1)
               enddo
               occ_out(k,ndetout(i),i)=-1
               cicoef_out(ndetout(i),i)=cicoef_out(j,i-1)
            else
              cicoef_out(idup,i)=cicoef_out(idup,i) + cicoef_out(j,i-1)
            endif
         endif
      enddo
   enddo
   dnorm=0.0d0
   do j=1,ndetout(i)
      dnorm=dnorm+cicoef_out(j,i)**2
   enddo
   do j=1,ndetout(i)
      cicoef_out(j,i)=cicoef_out(j,i)*(1.0d0/dsqrt(dnorm))
   enddo
enddo
! Normalize the MS(max) component
dnorm = 0.0
do j = 1, ndetout(1)
  dnorm = dnorm + cicoef_out(j,1)**2
end do
dnorm = 1.0/dsqrt(dnorm)
do j = 1, ndetout(1)
  cicoef_out(j,1) = cicoef_out(j,1)*dnorm
end do

return
end

! =============================================================================== !

subroutine ovlf_check_duplicate(list,det,korb,nact,ndetlist,dup,i_dup)

implicit none

integer :: nact,ndetlist,i,korb
integer :: list(nact,ndetlist),det(nact),newdet(nact)
integer :: i_dup
logical :: dup
do i=1,nact
   newdet(i)=det(i)
enddo
newdet(korb)=-1
dup=.false.
do i=1,ndetlist
   if(all(newdet.eq.list(:,i)))then
      dup=.true.
      i_dup=i
      goto 99
   else
      dup=.false.
   endif
enddo
99 continue
return
end

! =============================================================================== !

integer function ovlf_pascal(i,j)
implicit none
integer                :: i,j,ii,jj
integer, allocatable   :: table(:,:)

allocate ( table(i,i) )

table = 0
table(:,1) = 1
do ii = 2, i
  do jj = 2,ii
    table(ii,jj) = table(ii-1,jj-1)+table(ii-1,jj)
  end do
end do
ovlf_pascal = table(i,j)
deallocate(table)
return
end function ovlf_pascal

! =============================================================================== !

real(kind=8) function ovlf_clebsch_gordan(s1,s2,m1,m2,s,m)
!     Returns the requested Clebsch-Gordan coefficient
implicit none

external :: factorial,recurrence_factor,A0

integer           :: factorial
real(kind=8)      :: s1,s2,s
real(kind=8)      :: m1,m2,m
real(kind=8)      :: ispin,k,j
real(kind=8)      :: recurrence_factor,A0,cg_p1,cg_p2,cg

if((abs(m1).gt.s1).or.(abs(m2).gt.s2))then
   stop 'CLEBSCH-GORDAN:: Input data error'
else if( (m1+m2) .ne. m )then
   cg = 0.0d0
else if ((s.lt.abs(s1-s2)).or.(s.gt.(s1+s2)))then
   cg = 0.0d0
else
   k = factorial(int(2*s1))*factorial(int(2*s2))*factorial(int(s1+s2+m))*factorial(int(s1+s2-m))
   j = factorial(int(2*s1+2*s2)) * factorial(int(s1+m1))*factorial(int(s1-m1))*factorial(int(s2+m2))*factorial(int(s2-m2))
   cg = dsqrt(k/j)         ! starting point
   cg_p2 = cg
   ispin = s1 + s2         ! initial spin
   if(ispin.ne.s)then     ! s-1
      cg = cg*A0(ispin,m,s1,s2,m1,m2)/recurrence_factor(ispin,m,s1,s2)
      cg_p1 = cg
      ispin = ispin - 1
   endif
   do while (ispin .ne. s) ! start loop with s-2, need s-1 and s terms
      cg = (cg_p1*A0(ispin,m,s1,s2,m1,m2)-cg_p2*recurrence_factor(ispin+1,m,s1,s2))/recurrence_factor(ispin,m,s1,s2)
      ispin = ispin - 1
      cg_p2 = cg_p1
      cg_p1 = cg
   enddo
 endif
 ovlf_clebsch_gordan=cg
return
end

! =============================================================================== !

real(kind=8) function recurrence_factor(s,m,s1,s2)
implicit none

real(kind=8)   :: s,m,s1,s2

recurrence_factor = (s*s-m*m)*((s1+s2+1)**2-s*s)*(s*s-(s1-s2)**2)
recurrence_factor = recurrence_factor/(s*s*(4*s*s-1))
recurrence_factor = dsqrt(recurrence_factor)
return
end

! =============================================================================== !

real(kind=8) function A0(s,m,s1,s2,m1,m2)
implicit none

real(kind=8)   :: s,m,s1,s2,m1,m2

A0  = m1-m2+m*(s2*(s2+1)-s1*(s1+1))/(s*(s+1))
return
end

! =============================================================================== !

integer function factorial(n)
implicit none

integer :: n, j
if(n==0)then
   factorial = 1
else
   factorial = 1
   do j = 1, n
      factorial = factorial*j
   enddo
endif
return
end

! =============================================================================== !

subroutine ovlf_readin
use ovlfdets_data
implicit none

external :: ovlf_capitalize,ovlf_locate

integer, parameter                   :: nKeys = 5
integer                              :: jj,iKey,inorm

character (len=4)                    :: key
character (len=4), dimension(nKeys)  :: keyword
character (len=132)                  :: line

logical                              :: all_ok = .true.
logical, dimension(nkeys)            :: hit = .false.

data keyword /'SVDT','CITH','SPIN','NORM','FROZ'/

svdThreshold = 0.2
ciThreshold  = 1.0e-5
target_spin  = 1
normalize    = .true.
nFrozen1     = 0
nFrozen2     = 0

do while (all_ok)
  read(5,*,iostat=jj) line
  line = adjustl(line)
  key = line(1:4)
  call ovlf_capitalize(key)
  do iKey = 1, nKeys
    if ( key .eq. keyword(iKey) ) hit(iKey) = .true.
  end do
  if (  jj .lt. 0 ) all_ok = .false.
end do

do iKey = 1, nKeys
  if ( hit(iKey) ) then
    select case(iKey)
      case(1)
        call ovlf_locate('SVDT')
        read(*,*) svdThreshold
      case(2)
        call ovlf_locate('CITH')
        read(*,*) ciThreshold
      case(3)
        call ovlf_locate('SPIN')
        read(*,*) target_spin
      case(4)
        call ovlf_locate('NORM')
        inorm = 1
        read(*,*) inorm
        if ( inorm .eq. 0 ) normalize = .false.
      case(5)
        call ovlf_locate('FROZ')
        read(*,*) nFrozen1(:)
        read(*,*) nFrozen2(:)
    end select
  end if
end do
return
end subroutine ovlf_readin


! =============================================================================== !

subroutine ovlf_capitalize(string)
implicit none
integer      :: i
character(*) string

do i = 1, len(string)
  if (ichar(string(i:i)).gt.96) then
    string(i:i) = char(ichar(string(i:i))-32)
  endif
end do
return
end subroutine ovlf_capitalize

! =============================================================================== !

subroutine ovlf_locate(string)
implicit none
external :: ovlf_capitalize
character(4)   ::  string,string2
character(132) ::  line
rewind(5)
40 read(5,*) line
line=adjustl(line)
string2=line(1:4)
call ovlf_capitalize(string2)
if (string2.ne.string) goto 40
return
end subroutine ovlf_locate

! =============================================================================== !

recursive subroutine ovlf_quicksort(coef,occu,n)
implicit none

real (kind = 8)    :: coef(n)
character(len=255)  :: occu(n),pivot,aux2
real (kind = 8)    :: random
real (kind = 8)    :: aux1
integer            :: n,left,right,marker

if (n .gt. 1) then
  call random_number(random)
  pivot = occu(int(random*real(n-1))+1)
  left = 0
  right = n + 1
  do while (left .lt. right)
    right = right - 1
    do while ( occu(right) .lt. pivot )
      right = right - 1
    end do
    left = left + 1
    do while ( occu(left) .gt. pivot )
      left = left + 1
    end do
    if ( left .lt. right) then
      aux1 = coef(left)
      coef(left) = coef(right)
      coef(right) = aux1
      aux2 = occu(left)
      occu(left) = occu(right)
      occu(right) = aux2
    end if
  end do
  if (left .eq. right) then
    marker = left + 1
  else
    marker = left
  end if
  call ovlf_quicksort(coef(:marker-1),occu(:marker-1),marker-1)
  call ovlf_quicksort(coef(marker:),occu(marker:),n-marker+1)
end if
end subroutine ovlf_quicksort

! ===============================================================================
