! ===============================================================================
!  Program for rotating and translating a molecule to match as close as possible
!  a target molecule. The target molecule can be identical to the original
!  molecule; chemically the same, but with different internal coordinates:
!  and can even be chemically different. When available, the wave function of
!  the original molecule will also be rotated.
! 
!  grotate < input_file
!
!  input_file:
!  fragments
!    fragA  fragB  fragC
!  target
!    1  3   6
!    2  5   1
!    6  10  9
!  states
!   S0 S1 T1
!
!  fragA.xyz contains the coordinates of the original molecule, fragB.xyz and
!  fragC.xyz of the target molecule. Target atoms are given in pairs. Atom 1,
!  3 and 6 of A, B and C are translated to the origin; then atoms 2, 5 and 1 
!  are rotated on the positive part of the x-axis; next atoms 6, 10 and 9 are
!  rotated into the xy-plane (positive y-coordinate). Finally, fragment A is
!  rotated into the orientations of B and C using the negative of the angles
!  that were used to rotate B and C in their present orientation in the xy-plane.
!  When only the atoms should be modified and not the wave function, the keyword
!  'noinporb' must be added to the input
!
!  Required files: - input file
!                  - at least 2 xyz files (targets and original molecule, names
!                    constructed from the 'fragment' keyword.
!                  - orbital files (here fragA_S0.orb fragA_S1.orb fragA_T1.orb)
!                  - RUNFILE (additional info such as the basis set labels)
!
!  Output files:   - rotated wave functions, here fragB_S0.orb, fragB_S0.orb,
!                    fragB_T1.orb, fragC_S0.orb, fragC_S0.orb and fragC_T1.orb
!                  - fragB.xyz and fragC.xyz (coordinates of the rotated fragA
!                    molecule). Original xyzfiles can be saved using the 'KEEP'
!                    keyword
!
!  Coen de Graaf, URV/ICREA
!  November 2021
!
!    For now, only up to g-functions, but easily extendable to higher angular
!    moments in the subroutine rotate_vec
!    
! ===============================================================================

module grotate_files_data
implicit none
integer, dimension(2)          :: nAtoms,t1,t2,t3
integer, allocatable           :: tt1(:),tt2(:),tt3(:)
integer                        :: nStates,nFrags
real(kind=8),allocatable       :: x(:,:),y(:,:),z(:,:)
real(kind=8),allocatable       :: transx(:),transy(:),transz(:)
real(kind=8),allocatable       :: rotalpha(:),rotbeta(:),rotgamma(:)
real(kind=8)                   :: pi,zero,one
character(len=6),allocatable   :: label(:,:)
character(len=80),dimension(2) :: xyzfile
character(len=3),dimension(32) :: state
logical                        :: inporb,xyztarget,keep,rotate,translate
character(len=80),dimension(2) :: orbfile
character(len=80),allocatable  :: root(:),saveLabel(:)
end module grotate_files_data

module grotate_basis_set_data
implicit none
integer, allocatable              :: nCntr(:)
integer, parameter                :: LblL  = 6           !LenIN in OpenMolcas
integer, parameter                :: LblL8 = lblL + 8    !LenIN8 in OpenMolcas
character(len=LblL8), allocatable :: BasfnLbl(:)
end module grotate_basis_set_data

module grotate_orbital_data
implicit none
integer                         :: nBas,nSym
real (kind=8), allocatable      :: nOcc(:),vec(:,:),rotvec(:,:)
character (len=1), allocatable  :: orbLabel(:)
character(len=132)              :: title_line
end module grotate_orbital_data

module grotate_rotation_data
implicit none
integer                   :: lmax
real(kind=8),allocatable  :: rotmat(:,:,:)
real(kind=8)              :: phi,psi,theta
end module grotate_rotation_data

! ===============================================================================

program grotate
  use grotate_files_data
  use grotate_orbital_data
  implicit none

  external :: write_rotvec
  external :: Get_iScalar,nameRun,Get_iArray
  external :: grotate_readin,grotate_get_coordinates
  external :: grotate_readvec,grotate_get_basis_info
  external :: grotate_by_target,grotate_by_input
  
  real(kind=8),dimension(2) :: alpha, beta, delta
  integer                   :: j,iState,iFrag

  call SetMem('clear=off')
  
  zero  = 0.0
  one   = 1.0
  pi    = 4.0 * atan(one)
  alpha = 0.0
  beta  = 0.0
  delta = 0.0

  call grotate_readin

  if ( inporb ) then
    call NameRun('RUNFILE')
    call Get_iScalar('nSym',nSym)
    if ( nSym .ne. 1 ) then
      write(*,*) '  Symmetry is not implemented'
      write(*,*) 'Remove the symmetry elements and come back'
      stop
    endif
    call Get_iArray('nBas',nBas,nSym)              ! get the number of basis functions
    call grotate_get_basis_info
  endif
  
  write(xyzfile(1),'(a,a)') trim(root(1)),'.xyz'

  t1(1) = tt1(1)
  t2(1) = tt2(1)
  t3(1) = tt3(1)

  do iFrag=2,nFrags

    write(*,'(a,i2,a,i2,a,a)') 'Fragment ',iFrag-1,' of ',nFrags-1,': ',trim(root(iFrag))
    write(xyzfile(2),'(a,a)') trim(root(iFrag)),'.xyz'

    t1(2) = tt1(iFrag)
    t2(2) = tt2(iFrag)
    t3(2) = tt3(iFrag)

    do iState=1,nStates

      call grotate_get_coordinates(iFrag,iState)   

      if ( inporb ) then
        write(*,'(a,i2,a,i2,a,a)') 'State ',istate,' of ',nStates,': ',trim(state(iState))
      endif

      write(orbfile(1),'(a,a,a,a)') trim(root(1)),'_',trim(state(iState)),'.orb'
      write(orbfile(2),'(a,a,a,a)') trim(root(iFrag)),'_',trim(state(iState)),'.orb'

      if ( inporb ) call grotate_readvec

      if ( xyztarget ) then
        call grotate_by_target(iState)
      else 
        call grotate_by_input(iFrag)
      endif

      if ( inporb ) call write_rotvec

      if(allocated(orbLabel)) deallocate(orbLabel)
      if(allocated(vec)) deallocate(vec)
      if(allocated(rotvec)) deallocate(rotvec)
      if(allocated(nOcc)) deallocate(nOcc)

    enddo
    write(*,'(a,a)') "Writing new coordinates to ",trim(xyzfile(2))
    open(14,file=trim(xyzfile(2)))
    write(14,'(I5)') nAtoms(1)
    write(14,'(A)') 'Translated and rotated molecule'
    do j = 1, nAtoms(1)
      write(14,'(A6,3x,3F16.8)')label(1,j),x(1,j),y(1,j),z(1,j)
    end do
    close(14)
    if(allocated(x)) deallocate(x)
    if(allocated(y)) deallocate(y)
    if(allocated(z)) deallocate(z)
    if(allocated(label)) deallocate(label)
    write(*,*)
  enddo

end program

! ===============================================================================

subroutine grotate_readin
use grotate_files_data
implicit none

external :: grotate_locate,grotate_capitalize
integer, parameter                  :: nKeys = 7
integer                             :: i,iKey,jj,lc,bc

character(len=4), dimension(nKeys)  :: keyword
character(len=4)                    :: key
character(len=132)                  :: line,card

logical                             :: all_ok = .true.
logical, dimension(nkeys)           :: hit = .false.


data keyword /'FRAG','TARG','NOIN','STAT','TRAN','ROTA','KEEP' /
! FRAGments          - fragment labels, first comes the fragment that is to moved, followed by the target fragments
! TARGet atoms       - Atoms that are to be matched
! NOINporb           - No orbitals are provided, only coordinates will be processed
! STATes             - Labels of the electronic states that need to be rotated
! ROTAte             - Rotation angles are given in the input
! TRANslate          - Translation vector is given in the input
! KEEP               - Labels for saving a copy of the original target xyzfiles (Label.xyz), the original xyzfiles
!                      are overwritten by the rotated and/or translated coordinates
   
inporb = .true.
xyztarget = .true.
keep = .false.
nStates = 1
rotate = .false.
translate = .false.
do while (all_ok)
  read(5,*,iostat=jj) line
  line=adjustl(line)
  key = line(1:4)
  call grotate_capitalize(key)
  do iKey = 1, nKeys
    if ( key .eq. keyword(iKey) ) hit(iKey) = .true.
  end do
  if ( jj .lt. 0 ) all_ok = .false.
end do

if (hit(2) .and. (hit(5).or.hit(6))) then
  write(*,*) 'TARGet atoms is incompatible with TRANslation and/or ROTAtion'
  stop
endif

do iKey = 1, nKeys
  if ( hit(iKey) ) then
    select case(iKey)
      case(1)
        call grotate_locate('FRAG')
        read(*,'(a)') card
        nFrags=0
        card=adjustl(trim(card))
        lc=len(trim(card))
        bc=index(card,' ')
        do while(lc.ge.1)
          nFrags=nFrags+1
!          root=trim(card(1:bc-1))
          card=trim(card(bc+1:lc))
          lc=len(trim(card))
          bc=index(card,' ')
        enddo
        allocate(tt1(nFrags))
        allocate(tt2(nFrags))
        allocate(tt3(nFrags))
        allocate(root(nFrags))
        tt1 = 1     ! By default atoms 1, 2 and 3 are paired between original and target
        tt2 = 2
        tt3 = 3
        call grotate_locate('FRAG')
        read(*,'(a)') card
        nFrags=0
        card=adjustl(trim(card))
        lc=len(trim(card))
        bc=index(card,' ')
        do while(lc.ge.1)
          nFrags=nFrags+1
          root(nFrags)=trim(card(1:bc-1))
          card=trim(card(bc+1:lc))
          lc=len(trim(card))
          bc=index(card,' ')
        enddo
      case(2)
        call grotate_locate('TARG')
        read(*,*) (tt1(i),i=1,nFrags)
        read(*,*) (tt2(i),i=1,nFrags)
        read(*,*) (tt3(i),i=1,nFrags)
      case(3)
        call grotate_locate('NOIN')
        inporb = .false.
      case(4)
        call grotate_locate('STAT')
        read(*,'(a)') card
        nStates=0
        card=adjustl(trim(card))
        lc=len(trim(card))
        bc=index(card,' ')
        do while(lc.ge.1)
          nStates=nStates+1
          state(nStates)=trim(card(1:bc-1))
          card=trim(card(bc+1:lc))
          lc=len(trim(card))
          bc=index(card,' ')
        enddo
      case(5)
        call grotate_locate('TRAN')
        allocate(transx(nFrags))
        allocate(transy(nFrags))
        allocate(transz(nFrags))
        transx(1) = 0.0d0
        transy(1) = 0.0d0
        transz(1) = 0.0d0
        do i = 2, nFrags
          read(*,*) transx(i),transy(i),transz(i)
        enddo
        xyztarget = .false.
        translate = .true.
      case(6)
        call grotate_locate('ROTA')
        allocate(rotalpha(nFrags))
        allocate(rotbeta(nFrags))
        allocate(rotgamma(nFrags))
        rotalpha(1) = 0.0d0
        rotbeta(1)  = 0.0d0
        rotgamma(1) = 0.0d0
        do i = 2, nFrags
          read(*,*) rotalpha(i),rotbeta(i),rotgamma(i)
          rotalpha(i) = rotalpha(i) * pi / 180.0d0
          rotbeta(i)  = rotbeta(i)  * pi / 180.0d0
          rotgamma(i) = rotgamma(i) * pi / 180.0d0
          xyztarget = .false.
          rotate = .true.
        enddo
      case(7)
        call grotate_locate('KEEP')
        keep=.true.
        allocate(saveLabel(nFrags))
        read(*,'(a)') card
        nFrags=0
        card=adjustl(trim(card))
        lc=len(trim(card))
        bc=index(card,' ')
        do while(lc.ge.1)
          nFrags=nFrags+1
          saveLabel(nFrags)=trim(card(1:bc-1))
          card=trim(card(bc+1:lc))
          lc=len(trim(card))
          bc=index(card,' ')
        enddo
    end select
  endif
end do
end subroutine grotate_readin

! ===============================================================================

subroutine grotate_by_input(iFrag)
use grotate_files_data
use grotate_orbital_data
implicit none

external :: rotate_atoms,rotate_vec
integer  :: fragment,j,iFrag

if (rotate) then
  fragment = 1
  call rotate_atoms(fragment,rotalpha(iFrag),zero,zero)
  call rotate_atoms(fragment,zero,rotbeta(iFrag) ,zero)
  call rotate_atoms(fragment,zero,zero,rotgamma(iFrag))
  
  if ( inporb ) then
    call rotate_vec(rotalpha(iFrag),zero,zero)
    call rotate_vec(zero,rotbeta(iFrag) ,zero)
    call rotate_vec(zero,zero,rotgamma(iFrag))
  endif
endif

if (translate) then
  do j = 1, nAtoms(1)
    x(1,j) = x(1,j) + transx(iFrag)
    y(1,j) = y(1,j) + transy(iFrag)
    z(1,j) = z(1,j) + transz(iFrag)
  end do
  if ( inporb ) rotvec = vec
endif

end subroutine grotate_by_input

! ===============================================================================

subroutine grotate_by_target(iState)
use grotate_files_data
implicit none
external :: rotate_atoms,rotate_vec,grotate_put_on_x
external :: grotate_put_on_xy

real(kind=8),dimension(2) :: alpha, beta, delta
real(kind=8),dimension(2) :: xcenter,ycenter,zcenter
integer                   :: fragment,i,j,iState

do i = 1, 2   ! loop over fragments
  xcenter(i) = x(i,t1(i))
  ycenter(i) = y(i,t1(i))
  zcenter(i) = z(i,t1(i))
  do j = 1, nAtoms(i)
    x(i,j) = x(i,j) - xcenter(i)
    y(i,j) = y(i,j) - ycenter(i)
    z(i,j) = z(i,j) - zcenter(i)
  end do

  call grotate_put_on_x(i,alpha(i),beta(i))
  call grotate_put_on_xy(i,delta(i))

  if ( i .eq. 1 .and. inporb ) then
    call rotate_vec(alpha(i),zero,zero)
    call rotate_vec(zero,beta(i),zero)
    call rotate_vec(zero,zero,delta(i))
  endif
end do

if (iState .eq. 1) then
  write(*,*)'rotation angles'
  write(*,101)'alpha ',alpha(1),alpha(2),   &
      ' (',alpha(1)*180.0/pi,alpha(2)*180.0/pi,' )'
  write(*,101)'beta  ',beta(1),beta(2),     &
      ' (',beta(1)*180.0/pi,beta(2)*180.0/pi,' )'
  write(*,101)'gamma ',delta(1),delta(2),   &
      ' (',delta(1)*180.0/pi,delta(2)*180.0/pi,' )'
endif
101 format(4x,A,2F9.3,A,F6.1,3x,F6.1,A)

fragment = 1
call rotate_atoms(fragment,zero,zero,-delta(2))
call rotate_atoms(fragment,zero,-beta(2),zero)
call rotate_atoms(fragment,-alpha(2),zero,zero)

if ( inporb ) then
  call rotate_vec(zero,zero,-delta(2))
  call rotate_vec(zero,-beta(2),zero)
  call rotate_vec(-alpha(2),zero,zero)
endif

do j = 1, nAtoms(1)
  x(1,j) = x(1,j) + xcenter(2)
  y(1,j) = y(1,j) + ycenter(2)
  z(1,j) = z(1,j) + zcenter(2)
end do

end subroutine grotate_by_target

! ===============================================================================

subroutine grotate_get_coordinates(iFrag,iState)
use grotate_files_data
implicit none
integer             :: iAtom,fragment,iFrag,iState,nfragments
character(len=132)  :: filename

if (xyztarget) then
  nfragments = 2
else
  nfragments = 1
endif

do fragment = 1,nfragments 
  open(12,file=xyzfile(fragment))
  read(12,*) nAtoms(fragment)
  close(12)
end do
if(allocated(x)) deallocate(x)
if(allocated(y)) deallocate(y)
if(allocated(z)) deallocate(z)
if(allocated(label)) deallocate(label)
allocate( x(nfragments,maxval(nAtoms)) )
allocate( y(nfragments,maxval(nAtoms)) )
allocate( z(nfragments,maxval(nAtoms)) )
allocate( label(nfragments,maxval(nAtoms)) )
do fragment = 1, nfragments
  open(12,file=xyzfile(fragment))
  read(12,*)
  read(12,*)
  do iAtom = 1, nAtoms(fragment)
    read(12,*) label(fragment,iAtom),      &
        x(fragment,iAtom),y(fragment,iAtom),z(fragment,iAtom)
  end do
  close(12)
end do
if ( keep .and. iState .eq. 1) then
  write(filename,'(a,a)')trim(saveLabel(iFrag)),'.xyz'
  open(12,file=filename)
  write(12,'(i4)')nAtoms(2)
  write(12,*) 'original coordinates of ',trim(root(iFrag))
  do iAtom = 1, nAtoms(2)
    write(12,'(A6,3x,3F16.8)') label(2,iAtom),x(2,iAtom),y(2,iAtom),z(2,iAtom)
  end do
  close(12)
endif

end subroutine grotate_get_coordinates

! ===============================================================================

subroutine grotate_readvec
use grotate_files_data
use grotate_orbital_data
implicit none

integer                :: j,k
character(len=6)       :: mark
character(len=132)     :: line

if(.not.allocated(orbLabel)) allocate( orbLabel(nBas))
if(.not.allocated(vec)) allocate( vec(nBas,nBas))
if(.not.allocated(rotvec)) allocate( rotvec(nBas,nBas))
if(.not.allocated(nOcc)) allocate( nOcc(nBas))

write(*,'(2a)')' Opening:',trim(orbfile(1))
open(35,file=trim(orbfile(1)))

   mark = '#INDEX'
46 read(35,'(A132)') line
   if (line(1:6).ne.mark) goto 46
   read(35,'(A132)') line
   orbLabel = ' '
   read(35,'(2x,10A)')(orbLabel(j),j=1,nBas)

   rewind(35)
   mark = '#ORB'
47 read(35,'(A132)') line
   if (line(1:4).ne.mark) goto 47
   do j = 1, nBas
     read(35,'(A132)') line
     read(35,'(5E22.14)') (vec(j,k),k=1,nBas)
   end do

   rewind(35)
   mark = '#OCC'
48 read(35,'(A132)') line
   if (line(1:4).ne.mark) goto 48
   read(35,*) line
   read(35,'(5e22.14)')(nOcc(j),j=1,nBas)

   rewind(35)
   mark = '#INFO'
49 read(35,'(A132)') line
   if (line(1:5).ne.mark) goto 49
   read(35,'(A132)') title_line
   title_line=trim(title_line)//' (rotated with rotharm)'

close(35)
end subroutine grotate_readvec

! ===============================================================================

subroutine grotate_get_basis_info
use grotate_basis_set_data
use grotate_orbital_data, only : nBas

implicit none

external Get_cArray,Get_iScalar

integer                           :: i,j,jj,nAtoms
character(len=LblL8), allocatable :: AtomLbl(:)


if(.not.allocated(AtomLbl)) allocate( AtomLbl(nBas)  )
if(.not.allocated(BasfnLbl)) allocate( BasfnLbl(nBas) )
if(.not.allocated(nCntr)) allocate( nCntr(nBas)   )

call Get_iScalar('Unique atoms',nAtoms)
Call Get_cArray('Unique Atom Names',AtomLbl,LblL*nAtoms)
Call Get_cArray('Unique Basis Name',BasfnLbl,nBas*LblL8)
nCntr = 0
jj = 1
nCntr(jj)=1
do j = 2, nBas
  if ( BasfnLbl(j)(9:11) .ne. BasfnLbl(j-1)(9:11)  .or. &
        BasfnLbl(j)(1:6) .ne. BasfnLbl(j-1)(1:6) ) then
    jj = j
    nCntr(jj) = 1
  else
    nCntr(jj) = nCntr(jj) + 1
  endif
end do

do i = 1, nBas
  if ( nCntr(i) .eq. 0 ) nCntr(i) = nCntr(i-1)
end do

deallocate(AtomLbl)
end subroutine grotate_get_basis_info

! ===============================================================================

subroutine grotate_put_on_x(i,alpha,beta)
use grotate_files_data
implicit none

integer        :: i,j
real(kind=8)   :: alpha,beta,distance
real(kind=8)   :: xnew,ynew,znew

! rotation around z
if (abs(y(i,t2(i))) .lt. 1e-14          &
            .and. abs(x(i,t2(i))) .lt. 1e-14) then
  alpha = 0.0
else
  alpha = atan(y(i,t2(i))/x(i,t2(i)))
endif
do j = 1, nAtoms(i)
  xnew = cos(alpha)*x(i,j) + sin(alpha)*y(i,j)
  ynew = cos(alpha)*y(i,j) - sin(alpha)*x(i,j)
  x(i,j) = xnew
  y(i,j) = ynew
end do
! rotation around y
if (abs(z(i,t2(i))) .lt. 1e-14          &
            .and. abs(x(i,t2(i))) .lt. 1e-14) then
  beta = 0.0
else
  beta = atan(z(i,t2(i))/x(i,t2(i)))
endif
do j = 1, nAtoms(i)
  xnew = cos(beta)*x(i,j) + sin(beta)*z(i,j)
  znew = cos(beta)*z(i,j) - sin(beta)*x(i,j)
  x(i,j) = xnew
  z(i,j) = znew
end do
if ( x(i,t2(i)) .lt. 0 ) then     ! ensure target atom two is on the positive x-axis
  if ( beta .lt. 0 ) then
    beta =  pi + beta
  else
    beta = -pi + beta
  endif
  do j = 1, nAtoms(i)
    x(i,j) = -x(i,j)
    z(i,j) = -z(i,j)
  end do
endif
! check the distance between the second target atoms
! should be small
if ( i .eq. 2 ) then
  distance = sqrt ( (x(1,t2(1))-x(2,t2(2)))**2 +      & 
                    (y(1,t2(1))-y(2,t2(2)))**2 +      &
                    (z(1,t2(1))-z(2,t2(2)))**2 )
  if (distance .gt. 0.3) then
    write(*,*) 'warning'
    write(*,*) 'distance between third target atoms ',distance
    write(*,*) 'please, check the results'
  endif
endif
end subroutine grotate_put_on_x

! ===============================================================================

subroutine grotate_put_on_xy(i,delta)
  use grotate_files_data
  implicit none

  integer        :: i,j
  real(kind=8)   :: delta,distance
  real(kind=8)   :: ynew,znew

! rotation around x
if (abs(z(i,t3(i))) .lt. 1e-14        &
          .and. abs(y(i,t3(i))) .lt. 1e-14) then
  delta = 0.0
else
  delta = atan(z(i,t3(i))/y(i,t3(i)))
endif
do j = 1, nAtoms(i)
  ynew = cos(delta)*y(i,j) + sin(delta)*z(i,j)
  znew = cos(delta)*z(i,j) - sin(delta)*y(i,j)
  y(i,j) = ynew
  z(i,j) = znew
end do
if ( y(i,t3(i)) .lt. 0 ) then    ! ensure target atom 3 has a positive y-coordinate
  if ( delta .lt. 0) then
    delta =  pi + delta
  else
    delta = -pi + delta
  endif
  do j = 1, nAtoms(i)
    y(i,j) = -y(i,j)
  end do
endif
! check the distance between the third target atoms
! should be small
if ( i .eq. 2 ) then
  distance = sqrt ( (x(1,t3(1))-x(2,t3(2)))**2 +    &
                    (y(1,t3(1))-y(2,t3(2)))**2 +    &
                    (z(1,t3(1))-z(2,t3(2)))**2 )
  if (distance .gt. 0.3) then
    write(*,*) 'warning'
    write(*,*) 'distance between third target atoms ',distance
    write(*,*) 'please, check the results'
  endif
endif
end subroutine grotate_put_on_xy

! ===============================================================================

subroutine grotate_capitalize(string)
implicit none
integer      :: i
character(*) string
do i = 1, len(string)
  if (ichar(string(i:i)).gt.96) then
    string(i:i) = char(ichar(string(i:i))-32)
  endif
end do
end subroutine grotate_capitalize

! ===============================================================================

subroutine grotate_locate(string)
implicit none
external :: grotate_capitalize
character(4)   ::  string,string2
character(132) ::  line
rewind(5)
40 read(5,*) line
line=adjustl(line)
string2=line(1:4)
call grotate_capitalize(string2)
if (string2.ne.string) goto 40
end subroutine grotate_locate

! ===============================================================================

subroutine rotate_harmonics
use grotate_rotation_data
implicit none

external :: u,v,w,capU,capV,capW
external :: init_rotmat

integer       :: i,j,l,m1,m2
real(kind=8)  :: r(3,3)
real(kind=8)  :: u,v,w,capU,capV,capW

allocate( rotmat(lmax,2*lmax+1,2*lmax+1) )
rotmat = 0.0

call init_rotmat(phi,psi,theta,r)
if ( lmax .eq. 1 ) then
  do i = 1, 3
    do j = 1, 3
      rotmat(1,i+lmax-1,j+lmax-1) = r(i,j)
    end do
  end do
else
! due to the definition of x,y,z in Ivanic and Ruedenberg, JCP 100, 6342 (1996).
  rotmat(1,lmax  ,lmax  ) = r(2,2)
  rotmat(1,lmax  ,lmax+1) = r(2,3)
  rotmat(1,lmax  ,lmax+2) = r(2,1)
  rotmat(1,lmax+1,lmax  ) = r(3,2)
  rotmat(1,lmax+1,lmax+1) = r(3,3)
  rotmat(1,lmax+1,lmax+2) = r(3,1)
  rotmat(1,lmax+2,lmax  ) = r(1,2)
  rotmat(1,lmax+2,lmax+1) = r(1,3)
  rotmat(1,lmax+2,lmax+2) = r(1,1)
endif

do l = 2,lmax
  do i = 1,2*l+1
    m1 = i - (l+1)
    do j = 1,2*l+1
      m2 = j - (l+1)
      rotmat(l,i+lmax-l,j+lmax-l) = u(l,m1,m2) * capU(l,m1,m2) + v(l,m1,m2) * capV(l,m1,m2) + &
                   w(l,m1,m2) * capW(l,m1,m2)
    end do
  end do
end do

!do l = 1, lmax
!  write(*,'(a,i4)')'Rotation matrix for l = ',l
!  do i = 1,2*lmax+1
!    write(*,'(11F14.7)')(rotmat(l,i,j),j=1,2*lmax+1)
!  end do
!  write(*,*)
!end do

end subroutine rotate_harmonics

! ===============================================================================

function capU(l,m1,m2) result(x)
implicit none

external :: P

integer       :: l,m1,m2
real(kind=8)  :: x,P

x = P(0,l,m1,m2)
!   write(*,'(A,2I4,F15.8)')'U',m1,m2,x
end function capU

! ===============================================================================

function capV(l,m1,m2) result(x)
implicit none

external :: delta,P

integer       :: l,m1,m2,delta
real(kind=8)  :: x,a,b,c,P

x=0
if ( m1 .eq. 0 ) then
  x = P(1,l,1,m2) + P(-1,l,-1,m2)
endif
if ( m1 .gt. 0 ) then
  a = 1+delta(m1,1)
  b = sqrt(a)
  c = 1-delta(m1,1)
  x = P(1,l,m1-1,m2) * b - P(-1,l,-m1+1,m2) * c
endif
if ( m1 .lt. 0 ) then
  a = 1-delta(m1,-1)
  b = 1+delta(m1,-1)
  c = sqrt(b)
  x = P(1,l,m1+1,m2) * a + P(-1,l,-m1-1,m2) * c
endif
!   write(*,'(A,2I4,F15.8)')'V',m1,m2,x
end function capV
       
! ===============================================================================
   
function capW(l,m1,m2) result(x)
implicit none

external :: P

integer       :: l,m1,m2
real(kind=8)  :: x,P

x=0
if ( m1 .gt. 0 ) then
  x = P(1,l,m1+1,m2) + P(-1,l,-m1-1,m2)
endif
if ( m1 .lt. 0 ) then
  x = P(1,l,m1-1,m2) - P(-1,l,-m1+1,m2)
endif
!   write(*,'(A,2I4,F15.8)')'W',m1,m2,x
end function capW

! ===============================================================================

function P(i,l,mu,m) result(x)
use grotate_rotation_data
implicit none

integer       :: i,l,mu,m
integer       :: ii,imu,im
real(kind=8)  :: x

x=0
ii  = i  + lmax+1
imu = mu + lmax+1
im  = m  + lmax+1
if ( abs(m) .lt. l) then
! x = rotmat(1,i,0) * rotmat(l-1,mu,m)
  x = rotmat(1,ii,lmax+1) * rotmat(l-1,imu,im)
endif
if ( m .eq. l ) then
! x = rotmat(1,i,1) * rotmat(l-1,mu,m-1) - rotmat(1,i,-1) * rotmat(l-1,mu,-m+1)
  x = rotmat(1,ii,lmax+2) * rotmat(l-1,imu,im-1) - rotmat(1,ii,lmax) *  &
                                             rotmat(l-1,imu,im-2*m+1)
endif
if ( m .eq. -l ) then
! x = rotmat(1,i,1) * rotmat(l-1,mu,m+1) + rotmat(1,i,-1) * rotmat(l-1,mu,-m-1)
  x = rotmat(1,ii,lmax+2) * rotmat(l-1,imu,im+1) + rotmat(1,ii,lmax) *  &
                                             rotmat(l-1,imu,im-2*m-1)
endif
end function P

! ===============================================================================

function u(l,m1,m2) result(x2)
implicit none

integer,intent(in)   :: l,m1,m2
real(kind=8)         :: num,denom,x1,x2

if (l .eq. abs(m2) ) then
  num   = (l+m1)*(l-m1)
  denom = 2*l*(2*l-1) 
  x1 = num / denom
  x2 = sqrt(x1)
else
  num   = (l+m1)*(l-m1)
  denom = (l+m2)*(l-m2)
  x1 = num / denom
  x2 = sqrt(x1)
endif
!   write(*,'(2I4,F25.15)')m1,m2,x2
end function u

! ===============================================================================

function v(l,m1,m2) result(x2)
implicit none

external :: delta

integer, intent(in)   :: l,m1,m2
real(kind=8)          :: num,denom,x1,x2
integer               :: delta
if (l .eq. abs(m2) ) then
  num   = (1+delta(0,m1))*(l+abs(m1)-1)*(l+abs(m1))
  denom = 2*l*(2*l-1) 
  x1 = num / denom
  x2 = 0.5 * sqrt(x1) * (1-2*delta(0,m1))
else
  num   = (1+delta(0,m1))*(l+abs(m1)-1)*(l+abs(m1))
  denom = (l+m2)*(l-m2) 
  x1 = num / denom
  x2 = 0.5 * sqrt(x1) * (1-2*delta(0,m1))
endif
!   write(*,'(2I4,F25.15)')m1,m2,x2
end function v

! ===============================================================================

function w(l,m1,m2) result(x2)
implicit none

external :: delta

integer, intent(in)   :: l,m1,m2
real(kind=8)          :: num,denom,x1,x2
integer               :: delta

if (l .eq. abs(m2) ) then
  num   = (l-abs(m1)-1)*(l-abs(m1))
  denom = 2*l*(2*l-1)
  x1 = num / denom
  x2 = -0.5 * sqrt(x1) * (1-delta(0,m1))
else
  num   = (l-abs(m1)-1)*(l-abs(m1))
  denom = (l+m2)*(l-m2)
  x1 = num / denom
  x2 = -0.5 * sqrt(x1) * (1-delta(0,m1))
endif
!   write(*,'(2I4,F25.15)')m1,m2,x2
end function w

! ===============================================================================

function delta(i,j) result(d)
implicit none

integer, intent(in)   :: i,j
integer               :: d

if ( i .eq. j ) then
  d = 1
else
  d = 0
endif
end function delta
   
! ===============================================================================

subroutine init_rotmat(a,b,c,r)
implicit none
real(kind=8)       :: a,b,c
real(kind=8)       :: r(3,3)

r(1,1) =  cos(a) * cos(b)
r(1,2) =  cos(a) * sin(b) * sin(c) - sin(a) * cos(c)
r(1,3) = -sin(b)
r(2,1) =  sin(a) * cos(b)
r(2,2) =  sin(a) * sin(b) * sin(c) + cos(a) * cos(c)
r(2,3) =  sin(a) * sin(b) * cos(c) - cos(a) * sin(c)
r(3,1) =  cos(a) * sin(b) * cos(c) + sin(a) * sin(c)
r(3,2) =  cos(b) * sin(c)
r(3,3) =  cos(b) * cos(c)

end subroutine init_rotmat

! ===============================================================================

subroutine rotate_atoms(i,a,b,c)
use grotate_files_data
implicit none

external :: init_rotmat

integer        :: i,j
real(kind=8)   :: a,b,c,xnew,ynew,znew
real(kind=8)   :: r(3,3)

call init_rotmat(a,b,c,r)

do j = 1, nAtoms(i)
  xnew = x(i,j) * r(1,1) + y(i,j) * r(2,1) + z(i,j) * r(3,1)
  ynew = x(i,j) * r(1,2) + y(i,j) * r(2,2) + z(i,j) * r(3,2)
  znew = x(i,j) * r(1,3) + y(i,j) * r(2,3) + z(i,j) * r(3,3)
  x(i,j) = xnew
  y(i,j) = ynew
  z(i,j) = znew
end do

end subroutine rotate_atoms

! ===============================================================================

subroutine write_rotvec
use grotate_files_data
use grotate_orbital_data
implicit none
integer   :: j,k

write(*,'(a,a)') ' Write ',trim(orbfile(2))
open(12,file=trim(orbfile(2)))
write(12,'(A11)') '#INPORB 2.2'
write(12,'(A5)')  '#INFO'
write(12,'(A)') title_line
write(12,'(3I8)')0,nSym,0
write(12,'(I8)') nBas
write(12,'(I8)') nBas
write(12,'(A4)') '#ORB'
do j = 1, nBas
  write(12,'(A12,I4)')'* ORBITAL  1',j
  write(12,'(5E22.14)') (rotvec(j,k),k=1,nBas)
end do
write(12,'(A4)') '#OCC'
write(12,'(A20)') '* Occupation Numbers'
write(12,'(5e22.14)')(nOcc(j),j=1,nBas)
write(12,'(A6)')'#INDEX'
write(12,'(A12)')'* 1234567890'
do j = 1, nBas, 10
  if ( j + 9 .le. nBas) then
    write(12,602)mod(int(j/10),10),(orbLabel(k),k=j,j+9)
  else
    write(12,602)mod(int(j/10),10),(orbLabel(k),k=j,nBas)
  endif
end do
602  format(I1,x,10A1)
deallocate(orbLabel)

end subroutine write_rotvec

! ===============================================================================

subroutine rotate_vec(a,b,c)
use grotate_orbital_data
use grotate_basis_set_data
use grotate_rotation_data
implicit none

external :: r_times_v,rotate_harmonics

integer        :: i,j
real(kind=8)   :: a,b,c

phi = a
psi = b
theta = c
do i = 1, nBas           ! loop over the orbitals
  do j = 1, nBas         ! loop over the basis functions
    if ( BasFnLbl(j)(9:9) .eq. 's' ) rotvec(i,j) = vec(i,j)
    if ( BasFnLbl(j)(9:10) .eq. 'px' ) then
      lmax = 1
      call rotate_harmonics
      call r_times_v(i,j)
      deallocate(rotmat)
    endif
    if ( BasFnLbl(j)(9:12) .eq. 'd02-' ) then
      lmax = 2
      call rotate_harmonics
      call r_times_v(i,j)
      deallocate(rotmat)
    endif
    if ( BasFnLbl(j)(9:12) .eq. 'f03-' ) then
      lmax = 3
      call rotate_harmonics
      call r_times_v(i,j)
      deallocate(rotmat)
    endif
    if ( BasFnLbl(j)(9:12) .eq. 'g04-' ) then
      lmax = 4
      call rotate_harmonics
      call r_times_v(i,j)
      deallocate(rotmat)
    endif
  end do
end do

end subroutine rotate_vec

! ===============================================================================

subroutine r_times_v(i,j)
use grotate_orbital_data
use grotate_basis_set_data
use grotate_rotation_data
implicit none

integer                    :: i,j,ll,mm
real(kind=8)               :: rv
real(kind=8),allocatable   :: v(:)


allocate(v(2*lmax+1))
v=0
do ll = 1, 2*lmax + 1
  v(ll) = vec( i,j+(ll-1)*nCntr(j) )
end do
do ll = 1, 2*lmax + 1
  rv = 0.0
  do mm = 1, 2*lmax + 1
    rv = rv + v(mm) * rotmat(lmax,mm,ll)
  end do
  rotvec(i,j+(ll-1)*nCntr(j)) = rv
  vec(i,j+(ll-1)*nCntr(j)) = rv
end do
deallocate(v)

end subroutine r_times_v

! ===============================================================================
