! ===============================================================================
!  Program for rotating and translating a molecule to match as close as possible
!  a target molecule. The target molecule can be identical to the original
!  molecule; chemically the same, but with different internal coordinates:
!  and can even be chemically different. When available, the wave function of
!  the original molecule will also be rotated.
! 
!  rotharm < input_file
!
!  input_file:
!  xyzfiles
!    fragA.xyz
!    fragB.xyz
!  target
!    1  3
!    2  5
!    6  10
!
!  fragA.xyz contains the coordinates of the target molecule, fragB.xyz of
!  the original molecule. Target atoms are given in pairs. Atom 1 and 3 of
!  A and B are translated to the origin; then atoms 2 and 5 are rotated on
!  the positive part of the x-axis; next atoms 6 and 10 are rotated into
!  the xy-plane (positive y-coordinate). Finally, fragment B is rotated
!  into the orientation of A using the negative of the angles that were 
!  used to rotate A in its present orientation.
!  When only the atoms should be modified and not the wave function, the
!  keyword 'noinporb' must be added to the input
!
!  Required files: - input file
!                  - 2 xyz files (target and original molecule, names
!                    specified in the input by the 'xyzfiles' keyword)
!                  - INPORB (wave function file)
!                  - RUNFILE (additional info such as the basis set labels)
!
!  Output files:   - ROTORB (rotated wave function)
!                  - transrot.xyz (coordinates of the translated and rotated
!                                      fragment B)
!
!  Coen de Graaf, URV/ICREA
!  November 2021
!
!    For now, only up to g-functions, but easily extendable to higher angular
!    moments in the subroutine rotate_vec
!    
! ===============================================================================

module useful_data
implicit none
integer, dimension(2)          :: nAtoms,t1,t2,t3
real(kind=8),allocatable       :: x(:,:),y(:,:),z(:,:)
real(kind=8)                   :: pi
character(len=6),allocatable   :: label(:,:)
character(len=80),dimension(2) :: xyzfile
logical                        :: inporb
end module useful_data

module basis_set_data
implicit none
integer, allocatable              :: nCntr(:)
integer, parameter                :: LblL  = 6           !LenIN in OpenMolcas
integer, parameter                :: LblL8 = lblL + 8    !LenIN8 in OpenMolcas
character(len=LblL8), allocatable :: BasfnLbl(:)
end module basis_set_data

module orbital_data
implicit none
integer                         :: nBas,nSym
real (kind=8), allocatable      :: nOcc(:),vec(:,:),rotvec(:,:)
character (len=1), allocatable  :: orbLabel(:)
end module orbital_data

module rotation_data
implicit none
integer                   :: lmax
real(kind=8),allocatable  :: rotmat(:,:,:)
real(kind=8)              :: phi,psi,theta
end module rotation_data

! ===============================================================================

program rotate
use useful_data
implicit none

real(kind=8),dimension(2) :: alpha, beta, delta
real(kind=8),dimension(2) :: xcenter,ycenter,zcenter
real(kind=8)              :: zero,one
integer                   :: i,j,fragment

zero  = 0.0
one   = 1.0
pi    = 4.0 * atan(one)
alpha = 0.0
beta  = 0.0
delta = 0.0
call readin
call get_coordinates
if ( inporb ) then
  call readvec
  call get_basis_info
end if
do i = 1, 2   ! loop over fragments
  xcenter(i) = x(i,t1(i))
  ycenter(i) = y(i,t1(i))
  zcenter(i) = z(i,t1(i))
  do j = 1, nAtoms(i)
    x(i,j) = x(i,j) - xcenter(i)
    y(i,j) = y(i,j) - ycenter(i)
    z(i,j) = z(i,j) - zcenter(i)
  end do

  call put_on_x(i,alpha(i),beta(i))
  call put_on_xy(i,delta(i))

  if ( i .eq. 2 .and. inporb ) then
    call rotate_vec(alpha(i),zero,zero)
    call rotate_vec(zero,beta(i),zero)
    call rotate_vec(zero,zero,delta(i))
  end if
end do

write(*,*)'rotation angles'
write(*,101)'alpha ',alpha(1),alpha(2),   &
    ' (',alpha(1)*180.0/pi,alpha(2)*180.0/pi,' )'
write(*,101)'beta  ',beta(1),beta(2),     &
    ' (',beta(1)*180.0/pi,beta(2)*180.0/pi,' )'
write(*,101)'gamma ',delta(1),delta(2),   &
    ' (',delta(1)*180.0/pi,delta(2)*180.0/pi,' )'
101 format(A,2F9.3,A,F6.1,3x,F6.1,A)

fragment = 2
call rotate_atoms(fragment,zero,zero,-delta(1))
call rotate_atoms(fragment,zero,-beta(1),zero)
call rotate_atoms(fragment,-alpha(1),zero,zero)

if ( inporb ) then
  call rotate_vec(zero,zero,-delta(1))
  call rotate_vec(zero,-beta(1),zero)
  call rotate_vec(-alpha(1),zero,zero)
endif

do j = 1, nAtoms(2)
  x(2,j) = x(2,j) + xcenter(1)
  y(2,j) = y(2,j) + ycenter(1)
  z(2,j) = z(2,j) + zcenter(1)
end do
open(14,file='transrot.xyz')
write(14,'(I5)') nAtoms(2)
write(14,'(A)') 'Translated and rotated molecule'
do j = 1, nAtoms(2)
  write(14,'(A6,3x,3F16.8)')label(2,j),x(2,j),y(2,j),z(2,j)
end do
close(14)
deallocate(x,y,z,label)
if ( inporb ) call write_rotvec

end program rotate

! ===============================================================================

subroutine readin
use useful_data
implicit none

integer, parameter                  :: nKeys = 3
integer                             :: iKey,jj

character(len=4), dimension(nKeys)  :: keyword
character(len=4)                    :: key
character(len=132)                  :: line

logical                             :: all_ok = .true.
logical, dimension(nkeys)           :: hit = .false.

data keyword /'XYZF','TARG','NOIN' /
! XYZFiles           - Two xyz files (one per line) with the original and target coordinates
! TARGet atoms       - Atoms that to be matched
! CHANge orientation - Change the orientation of the original molecule 
! NOINporb           - No orbitals are provided, only coordinates will be processed
   
t1(1) = 1     ! target atom of molecule 1 for centering 
t1(2) = 1     ! same for molecule 2
t2(1) = 2     ! target atom of molecule 1 to rotate on the x-axis
t2(2) = 2     ! same for molecule 2
t3(1) = 3     ! target atom of molecule 1 to rotate in the yz plane
t3(2) = 3     ! same for molecule 2
inporb = .true.
do while (all_ok)
  read(5,*,iostat=jj) line
  key = adjustl(line)
  call capitalize(key)
  do iKey = 1, nKeys
    if ( key .eq. keyword(iKey) ) hit(iKey) = .true.
  end do
  if ( jj .lt. 0 ) all_ok = .false.
end do

do iKey = 1, nKeys
  if ( hit(iKey) ) then
    select case(iKey)
      case(1)
        call locate('XYZF')
        read(*,*) xyzfile(1)
        read(*,*) xyzfile(2)
      case(2)
        call locate('TARG')
        read(*,*) t1(1),t1(2)
        read(*,*) t2(1),t2(2)
        read(*,*) t3(1),t3(2)
      case(3)
        call locate('NOIN')
        inporb = .false.
    end select
  end if
end do
end subroutine readin

! ===============================================================================

subroutine get_coordinates
use useful_data
implicit none
integer           :: iAtom,fragment

do fragment = 1, 2
  open(12,file=xyzfile(fragment))
  read(12,*) nAtoms(fragment)
  close(12)
end do
allocate( x(2,maxval(nAtoms)) )
allocate( y(2,maxval(nAtoms)) )
allocate( z(2,maxval(nAtoms)) )
allocate( label(2,maxval(nAtoms)) )
do fragment = 1, 2
  open(12,file=xyzfile(fragment))
  read(12,*)
  read(12,*)
  do iAtom = 1, nAtoms(fragment)
    read(12,*) label(fragment,iAtom),      &
        x(fragment,iAtom),y(fragment,iAtom),z(fragment,iAtom)
  end do
  close(12)
end do
end subroutine get_coordinates

! ===============================================================================

subroutine readvec
use orbital_data
implicit none

integer                :: j,k
character(len=6)       :: mark
character(len=132)     :: line

call NameRun('RUNFILE')
call Get_iScalar('nSym',nSym)
if ( nSym .ne. 1 ) then
  write(*,*) '  Symmetry is not implemented'
  write(*,*) 'Remove the symmetry elements and come back'
  stop
end if
Call Get_iArray('nBas',nBas,nSym)              ! get the number of basis functions

allocate( orbLabel(nBas)    )
allocate( vec(nBas,nBas)    )
allocate( rotvec(nBas,nBas) )
allocate( nOcc(nBas)        )

open(35,file = 'INPORB')

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

close(35)
end subroutine readvec

! ===============================================================================

subroutine get_basis_info
use basis_set_data
use orbital_data, only : nBas
use useful_data , only : nAtoms

integer                           :: j
character(len=LblL8), allocatable :: AtomLbl(:)


allocate( AtomLbl(nBas)  )
allocate( BasfnLbl(nBas) )
allocate( nCntr(nBas)   )

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
  end if
end do

do i = 1, nBas
  if ( nCntr(i) .eq. 0 ) nCntr(i) = nCntr(i-1)
end do

deallocate(AtomLbl)
end subroutine get_basis_info

! ===============================================================================

subroutine put_on_x(i,alpha,beta)
use useful_data
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
end if
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
end if
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
  end if
  do j = 1, nAtoms(i)
    x(i,j) = -x(i,j)
    z(i,j) = -z(i,j)
  end do
end if
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
  end if
end if
end subroutine put_on_x

! ===============================================================================

subroutine put_on_xy(i,delta)
  use useful_data
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
end if
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
  end if
end if
end subroutine put_on_xy

! ===============================================================================

subroutine capitalize(string)
implicit none
integer      :: i
character(*) string
do i = 1, len(string)
  if (ichar(string(i:i)).gt.96) then
    string(i:i) = char(ichar(string(i:i))-32)
  endif
end do
end subroutine capitalize

! ===============================================================================

subroutine locate(string)
implicit none
character(4)   ::  string,string2
character(132) ::  line
rewind(5)
40 read(5,*) line
string2=adjustl(line)
call capitalize(string2)
if (string2.ne.string) goto 40
end subroutine locate

! ===============================================================================

subroutine rotate_harmonics
use rotation_data
implicit none

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
end if

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
integer       :: l,m1,m2
real(kind=8)  :: x,P
x = P(0,l,m1,m2)
!   write(*,'(A,2I4,F15.8)')'U',m1,m2,x
end function capU

function capV(l,m1,m2) result(x)
implicit none
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
end if
!   write(*,'(A,2I4,F15.8)')'V',m1,m2,x
end function capV
       
   
function capW(l,m1,m2) result(x)
implicit none
integer       :: l,m1,m2
real(kind=8)  :: x,P
x=0
if ( m1 .gt. 0 ) then
  x = P(1,l,m1+1,m2) + P(-1,l,-m1-1,m2)
end if
if ( m1 .lt. 0 ) then
  x = P(1,l,m1-1,m2) - P(-1,l,-m1+1,m2)
end if
!   write(*,'(A,2I4,F15.8)')'W',m1,m2,x
end function capW


function P(i,l,mu,m) result(x)
use rotation_data
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
end if
if ( m .eq. l ) then
! x = rotmat(1,i,1) * rotmat(l-1,mu,m-1) - rotmat(1,i,-1) * rotmat(l-1,mu,-m+1)
  x = rotmat(1,ii,lmax+2) * rotmat(l-1,imu,im-1) - rotmat(1,ii,lmax) *  &
                                             rotmat(l-1,imu,im-2*m+1)
end if
if ( m .eq. -l ) then
! x = rotmat(1,i,1) * rotmat(l-1,mu,m+1) + rotmat(1,i,-1) * rotmat(l-1,mu,-m-1)
  x = rotmat(1,ii,lmax+2) * rotmat(l-1,imu,im+1) + rotmat(1,ii,lmax) *  &
                                             rotmat(l-1,imu,im-2*m-1)
end if
end function P


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
end if
!   write(*,'(2I4,F25.15)')m1,m2,x2
end function u



function v(l,m1,m2) result(x2)
implicit none
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
end if
!   write(*,'(2I4,F25.15)')m1,m2,x2
end function v


function w(l,m1,m2) result(x2)
implicit none
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
end if
!   write(*,'(2I4,F25.15)')m1,m2,x2
end function w

function delta(i,j) result(d)
implicit none
integer, intent(in)   :: i,j
integer               :: d
if ( i .eq. j ) then
  d = 1
else
  d = 0
end if
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
use useful_data
implicit none

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
use orbital_data
implicit none
integer   :: j,k

open(12,file = 'ROTORB')
write(12,'(A11)') '#INPORB 2.2'
write(12,'(A5)')  '#INFO'
write(12,'(A28)') '* Rotated molecular orbitals'
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
  end if
end do
602  format(I1,x,10A1)
deallocate(orbLabel)

end subroutine write_rotvec

! ===============================================================================

subroutine rotate_vec(a,b,c)
use orbital_data
use basis_set_data
use rotation_data
implicit none

integer        :: i,j
real(kind=8)   :: a,b,c

!write(*,'(A,3F12.6)') 'rotation angles',a,b,c
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
    end if
    if ( BasFnLbl(j)(9:12) .eq. 'd02-' ) then
      lmax = 2
      call rotate_harmonics
      call r_times_v(i,j)
      deallocate(rotmat)
    end if
    if ( BasFnLbl(j)(9:12) .eq. 'f03-' ) then
      lmax = 3
      call rotate_harmonics
      call r_times_v(i,j)
      deallocate(rotmat)
    end if
    if ( BasFnLbl(j)(9:12) .eq. 'g04-' ) then
      lmax = 4
      call rotate_harmonics
      call r_times_v(i,j)
      deallocate(rotmat)
    end if
  end do
end do

end subroutine rotate_vec

! ===============================================================================

subroutine r_times_v(i,j)
use orbital_data
use basis_set_data
use rotation_data
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
