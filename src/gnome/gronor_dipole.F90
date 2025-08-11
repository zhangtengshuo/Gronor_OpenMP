!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ bases, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!> @brief
!! Dipole calculation
!!
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!

subroutine gronor_dipole(lfndbg)
  use cidist
  use gnome_parameters
  use gnome_data
#ifdef DEBUG_HDF5
  use debug_hdf5
#endif
  implicit none
  integer :: lfndbg

  integer :: i,j
  real (kind=8) :: x,y,z,aij,dx,dy,dz,z1,xc,yc,zc,dtot

600 format(' The dipole matrix elements will be calculated')

601 format(4x,'i',3x,'j',13x,'dx',20x,'dy',20x,'dz',//)

  x=0.0d0
  y=0.0d0
  z=0.0d0
  dx=0.0d0
  dy=0.0d0
  dz=0.0d0
  do i=1,nbas
    do j=1,nbas
      if(ising.eq.0) aij=ta(i,j)*deta*2
      if(ising .gt. 0) aij=diag(i)*sdiag(j)
      if(abs(aij).ge.1.0e-12) then
        dx=-1.0d0*dx*aij
        dy=-1.0d0*dy*aij
        dz=-1.0d0*dz*aij
        x=x+dx
        y=y+dy
        z=z+dz
110     format(1x,2i4,3(2x,e20.10))
      endif
    enddo
  enddo

  if(idipole.eq.1) then
    if(me.eq.mstr) then
      return
    endif
  else
#ifdef DEBUG_HDF5
    if(idbg.ge.13) then
      call dbg_log_msg('dipole','Dipole moment: electronic contribution')
    endif
#endif
  endif

  !     nuclear contribution to the dipole moment

#ifdef DEBUG_HDF5
  if(idbg.ge.13) then
    call dbg_write_scalar('dipole','x_elec',x)
    call dbg_write_scalar('dipole','y_elec',y)
    call dbg_write_scalar('dipole','z_elec',z)
  endif
#endif

  dx=0.0d0
  dy=0.0d0
  dz=0.0d0

  do i=1,nnucl
    z1=znuc(i)
    xc=xcord(i)
    yc=ycord(i)
    zc=zcord(i)
    dx=dx+z1*xc
    dy=dy+z1*yc
    dz=dz+z1*zc
  enddo

  !     nuclear contribution to the dipole moment
#ifdef DEBUG_HDF5
  if(idbg.ge.12) call dbg_log_msg('dipole','nuclear contribution')
#endif

  !     calculation of the total dipole moment

  dx=dx+x
  dy=dy+y
  dz=dz+z
  dtot=sqrt(dx*dx+dy*dy+dz*dz)

#ifdef DEBUG_HDF5
  if(idbg.ge.12) then
    call dbg_log_msg('dipole','Dipole moments')
    call dbg_write_scalar('dipole','dx',dx)
    call dbg_write_scalar('dipole','dy',dy)
    call dbg_write_scalar('dipole','dz',dz)
    call dbg_write_scalar('dipole','dtot',dtot)
  endif
#endif

  return
end subroutine gronor_dipole
