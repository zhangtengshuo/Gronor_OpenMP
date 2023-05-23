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
!! Cororb
!!
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!

subroutine gronor_cororb()
  use gnome_parameters
  use gnome_data
  implicit none
  integer :: i,j,k

  do i=1,nelecs
    if(ev(i).lt.0.99d0) then
      do j=1,mbasel
        veca(j)=0.0d0
        vecb(j)=0.0d0
      enddo
      do j=1,nveca
        do k=1,mbasel
          veca(k)=veca(k)+u(j,i)*va(j,k)
          vecb(k)=vecb(k)+w(j,i)*vb(j,k)
        enddo
      enddo
    endif
  enddo

  return
end subroutine gronor_cororb
