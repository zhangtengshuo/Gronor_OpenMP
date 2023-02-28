!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or
!     modified under
!     the Apache License version 2.0
!     (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license.
!     Unless required
!     by applicable law or agreed to in writing, software distributed
!     under the license
!     is distributed on an ‘as is’ basis, without warranties or
!     conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions
!     and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!>  @brief
!>    Get a rough estimate of the number of microstates for a fragment
!wave
!>    function, not taking into the duplicates to stay on the save side.

subroutine determine_nci(ndets,nci)
use makebasedata, only  : spinFrag,occ
implicit none

integer,intent(in)  :: ndets
integer,intent(out) :: nci
integer             :: idet,alpha,ms,pascal,iact,nOrb
character(len=255)  :: dumstr

nOrb = len_trim(occ(1))
nci = 0
do idet = 1, ndets
  alpha = 0
  dumstr=occ(idet)
  do iact = 1, nOrb
      if (dumstr(iact:iact) .eq. 'a') alpha = alpha + 1
  end do
  do ms = 1, spinFrag
    nci = nci + pascal(alpha+1,ms)
  end do
end do
return
end subroutine determine_nci

! ===============================================================================

integer function pascal(i,j)
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
pascal = table(i,j)
deallocate(table)
return
end function
