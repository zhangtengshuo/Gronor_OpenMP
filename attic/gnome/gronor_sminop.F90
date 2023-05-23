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
!>    Application of the S- operator on a determinant

subroutine gronor_sminop(coef_ms,occ_ms,new)
use makebasedata, only : occ_new,coef_new
implicit none

real(kind=8),intent(in)       :: coef_ms
character(len=255),intent(in) :: occ_ms
integer, intent(out)          :: new
integer                       :: iact,inew,nOrb
character(len=255)            :: dumstr

new = 0
nOrb = len_trim(occ_ms)
do iact = 1, nOrb
  if (occ_ms(iact:iact) .eq. 'a') new = new + 1
end do
if (new .ne. 0) then
  allocate(occ_new(new))
  allocate(coef_new(new))
  inew = 0
  do iact = 1, nOrb
    if (occ_ms(iact:iact) .eq. 'a') then
      inew = inew + 1
      dumstr = occ_ms
      dumstr(iact:iact) = 'b'
      occ_new(inew) = dumstr
      coef_new(inew) = coef_ms
    endif
  end do
endif

return
end subroutine gronor_sminop
