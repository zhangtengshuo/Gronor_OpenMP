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
!>    Generate all M_s components (micro states) of a determinant with
!>    maximum M_s

subroutine gronor_generate_microstates(ndets,microdets)
use makebasedata
use cidef, only : lfndbg

implicit none

external :: gronor_quicksort_string,gronor_sminop,gronor_determine_nci

integer, intent(in)            :: ndets
integer, intent(out)           :: microdets
integer                        :: idet,jdet,ms,nci
integer                        :: first,last,new,all_new
real(kind=8)                   :: ci_seed,norm,invsqnorm
real(kind=8),allocatable       :: coef_tmp(:)
character(len=255)             :: occ_seed
character(len=255),allocatable :: occ_tmp(:)

real(kind=8),external           :: timer_wall_total

call gronor_determine_nci(ndets,nci)
allocate(micro_coef(nci))
allocate(micro_occ(nci))
allocate(micro_ndets(spinFrag))
allocate(coef_tmp(nci))
allocate(occ_tmp(nci))

do idet = 1, nci
  micro_coef(idet) = 0.0
  micro_occ(idet) = ''
end do
do ms = 1, spinFrag
  micro_ndets(ms) = 0
end do

micro_ndets(1) = ndets
do idet = 1, ndets
  micro_coef(idet) = coef(idet)
  micro_occ(idet) = occ(idet)
end do
first = 1
last = ndets
microdets = last
do ms = 2, spinFrag
  norm = 0.0d0
  all_new = 0
  do idet = 1, nci
    coef_tmp(idet) = 0.0
    occ_tmp(idet) = ''
  end do
  do idet = first,last
    ci_seed = micro_coef(idet)
    occ_seed = micro_occ(idet)
    call gronor_sminop(ci_seed,occ_seed,new)
    do jdet = 1, new
      coef_tmp(jdet + all_new) = coef_new(jdet)
      occ_tmp(jdet + all_new) = occ_new(jdet)
    end do
    all_new = all_new + new
    deallocate(occ_new,coef_new)
! occ_new and coef_new are allocated in gronor_sminop
  end do
  call gronor_quicksort_string(coef_tmp,occ_tmp,all_new)
  micro_ndets(ms) = 1
  micro_coef(microdets + micro_ndets(ms)) = coef_tmp(1)
  micro_occ(microdets + micro_ndets(ms)) = occ_tmp(1)
  do idet = 2, all_new
    if (occ_tmp(idet) .ne. occ_tmp(idet-1)) then
      micro_ndets(ms) = micro_ndets(ms) + 1
      micro_coef(microdets + micro_ndets(ms)) = coef_tmp(idet)
      micro_occ(microdets + micro_ndets(ms)) = occ_tmp(idet)
    else
      micro_coef(microdets + micro_ndets(ms)) = micro_coef(microdets +   &
                                 micro_ndets(ms)) + coef_tmp(idet)
    endif
  end do
  do idet = 1, micro_ndets(ms)
    norm = norm + micro_coef(idet + microdets)**2
  end do
  invsqnorm = 1.0d0/sqrt(norm)
  do idet = 1, micro_ndets(ms)
    micro_coef(idet + microdets) = micro_coef(idet + microdets)*invsqnorm
  end do
  write(lfndbg,'(a,i4,a,f8.4)')'norm of ms-component ',ms,' : ',norm
  first = last + 1
  last = first-1 + micro_ndets(ms)
  microdets = microdets + micro_ndets(ms)
end do


write(lfndbg,*)'micro states'
do idet = 1, microdets
  write(lfndbg,'(i4,f8.4,4x,a)') idet,micro_coef(idet),trim(micro_occ(idet))
end do

deallocate(coef_tmp,occ_tmp)

return
end subroutine gronor_generate_microstates
