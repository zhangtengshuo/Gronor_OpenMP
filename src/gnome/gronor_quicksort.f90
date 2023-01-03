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

      recursive subroutine quicksort(coef,occu,n)
      implicit none

      real (kind = 8) :: coef(n)
      integer         :: occu(n)
      real (kind = 8) :: random
      real (kind = 8) :: pivot,aux1
      integer         :: n,left,right,marker,aux2

      if (n .gt. 1) then
        call random_number(random)
        pivot = coef(int(random*real(n-1))+1)
        left = 0
        right = n + 1
        do while (left .lt. right)
          right = right - 1
          do while (abs(coef(right)) .lt. abs(pivot) )
            right = right - 1
          end do
          left = left + 1
          do while ( abs(coef(left)) .gt. abs(pivot) )
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
        call quicksort(coef(:marker-1),occu(:marker-1),marker-1)
        call quicksort(coef(marker:),occu(marker:),n-marker+1)
      end if
      end subroutine quicksort

! ===============================================================================
! sorting on the occupation string

recursive subroutine quicksort_string(coef,occu,n)
implicit none

real (kind = 8)    :: coef(n)
character(len=255) :: occu(n),pivot,aux2
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
    endif
  end do
  if (left .eq. right) then
    marker = left + 1
  else
    marker = left
  endif
  call quicksort_string(coef(:marker-1),occu(:marker-1),marker-1)
  call quicksort_string(coef(marker:),occu(marker:),n-marker+1)
endif
end subroutine quicksort_string


! ===============================================================================
! sorting on the absolute value of the coefficients

recursive subroutine quicksort_number(coef,occu,n)
implicit none

real (kind = 8)    :: coef(n),pivot
character(len=255) :: occu(n),aux2
real (kind = 8)    :: random
real (kind = 8)    :: aux1
integer            :: n,left,right,marker

if (n .gt. 1) then
  call random_number(random)
  pivot = abs(coef(int(random*real(n-1))+1))
  left = 0
  right = n + 1
  do while (left .lt. right)
    right = right - 1
    do while ( abs(coef(right)) .lt. pivot )
      right = right - 1
    end do
    left = left + 1
    do while ( abs(coef(left)) .gt. pivot )
      left = left + 1
    end do
    if ( left .lt. right) then
      aux1 = coef(left)
      coef(left) = coef(right)
      coef(right) = aux1
      aux2 = occu(left)
      occu(left) = occu(right)
      occu(right) = aux2
    endif
  end do
  if (left .eq. right) then
    marker = left + 1
  else
    marker = left
  endif
  call quicksort_number(coef(:marker-1),occu(:marker-1),marker-1)
  call quicksort_number(coef(marker:),occu(marker:),n-marker+1)
endif
end subroutine quicksort_number

