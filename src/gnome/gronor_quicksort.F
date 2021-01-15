      recursive subroutine quicksort(coef,occu,n)
      implicit none

      real (kind = 8) :: coef(n)
      integer         :: occu(n)
      real (kind = 8) :: random
      real (kind = 8) :: pivot,aux1
      integer         :: n,i,left,right,marker,aux2

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
