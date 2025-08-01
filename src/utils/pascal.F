!> @brief   Calculate an element of Pascal's triangle
!!
!! @detail  used to determine the number of MS microstates generated by
!one particluar MS(max) determinant
!! @author  Coen de Graaf, URV
!! @date    November 2020


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

