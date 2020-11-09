
!> @brief
!! Read records from integral file
!!     
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  
      subroutine gronor_read0(lfn,n)

      implicit none
      integer :: lfn,n
      integer :: i
      real (kind=8) :: b(510)

      do i=1,n,510
       read(lfn) b
      enddo

      return
      end

      subroutine gronor_read1(lfn,a,n)

      implicit none
      integer :: lfn,n
      real (kind=8) :: a(n)
      integer :: i,j,mxj
      real (kind=8) :: b(510)

      do i=1,n,510
       read(lfn) b
       mxj=min0(510,(n-i+1))
       do j=1,mxj
        a(i+j-1)=b(j)
       enddo
      enddo

      return
      end
