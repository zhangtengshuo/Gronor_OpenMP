      integer function ipopcnt(ipack,n)
      use cidef
      implicit none
      integer :: ipack,n,iabs,i
      ipopcnt=0
      call unpack(ipack,iocch,n)
      do i=1,n
       ipopcnt=ipopcnt+iabs(iocch(i))
      enddo
      return
      end
