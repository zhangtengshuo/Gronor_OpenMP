#ifdef GCC4
#else
      integer function ishft(i,j)
!
! $Id: linux_shift.F 19707 2010-10-29 17:59:36Z d3y133 $
!
      implicit none
      integer i,j
      integer lshift
      real*8 rshift
      external lshift,rshift
!
!     routine to return i shifted by j bits
!
      if(j.eq.0) ishft=i
      if(j.gt.0) ishft=lshift(i,j)
      if(j.lt.0) ishft=rshift(i,-j)
      return
      end
#endif
