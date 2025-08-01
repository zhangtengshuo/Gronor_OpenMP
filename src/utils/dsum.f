      double precision function dsum(n,sx,incx)
!
! $Id: dsum.f 19707 2010-10-29 17:59:36Z d3y133 $
!
!
!     takes the sum of the array elements.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!
      double precision sx(*),stemp
      integer i,incx,m,mp1,n,nincx
!
      dsum = 0.0d0
      stemp = 0.0d0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        stemp = stemp + sx(i)
   10 continue
      dsum = stemp
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        stemp = stemp + sx(i) + sx(i + 1) + sx(i + 2)                           &
     &  + sx(i + 3) + sx(i + 4) + sx(i + 5)
   50 continue
   60 dsum = stemp
      return
      end
