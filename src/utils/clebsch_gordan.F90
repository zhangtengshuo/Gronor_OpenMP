      real(kind=8) function clebsch_gordon(s1,s2,m1,m2,s,m)
!     Returns the requested Clebsch-Gordan coefficient
      implicit none

      external :: factorial,A0,recurrence_factor
      
      integer           :: factorial
      real(kind=8)      :: s1,s2,s
      real(kind=8)      :: m1,m2,m
      real(kind=8)      :: ispin,k,j
      real(kind=8)      :: recurrence_factor,A0,cg_p1,cg_p2,cg

      if((abs(m1).gt.s1).or.(abs(m2).gt.s2))then
         stop 'CLEBSCH-GORDAN:: Input data error'
      else if( (m1+m2) .ne. m )then
         cg = 0.0d0
      else if ((s.lt.abs(s1-s2)).or.(s.gt.(s1+s2)))then
         cg = 0.0d0
      else
         k = factorial(int(2*s1))*factorial(int(2*s2))*                         &
     &        factorial(int(s1+s2+m))*factorial(int(s1+s2-m))
         j = factorial(int(2*s1+2*s2)) * factorial(int(s1+m1))*                 &
     &        factorial(int(s1-m1))*factorial(int(s2+m2))*                      &
     &        factorial(int(s2-m2))
         cg = dsqrt(k/j)         ! starting point
         cg_p2 = cg
         ispin = s1 + s2         ! initial spin
         if(ispin.ne.s)then     ! s-1
            cg = cg*A0(ispin,m,s1,s2,m1,m2)/                                    &
     &           recurrence_factor(ispin,m,s1,s2)
            cg_p1 = cg
            ispin = ispin - 1
         endif
         do while (ispin .ne. s) ! start loop with s-2, need s-1 and s terms
            cg = (cg_p1*A0(ispin,m,s1,s2,m1,m2)-                                &
     &           cg_p2*recurrence_factor(ispin+1,m,s1,s2))/                     &
     &           recurrence_factor(ispin,m,s1,s2)
            ispin = ispin - 1
            cg_p2 = cg_p1
            cg_p1 = cg
         enddo
       endif
       clebsch_gordon=cg
      return
      end function clebsch_gordon

      real(kind=8) function recurrence_factor(s,m,s1,s2)
      implicit none

      real(kind=8)   :: s,m,s1,s2

      recurrence_factor = (s*s-m*m)*((s1+s2+1)**2-s*s)*(s*s-(s1-s2)**2)
      recurrence_factor = recurrence_factor/(s*s*(4*s*s-1))
      recurrence_factor = dsqrt(recurrence_factor)
      return
      end


      real(kind=8) function A0(s,m,s1,s2,m1,m2)
      implicit none

      real(kind=8)   :: s,m,s1,s2,m1,m2

      A0  = m1-m2+m*(s2*(s2+1)-s1*(s1+1))/(s*(s+1))
      return
      end

      integer function factorial(n)
      implicit none

      integer :: n, j
      if(n==0)then
         factorial = 1
      else
         factorial = 1
         do j = 1, n
            factorial = factorial*j
         enddo
      endif
      return
      end
