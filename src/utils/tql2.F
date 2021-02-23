      subroutine tql2(nm, n, d,nd, e,ne, z,nz1,nz2, ierr)

      implicit none
      integer :: nm,n,ierr
      integer :: nd,ne,nz1,nz2
      real (kind=8) :: d(nd), e(ne), z(nz1,nz2)

      integer :: i,ii,j,k,mml,m,l
      real (kind=8) :: f,b,g,p,s,h,r,c

      real (kind=8) :: machep

      machep=2.0d0**(-40)

!     ********** machep is a machine dependent parameter specifying
!                the relative precision of floating point arithmetic.

      ierr=0
      if(n.eq.1) goto 160

      do 10 i=2,n
         e(i-1)=e(i)
10    continue

      f=0.0d0
      b=0.0d0
      e(n)=0.0d0
!
      do 110 l=1,n
         j=0
         h=machep*(abs(d(l))+abs(e(l)))
         if(b.lt.h) b=h

!     ********** look for small sub-diagonal element **********

         do 20 m=l,n
            if(abs(e(m)).le.b) goto 30

!     ********** e(n) is always zero, so there is no exit
!                through the bottom of the loop **********

20       continue
30       if(m.eq.l) goto 100
40       if(j.eq.30) goto 150
         j=j+1

!     ********** form shift **********

         p=(d(l+1)-d(l))/(2.0d0*e(l))
         r=sqrt(p*p+1.0d0)
         h=d(l)-e(l)/(p+sign(r,p))

         do 50 i=l,n
            d(i)=d(i)-h
50       continue

         f=f+h

!     ********** ql transformation **********

         p=d(m)
         c=1.0d0
         s=0.0d0
         mml=m-l

!     ********** for i=m-1 step -1 until l do -- **********

         do 90 ii=1,mml
            i=m-ii
            g=c*e(i)
            h=c*p
            if(abs(p).lt.abs(e(i))) goto 60
            c=e(i)/p
            r=sqrt(c*c+1.0d0)
            e(i+1)=s*p*r
            s=c/r
            c=1.0/r
            goto 70
60          c=p/e(i)
            r=sqrt(c*c+1.0d0)
            e(i+1)=s*e(i)*r
            s=1.0/r
            c=c*s
70          p=c*d(i)-s*g
            d(i+1)=h+s*(c*g+s*d(i))

!     ********** form vector **********

            do 80 k=1,n
               h=z(k,i+1)
               z(k,i+1)=s*z(k,i)+c*h
               z(k,i)=c*z(k,i)-s*h
80          continue

90       continue

         e(l)=s*p
         d(l)=c*p
         if(abs(e(l)).gt.b) goto 40
100      d(l)=d(l)+f
110   continue

!     ********** order eigenvalues and eigenvectors **********

      do 140 ii=2,n
         i=ii-1
         k=i
         p=d(i)

         do 120 j=ii,n
            if(d(j).ge.p) goto 120
            k=j
            p=d(j)
120      continue

         if(k.eq.i) goto 140
         d(k)=d(i)
         d(i)=p

         do 130 j=1,n
            p=z(j,i)
            z(j,i)=z(j,k)
            z(j,k)=p
130      continue

140   continue

      goto 160

!     ********** set error -- no convergence to an
!                eigenvalue after 30 iterations **********

150   ierr=l
160   return
      end
