
      subroutine tred2(nm,n,a,na1,na2,d,nd,e,ne,z,nz1,nz2)

      implicit none

      integer :: nm,n
      integer :: na1,na2,nd,ne,nz1,nz2
      real (kind=8) ::  a(na1,na2), d(nd), e(ne), z(nz1,nz2)

      integer :: i,j,k,l,jp1,ii
      real (kind=8) :: hh,h,scale,f,g

      do i=1,n
       do j=1,i
        z(i,j)=a(i,j)
       enddo
      enddo

      if(n.ne.1) then

!     ********** for i=n step -1 until 2 do -- **********

       do ii=2,n
        i=n+2-ii
        l=i-1
        h=0.0d0
        scale=0.0d0

        if(l.lt.2) goto 40

!     ********** scale row (algol tol then not needed) **********

         do k=1,l
          scale=scale+abs(z(i,k))
         enddo

         if(scale.ne.0.0) goto 50
40       e(i)=z(i,l)
         goto 140

50       continue

         do k=1,l
          z(i,k)=z(i,k)/scale
          h=h+z(i,k)*z(i,k)
         enddo

         f=z(i,l)
         g=-sign(sqrt(h),f)
         e(i)=scale*g
         h=h-f*g
         z(i,l)=f-g
         f=0.0

         do 100 j=1,l
            z(j,i)=z(i,j)/(scale*h)
            g=0.0

!      ********** form element of a*u **********

            do 70 k=1,j
               g=g+z(j,k)*z(i,k)
70          continue

            jp1=j+1
            if(l.lt.jp1) goto 90

            do 80 k=jp1,l
               g=g+z(k,j)*z(i,k)
80          continue

!     ********** form element of p **********

90          e(j)=g/h
            f=f+e(j)*z(i,j)
100      continue

         hh=f/(h+h)

!     ********** form reduced a **********

         do 120 j=1,l
            f=z(i,j)
            g=e(j)-hh*f
            e(j)=g

            do 110 k=1,j
               z(j,k)=z(j,k)-f*e(k)-g*z(i,k)
110         continue
120      continue

         do 130 k=1,l
            z(i,k)=scale*z(i,k)
130      continue

140      d(i)=h

       enddo

      endif

160   d(1)=0.0d0
      e(1)=0.0d0

!     ********** accumulation of transformation matrices **********

      do 220 i=1,n
         l=i-1
         if(d(i).eq.0.0d0) goto 200

         do 190 j=1,l
            g=0.0d0

            do 170 k=1,l
               g=g+z(i,k)*z(k,j)
170         continue

            do 180 k=1,l
               z(k,j)=z(k,j)-g*z(k,i)
180         continue
190      continue

200      d(i)=z(i,i)
         z(i,i)=1.0d0
         if(l.lt.1) goto 220

         do 210 j=1,l
            z(i,j)=0.0d0
            z(j,i)=0.0d0
210      continue

220   continue

      return
      end
