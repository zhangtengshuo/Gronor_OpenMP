      subroutine svd(nm,m,n,a,na1,na2,w,nw,matu,u,nu1,nu2,matv,
     & v,nv1,nv2,ierr,rv1,nrv)

      implicit none

      integer ::  i,j,k,l,m,n,ii,i1,kk,k1,ll,l1,mn,nm,its,ierr
      integer :: na1,na2,nw,nu1,nu2,nv1,nv2,nrv
      real (kind=8) :: a(na1,na2),w(nw),u(nu1,nu2),v(nv1,nv2),rv1(nrv)
      real (kind=8) :: c,f,g,h,s,x,y,z,scale,anorm
      real (kind=8) :: sqrt,abs,sign
      logical :: matu,matv,ltest

!     this subroutine is a translation of the algol procedure svd,
!     num. math. 14, 403-420(1970) by golub and reinsch.
!     handbook for auto. comp., vol ii-linear algebra, 134-151(1971).

!     this subroutine determines the singular value decomposition
!          t
!     a=usv  of a real m by n rectangular matrix.  householder
!     bidiagonalization and a variant of the qr algorithm are used.

!     on input.

!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.  note that nm must be at least
!          as large as the maximum of m and n.

!        m is the number of rows of a (and u).

!        n is the number of columns of a (and u) and the order of v.

!        a contains the rectangular input matrix to be decomposed.

!        matu should be set to .true. if the u matrix in the
!          decomposition is desired, and to .false. otherwise.

!        matv should be set to .true. if the v matrix in the
!          decomposition is desired, and to .false. otherwise.

!     on output.

!        a is unaltered (unless overwritten by u or v).

!        w contains the n (non-negative) singular values of a (the
!          diagonal elements of s).  they are unordered.  if an
!          error exit is made, the singular values should be correct
!          for indices ierr+1,ierr+2,...,n.

!        u contains the matrix u (orthogonal column  m.o.'s) of the
!          decomposition if matu has been set to .true.  otherwise
!          u is used as a temporary array.  u may coincide with a.
!          if an error exit is made, the columns of u corresponding
!          to indices of correct singular values should be correct.

!        v contains the matrix v (orthogonal) of the decomposition if
!          matv has been set to .true.  otherwise v is not referenced.
!          v may also coincide with a if u is not needed.  if an error
!          exit is made, the columns of v corresponding to indices of
!          correct singular values should be correct.

!        ierr is set to
!          zero       for normal return,
!          k          if the k-th singular value has not been
!                     determined after 30 iterations.
!
!        rv1 is a temporary storage array.
!
!     this is a modified version of a routine from the eispack
!     collection by the nats project

!     modified to eliminate machep

      ierr=0

      do i=1, m
         do j=1, n
            u(i,j)=a(i,j)
         enddo
      enddo

!     householder reduction to bidiagonal form

      g=0.0d0
      scale=0.0d0
      anorm=0.0d0

      do i=1, n
         l=i+1
         rv1(i)=scale*g
         g=0.0d0
         s=0.0d0
         scale=0.0d0
         if(i .le. m) then

            do k=i, m
               scale=scale+abs(u(k,i))
            enddo

            if(scale.ne.0.0d0) then

               do k=i, m
                  u(k,i)=u(k,i)/scale
                  s=s+u(k,i)**2
               enddo

               f=u(i,i)
               g=-sign(sqrt(s),f)
               h=f*g-s
               u(i,i)=f-g
               
               if(i.ne.n) then
                  do j=l, n
                     s=0.0d0

                     do k=i, m
                        s=s+u(k,i)*u(k,j)
                     enddo

                     f=s/h

                     do k=i, m
                        u(k,j)=u(k,j)+f*u(k,i)
                     enddo
                  enddo
               endif
               
               do k=i, m
                  u(k,i)=scale*u(k,i)
               enddo

            endif
         endif
         
         w(i)=scale*g
         g=0.0d0
         s=0.0d0
         scale=0.0d0
         if(i .le. m .and. i.ne.n) then

            do k=l, n
               scale=scale+abs(u(i,k))
            enddo

            if(scale.ne.0.0d0) then

               do k=l, n
                  u(i,k)=u(i,k)/scale
                  s=s+u(i,k)**2
               enddo

               f=u(i,l)
               g=-sign(sqrt(s),f)
               h=f*g-s
               u(i,l)=f-g

               do k=l, n
                  rv1(k)=u(i,k)/h
               enddo

               if(i.ne.m) then

                  do j=l, m
                     s=0.0d0

                     do k=l, n
                        s=s+u(j,k)*u(i,k)
                     enddo

                     do k=l, n
                        u(j,k)=u(j,k)+s*rv1(k)
                     enddo
                  enddo

               endif

               do k=l, n
                  u(i,k)=scale*u(i,k)
               enddo

            endif
         endif
         
         anorm=max(anorm,abs(w(i))+abs(rv1(i)))

      enddo

!     accumulation of right-hand transformations ..........

      if(matv) then

!     for i=n step -1 until 1 do -- ..........

         do ii=1, n
            i=n+1-ii
            if(i.ne.n) then
               if(g.ne.0.0d0) then

                  do j=l, n

!     double division avoids possible underflow ..........

                     v(j,i)=(u(i,j)/u(i,l))/g
                  enddo

                  do j=l, n
                     s=0.0d0

                     do k=l, n
                        s=s+u(i,k)*v(k,j)
                     enddo

                     do k=l, n
                        v(k,j)=v(k,j)+s*v(k,i)
                     enddo
                  enddo
               endif

               do j=l, n
                  v(i,j)=0.0d0
                  v(j,i)=0.0d0
               enddo
            endif

            v(i,i)=1.0d0
            g=rv1(i)
            l=i

         enddo
      endif

!     accumulation of left-hand transformations

      if(matu) then

!     for i=min(m,n) step -1 until 1 do --

         mn=n
         if(m .lt. n) mn=m

         do ii=1, mn
            i=mn+1-ii
            l=i+1
            g=w(i)
            
            if(i.ne.n) then
               do j=l, n
                  u(i,j)=0.0d0
               enddo
            endif

            if(g.ne.0.0d0) then
               if(i.ne.mn) then

                  do j=l, n
                     s=0.0d0

                     do k=l, m
                        s=s+u(k,i)*u(k,j)
                     enddo

!     double division avoids possible underflow

                     f=(s/u(i,i))/g

                     do k=i, m
                        u(k,j)=u(k,j)+f*u(k,i)
                        
                     enddo
                  enddo
               endif
               
               do j=i, m
                  u(j,i)=u(j,i)/g
               enddo

            else

               
               do j=i, m
                  u(j,i)=0.0d0
               enddo

            endif

            u(i,i)=u(i,i)+1.0d0
         enddo

      endif

!     .......... diagonalization of the bidiagonal form ..........
!     .......... for k=n step -1 until 1 do -- ..........

      do kk=1,n
         
         k1=n-kk
         k=k1+1
         its=0
!     .......... test for splitting.
!     for l=k step -1 until 1 do -- ..........

 520     continue

         ltest=.false.
         do ll=1,k
            l1=k-ll
            l=l1+1
            if(abs(rv1(l))+anorm.eq.anorm) then
               ltest=.true.
               exit
            endif
!     .......... rv1(1) is always zero, so there is no exit
!     through the bottom of the loop ..........
            if(abs(w(l1))+anorm.eq.anorm) exit

         enddo

         if(.not.ltest) then
!     .......... cancellation of rv1(l) if l greater than 1 ..........
            c=0.0d0
            s=1.0d0

            do i=l,k
               f=s*rv1(i)
               rv1(i)=c*rv1(i)
               if(abs(f)+anorm.eq.anorm) exit
               g=w(i)
               h=sqrt(f*f+g*g)
               w(i)=h
               c=g/h
               s=-f/h
               if(matu) then

                  do j=1, m
                     y=u(j,l1)
                     z=u(j,i)
                     u(j,l1)=y*c+z*s
                     u(j,i)=-y*s+z*c
                  enddo
               endif
            enddo
         endif
!     .......... test for convergence ..........
         z=w(k)
         if(l.ne.k) then
!     .......... shift from bottom 2 by 2 minor ..........
            if(its.eq.30) then
!     .......... set error -- no convergence to a
!     singular value after 30 iterations ..........
               ierr=k
               return
            endif
            its=its+1
            x=w(l)
            y=w(k1)
            g=rv1(k1)
            h=rv1(k)
            f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
            g=sqrt(f*f+1.0d0)
            f=((x-z)*(x+z)+h*(y/(f+sign(g,f))-h))/x
!     .......... next qr transformation ..........
            c=1.0d0
            s=1.0d0

            do i1=l, k1
               i=i1+1
               g=rv1(i)
               y=w(i)
               h=s*g
               g=c*g
               z=sqrt(f*f+h*h)
               rv1(i1)=z
               c=f/z
               s=h/z
               f=x*c+g*s
               g=-x*s+g*c
               h=y*s
               y=y*c
               
               if(matv) then
                  do j=1, n
                     x=v(j,i1)
                     z=v(j,i)
                     v(j,i1)=x*c+z*s
                     v(j,i)=-x*s+z*c
                  enddo
               endif

               z=sqrt(f*f+h*h)
               w(i1)=z
!     .......... rotation can be arbitrary if z is zero ..........
               if(z.ne.0.0d0) then
                  c=f/z
                  s=h/z
               endif
               f=c*g+s*y
               x=-s*g+c*y

               if(matu) then

                  do j=1, m
                     y=u(j,i1)
                     z=u(j,i)
                     u(j,i1)=y*c+z*s
                     u(j,i)=-y*s+z*c
                  enddo
               endif

            enddo

            rv1(l)=0.0d0
            rv1(k)=f
            w(k)=x
            goto 520
!     .......... convergence ..........
         endif

         if(z.ne.0.0d0) then
!     .......... w(k) is made non-negative ..........
            w(k)=-z
            
            if(matv) then
               do j=1, n
                  v(j,k)=-v(j,k)
               enddo
            endif
         endif

      enddo

      return
      end
