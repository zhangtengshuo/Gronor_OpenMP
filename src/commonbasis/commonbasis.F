        program commonbasis
        implicit none

        call transvecmain()

        end
        
        subroutine transvecmain()
        implicit none

        integer :: i,j,k,l,m,n,nbas,l1,l2,kbas,ntot,ind3
        integer, parameter :: maxvec = 20, maxorb=4096
        integer :: nstates,lfnov,lfnvec

        character (len=255) :: item,root,ovname
        character (len=255),dimension(maxvec) :: vecname
        integer, dimension(maxvec) :: nocc,vecsec
        character (len=255) :: filov,filvec
        real(kind=8), allocatable :: s(:), vec(:,:)

        integer :: ipos,jtype,ilen,ipp,len,lensec,ispace
        real(kind=8) :: apos,deig,potnuc,thresh
!       common blocks
        common /blockinfo/ istart,apos(512)
        common /info/ ispace(2704)

        real(kind=8) :: q(512),pop(maxorb),potn,core
        integer :: istart,isec1,iseca,isecg,isecg2,ib,nrecl
        integer(kind=4) :: nrec

        logical :: inp_read,inp_compare,inp_a,inp_i,inp_strlen,inp_f

! thresh must be input
      data thresh/1.0d-6/

!       default values 
        nstates=0
        thresh=0.0d0

  1     continue

        if(inp_read()) then

  2      continue

         if(.not.inp_a(item)) goto 1

  3      continue

         if(inp_compare(.false.,'rootname',item)) then
          if(.not.inp_a(root)) root='rootname'
          goto 2
         endif

         if(inp_compare(.false.,'overlap',item)) then
          if(.not.inp_a(ovname)) ovname='ftn058'
         goto 2
         endif

!       read vectors from ftn file
         if(inp_compare(.false.,'ftnfile',item)) then
          nstates=nstates+1
!       If this is greater than the highest allocated (maxvec)
          if (nstates .gt. maxvec) then
           write(6,*) 'Too many vector files'
           stop
          end if
      if(.not.(inp_a(vecname(nstates)).and.inp_i(nocc(nstates))))then
           vecname(nstates)='ftn058'
           nocc(nstates)=0
          end if
         goto 2
         endif

!       read from ed3 file
         if(inp_compare(.false.,'ed3file',item)) then
          nstates=nstates+1
!       If this is greater than the highest allocated (maxvec)
          if (nstates .gt. maxvec) then
           write(6,*) 'Too many vector files'
           stop
          end if
          if(.not.(inp_a(vecname(nstates)).and.inp_i(vecsec(nstates))
     &    .and.inp_i(nocc(nstates))))then
           vecname(nstates)='ed3'
           vecsec(nstates)=1
           nocc(nstates)=0
          end if
         goto 2
         endif


         if(inp_compare(.false.,'thresh',item)) then
         if(.not.inp_f(thresh)) thresh=1.0d-6
         goto 2
         endif

         goto 1
        endif

        filov=ovname(1:inp_strlen(ovname))
        lfnov=10
!       section parameters, data from input
!       starting block
        istart=1
!       1-ints AO basis
        isec1=492
!       basis & geometry info
        iseca=501

        inquire(iolength=nrec)apos(1)
        nrecl=512*nrec
        open(unit=lfnov,file=filov,access='direct',recl=nrecl)
        ib=istart
        read(lfnov,rec=ib) apos

!       read geometry
        call upack3(apos(iseca),ipos,jtype,ilen)
        len=50
        ipp=ipos+istart+lensec(len)
        len=2704
        call reads(ispace,ipp,len,lfnov)
!       nba=number of cartesian basis functions
        nbas=ispace(2700)

        l1=nbas*(nbas+1)/2
!       allocate (one(l1,1))
        allocate(s(l1))
        call upack3(apos(isec1),ipos,jtype,ilen)
        ipp=ipos+istart
        call reads(q,ipp,7,lfnov)
        potnuc=q(1)
        write(6,'(A,F15.7)')'Nuclear repulsion energy: ',potnuc
        len=lensec(l1)
        ipp=ipos+istart+1+(i-1)*len
!       call reads(one(1,1),ipp,l1,lfnov)
        call reads(s(1),ipp,l1,lfnov)
        call prtri(s(1),nbas)
        close(lfnov)

        print*, 'nbas, nstates',nbas, nstates
!       read the vector sets
        kbas=0
        do i=1,nstates
         kbas=kbas+nocc(nstates)
        enddo

! Each set of vectors may have different number of occupied orbitals(AO
! basis),hence we need to total number of orbitals for the combination,
! (in MO basis). kbas is the total. ntot is used for filling vec(:,:).   

        allocate(vec(nbas,kbas))
        ntot=0

        if (vecname(1) .eq. 'ftn059') then
!       read ftn file vector info
        do i=1,nstates
         call readbasis(nbas,kbas,ntot,nocc(i),vec,vecname(i))   
         ntot=ntot+nocc(i)
        end do
        else
!       read ed3 file vector info 
        do i=1,nstates
         call read_write_vec(nbas,kbas,ntot,nocc(i),
     &   vec,vecname(i),vecsec(i))
         ntot=ntot+nocc(i)
        end do
        endif

        print *,'Vectors'
        call prev(vec,vec,kbas,nbas,nbas)

        call commonbasis(nocc,nstates,s,vec,thresh,nbas,kbas,l1,root)   
         
        end subroutine

        subroutine readbasis(nbas,kbas,ntot,nocc,vec,vecfile)
        implicit none

        integer :: i,j,nbas,l1,l2,kbas,ntot,nocc

        character (len=255) :: vecfile,filvinp
        integer :: lfnvinp
        dimension :: vec(nbas,kbas)            
        real(kind=8) :: vec
        logical :: inp_strlen

        filvinp=vecfile(1:inp_strlen(vecfile))
        lfnvinp=11
        open(unit=lfnvinp,file=filvinp,form='formatted')
        read(lfnvinp,*)
        do i=1,nocc
!        read(lfnvinp,'(6F20.15)')(vec(j,(n-1)*nocc+i),j=1,nbas)
         read(lfnvinp,'(6F20.15)')(vec(j,ntot+i),j=1,nbas)
        enddo
        close(lfnvinp)

        end subroutine

        subroutine commonbasis(nocc,nstates,s,vec,thresh,
     &  nbas,kbas,l1,root)
!       implicit real*8 (a-h,o-z)
        implicit none

        integer :: i,j,k,l,m,n,l1,l2,nstates,nbas,kbas,nmax,noff,ind3,
     &  newbas,ist,nocc(nstates) 
        real(kind=8) :: thresh,s(l1),vec(nbas,kbas)
        real(kind=8), allocatable :: scr(:),is(:),dummy(:)
        real(kind=8), allocatable :: ov(:),vv(:,:),ev(:),ov2(:)
        real(kind=8), allocatable :: vecbas(:,:),su(:,:),vecnew(:,:,:)
        character (len=255) :: vecfile,root,filvec 
        integer :: lfnvec 
        logical :: inp_strlen

        l2=kbas*(kbas+1)/2
        allocate(ov(l2),vv(kbas,kbas),ev(kbas),scr(2*kbas),is(kbas+1),
     &   ov2(l2))
!       calculate the overlap between vector sets
        do i=1,kbas
          do j=1,i
           ov(ind3(i,j))=0.0d0
           do k=1,nbas
            do l=1,nbas
             ov(ind3(i,j))=ov(ind3(i,j))+vec(k,i)*vec(l,j)*s(ind3(k,l))
            enddo
           enddo
           ov2(ind3(i,j))=ov(ind3(i,j))
          enddo
        enddo
        print *,'Overlap'
        call prtri(ov,kbas)
!       diagonalize the initial overlap matrix 
        call jacodiag(ov,kbas,vv,ev,scr,is)
        print *,'Eigenvalues'
        do i=1,kbas
           write(6,'(I5,E18.7)')i,ev(i)
        enddo
        call prev(vv,ev,kbas,kbas,kbas)
        newbas=0
        do i=1,kbas
          if (ev(i).gt.thresh) newbas=newbas+1
        enddo
        if (newbas.gt.nbas) stop 'onzin'
        allocate(vecbas(nbas,nbas))
        do i=1,nbas
          do j=1,nbas
            vecbas(i,j)=0.0d0
          enddo
        enddo
!       get new set of vectors from: old vector sets and 
!       eigenvalues and eigenvectors of the overlap matrix 
        do i=1,newbas
         do k=1,kbas
          do j=1,nbas
           vecbas(j,i)=vecbas(j,i)+vv(k,i)*vec(j,k)/sqrt(ev(i))
          enddo
         enddo
        enddo
        print *,'New vectors'
        open(20,file='common.basis',form='formatted')
        write(20,'(A11,2I5)')'read input ',nbas,nbas
        do i=1,nbas
           write(20,'(6F20.15)')(vecbas(j,i),j=1,nbas)
        enddo
        close(20)
!       nmax-largest number of occupied orbitals among the vector sets 
        nmax=0
!       print*, 'nstates:',nstates
        do i=1,nstates
          nmax=max(nmax,nocc(i))
        enddo
        allocate(vecnew(newbas,nmax,nstates))
        allocate(su(nbas,nmax),dummy(newbas))
        do i=1,newbas
           dummy(i)=0.0d0
        enddo
!       get the intermediate matrix SU of the transformation V'SU
        noff=0
        do ist=1,nstates
         do i=1,nbas
          do j=1,nocc(ist)
           su(i,j)=0.0d0
           do k=1,nbas
            su(i,j)=su(i,j)+s(ind3(i,k))*vec(k,noff+j)
           enddo
          enddo
         enddo
         noff=noff+nocc(ist)
!       get the final vectors after the transformation V'SU
         do i=1,newbas
          do j=1,nocc(ist)
           vecnew(i,j,ist)=0.0d0
           do k=1,nbas
            vecnew(i,j,ist)=vecnew(i,j,ist)+vecbas(k,i)*su(k,j)
           enddo
          enddo
         enddo
         print *,'Vectors of state ',ist,' in common MO basis'

         write(filvec,'(A,A,I3.3,A)') root(1:inp_strlen(root)),
     &   '_',ist,'.vec'
         lfnvec=12

         open(unit=lfnvec,file=filvec,action='readwrite')
         write(lfnvec,'(I)') newbas
         do i=1,nocc(ist)
          write(lfnvec,'(4F18.14)') (vecnew(j,i,ist),j=1,newbas)
         enddo
         do i=nocc(ist)+1,newbas
          write(lfnvec,'(4F18.14)')(dummy(j),j=1,newbas)
         enddo
         close(lfnvec)
        enddo
        end subroutine

      subroutine prtri(d,n)
c
c     ----- print out a triangular matrix -----
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension d(*),dd(12)
      iwr=6
      max = 12
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. n) imax = n
      write (iwr,9008)
      write (iwr,8028) (i,i = imin,imax)
      do 160 j = 1,n
      k = 0
      do 140 i = imin,imax
      k = k+1
      m = max0(i,j)*(max0(i,j)-1)/2 + min0(i,j)
  140 dd(k) = d(m)
      write (iwr,8048) j,(dd(i),i = 1,k)
  160 continue
      if (imax .lt. n) go to 100
      return
 9008 format(/)
 9028 format(6x,7(6x,i3,6x))
 9048 format(i5,1x,7f15.10)
 8028 format(6x,12(3x,i3,3x))
 8048 format(i5,1x,12f9.4)
      end

      subroutine prev(v,e,m,n,ndim)
c
c     ----- print out e and v-matrices
c
      implicit real*8  (a-h,p-w),integer   (i-n),logical    (o)
      implicit character *8 (z),character *1 (x)
      implicit character *4 (y)
      dimension v(ndim,*),e(*)
      iwr=6
      max = 10
      imax = 0
  100 imin = imax+1
      imax = imax+max
      if (imax .gt. m) imax = m
      write (iwr,9008)
      write (iwr,8068) (e(i),i = imin,imax)
      write (iwr,9008)
      write (iwr,8028) (i,i = imin,imax)
      write (iwr,9008)
      do 150 j = 1,n
  150 write (iwr,8048) j,(v(j,i),i = imin,imax)
 140  if (imax .lt. m) go to 100
      return
 9008 format(/)
 8028 format(17x,10(3x,i3,3x))
 8048 format(i5,12x,10f9.4)
 8068 format(17x,10f9.4)
 9028 format(17x,12(6x,i3,6x))
 9048 format(i5,2x,10x,7f15.10)
 9068 format(17x,7f15.10)
      end

        function ind3(y,z)
!       index on 2-e list when e-labels are given
        integer :: y,z,ind3
        ind3=max0(y,z)*(max0(y,z)-1)/2 + min0(y,z)
        end function

!
! subroutines for jacobi diagonalization
!

      subroutine jacodiag(a,n,evecs,evals,si,is)
c
      implicit real*8   (a-h,o-z) , integer   (i-n)
c....
c.... diagonalises triangular matrix with dimension n eigenvectors
c.... are returned in evecs, eigenvalues in evals, 
c.... si is een scratch array 2*n
c
      dimension a(*),evecs(*),evals(*),si(*),is(*)
c
      small = 1.0d-12
      do i=1,n+1
         is(i)=i*(i-1)/2
      enddo
c
      call jacobi(a,n,is,evecs,n,evals,0,3,small,si)
      return
      end

      subroutine jacobi(h,nbasis,iky,v,nrow,e,init,iorder,small,y)
c
      implicit real*8   (a-h,o-z) , integer   (i-n)
c
c.....
c.....standard jacobi diagonaliser. h is diagonalised, v will contain
c.....the eigenvectors, e the eigenvalues. iky should contain i*(i-1)/2
c.....at i. y is a scratch array of length 2*nbasis.
c.....init   : = 1 ? => v is assumed to be initialised elsewhere
c.....  ,,   : other values will cause it to be a unit matrix at first
c.....iorder : = 0 ? => no ordering in the vectors, no e.val. at all
c.....  ,,   : = 1 ? => no ordering in the vectors, e.val's in e.
c.....  ,,   : = 2 ? => increasing eigenvalues/vectors
c.....  ,,   : = 3 ? => decreasing eigenvalues/vectors
c.....  ,,   : = 4 ? => lock mode (i.e. try to change the initial
c.....                             vectors the least possible)
c.....
      dimension h(*),v(nrow,nbasis),e(*),y(*),iky(*)     
      if (init.ne.1) then
         call zero(v,nrow*nbasis)
         do 10 i=1,nbasis
            v(i,i) = 1.0d0
10       continue
      end if
      if (nbasis.eq.1) then
         e(1)   = h(1)
         return
      end if
11    rlarge = 0.0d0
      do 30 i=2,nbasis
         ii = iky(i)
         do 20 j=1,i-1
           if (dabs(h(ii+j)).gt.rlarge) rlarge = dabs(h(ii+j))
20       continue
30    continue
      if (rlarge.gt.small) then
         reason = rlarge * .1d0
         do 60 i=1,nbasis - 1
            call fmove( h(iky(i)+1),y,i)
            call gather(nbasis-i,y(i+1),h(i+1),iky(i+1))
            hii = y(i)
            do 50 j=i+1,nbasis
               if (dabs(y(j)).ge.reason) then
                  vii  = hii * 0.5d0
                  vij  = y(j)
                  vjj  = h(iky(j)+j) * 0.5d0
                  diff = vii - vjj
                  root = dsqrt( diff*diff + vij*vij )
                  if (diff.lt.0.0d0) root = -root
                  cosr = (diff + root)/vij
                  sinr = dsqrt( 1.0 d0/ (1.0d0+cosr*cosr) )
                  cosr = cosr * sinr
                  diff = vii + vjj
                  call drotj(j-1,y,1,h(iky(j)+1),1,cosr,sinr)
                  hii  = diff + root
                  h( iky(j)+j ) = diff - root
                  y(j) = 0.0d0
                  do 40 k=j+1,nbasis
                     vij   = y(k)
                     kj    = iky(k) + j
                     y(k)  = vij   * cosr + h(kj) * sinr
                     h(kj) = h(kj) * cosr - vij   * sinr
40                continue
                  call drotj(nbasis,v(1,i),1,v(1,j),1,cosr,sinr)
               end if
50          continue
            call fmove(y,h(iky(i)+1),i-1)
            call scatter(nbasis-i,h(i+1),iky(i+1),y(i+1))
            h(iky(i)+i) = hii
60       continue
         goto 11
      end if
      if (iorder.eq.0) then
         return
      else if (iorder.eq.1) then
         call gather(nbasis,e,h,iky(2))
         return
      else if (iorder.eq.2) then
         call gather(nbasis,e,h,iky(2))
         call order(e,v,nbasis,nrow,y,y(nbasis+1), 1)
         call scatter(nbasis,h,iky(2),e)
      else if (iorder.eq.3) then
         call gather(nbasis,e,h,iky(2))
         call order(e,v,nbasis,nrow,y,y(nbasis+1),-1)
         call scatter(nbasis,h,iky(2),e)
      end if
      return
      end
      subroutine gather(n,r,a,map)
c
      implicit real*8 (a-h,o-z), integer (i-n)
      dimension r(n),a(*),map(n)
c
      do 10 loop=1,n
   10 r(loop) = a(map(loop))
c
      return
      end

      subroutine igather(n,r,a,map)
      implicit integer (a-z)
      dimension a(*),r(n),map(n)
c
      do 10 loop=1,n
   10 r(loop) = a(map(loop))
c
      return
      end

      subroutine scatter(n,r,index,a)
c
c...  cray scilib imitation/ but not as in my manual (jvl 1986)
c...  arguments as in cyber205 atmol-library
c
      implicit real*8 (a-h,o-z), integer (i-n)
      dimension r(*),index(n),a(n)
c
      do 10 i=1,n
10    r(index(i)) = a(i)
c
      return
      end
      subroutine order(e,v,n,nrow,y,imap,ipar)
c
      implicit real*8   (a-h,o-z) , integer   (i-n)
c
      dimension e(*),v(nrow,n),y(*),imap(*)
      do 10 i=1,n
         imap(i) = i
10    continue
      call bubble(e,n,ipar,imap)
      do 20 iold=1,n-1
         inew = imap(iold)
         imap(iold) = 0
         if (iold.eq.inew) goto 20
         call fmove(v(1,iold),y,n)
         call fmove(v(1,inew),v(1,iold),n)
         call fmove(y,v(1,inew),n)
         imap(locati(iold,imap,n)) = inew
20    continue
      return
      end
      subroutine bubble(a,n,ipar,imap)
c
      implicit real*8   (a-h,o-z) , integer   (i-n)
c
c.....
c.....simple bubble sort. imap keeps track of the permutations
c.....ipar enables one to use this routine for sorting increasingly
c.....(ipar > 0 => first element smallest) or the other way around
c.....
      dimension a(n),imap(n)
      logical ready
10    ready = .true.
      do 20 i=2,n
         if (a(i-1)*ipar.gt.a(i)*ipar) then
            r         = a(i-1)
            a(i-1)    = a(i  )
            a(i  )    = r
            ready     = .false.
            ii        = imap(i-1)
            imap(i-1) = imap(i  )
            imap(i  ) = ii
         end if
20    continue
      if (.not.ready) goto 10
      return
      end
      subroutine zero(v,n)
c
c...  zero vector
c
      implicit double precision (a-h,o-z), integer(i-n)
      dimension v(*)
c
      do 10 i=1,n
         v(i)=0.0d0
10    continue
c
      return
      end
      subroutine fmove(a,b,n)
c
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(n),b(n)
c
c...  b = a
c
      do 10 i=1,n
         b(i)=a(i)
10    continue
c
      return
      end
      function locati(i,ii,n)
      implicit real*8(a-h,o-z)
      dimension ii(*)
      do 10 j=1,n
         if (ii(j).eq.i) then
            locati = j
            return
         end if
10    continue
      locati = 0
      return
      end

      subroutine drotj(n,sx,incx,sy,incy,sc,ss)
      implicit real*8  (a-h,o-z)
c
c                b l a s  subprogram
c    description of parameters
c
c     --input--
c        n  number of elements in input vector(s)
c       sx  double precision vector with n elements
c     incx  storage spacing between elements of sx
c       sy  single precision vector with n elements
c     incy  storage spacing between elements of sy
c       sc  element of rotation matrix
c       ss  element of rotation matrix
c
c     --output--
c       sx  rotated vector sx (unchanged if n .le. 0)
c       sy  rotated vector sy (unchanged if n .le. 0)
c
c     multiply the 2 x 2 matrix  ( sc ss) times the 2 x n matrix (sx**t)
c                                (-ss sc)                        (sy**t)
c     where **t indicates transpose.  the elements of sx are in
c     sx(lx+i*incx), i = 0 to n-1, where lx = 1 if incx .ge. 0, else
c     lx = (-incx)*n, and similarly for sy using ly and incy.
      dimension sx(*),sy(*)
      data zero,one/0.d0,1.d0/
      if(n .le. 0 .or. (ss .eq. zero .and. sc .eq. one)) go to 40
      if(.not. (incx .eq. incy .and. incx .gt. 0)) go to 20
c
           nsteps=incx*n
           do 10 i=1,nsteps,incx
                w=sx(i)
                z=sy(i)
                sx(i)=sc*w+ss*z
                sy(i)=-ss*w+sc*z
   10           continue
           go to 40
c
   20 continue
           kx=1
           ky=1
c
           if(incx .lt. 0) kx=1-(n-1)*incx
           if(incy .lt. 0) ky=1-(n-1)*incy
c
           do 30 i=1,n
                w=sx(kx)
                z=sy(ky)
                sx(kx)=sc*w+ss*z
                sy(ky)=-ss*w+sc*z
                kx=kx+incx
                ky=ky+incy
   30           continue
   40 continue
c
      return
      end
        subroutine reads(a,iblk,len,inpfile)
!       read the record, gamess-uk has a record length of 512 words
        implicit none

        integer iread,ib,iblk,len,i,jj,nword
        real*8 a(len), q(512)
        integer :: inpfile

        iread=0
        ib=iblk

   10   read(inpfile,rec=ib) q
        call upack2(q(512),jj,nword)
        do i=1,nword
         iread=iread+1
         a(iread)=q(i)
        enddo
        if (iread.lt.len) then
         ib=ib+1
         goto 10
        endif
        return
        end

        subroutine upack2(i4,left,iright)
        implicit none

        integer left,iright
        integer*4 i4(2)
        left=i4(2)
        iright=i4(1)
        return
        end

        subroutine upack3(packed,i,j,k)
        implicit none
        integer i,j,k
        integer packed

        i=ishft(packed,-40)
        j=iand(ishft(packed,-20),z'FFFFF')
        k=iand(packed,z'FFFFF')
        return
        end

        function lensec(nwords)
        implicit none

        integer nwords,lensec
        lensec=(nwords-1)/511+1
        return
        end

        subroutine unpack(ip,nbits,iu,nw)
        implicit none

c     Unpack performs the inverse of the operation performed by pack.
c     Note that the number of integers to unpack does not have to be a
c     multiple of the integer size divided by nbits. See the subroutine
c     pack for further details.

c     Input:

        integer nw     ! the number of integer to unpack
        integer nbits  ! the number of bits to unpack for every integer
        integer ip(*)  ! the closely packed integers

c     Output:

        integer iu(nw) ! the unpacked integers

c        ((nw/nbits)*nbits.ne.nw)  does happen

        if (nbits.eq.8) then
        call viup8(ip,iu,nw)
        else if (nbits.eq.16) then
        call viup16(ip,iu,nw)
        else
         write(6,*) '*** unpack called with nbits ne 8,16,32 ***'
        end if

        return
        end

        subroutine viup8(ip,iu,nw)
        implicit none

c     Unpack the 8 rightmost bits of every integer for iu from the
c     closely packed array ip. See the subroutine pack for the
c     definitions of the bit operations. The implementation of this
c     subroutine is derived of the implementation of vipk8.

c     Parameters:

        integer mx               ! the size of integer divide by 8
!_IF(i8)
        parameter( mx = 8 )
!_ELSE
!      parameter( mx = 4 )
!_ENDIF

c     Input:

        integer nw               ! the number of integer to unpack into iu
        integer ip((nw+mx-1)/mx) ! the packed integers

c     Output:

        integer iu(nw)           ! the unpacked integers

c     Local:

        integer ibits8           ! all rightmost 8 bits set
        parameter(ibits8 = z'ff')
        integer indu             ! the number of unpacked integers
        integer indp             ! counter over elements in ip
        integer i                ! counter
        integer j                ! limit for the remaining integers
        integer ishifted         ! the bits of interest shifted to the
                               ! right position
        integer ibits            ! the bits of interest only
        integer nshift(mx)       ! the number of bits to shift
                               ! required for explicit unpacking
                               ! see gotcha

!IF(i8)
!_IF(littleendian)
          data nshift/56,0,8,16,24,32,40,48/
!_ELSE
!      data nshift/0,56,48,40,32,24,16,8/
!_ENDIF
!_ELSE
!_IF(littleendian)
!      data nshift/24,0,8,16/
!_ELSE
!      data nshift/0,24,16,8/
!_ENDIF
!_ENDIF

        indu = 0
        do indp = 1, nw/mx
          do i = 1, mx
            indu = indu + 1
            ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
            ibits    = iand( ishifted, ibits8 )

            iu(indu) = ibits
          enddo
        enddo

c     if nw is not a multiple of mx then unpack the remaining integers

        if (mod(nw,mx) .gt. 0) then
          indp = nw/mx+1
          do i = 1, mod(nw,mx)
            indu = indu + 1
            ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
            ibits    = iand( ishifted, ibits8 )
            iu(indu) = ibits
          enddo
        endif
        end subroutine

        subroutine viup16(ip,iu,nw)
        implicit none

        integer mx
        parameter( mx = 4 )
        integer nw
        integer ip((nw+mx-1)/mx)
        integer iu(nw)
        integer ibits16
        parameter(ibits16 = z'ffff')
        integer indu
        integer indp
        integer i
        integer j
        integer ishifted

        integer ibits
        integer nshift(mx)


        data nshift/48,0,16,32/

        indu = 0
        do indp = 1, nw/mx
          do i = 1, mx
            indu = indu + 1
            ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
            ibits    = iand( ishifted, ibits16 )

            iu(indu) = ibits
          enddo
        enddo

c     if nw is not a multiple of mx then unpack the remaining integers
        if (mod(nw,mx) .gt. 0) then
          indp = nw/mx+1
          do i = 1, mod(nw,mx)
            indu = indu + 1
            ishifted = ishft( ip(indp), -nshift(mod(indu,mx)+1) )
            ibits    = iand( ishifted, ibits16 )
            iu(indu) = ibits
          enddo
        endif
        end subroutine

        subroutine tdown(qnew,ilifn,q,ilifq,nnn,num,
     &  ilifc,ntran,itran,ctran,otran)
!       transforms vectors to AO basis
        implicit none

        integer,parameter :: maxorb=4096,mxorb1=maxorb+1,mxorb3=maxorb*3
        integer nnn,num,i,j,k,l,m,n
        integer ilifn(num),ilifq(num),ilifc(maxorb),ntran(maxorb),
     &  itran(mxorb3)
        real(kind=8)::ctran(mxorb3),qnew(num*num),q(num*num),dum(maxorb)
        logical otran

        do i=1,num
         ilifq(i) = (i-1)*num
         ilifn(i) = (i-1)*num
        end do

        if(otran)go to 60008
        do 731 i=1,nnn
        m=ilifq(i)
        call vclr(dum,1,num)
        do 733 j=1,num
        n=ntran(j)
        do 733 k=1,n
        l=ilifc(j)+k
  733   dum(itran(l))=ctran(l)*q(m+j)+dum(itran(l))
        call dcopy(num,dum(1),1,qnew(ilifn(i)+1),1)
  731   continue

60010   if (nnn.lt.num) then
!       clear vectors that are extra
         do i=nnn+1,num
!       fills an array with zeros
          call vclr(qnew(ilifn(i)+1),1,num)
         end do
        end if

        return
60008   do 60009 i=1,nnn

        call dcopy(num,q(ilifq(i)+1),1,qnew(ilifn(i)+1),1)
60009   continue
        return
        end subroutine

        subroutine vclr(a,incr,n)
        implicit none

        integer incr,n,loop,loopi
        real*8 a(*)

        if(incr.eq.1) then
         do 10 loop=1,n
  10     a(loop) = 0.0d0
        else
         loopi=1
         do 20 loop=1,n
         a(loopi)=0.0d0
         loopi=loopi+incr
  20    continue
        endif
        return
        end subroutine

        
        subroutine read_write_vec(nbas,kbas,ntot,nocc,
     &  vec_tran,vecfile,isecv)
        implicit none

!       gamess parameters
        integer,parameter :: maxorb=4096,mxorb1=maxorb+1,mxorb3=maxorb*3
        integer,parameter :: maxat=750,mxshel=2048,mxprim=8192

!       common blocks
        common/blockinfo/istart,apos(512)
        common /runlab/ zcom(19),ztitle(10),zaname(maxat),
     &  zbflab(maxorb),ztag(maxat),zsymm,zgroup,zscftp,zruntp,zguess,
     &  zconf,zstate,zorb(maxorb),zpseud(maxat)
!       zaname : names of atoms during calculation (after reorder)
!       ztag   : names of atoms as read in
        common /tran/ ilifc(maxorb),ntran(maxorb),itran(mxorb3),
     &  ctran(mxorb3),otran,otri
        common /blkqio/ deig(maxorb),dpop(mxorb1),nba,new,ncol,jeig,jpop
        common /info/ ispace(2704)

        integer :: nvec,nints,iv,ib,istart,isecv,nocc,ntot,nbas,kbas
        real(kind=8) :: apos
        logical :: inp_strlen
        character (len=255) :: item,root,vecfile
        integer, parameter :: maxvec = 10
        character (len=255),dimension(maxvec) :: vector
        integer, dimension(maxvec) :: isecv_arr
        character (len=255) :: filvinp,filvec
        integer :: lfnvinp,lfnvec
        integer :: m29,ipos,jtype,ilen,ipp,len,lensec
        real(kind=8), allocatable :: vec(:,:)
        integer :: ilifc,ntran,itran,ctran,ispace
        logical otran,otri

        integer :: ns(8),ideg(3),iddeg(6),inv(8),imp(8),ndri(8),
     &  nt(3),ilifq(maxorb),ilifn(maxorb)
        real(kind=8) :: q(512),buf(510),qnew(maxorb*maxorb),
     &  vec_tran(nbas,kbas)

        character *8 zcom,ztitle,zaname,ztag,zsymm,zgroup
        character *10 zbflab
        character *8 zscftp,zruntp,zguess,zconf,zstate,zorb,zpseud

        integer :: nat,ich,mul,num,nx,ne,na,nb,ix
        real(kind=8) :: z,ztitl
        real(kind=8) :: dpop
        integer :: new,ncol,jeig,jpop
        integer :: non,norb,nspace
        integer :: nrecl
        integer(kind=4) :: nrec

        real(kind=8) :: deig
        integer :: nba,l1,l3,mdum,nstep,lz,ioff
        integer :: i,j,k,l,m,n
        data m29/29/

        filvinp=vecfile(1:inp_strlen(vecfile))
        lfnvinp=12

        istart=1
        inquire(iolength=nrec)apos(1)
        nrecl=512*nrec
        open(unit=lfnvinp,file=filvinp,access='direct',recl=nrecl)
        ib=istart
        read(lfnvinp,rec=ib) apos

!       read title
        call upack3(apos(isecv),ipos,jtype,ilen)
        ipp=ipos+istart
        call reads(zcom,ipp,m29+2*maxat + maxorb,lfnvinp)
        len=2*maxorb+7
        ipp=ipos+1+istart
        call reads(deig,ipp,len,lfnvinp)
        len=2*maxorb+7
        ipp=ipos+istart+lensec(len)+1

!       read in ilifc,ntran,itran,ctran,otran,otri
        len=2*maxorb+2*mxorb3
        call reads(ilifc,ipp,len,lfnvinp)

!       read vectors
        allocate(vec(nba,nba))
        len=mxorb3+(2*maxorb+mxorb3+1)
        ipp=ipp+lensec(len)
        call reads(vec,ipp,nba*nba,lfnvinp)

!       transform vectors from symmetry adapted basis back to AO basis
        call tdown(qnew,ilifn,vec,ilifq,nba,nba,
     &  ilifc,ntran,itran,ctran,otran)
!       qnew contains transformed vectors
!       open(unit=lfnvec,file=filvec,action='readwrite')
!       write(lfnvec,'(I)') nba

!       do j=1,nba
!        write(lfnvec,'(4F18.14)') (qnew((j-1)*nba+i), i=1,nba)
!       end do
!       close(lfnvec)

!       do i=1,nocc
!!       read(lfnvinp,'(6F20.15)')(vec(j,(n-1)*nocc+i),j=1,nbas)
!        read(lfnvinp,'(6F20.15)')(vec(j,ntot+i),j=1,nbas)
!       enddo


        do i=1,nocc
         do j=1,nbas
          vec_tran(j,ntot+i)=qnew(j+(i-1)*nbas)
          print*, vec_tran(j,ntot+i)
         enddo
        enddo

        close(lfnvinp)

        end subroutine

