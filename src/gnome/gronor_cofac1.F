      subroutine gronor_cofac1(lfndbg)
      use cidist
      use gnome_parameters
      use gnome_data

#ifdef CUSOLVER
      use cusolverDn
c      use iso_c_binding
      use cuda_cusolver
c      use cuda_functions
#endif
      implicit none
      integer :: lfndbg
      integer :: i,ierr,j,idetuw,k
      real (kind=8) :: coef
      real (kind=8) :: cmax, cnorm, coefu
      integer :: nz1, nz2
      character (len=1), target :: jobu, jobvt
      real (kind=8), allocatable :: temp(:,:), ttt(:,:)
      
#ifdef CUSOLVER
c      integer (c_int) :: ndim
c      integer (c_int) :: ldim
#endif


      if(idbg.ge.13) write(lfndbg,600)
 600  format(/,' Cofactor matrix will be calculated')

      call timer_start(41)

#ifdef CUSOLVER

      if(itest.eq.1.or.itest.eq.3) then

       jobu = 'A'  ! all m columns of U
       jobvt= 'A'  ! all m columns of VT

       ndim=nelecs
       mdim=mbasel
c       print*,me,' NDIM,MDIM',ndim,mdim

       allocate(temp(ndim,ndim))

       do j=1,ndim
        do i=1,ndim
         temp(i,j)=ta(i,j)
        enddo
       enddo

c       write(*,'(a)') 'ta'
c       do j=1,nelecs
c        write(*,'(12f12.6)') (ta(i,j),i=1,nelecs)
c       enddo

!$acc data copy(temp,ev,u,w,dev_info_d) create(workspace_d,rwork)
!$acc host_data use_device(temp,ev,u,w,dev_info_d,workspace_d,rwork)
       cusolver_status=cusolverDnDgesvd(cusolver_handle,jobu,jobvt,
     &  ndim,ndim,temp,ndim,ev,u,ndim,w,ndim,workspace_d,
     &      worksize,rwork,dev_info_d)
!$acc end host_data
!$acc wait

       if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &   write(*,*) 'cusolverDnDgesvd failed',cusolver_status

!$acc update self(temp,u,w,ev,dev_info_d)
!$acc end data

c       do j=1,ndim
c        do i=1,ndim
c         ta(i,j)=temp(i,j)
c        enddo
c       enddo

       do i=1,nelecs
        do j=1,nelecs
         temp(i,j)=w(j,i)
        enddo
       enddo
       do i=1,nelecs
        do j=1,nelecs
         w(i,j)=temp(i,j)
        enddo
      enddo
      
c       write(*,'(a)') 'ev_c'
c       write(*,'(12f12.6)') (ev(i),i=1,nelecs)
c       write(*,'(a)') 'u_c'
c       do j=1,nelecs
c        write(*,'(12f12.6)') (u(i,j),i=1,nelecs)
c       enddo
c       write(*,'(a)') 'w_c'
c       do j=1,nelecs
c        write(*,'(12f12.6)') (w(i,j),i=1,nelecs)
c       enddo

c       do j=1,ndim
c        do i=1,mdim
c         ta(i,j)=temp(i,j)
c        enddo
c     enddo
      
       deallocate(temp)
      else
      call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     & ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     & ierr,sdiag,nelecs)
      endif

c      write(lfndbg,'(a)') 'ta'
c      do j=1,4
c       write(lfndbg,'(12f10.5)') (ta(i,j),i=1,nelecs)
c      enddo
c      write(lfndbg,'(a)') 'u'
c      do j=1,4
c       write(lfndbg,'(12f10.5)') (u(i,j),i=1,nelecs)
c      enddo
c      write(lfndbg,'(a)') 'w'
c      do j=1,4
c       write(lfndbg,'(12f10.5)') (w(i,j),i=1,nelecs)
c      enddo

c       write(*,'(a)') 'ev_q'
c       write(*,'(12f12.6)') (ev(i),i=1,nelecs)
c       write(*,'(a)') 'u_q'
c       do j=1,nelecs
c        write(*,'(12f12.6)') (u(i,j),i=1,nelecs)
c       enddo
c       write(*,'(a)') 'w_q'
c       do j=1,nelecs
c        write(*,'(12f12.6)') (w(i,j),i=1,nelecs)
c       enddo
c       write(*,'(a)') 'sdiag_q'
c       write(*,'(12f12.6)') (sdiag(i),i=1,nelecs)

c       stop
#else
      call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     & ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     & ierr,sdiag,nelecs)

#endif


      
#ifdef CUSOLVEROLD

      if(itest.gt.0) then

      jobu='A'
      jobvt='A'

      do j=1,nelecs
       do i=1,nelecs
        acu((j-1)*nelecs+i)=ta(i,j)
       enddo
      enddo

      c_A=c_loc(acu)
      cusolver_stat=cudaMalloc(d_A,sizeof(acu))
      write(*,997) cusolver_stat
  997 format('allstat',i5)
      cusolver_stat=cudaMemcpy(d_A,c_A,sizeof(acu),1)
      write(*,999) cusolver_stat
      c_U=c_loc(ucu)
      cusolver_stat=cudaMalloc(d_U,sizeof(ucu))
      write(*,997) cusolver_stat
      c_S=c_loc(scu)
      cusolver_stat=cudaMalloc(d_S,sizeof(scu))
      write(*,997) cusolver_stat
      c_W=c_loc(wcu)
      cusolver_stat=cudaMalloc(d_W,sizeof(wcu))
      write(*,997) cusolver_stat
      cusolver_stat=cudaMalloc(d_WK,lworkc)
      write(*,997) cusolver_stat
      cusolver_stat=cudaMalloc(d_RW,nelecs)
      write(*,997) cusolver_stat
      cusolver_stat=cudaMalloc(d_devinfo,sizeof(devinfo))
      write(*,997) cusolver_stat
      ndim=nelecs
      ldim=lworkc
      write(*,995) ndim,ldim,sizeof(acu),sizeof(ucu),sizeof(scu),
     & sizeof(wcu)
  995 format('dims',6i8)

      write(*,1000) ((acu((j-1)*nelecs+i),i=1,5),j=1,5)

      cusolver_stat=cusolverDnDgesvd(cusolver_Hndl,c_loc(jobu),
     & c_loc(jobvt),ndim,ndim,
     & d_A,ndim,d_S,d_U,ndim,d_W,ndim,d_WK,ldim,
     & d_RW,d_devinfo)

      write(*,996) cusolver_stat
c,$CUSOLVER_STATUS_SUCCESS
  996 format('svdstat',2i5)

      write(*,998) devinfo
  998 format('devinfo',i5)

      cusolver_stat=cudaMemcpy(c_A,d_A,sizeof(acu),2)

      write(*,999) cusolver_stat
  999 format('cpystat',i5)

      cusolver_stat=cudaMemcpy(c_U,d_U,sizeof(ucu),2)

      write(*,999) cusolver_stat

      cusolver_stat=cudaMemcpy(c_W,d_W,sizeof(wcu),2)

      write(*,999) cusolver_stat

      cusolver_stat=cudaMemcpy(c_S,d_S,sizeof(scu),2)

      write(*,999) cusolver_stat



      write(*,1000) ((acu((j-1)*nelecs+i),i=1,5),j=1,5)
      write(*,1001) ((ucu((j-1)*nelecs+i),i=1,5),j=1,5)
      write(*,1002) ((wcu((j-1)*nelecs+i),i=1,5),j=1,5)
      write(*,1003) (scu(i),i=1,5)
 1000 format('A',/,(5f12.6))
 1001 format('U',/,(5f12.6))
 1002 format('W',/,(5f12.6))
 1003 format('S',/,(5f12.6))

      endif

      call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     & ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     & ierr,sdiag,nelecs)


      if(itest.gt.0) then

       write(*,1000) ((ta(i,j),i=1,5),j=1,5)
       write(*,1001) ((u(i,j),i=1,5),j=1,5)
       write(*,1002) ((w(i,j),i=1,5),j=1,5)
       write(*,1003) (sdiag(i),i=1,5)

c       stop
       itest=itest-1
       if(itest.eq.0) stop
      endif

#else

      call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     & ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     & ierr,sdiag,nelecs)

#endif
      call timer_stop(41)

      if(idbg.ge.12) write(lfndbg,601) (ev(i),i=1,nelecs)
 601  format(//,' Eigenvalues of diagonalized overlap matrix:',
     & //,(3x,6e20.12))

!     Calculation of det(uw) , by determination of the number
!     of eigenvalues -2 of " a=uw+tranposed(uw) "

      call timer_start(42)
      idetuw=1
c!$acc kernels copyin(u,w) copyout(ta)
      do i=1,nelecs
c!$acc loop vector
       do j=1,i
        coef=0.0d0
c!$acc loop seq
        do k=1,nelecs
         coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
        enddo
        ta(i,j)=coef
       enddo
      enddo
c!$acc end kernels
      call timer_stop(42)

      call timer_start(43)
      call tred2(nelecs,nelecs,ta,mbasel,max(mbasel,nveca),diag,nelecs,
     & sdiag,nelecs,ta,mbasel,max(mbasel,nveca))
      call timer_stop(43)

      call timer_start(44)
      call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,ta,mbasel,
     & max(mbasel,nveca),ierr)
      call timer_stop(44)

      call timer_start(45)
      do i=1,nelecs
       if(diag(i).lt.-1.999999d0) idetuw=-idetuw
      enddo

      do i=1,nelecs
       cdiag(i)=diag(i)
       csdiag(i)=sdiag(i)
      enddo

!      calculation of det(a) and x and y

      cmax=0.00

      do i=1,nelecs
       if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
      enddo

      if(cmax.le.0.01) stop 'no overlap between m.o.s'

      cnorm=thresh_SIN*cmax
      coef=idetuw
      nz1=0
      nz2=0
      deta=0.0

      do i=1,nelecs
       coefu=ev(i)
       if(abs(coefu) .le. cnorm) then
        if(nz1.gt.0.and.nz2.gt.0) then
         ising=3
         call timer_stop(45)
         return
        endif
        if(nz1 .gt. 0) nz2=i
        if(nz1.eq.0) nz1=i
       else
        coef=coef*coefu
       endif
      enddo

      ising=2
      if(nz2.le.0) then
       ising=1
       if(nz1.le.0) then
       deta=coef
       ising=0

cc !$acc kernels copyout(ta) copyin(u,w,ev)
       do i=1,nelecs
        do j=1,nelecs
         coefu=0.0d0
         do k=1,nelecs
          coefu=coefu+u(i,k)*w(j,k)/ev(k)
         enddo
         ta(i,j)=coefu*0.5d0
        enddo
       enddo
cc !$acc end kernels

       call timer_stop(45)
       return
       ising=1
      endif

cc !$acc kernels copyout(diag,sdiag,ta) copyin(u,w,ev)
      do i=1,nelecs
       diag(i)=u(i,nz1)*coef
       sdiag(i)=w(i,nz1)
       do j=1,nelecs
        coefu=0.0d0
        do k=1,nelecs
         if(k.ne.nz1) coefu=coefu+u(i,k)*w(j,k)/ev(k)
        enddo
        ta(i,j)=coefu
       enddo
      enddo
cc !$acc end kernels

       do i=1,nelecs
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
       enddo
       call timer_stop(45)
      return
      endif

c      if(abs(coef) .lt. 1.0d-12) then
      if(abs(coef).lt.thresh_SIN) then
       ising=3
      else
       do i=1,nelecs
        diag(i)=u(i,nz1)*coef
        sdiag(i)=w(i,nz1)
        coefu=u(i,nz2)
        do j=1,nelecs
         ta(i,j)=coefu*w(j,nz2)
        enddo
       enddo
       do i=1,nelecs
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
       enddo
       call timer_stop(45)
       return
      endif

      call timer_stop(45)
      return
      end
