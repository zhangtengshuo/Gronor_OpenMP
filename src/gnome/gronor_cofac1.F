      subroutine gronor_cofac1(lfndbg)
      use cidist
      use gnome_parameters
      use gnome_data

#ifdef CUSOLVER
      use cusolverDn
      use cuda_cusolver
#endif
      implicit none
      integer :: lfndbg
      integer :: i,ierr,j,idetuw,k
      real (kind=8) :: coef
      real (kind=8) :: cmax, cnorm, coefu
      integer :: nz1, nz2
      character (len=1), target :: jobu, jobvt
      
      if(idbg.ge.13) write(lfndbg,600)
 600  format(/,' Cofactor matrix will be calculated')

      call timer_start(41)
      
#ifdef CUSOLVER

!$acc data copy(dev_info_d)
!$acc& create(temp,workspace_d,rwork)
      
      if(isolver.eq.1) then

       jobu = 'A'  ! all m columns of U
       jobvt= 'A'  ! all m columns of VT

       ndim=nelecs
       mdim=mbasel

!$acc kernels present(temp,ta)
       do j=1,ndim
         do i=1,ndim
           temp(i,j)=ta(i,j)
         enddo
       enddo
!$acc end kernels
       
!$acc host_data use_device(temp,ev,u,w,dev_info_d,workspace_d,rwork)
       cusolver_status=cusolverDnDgesvd(cusolver_handle,jobu,jobvt,
     &      ndim,ndim,temp,ndim,ev,u,ndim,w,ndim,workspace_d,
     &      worksize,rwork,dev_info_d)
!$acc end host_data
!$acc wait

       if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &      write(*,*) 'cusolverDnDgesvd failed',cusolver_status

      call timer_start(42)
      
!$acc kernels present(temp,w)
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
!$acc end kernels
     
!$acc kernels present(u,w,ta)
      do i=1,nelecs
!$acc loop vector
        do j=1,i
          coef=0.0d0
!$acc loop seq
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          ta(i,j)=coef
        enddo
      enddo
!$acc end kernels

      call timer_stop(42)
      
      elseif(isolver.eq.2) then

!     here the Jacobi svd from cusolver should be implemented
        
      else

!$acc update host (ta)

        call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     &       ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &       ierr,sdiag,nelecs)
        
!$acc update device (ev,u,w)
       
       call timer_start(42)
       
!$acc kernels present(u,w,ta)
       do i=1,nelecs
         do j=1,i
           coef=0.0d0
           do k=1,nelecs
             coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
           enddo
           ta(i,j)=coef
         enddo
       enddo
!$acc end kernels

       call timer_stop(42)
       
      endif

!$acc end data
      
#else
      
!$acc update host (ta)
      
      call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     &     ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &     ierr,sdiag,nelecs)
      
!$acc update device (ev,u,w)
      
      call timer_start(42)
      
!$acc kernels present(u,w,ta)
      do i=1,nelecs
        do j=1,i
          coef=0.0d0
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          ta(i,j)=coef
        enddo
      enddo
!$acc end kernels
      
      call timer_stop(42)

#endif
      
      call timer_stop(41)

      if(idbg.ge.12) write(lfndbg,601) (ev(i),i=1,nelecs)
 601  format(//,' Eigenvalues of diagonalized overlap matrix:',
     & //,(3x,6e20.12))

!     Calculation of det(uw) , by determination of the number
!     of eigenvalues -2 of " a=uw+tranposed(uw) "

      idetuw=1
      
#ifdef CUSOLVER

!$acc data copy(dev_info_d)
!$acc& create(temp,workspace_d,rwork)
      
      if(jsolver.eq.1) then
      elseif(jsolver.eq.2) then
        
       ndim=nelecs
       mdim=mbasel

!$acc kernels present(temp,ta)
       do j=1,ndim
         do i=1,ndim
           temp(i,j)=ta(i,j)
         enddo
       enddo
!$acc end kernels
       
!$acc host_data use_device(temp,ta,w,diag,dev_info_d,workspace_d,rwork)
       cusolver_status = cusolverDnDsyevd(cusolver_handle,
     &       CUSOLVER_EIG_MODE_VECTOR,CUBLAS_FILL_MODE_LOWER,
     &   ndim,temp,mdim,diag,workspace_d,worksize2,dev_info_d)
!$acc end host_data
       
      else

!$acc update host (ta)

      call timer_start(43)
      call tred2(nelecs,nelecs,ta,mbasel,max(mbasel,nveca),diag,nelecs,
     &     sdiag,nelecs,ta,mbasel,max(mbasel,nveca))
      call timer_stop(43)

      call timer_start(44)
      call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,ta,mbasel,
     &     max(mbasel,nveca),ierr)
      call timer_stop(44)

!$acc update device (ta,diag,sdiag)

      endif
       
!$acc end data
      
#else
      
!$acc update host (ta)

      call timer_start(43)
      call tred2(nelecs,nelecs,ta,mbasel,max(mbasel,nveca),diag,nelecs,
     &     sdiag,nelecs,ta,mbasel,max(mbasel,nveca))
      call timer_stop(43)

      call timer_start(44)
      call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,ta,mbasel,
     &     max(mbasel,nveca),ierr)
      call timer_stop(44)

!$acc update device (ta,diag,sdiag)

#endif
      
      call timer_start(45)
      do i=1,nelecs
        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
      enddo

!$acc kernels present(cdiag,diag,csdiag,sdiag)
      do i=1,nelecs
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
      enddo
!$acc end kernels

!      calculation of det(a) and x and y

      cmax=0.00

!$acc kernels present(ev)
      do i=1,nelecs
        if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
      enddo
!$acc end kernels

      if(cmax.le.0.01) stop 'no overlap between m.o.s'

      cnorm=thresh_SIN*cmax
      coef=idetuw
      nz1=0
      nz2=0
      deta=0.0

!$acc update host(ev)
      
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

!$acc kernels present(ta,u,w,ev)
          do i=1,nelecs
            do j=1,nelecs
              coefu=0.0d0
              do k=1,nelecs
                coefu=coefu+u(i,k)*w(j,k)/ev(k)
              enddo
              ta(i,j)=coefu*0.5d0
            enddo
          enddo
!$acc end kernels

          call timer_stop(45)
          return
          ising=1
        endif

!$acc kernels present(diag,sdiag,ta,u,w,ev)
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
!$acc end kernels

!$acc kernels present(cdiag,diag,csdiag,sdiag)
        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
!$acc end kernels
      
        call timer_stop(45)
        return
      endif

      if(abs(coef).lt.thresh_SIN) then
        ising=3
      else
!$acc kernels present(cdiag,diag,csdiag,sdiag,u,w,ta)
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
!$acc end kernels
      
        call timer_stop(45)
        return
      endif

      call timer_stop(45)
      
      return
      end

      
      subroutine gronor_cofac1_omp(lfndbg)
      use cidist
      use gnome_parameters
      use gnome_data

      implicit none
      integer :: lfndbg
      integer :: i,ierr,j,idetuw,k
      real (kind=8) :: coef
      real (kind=8) :: cmax, cnorm, coefu
      integer :: nz1, nz2
      character (len=1), target :: jobu, jobvt
      
      if(idbg.ge.13) write(lfndbg,600)
 600  format(/,' Cofactor matrix will be calculated')

      call timer_start(41)

cc!$acc data copy(dev_info_d)
cc!$acc& create(temp,workspace_d,rwork)
      

      
cc!$acc update host (ta)
      
      call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     &     ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &     ierr,sdiag,nelecs)
      
cc!$acc update device (u,w)
      
      call timer_start(42)
      
cc!$acc kernels present(u,w,ta)
      do i=1,nelecs
        do j=1,i
          coef=0.0d0
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          ta(i,j)=coef
        enddo
      enddo
cc!$acc end kernels
      
      call timer_stop(42)

cc!$acc end data
      
      call timer_stop(41)

      if(idbg.ge.12) write(lfndbg,601) (ev(i),i=1,nelecs)
 601  format(//,' Eigenvalues of diagonalized overlap matrix:',
     & //,(3x,6e20.12))

!     Calculation of det(uw) , by determination of the number
!     of eigenvalues -2 of " a=uw+tranposed(uw) "

      idetuw=1
      
cc!$acc update host (ta)

      call timer_start(43)
      call tred2(nelecs,nelecs,ta,mbasel,max(mbasel,nveca),diag,nelecs,
     &     sdiag,nelecs,ta,mbasel,max(mbasel,nveca))
      call timer_stop(43)

      call timer_start(44)
      call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,ta,mbasel,
     &     max(mbasel,nveca),ierr)
      call timer_stop(44)

cc!$acc update device (ta,diag,sdiag)
      
      call timer_start(45)
      do i=1,nelecs
        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
      enddo

cc!$acc kernels present(cdiag,diag,csdiag,sdiag)
      do i=1,nelecs
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
      enddo
cc!$acc end kernels

!      calculation of det(a) and x and y

      cmax=0.00

cc!$acc kernels present(ev)
      do i=1,nelecs
        if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
      enddo
cc!$acc end kernels

      if(cmax.le.0.01) stop 'no overlap between m.o.s'

      cnorm=thresh_SIN*cmax
      coef=idetuw
      nz1=0
      nz2=0
      deta=0.0

cc!$acc update host(ev)
      
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

cc!$acc kernels present(ta,u,w,ev)
          do i=1,nelecs
            do j=1,nelecs
              coefu=0.0d0
              do k=1,nelecs
                coefu=coefu+u(i,k)*w(j,k)/ev(k)
              enddo
              ta(i,j)=coefu*0.5d0
            enddo
          enddo
cc!$acc end kernels

          call timer_stop(45)
          return
          ising=1
        endif

cc!$acc kernels present(diag,sdiag,ta,u,w,ev)
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
cc!$acc end kernels

cc!$acc kernels present(cdiag,diag,csdiag,sdiag)
        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
cc!$acc end kernels
      
        call timer_stop(45)
        return
      endif

      if(abs(coef).lt.thresh_SIN) then
        ising=3
      else
cc!$acc kernels present(cdiag,diag,csdiag,sdiag,u,w,ta)
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
cc!$acc end kernels
      
        call timer_stop(45)
        return
      endif

      call timer_stop(45)
      
      return
      end
