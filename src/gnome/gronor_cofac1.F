
!> @brief
!! Cofactor matrix evaluation and factorization
!!     
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!

#include "tgt_guard.fh"
      
      subroutine gronor_cofac1(lfndbg)
      use cidist
      use gnome_parameters
      use gnome_data

#ifdef CUSOLVER
      use cusolverDn
      use cuda_cusolver
      use cudafor
#endif
      implicit none
      integer :: lfndbg
      integer :: i,ierr,j,idetuw,k
      real (kind=8) :: coef,td
      real (kind=8) :: cmax, cnorm, coefu
      integer :: nz1, nz2
      character (len=1), target :: jobu, jobvt
      
      if(idbg.ge.13) write(lfndbg,600)
 600  format(/,' Cofactor matrix will be calculated')

      call timer_start(41)

#ifdef CUSOLVER

       ndim=nelecs
       mdim=mbasel

       if(isolver.eq.1) then

       jobu = 'A'  ! all m columns of U
       jobvt= 'A'  ! all m columns of VT

!_ACCTGT_($acc data copy(dev_info_d) create(temp,workspace_d,rwork))
       
!_ACCTGT_($acc host_data use_device(ta,ev,u,w,dev_info_d,workspace_d,rwork))
       cusolver_status=cusolverDnDgesvd(cusolver_handle,jobu,jobvt,
     &      ndim,ndim,ta,mdim,ev,u,ndim,w,ndim,workspace_d,
     &      worksize,rwork,dev_info_d)
!_ACCTGT_($acc end host_data)
!_ACCTGT_($acc wait)

       if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &      write(*,*) 'cusolverDnDgesvd failed',cusolver_status

      call timer_start(42)
      
!_ACCTGT_($acc kernels present(temp,u,w,ta))
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
      do i=1,nelecs
!_ACCTGT_($acc loop vector)
        do j=1,i
          coef=0.0d0
!_ACCTGT_($acc loop seq)
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          ta(i,j)=coef
          ta(j,i)=coef
        enddo
      enddo
!_ACCTGT_($acc end kernels)

      call timer_stop(42)
      
!_ACCTGT_($acc end data)
#ifdef CUSOLVERJ
      elseif(isolver.eq.2) then

        jobz=CUSOLVER_EIG_MODE_VECTOR
        
!     here the Jacobi svd from cusolver should be implemented
        
!_ACCTGT_($acc data copy(dev_info_d,gesvdj_params) create(workspace_d,temp))

!_ACCTGT_($acc host_data)
!_ACCTGT_($acc& use_device(ta,ev,u,w,dev_info_d,workspace_d))

        cusolver_status=cusolverDnDgesvdj
     &       (cusolver_handle,jobz,econ,
     &       ndim,ndim,ta,mdim,ev,u,ndim,w,ndim,workspace_d,
     &       worksize,dev_info_d,gesvdj_params)
        
!_ACCTGT_($acc end host_data)

        cusolver_status=cudaDeviceSynchronize()
        
        if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &       write(*,*) 'cusolverDnDgesvdj failed',cusolver_status
        
        cusolver_status = cusolverDnXgesvdjGetSweeps
     &       (cusolver_handle, gesvdj_params, exec_sweeps)
        cusolver_status = cusolverDnXgesvdjGetResidual
     &       (cusolver_handle, gesvdj_params, residual)

        call timer_start(42)
     
!_ACCTGT_($acc kernels present(u,w,ta))
        do i=1,nelecs
!_ACCTGT_($acc loop vector)
          do j=1,i
            coef=0.0d0
!_ACCTGT_($acc loop seq)
            do k=1,nelecs
              coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
            enddo
            ta(i,j)=coef
            ta(j,i)=coef
          enddo
        enddo
!_ACCTGT_($acc end kernels)

        call timer_stop(42)
        
!_ACCTGT_($acc end data)

      else
#endif

!_ACCTGT_($acc update host (ta))

        call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     &       ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &       ierr,sdiag,nelecs)
        
!_ACCTGT_($acc update device (ev,u,w))
       
       call timer_start(42)
       
!_ACCTGT_($acc kernels present(u,w,ta))
       do i=1,nelecs
         do j=1,i
           coef=0.0d0
           do k=1,nelecs
             coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
           enddo
           ta(i,j)=coef
           ta(j,i)=coef
         enddo
       enddo
!_ACCTGT_($acc end kernels)

       call timer_stop(42)
       
      endif

#else
      
!_ACCTGT_($acc update host (ta))
      
      call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     &     ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &     ierr,sdiag,nelecs)
      
!_ACCTGT_($acc update device (ev,u,w))
      
      call timer_start(42)
      
!_ACCTGT_($acc kernels present(u,w,ta))
      do i=1,nelecs
        do j=1,i
          coef=0.0d0
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          ta(i,j)=coef
          ta(j,i)=coef
        enddo
      enddo
!_ACCTGT_($acc end kernels)
      
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

!_ACCTGT_($acc data copy(dev_info_d))
!_ACCTGT_($acc& create(workspace_d))
      
      if(jsolver.eq.1) then
        
        ndim=nelecs
        mdim=mbasel
        
        call timer_start(43)

!_ACCTGT_($acc wait)
!_ACCTGT_($acc host_data use_device(ta,diag,dev_info_d,workspace_d))
        cusolver_status = cusolverDnDsyevd(cusolver_handle,
     &       CUSOLVER_EIG_MODE_NOVECTOR,CUBLAS_FILL_MODE_LOWER,
     &       ndim,ta,mdim,diag,workspace_d,worksize2,dev_info_d)
!_ACCTGT_($acc end host_data)
!_ACCTGT_($acc wait)

       if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &       write(*,*) 'cusolverDnDsyevd failed',cusolver_status
       
        call timer_stop(43)

!_ACCTGT_($acc kernels present(sdiag))
        do i=1,ndim
          sdiag(i)=0.0d0
        enddo
!_ACCTGT_($acc end kernels)

#ifdef CUSOLVERJ
      elseif(jsolver.eq.2) then

        ndim=nelecs
        mdim=mbasel
        
        jobz = CUSOLVER_EIG_MODE_NOVECTOR
        uplo = CUBLAS_FILL_MODE_LOWER
          
        call timer_start(43)

!_ACCTGT_($acc data copy(dev_info_d,syevj_params) create(workspace_d))

!_ACCTGT_($acc host_data use_device(ta,diag,dev_info_d,workspace_d))
        cusolver_status = cusolverDnDsyevj
     &       (cusolver_handle, jobz, uplo, ndim,ta,mdim,diag,
     &       workspace_d,worksize2,dev_info_d,syevj_params)
!_ACCTGT_($acc end host_data)

        cusolver_status=cudaDeviceSynchronize()
        
       if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &       write(*,*) 'cusolverDnDsyevj failed',cusolver_status
       
        call timer_stop(43)

!_ACCTGT_($acc kernels present(sdiag))
        do i=1,ndim
          sdiag(i)=0.0d0
        enddo
!_ACCTGT_($acc end kernels)
        
!_ACCTGT_($acc end data)
        
#endif
      else

!_ACCTGT_($acc update host (ta))

        call timer_start(43)
        call tred2(nelecs,nelecs,ta,mbasel,max(mbasel,nveca),diag,
     &       nelecs,sdiag,nelecs,ta,mbasel,max(mbasel,nveca))
        
        call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,ta,mbasel,
     &       max(mbasel,nveca),ierr)
        call timer_stop(43)
      
!_ACCTGT_($acc update device (ta,diag,sdiag))

      endif
       
!_ACCTGT_($acc end data)

#else
      
!_ACCTGT_($acc update host (ta))

      call timer_start(43)
      call tred2(nelecs,nelecs,ta,mbasel,max(mbasel,nveca),diag,nelecs,
     &     sdiag,nelecs,ta,mbasel,max(mbasel,nveca))

      call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,ta,mbasel,
     &     max(mbasel,nveca),ierr)
      call timer_stop(43)

!_ACCTGT_($acc update device (ta,diag,sdiag))

#endif
      
      call timer_start(44)
!_ACCTGT_($acc kernels present(diag))
      do i=1,nelecs
        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
      enddo
!_ACCTGT_($acc end kernels)

      cmax=0.0d0
      
!_ACCTGT_($acc kernels present(cdiag,diag,csdiag,sdiag,ev))
!     calculation of det(a) and x and y
      do i=1,nelecs
        if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
      enddo
      do i=1,nelecs
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
      enddo
!_ACCTGT_($acc end kernels)

      if(cmax.le.0.01) stop 'no overlap between m.o.s'

      cnorm=thresh_SIN*cmax
      coef=idetuw
      nz1=0
      nz2=0
      deta=0.0d0

!_ACCTGT_($acc update host(ev))
      
      do i=1,nelecs
        coefu=ev(i)
        if(abs(coefu) .le. cnorm) then
          if(nz1.gt.0.and.nz2.gt.0) then
            ising=3
            call timer_stop(44)
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

!_ACCTGT_($acc kernels present(ta,u,w,ev))
          do i=1,nelecs
            do j=1,nelecs
              coefu=0.0d0
              do k=1,nelecs
                coefu=coefu+u(i,k)*w(j,k)/ev(k)
              enddo
              ta(i,j)=coefu*0.5d0
            enddo
          enddo
!_ACCTGT_($acc end kernels)

          call timer_stop(44)
          return
          ising=1
        endif

!_ACCTGT_($acc kernels present(ta,u,w,ev,cdiag,diag,csdiag,sdiag))
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
        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
!_ACCTGT_($acc end kernels)
      
        call timer_stop(44)
        return
      endif

      if(abs(coef).lt.thresh_SIN) then
        ising=3
      else
!_ACCTGT_($acc kernels present(cdiag,diag,csdiag,sdiag,u,w,ta))
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
!_ACCTGT_($acc end kernels)
      
        call timer_stop(44)
        return
      endif

      call timer_stop(44)
      
      return
      end
      
      subroutine gronor_cofac1_omp(lfndbg)
      use cidist
      use gnome_parameters
      use gnome_data
#ifdef MKL
      use mkl_solver
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

#ifdef MKL

      if(isolver.eq.1) then

        ndim=nelecs
        mdim=mbasel

        call dgesvd('All','All',ndim,ndim,ta,mdim,ev,u,ndim,w,ndim,
     &       workspace_d,worksize,ierr)
        
        do i=1,nelecs
          do j=1,nelecs
            a(i,j)=w(i,j)
          enddo
        enddo
        do i=1,nelecs
          do j=1,nelecs
            w(j,i)=a(i,j)
          enddo
        enddo
      else
      
        call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     &       ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &       ierr,sdiag,nelecs)
        
      endif

#else
      
      call svd(nelecs,nelecs,nelecs,ta,mbasel,max(mbasel,nveca),
     &     ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &     ierr,sdiag,nelecs)
      
#endif

      call timer_start(42)
      
      do i=1,nelecs
        do j=1,i
          coef=0.0d0
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          ta(i,j)=coef
        enddo
      enddo
      
      call timer_stop(42)

      call timer_stop(41)

      if(idbg.ge.12) write(lfndbg,601) (ev(i),i=1,nelecs)
 601  format(//,' Eigenvalues of diagonalized overlap matrix:',
     & //,(3x,6e20.12))

!     Calculation of det(uw) , by determination of the number
!     of eigenvalues -2 of " a=uw+tranposed(uw) "

      idetuw=1
      
      call timer_start(43)
#ifdef MKL
      if(jsolver.eq.1) then

        ndim=nelecs
        mdim=mbasel

        call dsyevd('N','U',ndim,ta,mdim,diag,
     &       workspace_d,worksize,
     &       workspace_i,worksizei,ierr)

!        stop 'No correctly working MKL evd yet'

      else

        call tred2(nelecs,nelecs,ta,mbasel,max(mbasel,nveca),diag,
     &       nelecs,sdiag,nelecs,ta,mbasel,max(mbasel,nveca))
        
        call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,ta,mbasel,
     &       max(mbasel,nveca),ierr)
      endif
#else
      call tred2(nelecs,nelecs,ta,mbasel,max(mbasel,nveca),diag,nelecs,
     &     sdiag,nelecs,ta,mbasel,max(mbasel,nveca))

      call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,ta,mbasel,
     &     max(mbasel,nveca),ierr)
#endif
      
      call timer_stop(43)

      call timer_start(44)
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
            call timer_stop(44)
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

          do i=1,nelecs
            do j=1,nelecs
              coefu=0.0d0
              do k=1,nelecs
                coefu=coefu+u(i,k)*w(j,k)/ev(k)
              enddo
              ta(i,j)=coefu*0.5d0
            enddo
          enddo

          call timer_stop(44)
          return
          ising=1
        endif

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

        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
      
        call timer_stop(44)
        return
      endif

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
      
        call timer_stop(44)
        return
      endif

      call timer_stop(44)
      
      return
      end
