!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ bases, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen
      
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
      real (kind=8) :: coef
      real (kind=8) :: cmax, cnorm, coefu
      integer :: nz1, nz2
      
#ifdef CUSOLVER
      character (len=1), target :: jobu, jobvt
#endif
      
      if(idbg.ge.30) write(lfndbg,600)
 600  format(/,' Cofactor matrix will be calculated')

      call timer_start(41)

!_ACCTGT_($acc kernels present(sdiag))
      do i=1,nelecs
        sdiag(i)=0.0d0
      enddo
!_ACCTGT_($acc end kernels)


      if(idbg.ge.90) then
!_ACCTGT_($acc update host (a))
        write(lfndbg,1601) nelecs,nelecs,mbasel
 1601   format(//,' SVD input matrix:',3i6,/)
        do j=1,nelecs
          write(lfndbg,1602) (a(i,j),i=1,nelecs)
 1602     format((3x,6e20.12))
        enddo
        flush(lfndbg)
      endif
      
#ifdef CUSOLVER

       ndim=nelecs
       mdim=mbasel

       if(isolver.eq.1) then

       jobu = 'A'  ! all m columns of U
       jobvt= 'A'  ! all m columns of VT

!_ACCTGT_($acc data copy(dev_info_d) create(temp,workspace_d,rwork))
       
!_ACCTGT_($acc host_data use_device(a,ev,u,w,dev_info_d,workspace_d,rwork))
       cusolver_status=cusolverDnDgesvd(cusolver_handle,jobu,jobvt,
     &      ndim,ndim,a,ndim,ev,u,ndim,w,ndim,workspace_d,
     &      lwork1,rwork,dev_info_d)
!_ACCTGT_($acc end host_data)
!_ACCTGT_($acc wait)

       if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &      write(*,*) 'cusolverDnDgesvd failed',cusolver_status

      call timer_start(42)
      
!_ACCTGT_($acc kernels present(temp,u,w,ta))
!_ACCTGT_($acc loop collapse(2))
      do i=1,nelecs
        do j=1,nelecs
          temp(i,j)=w(j,i)
        enddo
      enddo
!_ACCTGT_($acc loop collapse(2))
      do i=1,nelecs
        do j=1,nelecs
          w(i,j)=temp(i,j)
        enddo
      enddo
!_ACCTGT_($acc end kernels)
!_ACCTGT_($acc parallel loop gang present(temp,u,w,a))
      do i=1,nelecs
!_ACCTGT_($acc loop vector)
        do j=1,i
          coef=0.0d0
!_ACCTGT_($acc loop seq)
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          a(i,j)=coef
          a(j,i)=coef
        enddo
      enddo
!_ACCTGT_($acc end parallel loop)

      call timer_stop(42)
      
!_ACCTGT_($acc end data)
#ifdef CUSOLVERJ
      elseif(isolver.eq.2) then

        jobz=CUSOLVER_EIG_MODE_VECTOR
        
!     here the Jacobi svd from cusolver should be implemented
        
!_ACCTGT_($acc data copy(dev_info_d,gesvdj_params) create(workspace_d,temp))

!_ACCTGT_($acc host_data)
!_ACCTGT_($acc& use_device(a,ev,u,w,dev_info_d,workspace_d))

        cusolver_status=cusolverDnDgesvdj
     &       (cusolver_handle,jobz,econ,
     &       ndim,ndim,a,ndim,ev,u,ndim,w,ndim,workspace_d,
     &       lwork1,dev_info_d,gesvdj_params)
        
!_ACCTGT_($acc end host_data)

        cusolver_status=cudaDeviceSynchronize()
        
        if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &       write(*,*) 'cusolverDnDgesvdj failed',cusolver_status
        
        cusolver_status = cusolverDnXgesvdjGetSweeps
     &       (cusolver_handle, gesvdj_params, exec_sweeps)
        cusolver_status = cusolverDnXgesvdjGetResidual
     &       (cusolver_handle, gesvdj_params, residual)

        call timer_start(42)
     
!_ACCTGT_($acc parallel loop gang present(temp,u,w,a))
        do i=1,nelecs
!_ACCTGT_($acc loop vector)
          do j=1,i
            coef=0.0d0
!_ACCTGT_($acc loop seq)
            do k=1,nelecs
              coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
            enddo
            a(i,j)=coef
            a(j,i)=coef
          enddo
        enddo
!_ACCTGT_($acc end parallel loop)

        call timer_stop(42)
        
!_ACCTGT_($acc end data)

      else
#endif

!_ACCTGT_($acc update host (a))
!_OMPTGT_($omp target update from(a))
        call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,
     &       ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &       ierr,sdiag,nelecs)
!_ACCTGT_($acc update device (ev,u,w))
!_OMPTGT_($omp target update to(ev,u,w))
       
       call timer_start(42)
       
!_ACCTGT_($acc parallel loop gang present(u,w,a))
!_OMPTGT_($omp target teams distribute parallel do private(coef))
       do i=1,nelecs
         do j=1,i
           coef=0.0d0
           do k=1,nelecs
             coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
           enddo
           a(i,j)=coef
           a(j,i)=coef
         enddo
       enddo
!_ACCTGT_($acc end parallel loop)
!_OMPTGT_($omp end target)
       call timer_stop(42)
       
      endif

#else
      
!_ACCTGT_($acc update host (a))
!_OMPTGT_($omp target update from(a))      
      call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,
     &     ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &     ierr,sdiag,nelecs)      
!_ACCTGT_($acc update device (ev,u,w))
!_OMPTGT_($omp target update to(ev,u,w))     
      call timer_start(42)
      
!_ACCTGT_($acc kernels present(u,w,a))
!_OMPTGT_($omp target teams distribute parallel do private(coef))
      do i=1,nelecs
        do j=1,i
          coef=0.0d0
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          a(i,j)=coef
          a(j,i)=coef
        enddo
      enddo
!_ACCTGT_($acc end kernels)
      call timer_stop(42)

#endif

      call timer_stop(41)

      if(idbg.ge.90) then
!_ACCTGT_($acc update host (u,w,a,ev))
        write(lfndbg,601) (ev(i),i=1,nelecs)
 601    format(//,' Eigenvalues of diagonalized overlap matrix:',
     &       //,(3x,6e20.12))
        write(lfndbg,1603) nelecs,nelecs,mbasel
 1603   format(//,' Matrix u:',3i6,/)
        do j=1,nelecs
          write(lfndbg,1602) (u(i,j),i=1,nelecs)
        enddo
        write(lfndbg,1604) nelecs,nelecs,mbasel
 1604   format(//,' Matrix w:',3i6,/)
        do j=1,nelecs
          write(lfndbg,1602) (w(i,j),i=1,nelecs)
        enddo
        write(lfndbg,1605) nelecs,nelecs,mbasel
 1605   format(//,' Coefficient matrix:',3i6,/)
        do j=1,nelecs
          write(lfndbg,1602) (a(i,j),i=1,nelecs)
        enddo
        flush(lfndbg)
      endif

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
!_ACCTGT_($acc host_data use_device(a,diag,dev_info_d,workspace_d))
        cusolver_status = cusolverDnDsyevd(cusolver_handle,
     &       CUSOLVER_EIG_MODE_NOVECTOR,CUBLAS_FILL_MODE_LOWER,
     &       ndim,a,ndim,diag,workspace_d,lwork2,dev_info_d)
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

!_ACCTGT_($acc host_data use_device(a,diag,dev_info_d,workspace_d))
        cusolver_status = cusolverDnDsyevj
     &       (cusolver_handle, jobz, uplo, ndim,a,ndim,diag,
     &       workspace_d,lwork2,dev_info_d,syevj_params)
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

!_ACCTGT_($acc update host (a))
!_OMPTGT_($omp target update from(a))
        call timer_start(43)
        call tred2(nelecs,nelecs,a,nelecs,nelecs,diag,
     &       nelecs,sdiag,nelecs,a,nelecs,nelecs)
        call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,
     &       a,nelecs,nelecs,ierr)
        call timer_stop(43)
!_ACCTGT_($acc update device (a,diag,sdiag))
!_OMPTGT_($omp target update to(a,diag,sdiag))
      endif
       
!_ACCTGT_($acc end data)

#else
      
!_ACCTGT_($acc update host (a))
!_OMPTGT_($omp target update from(a))
      call timer_start(43)
        call tred2(nelecs,nelecs,a,nelecs,nelecs,diag,
     &       nelecs,sdiag,nelecs,a,nelecs,nelecs)
        call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,
     &       a,nelecs,nelecs,ierr)
      call timer_stop(43)
!_ACCTGT_($acc update device (a,diag,sdiag))
!_OMPTGT_($omp target update to(a,diag,sdiag))

#endif
      
      
      call timer_start(44)
!_ACCTGT_($acc kernels present(diag))
!_OMPTGT_($omp target teams distribute parallel do)
      do i=1,nelecs
        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
      enddo
!_ACCTGT_($acc end kernels)
c!_OMPTGT_($omp end target)
      cmax=0.0d0
      
!_ACCTGT_($acc kernels present(cdiag,diag,csdiag,sdiag,ev))
!_OMPTGT_($omp target teams distribute parallel do)
!     calculation of det(a) and x and y
      do i=1,nelecs
        if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
      enddo
c!_OMPTGT_($omp end target)
!_OMPTGT_($omp target teams distribute parallel do)
      do i=1,nelecs
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
      enddo
!_ACCTGT_($acc end kernels)
c!_OMPTGT_($omp end target)
      if(cmax.le.0.01) call errquit(me,8,"No overlap between m.o.s")

      if(idbg.ge.90) then
!_ACCTGT_($acc update host (diag,cdiag,csdiag,sdiag))
        write(lfndbg,604)
     &       (diag(i),cdiag(i),csdiag(i),sdiag(i),i=1,nelecs)
 604    format(//,' Diagonals:',
     &       //,(3x,4e20.12))
        flush(lfndbg)
      endif
      
      cnorm=tau_SIN*cmax
      coef=idetuw
      nz1=0
      nz2=0
      deta=0.0d0

!_ACCTGT_($acc update host(ev))
!_OMPTGT_($omp target update from(ev))

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
c!_OMPTGT_($omp end target)
      ising=2
      if(nz2.le.0) then
        ising=1
        if(nz1.le.0) then
          deta=coef
          ising=0

!_OMPTGT_($omp target teams distribute parallel do collapse(2) private(coefu))
!_ACCTGT_($acc parallel loop gang collapse(2) present(ta,u,w,ev))
          do i=1,nelecs
            do j=1,nelecs
              coefu=0.0d0
              do k=1,nelecs
                coefu=coefu+u(i,k)*w(j,k)/ev(k)
              enddo
              ta(i,j)=coefu*0.5d0
            enddo
          enddo
!_ACCTGT_($acc end parallel loop)
          call timer_stop(44)
          return
          ising=1
        endif

!_ACCTGT_($acc kernels present(ta,u,w,ev,cdiag,diag,csdiag,sdiag))
!_OMPTGT_($omp target teams distribute parallel do private(coefu))
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
!_OMPTGT_($omp target teams distribute parallel do)
        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
c!_OMPTGT_($omp end target)
!_ACCTGT_($acc end kernels)
      
        call timer_stop(44)
        return
      endif

      if(abs(coef).lt.tau_SIN) then
        ising=3
      else
!_ACCTGT_($acc kernels present(cdiag,diag,csdiag,sdiag,u,w,ta))
!_OMPTGT_($omp target teams distribute parallel do private(coefu))
        do i=1,nelecs
          diag(i)=u(i,nz1)*coef
          sdiag(i)=w(i,nz1)
          coefu=u(i,nz2)
          do j=1,nelecs
            ta(i,j)=coefu*w(j,nz2)
          enddo
        enddo
!_OMPTGT_($omp target teams distribute parallel do)
        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
c!_OMPTGT_($omp end target)
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
      
      if(idbg.ge.30) write(lfndbg,600)
 600  format(/,' Cofactor matrix will be calculated')

      call timer_start(41)

      if(idbg.ge.90) then
        write(lfndbg,1601) nelecs,nelecs,mbasel
 1601   format(//,' SVD input matrix:',3i6,/)
        do j=1,nelecs
          write(lfndbg,1602) (a(i,j),i=1,nelecs)
 1602     format((3x,6e20.12))
        enddo
      endif
      
#ifdef MKL

      if(isolver.eq.1) then

        ndim=nelecs

        call dgesvd('All','All',ndim,ndim,a,ndim,ev,u,ndim,w,ndim,
     &       workspace_d,lwork1,ierr)
        
!$omp parallel do
!$omp& shared(temp,w) collapse(2)
        do i=1,nelecs
          do j=1,nelecs
            temp(i,j)=w(i,j)
          enddo
        enddo
!$omp end parallel do
!$omp parallel do
!$omp& shared(temp,w) collapse(2)
        do i=1,nelecs
          do j=1,nelecs
            w(j,i)=temp(i,j)
          enddo
        enddo
!$omp end parallel do
      else
      
        call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,
     &       ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &       ierr,sdiag,nelecs)
        
      endif

#else
      
      call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,
     &     ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,
     &     ierr,sdiag,nelecs)
      
#endif

      call timer_start(42)
      
!$omp parallel do
!$omp& shared(a,u,w) private(i,j,k,coef)
      do i=1,nelecs
        do j=1,i
          coef=0.0d0
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          a(i,j)=coef
        enddo
      enddo
!$omp end parallel do
      
      call timer_stop(42)

      call timer_stop(41)

      if(idbg.ge.90) write(lfndbg,601) (ev(i),i=1,nelecs)
 601  format(//,' Eigenvalues of diagonalized overlap matrix:',
     & //,(3x,6e20.12))

!     Calculation of det(uw) , by determination of the number
!     of eigenvalues -2 of " a=uw+tranposed(uw) "

      idetuw=1
      
      call timer_start(43)
#ifdef MKL
      if(jsolver.eq.1) then

        ndim=nelecs

        call dsyevd('N','L',ndim,a,nelecs,diag,
     &       workspace_d,lwork1,
     &       workspace_i,lworki,ierr)

      else

        call tred2(nelecs,nelecs,a,nelecs,nelecs,diag,
     &       nelecs,sdiag,nelecs,a,nelecs,nelecs)
        
        call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,a,nelecs,
     &       nelecs,ierr)
      endif
#else
      call tred2(nelecs,nelecs,a,nelecs,nelecs,diag,nelecs,
     &     sdiag,nelecs,a,nelecs,nelecs)

      call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,a,nelecs,
     &     nelecs,ierr)
#endif
      
      call timer_stop(43)

      call timer_start(44)
      do i=1,nelecs
        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
      enddo

!$omp parallel do
!$omp& shared(cdiag,csdiag,diag,sdiag)
      do i=1,nelecs
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
      enddo
!$omp end parallel do

!      calculation of det(a) and x and y

      cmax=0.00

      do i=1,nelecs
        if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
      enddo

      if(cmax.le.0.01) call errquit(me,8,"No overlap between m.o.s")

      cnorm=tau_SIN*cmax
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

!$omp parallel do
!$omp& shared(u,w,ev,ta) private(coefu) collapse(2)
          do i=1,nelecs
            do j=1,nelecs
              coefu=0.0d0
              do k=1,nelecs
                coefu=coefu+u(i,k)*w(j,k)/ev(k)
              enddo
              ta(i,j)=coefu*0.5d0
            enddo
          enddo
!$omp end parallel do

          call timer_stop(44)
          return
          ising=1
        endif

!$omp parallel do
!$omp& shared(u,w,ev,ta,diag,sdiag) private(coefu)
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
!$omp end parallel do

!$omp parallel do
!$omp& shared(cdiag,csdiag,diag,sdiag)
        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
!$omp end parallel do
      
        call timer_stop(44)
        return
      endif

      if(abs(coef).lt.tau_SIN) then
        ising=3
      else

!$omp parallel do
!$omp& shared(u,w,ev,ta,diag,sdiag,cdiag,csdiag) private(coefu)
        do i=1,nelecs
          diag(i)=u(i,nz1)*coef
          sdiag(i)=w(i,nz1)
          coefu=u(i,nz2)
          do j=1,nelecs
            ta(i,j)=coefu*w(j,nz2)
          enddo
        enddo
!$omp end parallel do
 
!$omp parallel do
!$omp& shared(cdiag,csdiag,diag,sdiag)
        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
!$omp end parallel do
      
        call timer_stop(44)
        return
      endif

      call timer_stop(44)
      
      return
      end
