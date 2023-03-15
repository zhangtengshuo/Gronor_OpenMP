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


      subroutine gronor_cofac1(lfndbg)
      use cidist
      use gnome_parameters
      use gnome_data
#ifdef MKL
      use mkl_solver
#endif

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

#ifdef ACC
!$acc kernels present(sdiag)
#endif
      do i=1,nelecs
        sdiag(i)=0.0d0
      enddo
#ifdef ACC
!$acc end kernels
#endif
      if(idbg.ge.90) then
#ifdef ACC
!$acc update host (a)
#endif
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

#ifdef ACC
!$acc data copy(dev_info_d) create(temp,workspace_d,rwork)
!$acc host_data use_device(a,ev,u,w,dev_info_d,workspace_d,rwork)
#endif
       cusolver_status=cusolverDnDgesvd(cusolver_handle,jobu,jobvt,             &
     &      ndim,ndim,a,ndim,ev,u,ndim,w,ndim,workspace_d,                      &
     &      lwork1,rwork,dev_info_d)
#ifdef ACC
!$acc end host_data
!$acc wait
#endif
       if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)                           &
     &      write(*,*) 'cusolverDnDgesvd failed',cusolver_status

      call timer_start(42)

#ifdef ACC
!$acc kernels present(temp,u,w,ta)
!$acc loop collapse(2)
#endif
      do i=1,nelecs
        do j=1,nelecs
          temp(i,j)=w(j,i)
        enddo
      enddo
#ifdef ACC
!$acc loop collapse(2)
#endif
      do i=1,nelecs
        do j=1,nelecs
          w(i,j)=temp(i,j)
        enddo
      enddo
#ifdef ACC
!$acc end kernels
!$acc parallel loop gang present(temp,u,w,a)
#endif
      do i=1,nelecs
#ifdef ACC
!$acc loop vector
#endif
        do j=1,i
          coef=0.0d0
#ifdef ACC
!$acc loop seq
#endif
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          a(i,j)=coef
          a(j,i)=coef
        enddo
      enddo
#ifdef ACC
!$acc end parallel loop
#endif
      call timer_stop(42)

#ifdef ACC
!$acc end data
#endif
#ifdef CUSOLVERJ
      elseif(isolver.eq.2) then

        jobz=CUSOLVER_EIG_MODE_VECTOR

!     here the Jacobi svd from cusolver should be implemented

#ifdef ACC
!$acc data copy(dev_info_d,gesvdj_params) create(workspace_d,temp)
!$acc host_data
!$acc& use_device(a,ev,u,w,dev_info_d,workspace_d)
#endif
        cusolver_status=cusolverDnDgesvdj                                       &
     &       (cusolver_handle,jobz,econ,                                        &
     &       ndim,ndim,a,ndim,ev,u,ndim,w,ndim,workspace_d,                     &
     &       lwork1,dev_info_d,gesvdj_params)

#ifdef ACC
!$acc end host_data
#endif
        cusolver_status=cudaDeviceSynchronize()

        if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)                          &
     &       write(*,*) 'cusolverDnDgesvdj failed',cusolver_status

        cusolver_status = cusolverDnXgesvdjGetSweeps                            &
     &       (cusolver_handle, gesvdj_params, exec_sweeps)
        cusolver_status = cusolverDnXgesvdjGetResidual                          &
     &       (cusolver_handle, gesvdj_params, residual)

        call timer_start(42)

#ifdef ACC
!$acc parallel loop gang present(temp,u,w,a)
#endif
        do i=1,nelecs
#ifdef ACC
!$acc loop vector
#endif
          do j=1,i
            coef=0.0d0
#ifdef ACC
!$acc loop seq
#endif
            do k=1,nelecs
              coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
            enddo
            a(i,j)=coef
            a(j,i)=coef
          enddo
        enddo
#ifdef ACC
!$acc end parallel loop
#endif
        call timer_stop(42)

#ifdef ACC
!$acc end data
#endif
      else
#endif

#ifdef ACC
!$acc update host (a)
#endif
#ifdef OMPTGT
!$omp target update from(a)
#endif
#ifdef MKL

        ndimm=nelecs

        call dgesvd('All','All',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,         &
     &       workspace_d,lwork1m,ierr)

        do i=1,nelecs
          do j=1,nelecs
            temp(i,j)=w(i,j)
          enddo
        enddo
        do i=1,nelecs
          do j=1,nelecs
            w(j,i)=temp(i,j)
          enddo
        enddo

#else
         call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,                         &
     &       ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,           &
     &        ierr,sdiag,nelecs)
#endif

#ifdef ACC
!$acc update device (ev,u,w)
#endif
#ifdef OMPTGT
!$omp target update to(ev,u,w)
#endif
       call timer_start(42)

#ifdef ACC
!$acc parallel loop gang present(u,w,a)
#endif
#ifdef OMPTGT
!$omp target teams loop private(coef)
#endif
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
#ifdef OMPTGT
!$omp end target teams loop
#endif
#ifdef ACC
!$acc end parallel loop
#endif
#ifdef OMPTGT
!$omp end target
#endif
       call timer_stop(42)

      endif

#else

#ifdef ACC
!$acc update host (a)
#endif
#ifdef OMPTGT
!$omp target update from(a)
#endif
      call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,                            &
     &     ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,             &
     &     ierr,sdiag,nelecs)
#ifdef ACC
!$acc update device (ev,u,w)
#endif
#ifdef OMPTGT
!$omp target update to(ev,u,w)
#endif
      call timer_start(42)

#ifdef ACC
!$acc kernels present(u,w,a)
#endif
#ifdef OMPTGT
!$omp target teams loop private(coef,k)
#endif
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
#ifdef OMPTGT
!$omp end target teams loop
#endif
#ifdef ACC
!$acc end kernels
#endif
      call timer_stop(42)

#endif

      call timer_stop(41)

      if(idbg.ge.90) then
#ifdef ACC
!$acc update host (u,w,a,ev)
#endif
        write(lfndbg,601) (ev(i),i=1,nelecs)
 601    format(//,' Eigenvalues of diagonalized overlap matrix:',               &
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

#ifdef ACC
!$acc data copy(dev_info_d)                                                    &
!$acc& create(workspace_d)
#endif
      if(jsolver.eq.1) then

        ndim=nelecs
        mdim=mbasel

        call timer_start(43)

#ifdef ACC
!$acc wait
!$acc host_data use_device(a,diag,dev_info_d,workspace_d)
#endif
        cusolver_status = cusolverDnDsyevd(cusolver_handle,                     &
     &       CUSOLVER_EIG_MODE_NOVECTOR,CUBLAS_FILL_MODE_LOWER,                 &
     &       ndim,a,ndim,diag,workspace_d,lwork2,dev_info_d)
#ifdef ACC
!$acc end host_data
!$acc wait
#endif
       if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)                           &
     &       write(*,*) 'cusolverDnDsyevd failed',cusolver_status

        call timer_stop(43)

#ifdef ACC
!$acc kernels present(sdiag)
#endif
        do i=1,ndim
          sdiag(i)=0.0d0
        enddo
#ifdef ACC
!$acc end kernels
#endif
#ifdef CUSOLVERJ
      elseif(jsolver.eq.2) then

        ndim=nelecs
        mdim=mbasel

        jobz = CUSOLVER_EIG_MODE_NOVECTOR
        uplo = CUBLAS_FILL_MODE_LOWER

        call timer_start(43)

#ifdef ACC
!$acc data copy(dev_info_d,syevj_params) create(workspace_d)
!$acc host_data use_device(a,diag,dev_info_d,workspace_d)
#endif
        cusolver_status = cusolverDnDsyevj                                      &
     &       (cusolver_handle, jobz, uplo, ndim,a,ndim,diag,                    &
     &       workspace_d,lwork2,dev_info_d,syevj_params)
#ifdef ACC
!$acc end host_data
#endif
        cusolver_status=cudaDeviceSynchronize()

       if(cusolver_status /= CUSOLVER_STATUS_SUCCESS)                           &
     &       write(*,*) 'cusolverDnDsyevj failed',cusolver_status

        call timer_stop(43)

#ifdef ACC
!$acc kernels present(sdiag)
#endif
        do i=1,ndim
          sdiag(i)=0.0d0
        enddo
#ifdef ACC
!$acc end kernels
!$acc end data
#endif
#endif
      else

#ifdef ACC
!$acc update host (a)
#endif
#ifdef OMPTGT
!$omp target update from(a)
#endif
        call timer_start(43)
#ifdef MKL
        ndimm=nelecs
        call dsyevd('N','L',ndimm,a,nelecs,diag,                                &
     &       workspace_d,lwork1m,                                               &
     &       workspace_i,lworki,ierr)
#else
        call tred2(nelecs,nelecs,a,nelecs,nelecs,diag,                          &
     &       nelecs,sdiag,nelecs,a,nelecs,nelecs)
        call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,                       &
     &       a,nelecs,nelecs,ierr)
#endif
        call timer_stop(43)
#ifdef ACC
!$acc update device (a,diag,sdiag)
#endif
#ifdef OMPTGT
!$omp target update to(a,diag,sdiag)
#endif
      endif

#ifdef ACC
!$acc end data
#endif
#else

#ifdef ACC
!$acc update host (a)
#endif
#ifdef OMPTGT
!$omp target update from(a)
#endif
      call timer_start(43)
        call tred2(nelecs,nelecs,a,nelecs,nelecs,diag,                          &
     &       nelecs,sdiag,nelecs,a,nelecs,nelecs)
        call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,                       &
     &       a,nelecs,nelecs,ierr)
      call timer_stop(43)
#ifdef ACC
!$acc update device (a,diag,sdiag)
#endif
#ifdef OMPTGT
!$omp target update to(a,diag,sdiag)
#endif
#endif


      call timer_start(44)
!c!_ACCTGT_($acc kernels present(diag))
!cc!_OMPTGT_($omp target teams loop)
!c      do i=1,nelecs
!c        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
!c      enddo
!c!_ACCTGT_($acc end kernels)
!cc!_OMPTGT_($omp end target)

      cmax=0.0d0

#ifdef ACC
!$acc kernels present(cdiag,diag,csdiag,sdiag,ev)
#endif
#ifdef OMPTGT
!$omp target teams loop
#endif
!     calculation of det(a) and x and y
      do i=1,nelecs
        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
        if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
      enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif

!cc!_OMPTGT_($omp end target)
!c!_OMPTGT_($omp target teams loop)
!c      do i=1,nelecs
!c        cdiag(i)=diag(i)
!c        csdiag(i)=sdiag(i)
!c     enddo

#ifdef ACC
!$acc end kernels
#endif
!cc   !_OMPTGT_($omp end target)
      if(cmax.le.0.01) call gronor_abort(310,"No overlap between m.o.s")

      if(idbg.ge.90) then
#ifdef ACC
!$acc update host (diag,cdiag,csdiag,sdiag)
#endif
        write(lfndbg,604)                                                       &
     &       (diag(i),cdiag(i),csdiag(i),sdiag(i),i=1,nelecs)
 604    format(//,' Diagonals:',                                                &
     &       //,(3x,4e20.12))
        flush(lfndbg)
      endif

      cnorm=tau_SIN*cmax
      coef=idetuw
      nz1=0
      nz2=0
      deta=0.0d0

#ifdef ACC
!$acc update host(ev)
#endif
#ifdef OMPTGT
!$omp target update from(ev)
#endif
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
!!_OMPTGT_($omp end target)
      ising=2
      if(nz2.le.0) then
        ising=1
        if(nz1.le.0) then
          deta=coef
          ising=0

#ifdef OMPTGT
!$omp target teams loop collapse(2) private(coefu)
#endif
#ifdef ACC
!$acc parallel loop gang collapse(2) present(ta,u,w,ev)
#endif
          do i=1,nelecs
            do j=1,nelecs
              coefu=0.0d0
              do k=1,nelecs
                coefu=coefu+u(i,k)*w(j,k)/ev(k)
              enddo
              ta(i,j)=coefu*0.5d0
            enddo
          enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif
#ifdef ACC
!$acc end parallel loop
#endif
          call timer_stop(44)
          return
          ising=1
        endif

#ifdef ACC
!$acc kernels present(ta,u,w,ev,cdiag,diag,csdiag,sdiag)
#endif
#ifdef OMPTGT
!$omp target teams loop private(coefu)
#endif
        do i=1,nelecs
          diag(i)=u(i,nz1)*coef
          sdiag(i)=w(i,nz1)
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
          do j=1,nelecs
            coefu=0.0d0
            do k=1,nelecs
              if(k.ne.nz1) coefu=coefu+u(i,k)*w(j,k)/ev(k)
            enddo
            ta(i,j)=coefu
          enddo
        enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif
!c!_OMPTGT_($omp target teams loop)
!c       do i=1,nelecs
!c          cdiag(i)=diag(i)
!c          csdiag(i)=sdiag(i)
!c        enddo
!cc!_OMPTGT_($omp end target)
#ifdef ACC
!$acc end kernels
#endif
        call timer_stop(44)
        return
      endif

      if(abs(coef).lt.tau_SIN) then
        ising=3
      else
#ifdef ACC
!$acc kernels present(cdiag,diag,csdiag,sdiag,u,w,ta)
#endif
#ifdef OMPTGT
!$omp target teams loop private(coefu)
#endif
        do i=1,nelecs
          diag(i)=u(i,nz1)*coef
          sdiag(i)=w(i,nz1)
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
          coefu=u(i,nz2)
          do j=1,nelecs
            ta(i,j)=coefu*w(j,nz2)
          enddo
        enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif
!c!_OMPTGT_($omp target teams distribute parallel do)
!c        do i=1,nelecs
!c          cdiag(i)=diag(i)
!c          csdiag(i)=sdiag(i)
!c        enddo
!cc!_OMPTGT_($omp end target)
#ifdef ACC
!$acc end kernels
#endif
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

        ndimm=nelecs

        call dgesvd('All','All',ndimm,ndimm,a,ndimm,ev,u,ndimm,w,ndimm,         &
     &       workspace_d,lwork1m,ierr)
#ifdef OMP
!$omp parallel do                                                              &
!$omp& shared(temp,w) collapse(2)
#endif
        do i=1,nelecs
          do j=1,nelecs
            temp(i,j)=w(i,j)
          enddo
        enddo
#ifdef OMP
!$omp end parallel do
!$omp parallel do
!$omp& shared(temp,w) collapse(2)
#endif
        do i=1,nelecs
          do j=1,nelecs
            w(j,i)=temp(i,j)
          enddo
        enddo
#ifdef OMP
!$omp end parallel do
#endif
      else

        call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,                          &
     &       ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,           &
     &       ierr,sdiag,nelecs)

      endif

#else

      call svd(nelecs,nelecs,nelecs,a,nelecs,nelecs,                            &
     &     ev,nelecs,.true.,u,nelecs,nelecs,.true.,w,nelecs,nelecs,             &
     &     ierr,sdiag,nelecs)

#endif

      call timer_start(42)

#ifdef OMP
!$omp parallel do shared(a,u,w) private(coef)
#endif
      do i=1,nelecs
        do j=1,i
          coef=0.0d0
          do k=1,nelecs
            coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
          enddo
          a(i,j)=coef
        enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
      call timer_stop(42)

      call timer_stop(41)

      if(idbg.ge.90) write(lfndbg,601) (ev(i),i=1,nelecs)
 601  format(//,' Eigenvalues of diagonalized overlap matrix:',                 &
     & //,(3x,6e20.12))

!     Calculation of det(uw) , by determination of the number
!     of eigenvalues -2 of " a=uw+tranposed(uw) "

      idetuw=1

      call timer_start(43)
#ifdef MKL
      if(jsolver.eq.1) then

        ndimm=nelecs

        call dsyevd('N','L',ndimm,a,nelecs,diag,                                &
     &       workspace_d,lwork1m,                                               &
     &       workspace_i,lworki,ierr)

      else

        call tred2(nelecs,nelecs,a,nelecs,nelecs,diag,                          &
     &       nelecs,sdiag,nelecs,a,nelecs,nelecs)

        call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,a,nelecs,              &
     &       nelecs,ierr)
      endif
#else
      call tred2(nelecs,nelecs,a,nelecs,nelecs,diag,nelecs,                     &
     &     sdiag,nelecs,a,nelecs,nelecs)

      call tql2(nelecs,nelecs,diag,nelecs,sdiag,nelecs,a,nelecs,                &
     &     nelecs,ierr)
#endif

      call timer_stop(43)

      call timer_start(44)
      do i=1,nelecs
        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
      enddo

#ifdef OMP
!$omp parallel do shared(cdiag,csdiag,diag,sdiag)
#endif
      do i=1,nelecs
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
      enddo
#ifdef OMP
!$omp end parallel do
#endif
!      calculation of det(a) and x and y

      cmax=0.00

      do i=1,nelecs
        if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
      enddo

      if(cmax.le.0.01) call gronor_abort(310,"No overlap between m.o.s")

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

#ifdef OMP
!$omp parallel do shared(u,w,ev,ta) private(coefu) collapse(2)
#endif
          do i=1,nelecs
            do j=1,nelecs
              coefu=0.0d0
              do k=1,nelecs
                coefu=coefu+u(i,k)*w(j,k)/ev(k)
              enddo
              ta(i,j)=coefu*0.5d0
            enddo
          enddo
#ifdef OMP
!$omp end parallel do
#endif
          call timer_stop(44)
          return
          ising=1
        endif

#ifdef OMP
!$omp parallel do shared(u,w,ev,ta,diag,sdiag) private(coefu)
#endif
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
#ifdef OMP
!$omp end parallel do
#endif
#ifdef OMP
!$omp parallel do shared(cdiag,csdiag,diag,sdiag)
#endif
        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
#ifdef OMP
!$omp end parallel do
#endif
        call timer_stop(44)
        return
      endif

      if(abs(coef).lt.tau_SIN) then
        ising=3
      else

#ifdef OMP
!$omp parallel do shared(u,w,ta,diag,sdiag,coef) private(coefu)
#endif
        do i=1,nelecs
          diag(i)=u(i,nz1)*coef
          sdiag(i)=w(i,nz1)
          coefu=u(i,nz2)
          do j=1,nelecs
            ta(i,j)=coefu*w(j,nz2)
          enddo
        enddo
#ifdef OMP
!$omp end parallel do
#endif
#ifdef OMP
!$omp parallel do shared(cdiag,csdiag,diag,sdiag)
#endif
        do i=1,nelecs
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
#ifdef OMP
!$omp end parallel do
#endif
        call timer_stop(44)
        return
      endif

      call timer_stop(44)

      return
      end
