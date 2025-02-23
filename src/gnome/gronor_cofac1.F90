
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
      use gnome_solvers
      
      implicit none

      external :: timer_start,timer_stop
      external :: gronor_abort
      external :: gronor_svd,gronor_evd
      
      integer :: lfndbg
      integer :: i,j,idetuw,k
      real (kind=8) :: coef
      real (kind=8) :: cmax, cnorm, coefu
      integer :: nz1, nz2


      if(idbg.ge.30) write(lfndbg,600)
 600  format(/,' Cofactor matrix will be calculated')


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

      call timer_start(41)
      
      call gronor_svd()

      call timer_stop(41)

      !  Calculation of det(uw) by determination of the number of eigenvalues -2 of a=uw+transpose(uw)

#ifdef OMPTGT
!$omp target update from(u,w,a)
#endif

!#ifdef ACC
!!$acc kernels present(u,w,a)
!#endif
!#ifdef OMPTGT
!#ifdef OMP5
!!$omp target teams loop private(coef,k)
!#else
!!$omp target teams distribute parallel do private(coef,k)
!#endif
!#endif
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
!#ifdef OMPTGT
!#ifdef OMP5
!!$omp end target teams loop
!#else
!!$omp end target teams distribute parallel do
!#endif
!#endif
!#ifdef ACC
!!$acc end kernels
!#endif
      
#ifdef OMPTGT
!$omp target update to(a)
#endif


      if(idbg.ge.90) then
#ifdef ACC
!$acc update host (u,w,a,ev)
#endif
#ifdef OMPTGT 
!$omp target update from(a,ev,u,w)
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

      idetuw=1

      call timer_start(43)
      
      call gronor_evd()

      call timer_stop(43)
      
      call timer_start(44)

      cmax=0.0d0

#ifdef OMPTGT
!$omp target update from(diag,sdiag)
#endif
!#ifdef ACC
!!$acc kernels present(cdiag,diag,csdiag,sdiag,ev)
!#endif
!#ifdef OMPTGT
!#ifdef OMP5
!!$omp target teams loop reduction(max:cmax)
!#else
!!$omp target teams distribute parallel do reduction(max:cmax)
!#endif
!#endif
      
      !  Calculation of det(a) and x and y
      
      do i=1,nelecs
        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
        if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
      enddo
!#ifdef OMPTGT
!#ifdef OMP5
!!$omp end target teams loop
!#else
!!$omp end target teams distribute parallel do
!#endif
!#endif

!#ifdef ACC
!!$acc end kernels
!#endif
#ifdef OMPTGT
!$omp target update to(cdiag,csdiag)
#endif

      if(cmax.le.0.01) call gronor_abort(310,"No overlap between m.o.s")

      if(idbg.ge.90) then
#ifdef ACC
!$acc update host (diag,cdiag,csdiag,sdiag)
#endif
        write(lfndbg,604) (diag(i),cdiag(i),csdiag(i),sdiag(i),i=1,nelecs)
 604    format(//,' Diagonals:',//,(3x,4e20.12))
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

      ising=2
      if(nz2.le.0) then
        ising=1
        if(nz1.le.0) then
          deta=coef
          ising=0

#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop collapse(2) private(coefu)
#else
!$omp target teams distribute parallel do collapse(2) private(coefu)
#endif
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
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
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
#ifdef OMP5
!$omp target teams loop
#else
!$omp target teams distribute parallel do
#endif
#endif
        do i=1,nelecs
          diag(i)=u(i,nz1)*coef
          sdiag(i)=w(i,nz1)
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif

#ifdef OMPTGT
!$omp target update from(u,w,ev)
#endif


!#ifdef OMPTGT
!#ifdef OMP5
!!$omp target teams loop collapse(2) private(coefu)
!#else
!!$omp target teams distribute parallel do collapse(2) private(coefu)
!#endif
!#endif


        do i=1,nelecs
          do j=1,nelecs
            coefu=0.0d0
            do k=1,nelecs
              if(k.ne.nz1) coefu=coefu+u(i,k)*w(j,k)/ev(k)
            enddo
            ta(i,j)=coefu
          enddo
        enddo
!#ifdef OMPTGT
!#ifdef OMP5
!!$omp end target teams loop
!#else
!!$omp end target teams distribute parallel do
!#endif
!#endif
#ifdef ACC
!$acc end kernels
#endif
#ifdef OMPTGT
!$omp target update to(ta)
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
#ifdef OMP5
!$omp target teams loop private(coefu)
#else
!$omp target teams distribute parallel do private(coefu)
#endif
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
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
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
      end subroutine gronor_cofac1

      subroutine gronor_cofac1_omp(lfndbg)
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_solvers
#ifdef MKL
      use mkl_solver
#endif

      implicit none
      
      external :: timer_start,timer_stop
      external :: gronor_abort
      external :: gronor_svd,gronor_evd
      
      integer :: lfndbg
      integer :: i,j,idetuw,k
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

      call gronor_svd()
            
      call timer_start(42)

#ifdef OMP
!$omp parallel do shared(a,u,w) private(coef) schedule(dynamic)
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

      call gronor_evd()

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
      end subroutine gronor_cofac1_omp
