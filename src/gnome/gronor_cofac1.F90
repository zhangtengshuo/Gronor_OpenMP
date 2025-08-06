
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

      module gronor_cofac1_mod
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_solvers
      use gronor_svd_mod, only: gronor_svd
      use gronor_evd_mod, only: gronor_evd
      use omp_lib
      implicit none
      contains
      subroutine gronor_cofac1(lfndbg,a,u,w,wt,ev,ta,diag,sdiag,cdiag,csdiag)

      implicit none

      real (kind=8), intent(inout) :: a(:,:),u(:,:),w(:,:),wt(:,:),ev(:),ta(:,:)
      real (kind=8), intent(inout) :: diag(:),sdiag(:),cdiag(:),csdiag(:)

      external :: timer_start,timer_stop
      external :: gronor_abort

      integer :: lfndbg
      integer :: i,j,idetuw,k
      integer :: imax,jmax,iout,jout
      real (kind=8) :: coef
      real (kind=8) :: cmax, cnorm, coefu
      integer :: nz1, nz2
      integer :: thread_id
      character(len=10) :: today, now


      thread_id = omp_get_thread_num()
      if(idbg.gt.10 .and. thread_id==0) then
        call swatch(today,now)
        write(lfndbg,'(a,1x,a,a)') today(1:8),now(1:8), &
             ' Array dimensions check in gronor_cofac1:'
        write(lfndbg,'(a,2i10)') ' a:     ', size(a,1), size(a,2)
        write(lfndbg,'(a,2i10)') ' u:     ', size(u,1), size(u,2)
        write(lfndbg,'(a,2i10)') ' w:     ', size(w,1), size(w,2)
        write(lfndbg,'(a,2i10)') ' wt:    ', size(wt,1), size(wt,2)
        write(lfndbg,'(a,2i10)') ' ta:    ', size(ta,1), size(ta,2)
        write(lfndbg,'(a,1i10)') ' ev:    ', size(ev)
        write(lfndbg,'(a,1i10)') ' diag:  ', size(diag)
        write(lfndbg,'(a,1i10)') ' sdiag: ', size(sdiag)
        write(lfndbg,'(a,1i10)') ' cdiag: ', size(cdiag)
        write(lfndbg,'(a,1i10)') ' csdiag:', size(csdiag)
        write(lfndbg,'(a,5i10)') ' scalars nelecs mbasel ntcla ntclb nveca:', &
             nelecs, mbasel, ntcla, ntclb, nveca
        write(lfndbg,'(a,1x,es12.4)') ' a(1,1)=', a(1,1)
        write(lfndbg,'(a,1x,es12.4)') ' u(1,1)=', u(1,1)
        write(lfndbg,'(a,1x,es12.4)') ' w(1,1)=', w(1,1)
        imax = min(size(a,1),500)
        jmax = min(size(a,2),500)
        write(lfndbg,'(a)') ' a contents:'
        do iout=1,imax
          write(lfndbg,'(1x,*(es12.4))') (a(iout,jout),jout=1,jmax)
        enddo
        imax = min(size(u,1),500)
        jmax = min(size(u,2),500)
        write(lfndbg,'(a)') ' u contents:'
        do iout=1,imax
          write(lfndbg,'(1x,*(es12.4))') (u(iout,jout),jout=1,jmax)
        enddo
        imax = min(size(w,1),500)
        jmax = min(size(w,2),500)
        write(lfndbg,'(a)') ' w contents:'
        do iout=1,imax
          write(lfndbg,'(1x,*(es12.4))') (w(iout,jout),jout=1,jmax)
        enddo
        imax = min(size(wt,1),500)
        jmax = min(size(wt,2),500)
        write(lfndbg,'(a)') ' wt contents:'
        do iout=1,imax
          write(lfndbg,'(1x,*(es12.4))') (wt(iout,jout),jout=1,jmax)
        enddo
        imax = min(size(ta,1),500)
        jmax = min(size(ta,2),500)
        write(lfndbg,'(a)') ' ta contents:'
        do iout=1,imax
          write(lfndbg,'(1x,*(es12.4))') (ta(iout,jout),jout=1,jmax)
        enddo
        jmax = min(size(ev),500)
        write(lfndbg,'(a)') ' ev contents:'
        write(lfndbg,'(1x,*(es12.4))') (ev(jout),jout=1,jmax)
        jmax = min(size(diag),500)
        write(lfndbg,'(a)') ' diag contents:'
        write(lfndbg,'(1x,*(es12.4))') (diag(jout),jout=1,jmax)
        jmax = min(size(sdiag),500)
        write(lfndbg,'(a)') ' sdiag contents:'
        write(lfndbg,'(1x,*(es12.4))') (sdiag(jout),jout=1,jmax)
        jmax = min(size(cdiag),500)
        write(lfndbg,'(a)') ' cdiag contents:'
        write(lfndbg,'(1x,*(es12.4))') (cdiag(jout),jout=1,jmax)
        jmax = min(size(csdiag),500)
        write(lfndbg,'(a)') ' csdiag contents:'
        write(lfndbg,'(1x,*(es12.4))') (csdiag(jout),jout=1,jmax)
        flush(lfndbg)
      end if

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
      call gronor_svd(a,ev,u,w,sdiag,wt)
      call timer_stop(41)

  !  Calculation of det(uw) by determination of the number of eigenvalues -2 of a=uw+transpose(uw)
 !!! the MOST TIME COMSUMING CONSUMING loop if this subroutine!!!
! !$acc kernels present(u,w,a)
!   do i=1,nelecs
!     do j=1,i
!       coef=0.0d0
!       do k=1,nelecs
!         coef=coef+u(i,k)*w(k,j)+u(j,k)*w(k,i)
!       enddo
!       a(i,j)=coef
!       a(j,i)=coef
!     enddo
!   enddo
! !$acc end kernels
      
  ! "do j=1,i" 导致 j 依赖 i , 编译器无法有效并行化。
  ! 直接写 “do j = 1, nelecs” ，计算量翻倍，但是有效并行，速度反而快。
  ! 因为多出来的计算量是放在原来就没占满的计算单元上！

  !$acc parallel loop gang vector collapse(2) present(u, w, a) 
  do i = 1, nelecs
    do j = 1, nelecs
      coef = 0.0d0
      !$acc loop seq reduction(+:coef)
      do k = 1, nelecs
        coef = coef + u(i,k) * w(k,j) + u(j,k) * w(k,i)
      enddo
      a(i,j) = coef
    enddo
  enddo

  !可能的更快代码，但没有必要了。
!   !$acc parallel loop gang vector collapse(2) present(u,w,a)
! do i = 1, nelecs
!   do j = 1, nelecs
!     if (j <= i) then  ! 仅计算下三角
!       coef = 0.0d0
!       !$acc loop seq reduction(+:coef)
!       do k = 1, nelecs
!         coef = coef + u(i,k)*w(k,j) + u(j,k)*w(k,i)
!       enddo
!       a(i,j) = coef
!       if (i /= j) a(j,i) = coef  ! 对称赋值
!     endif
!   enddo
! enddo
      
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

      idetuw=1

      call timer_start(43)
      
      call gronor_evd(a,diag,sdiag)

      call timer_stop(43)
      
      call timer_start(44)

      cmax=0.0d0

#ifdef ACC
!$acc kernels present(cdiag,diag,csdiag,sdiag,ev)
#endif
      
      !  Calculation of det(a) and x and y
      
      do i=1,nelecs
        if(diag(i).lt.-1.999999d0) idetuw=-idetuw
        if(abs(ev(i)).gt.cmax) cmax=abs(ev(i))
        cdiag(i)=diag(i)
        csdiag(i)=sdiag(i)
      enddo

#ifdef ACC
!$acc end kernels
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
        do i=1,nelecs
          diag(i)=u(i,nz1)*coef
          sdiag(i)=w(i,nz1)
          cdiag(i)=diag(i)
          csdiag(i)=sdiag(i)
        enddo
        do i=1,nelecs
          do j=1,nelecs
            coefu=0.0d0
            do k=1,nelecs
              if(k.ne.nz1) coefu=coefu+u(i,k)*w(j,k)/ev(k)
            enddo
            ta(i,j)=coefu
          enddo
        enddo
        
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

#ifdef ACC
!$acc end kernels
#endif
        call timer_stop(44)
        return
      endif

      call timer_stop(44)

      return
      end subroutine gronor_cofac1

      end module gronor_cofac1_mod
