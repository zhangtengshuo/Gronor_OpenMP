
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

!>    Driver routine for worker ranks
!!    @brief Driver for calculation Hamiltonian matrix elements on worker ranks
!!    @author T. P. Straatsma (ORNL)


      subroutine gronor_worker()

      use mpi
      use cidef
      use cidist
      use gnome_integrals
      use gnome_data
      use gnome_parameters
#ifdef CUSOLVER
      use cusolverDn
      use cuda_cusolver
      use cudafor
#endif
#ifdef MKL
      use mkl_solver
#endif
      implicit none

      real(kind=8), external :: timer_wall_total, timer_wall

      integer :: ibase,jbase,idet,jdet,nidet,njdet
      integer :: i,j,k,l2,n,iact
      integer (kind=4) :: ireq,ierr,ncount,mpitag,mpidest
      integer (kind=8) :: ibuf(4)
      integer (kind=4) :: status(MPI_STATUS_SIZE)
      real (kind=8) :: rbuf(17)

      logical (kind=4) :: flag

      flush(lfnout)

      if(iamacc.eq.0.and.ntask.eq.0) return

      otreq=.false.

      l2=0
      mnact=0
      mvec=0
      do ibase=1,nbase
       do jbase=1,ibase
        nidet=idetb(ibase)
        njdet=idetb(jbase)
        if(ibase.eq.jbase) then
         l2=max(l2,nidet*(nidet+1)/2)
        else
         l2=max(l2,nidet*njdet)
        endif
       enddo
       mnact=max(mnact,nactb(ibase))
       mvec=max(mvec,nactb(ibase)+inactb(ibase))
      enddo

      nelecs=0
      nveca=0
      n=0
      do ibase=1,nbase
       nveca=max(nveca,inactb(ibase)+nactb(ibase))
       n=2*inactb(ibase)
       do iact=1,nactb(ibase)
         n=n+iabs(int(iocc(1,ibase,iact),kind=kind(n)))
       enddo
       nelecs=max(nelecs,n)
      enddo

      nvecb=nveca
      nstdim=max(1,nelecs*nelecs,nbas*(nbas+1)/2)
      mbasel=max(nelecs,nbas)

      ibase0=0
      jbase0=0
      idet0=0
      jdet0=0

      if(iamacc.ne.0) then
#ifdef CUSOLVER

        ndim=nelecs
        mdim=mbasel
        lwork1=0
        lwork2=0

#ifdef CUSOLVERJ
        if(isolver.eq.1) then
#endif
#ifdef ACC
!$acc data copyin(w,ta) create(dev_info_d)
!$acc host_data use_device(ta)
#endif
          cusolver_status = cusolverDnDgesvd_bufferSize                         &
     &         (cusolver_handle,ndim,ndim,lwork1)
#ifdef ACC
!$acc end host_data
#endif
          if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)                       &
     &         write(*,*) 'cusolverDnDgesvd_bufferSize failed'

#ifdef ACC
!$acc end data
#endif
#ifdef CUSOLVERJ
        elseif(isolver.eq.2) then

          ndim=nelecs
          mdim=mbasel
          tol = tolsvj
          max_sweeps = iswsvj

          jobz=CUSOLVER_EIG_MODE_VECTOR

#ifdef ACC
!$acc data create(dev_info_d,u,w,ev,ta)
#endif
          cusolver_status = cudaStreamCreateWithFlags                           &
     &         (stream,cudaStreamNonBlocking)

          cusolver_status = cusolverDnSetStream(cusolver_handle,stream)

          cusolver_status = cusolverDnCreateGesvdjInfo(gesvdj_params)

          cusolver_status = cusolverDnXgesvdjSetTolerance                       &
     &         (gesvdj_params,tol)

          cusolver_status = cusolverDnXgesvdjSetMaxSweeps                       &
     &         (gesvdj_params,max_sweeps)

#ifdef ACC
!$acc host_data use_device(ta,ev,u,w)
#endif
          cusolver_status = cusolverDnDgesvdj_bufferSize                        &
     &         (cusolver_handle,jobz,econ,                                      &
     &         ndim,ndim,ta,mdim,ev,u,ndim,w,ndim,lwork1,                       &
     &         gesvdj_params)

#ifdef ACC
!$acc end host_data
#endif
          if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)                       &
     &         print *,"cusolverDnDgesvdj_bufferSize failed",                   &
     &         cusolver_status

#ifdef ACC
!$acc end data
#endif
        endif
#endif

#ifdef CUSOLVERJ
        if(jsolver.eq.1) then
#endif
#ifdef ACC
!$acc data copyin(w,ta) create(dev_info_d)
!$acc host_data use_device(ta,w)
#endif
          cusolver_status = cusolverDnDsyevd_bufferSize                         &
     &         (cusolver_handle,CUSOLVER_EIG_MODE_NOVECTOR,                     &
     &         CUBLAS_FILL_MODE_LOWER,                                          &
     &         ndim,ta,mdim,w,lwork2)
#ifdef ACC
!$acc end host_data
#endif
          if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)                       &
     &         print *,"cusolverDnDsyevd_bufferSize failed",                    &
     &         cusolver_status

#ifdef ACC
!$acc end data
#endif
#ifdef CUSOLVERJ
        elseif(jsolver.eq.2) then

! Jacobi EVD

          ndim=nelecs
          mdim=mbasel
          tol = tolevj
          max_sweeps = iswevj

          jobz = CUSOLVER_EIG_MODE_NOVECTOR
          uplo = CUBLAS_FILL_MODE_LOWER

#ifdef ACC
!$acc data create(dev_info_d,w,ta,ev)
#endif
          cusolver_status = cudaStreamCreateWithFlags                           &
     &         (stream,cudaStreamNonBlocking)

          cusolver_status = cusolverDnSetStream(cusolver_handle,stream)

          cusolver_status = cusolverDnCreateSyevjInfo(syevj_params)

          cusolver_status = cusolverDnXsyevjSetTolerance                        &
     &         (syevj_params,tol)

          cusolver_status = cusolverDnXsyevjSetMaxSweeps                        &
     &         (syevj_params,max_sweeps)

#ifdef ACC
!$acc host_data use_device(ta,w)
#endif
          cusolver_status = cusolverDnDsyevj_bufferSize                         &
     &         (cusolver_handle,jobz,uplo,                                      &
     &         ndim,ta,mdim,ev,lwork2,                                          &
     &         syevj_params)

#ifdef ACC
!$acc end host_data
#endif
          if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)                       &
     &         print *,"cusolverDnDsyevj_bufferSize failed",                    &
     &         cusolver_status

#ifdef ACC
!$acc end data
#endif
        endif
#endif

        lwork1=max(lwork1,lwork2)

        allocate(workspace_d(lwork1))
        allocate(rwork(nelecs))

#endif
      else
#ifdef MKL
        ndimm=nelecs
        mdimm=mbasel
        lwork1m=-1
        lwork2m=-1
        lworki=-1
        if(isolver.eq.1) then
          allocate(workspace_d(2448))
          call dgesvd('All','All',ndimm,ndimm,a,ndimm,                          &
     &         ev,u,ndimm,w,ndimm,                                              &
     &         workspace_d,lwork1m,ierr)
          lwork1m=int(workspace_d(1))
          deallocate(workspace_d)
        endif
        if(jsolver.eq.1) then
          allocate(workspace_d(2448))
          allocate(workspace_i(2448))
          call dsyevd('N','L',ndimm,a,ndimm,w,workspace_d,lwork2m,              &
     &         workspace_i,lworki,ierr)
          lwork2m=int(workspace_d(1))
          lworki=int(workspace_i(1))
          deallocate(workspace_d)
          deallocate(workspace_i)
        endif
        lwork1m=max(8,lwork1m,lwork2m)
        lworki=max(8,lworki)
        allocate(workspace_d(lwork1m))
        allocate(workspace_i(lworki))
#endif
      endif

      do i=1,17
       rbuf(i)=0.0d0
      enddo

      if(idbg.gt.0) then
       call swatch(date,time)
       write(lfndbg,'(a,1x,a,1x,a,5i5)') date(1:8),time(1:8),                   &
     & ' iamhead, numdev, master, mygroup =',                                   &
     & iamhead,numdev,master,mygroup
       call swatch(date,time)
       write(lfndbg,130) date(1:8),time(1:8),                                   &
     & ' thisgroup=',(thisgroup(i),i=1,mgr+1)
 130   format(a,1x,a,1x,a,t30,11i5,/,(t35,10i5))
       flush(lfndbg)
      endif

!     If head thread signal master thread to start sending tasks

      if(iamhead.eq.1) then
        ncount=17
        mpitag=1
       call MPI_iSend(rbuf,ncount,MPI_REAL8,master,mpitag,                      &
     & MPI_COMM_WORLD,ireq,ierr)
       if(idbg.gt.20) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),                      &
     &      ' Head signalled master'
        flush(lfndbg)
       endif
       if(idbg.gt.10) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,i5,a,4i7)') date(1:8),time(1:8),                  &
     &     me,' sent buffer   ',master
        flush(lfndbg)
       endif
      endif

      ibase=1

      do while(ibase.gt.0)

        call timer_start(39)

        if(iamhead.eq.1) then

!     Receive next task from master on head thread
          ncount=4
          mpitag=2
          call MPI_Recv(ibuf,ncount,MPI_INTEGER8,master,mpitag,                 &
     &         MPI_COMM_WORLD,status,ierr)

          if(idbg.gt.10) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,i5,a,7i7)') date(1:8),time(1:8),              &
     &           me,' received task ',master,mpitag,(ibuf(i),i=1,4),ierr
            flush(lfndbg)
          endif

!     Send task to other worker threads in the same group as current head thread

          if(mgr.gt.1) then

            do i=1,mgr-1
              ncount=4
              mpidest=thisgroup(i+2)
              mpitag=15
              call MPI_iSend(ibuf,ncount,MPI_INTEGER8,mpidest,mpitag,           &
     &             MPI_COMM_WORLD,ireq,ierr)
              if(idbg.gt.10) then
                call swatch(date,time)
                write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),          &
     &               me,' sent task to group rank ',thisgroup(i+2)
                flush(lfndbg)
              endif

            enddo
          endif

        else

!     Receive task from head thread

          if(idbg.gt.30) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),              &
     &           me,' waiting for task from head rank ',thisgroup(2)
            flush(lfndbg)
          endif
          ncount=4
          mpidest=thisgroup(2)
          mpitag=15
          call MPI_Recv(ibuf,ncount,MPI_INTEGER8,mpidest,mpitag,                &
     &         MPI_COMM_WORLD,status,ierr)
          if(idbg.gt.10) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),              &
     &           me,' received task from head rank ',thisgroup(2)
            flush(lfndbg)
          endif

        endif

!     Generate the ME list for ibase=ibuf(1) and jbase=ibuf(2)

        if((icur.ne.iabs(ibuf(1)).or.jcur.ne.ibuf(2)).and.                       &
     &       ibuf(2).gt.0) then
          if(icur.eq.0.and.jcur.eq.0) allocate(melist(memax,2))
          icur=iabs(ibuf(1))
          jcur=ibuf(2)
          if(icur.gt.0.and.jcur.gt.0) then
            ndeti=idetb(icur)
            ndetj=idetb(jcur)
            k=0
            if(icur.eq.jcur) then
              do i=1,ndeti
                do j=i,ndeti
                  if(dabs(civb(i,icur)*civb(j,jcur)).lt.tau_CI) exit
                  k=k+1
                  melist(k,1)=i
                  melist(k,2)=j
                enddo
              enddo
            else
              do i=1,ndeti
                do j=1,ndetj
                  if(dabs(civb(i,icur)*civb(j,jcur)).lt.tau_CI) exit
                  k=k+1
                  melist(k,1)=i
                  melist(k,2)=j
                enddo
              enddo
            endif
          endif
        endif

!     Check if this is a duplicate

        oterm=ibuf(2).lt.0.or.ibuf(3).lt.0.or.ibuf(4).lt.0
        if(oterm) then
          if(otreq) call MPI_Cancel(itreq,ierr)
          call timer_stop(39)
          return
        endif
        odupl=ibuf(1).lt.0
        ibuf(1)=iabs(ibuf(1))
        if(odupl) then
          ncount=4
          mpitag=99
          if(.not.otreq) then
            call MPI_iRecv(irbuf,ncount,MPI_INTEGER8,MPI_ANY_SOURCE,            &
     &           mpitag,MPI_COMM_WORLD,itreq,ierr)
            if(idbg.gt.10) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,a)')                                        &
     &             date(1:8),time(1:8),' Terminate iRecv posted '
            endif
            otreq=.true.
          endif
          call MPI_Test(itreq,flag,status,ierr)
          if(flag) then
            call MPI_Cancel(itreq,ierr)
            if(idbg.gt.10) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,a)') date(1:8),time(1:8),                   &
     &             ' Terminating in gronor_worker'
            endif
            call timer_stop(39)
            oterm=.true.
            return
          endif
        endif

        call timer_stop(39)

        if(idbg.gt.15) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,a,f12.6)') date(1:8),time(1:8),                 &
     &         ' Cumulative COMM1 Wait Time ',timer_wall_total(39)
          flush(lfndbg)
        endif

        ibase=ibuf(1)
        jbase=ibuf(2)
        idet=ibuf(3)
        jdet=ibuf(4)

        call timer_start(46)
        if(ibase.ne.0.and..not.oterm) then
          if(idbg.gt.30) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,i5,a,6i10)') date(1:8),time(1:8),             &
     &           me,' Entering gronor_calculate with ',                         &
     &           ibase,jbase,idet,jdet,ntask,nbatch
            flush(lfndbg)
          endif
          call timer_start(47)
          call gronor_calculate(ibase,jbase,idet,jdet)
          call timer_stop(47)

          if(oterm) return

          buffer(3)=timer_wall(47)
          if(idbg.gt.30) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,i5,a)') date(1:8),time(1:8),                  &
     &           me,' Returned from gronor_calculate '
            flush(lfndbg)
          endif
          if(idbg.ge.12) then
            write(lfndbg,*)'Multipoles after multiplying the coeffs',           &
     &           (buffer(i),i=9,17)          
          endif
          call timer_start(48)
          if(iamhead.eq.1) then
!     call MPI_Wait(ireq,status,ierr)
            do i=1,17
              rbuf(i)=buffer(i)
            enddo
            ncount=17
            mpitag=1
            call MPI_iSend(rbuf,ncount,MPI_REAL8,master,mpitag,                 &
     &           MPI_COMM_WORLD,ireq,ierr)
            if(idbg.gt.10) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,i5,a,7i7)') date(1:8),time(1:8),            &
     &             me,' sent results  ',master,(ibuf(i),i=1,4)
              flush(lfndbg)
            endif
          endif
          call timer_stop(48)
        endif
        call timer_stop(46)

      enddo

      if(iamacc.gt.0) then
#ifdef CUSOLVER
        cusolver_status = cusolverDnDestroy(cusolver_handle)
        if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)                         &
     &        write(*,*) 'cusolver_handle destruction failed'
#endif
      endif

      return
      end
