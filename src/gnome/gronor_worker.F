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
c      use cuda_functions
      use cudafor
#endif
      implicit none

      real(kind=8), external :: timer_wall_total, timer_wall
      
      integer :: ibase,jbase,idet,jdet,nidet,njdet
      integer :: i,j,l2,n,k,ict
      integer (kind=4) :: ireq,ierr
      integer :: ibuf(4)
      real (kind=8) :: sint
      integer :: status(MPI_STATUS_SIZE)
      integer :: nel2
      integer :: ipopcnt

      flush(lfnout)

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
       n=n+ipopcnt(ioccb(1,ibase),nactb(ibase))
       nelecs=max(nelecs,n)
      enddo
      nvecb=nveca
      nstdim=max(1,nelecs*nelecs,nbas*(nbas+1)/2)
      mbasel=max(nelecs,nbas)

      allocate(ioccup(mnact,2))
      allocate(vec(mvec,mbasel,2))
      allocate(vtemp(mvec,mbasel,2))
      allocate(itemp(21),ioccn(20,2))

      allocate(st(nstdim))

      allocate(va(nveca,mbasel))
      allocate(vb(nvecb,mbasel))
      allocate(vna(nsmf),vnb(nsmf))
      allocate(veca(mbasel))
      allocate(vecb(mbasel))
      allocate(ta(mbasel,max(mbasel,nveca)))
      allocate(tat(mbasel,max(mbasel,nveca)))
      allocate(taa(mbasel,max(mbasel,nveca)))
      allocate(aaa(mbasel,max(mbasel,nveca)))
      allocate(tb(mbasel,nvecb))
      allocate(s12d(mbasel,mbasel))
      allocate(tt(mbasel,max(mbasel,nveca)))
      allocate(aat(mbasel,max(mbasel,nveca)))
      allocate(sm(mbasel,max(mbasel,nveca)))

      allocate(a(nelecs,nelecs))
      allocate(u(nelecs,nelecs))
      allocate(w(nelecs,nelecs))
      allocate(ev(nelecs))
      allocate(sdiag(max(nelecs,nbas,mbasel)))
      allocate(diag(max(nelecs,nbas,mbasel)))
      allocate(bsdiag(max(nelecs,nbas,mbasel)))
      allocate(bdiag(max(nelecs,nbas,mbasel)))
      allocate(csdiag(max(nelecs,nbas,mbasel)))
      allocate(cdiag(max(nelecs,nbas,mbasel)))

      allocate(w1(max(nelecs,nbas,mbasel)))
      allocate(w2(max(nelecs,nbas,mbasel),max(nelecs,nbas,mbasel)))

      allocate(temp(nelecs,nelecs))

      ibase0=0
      jbase0=0
      idet0=0
      jdet0=0
      worksize=0
      worksize2=0

      if(iamacc.ne.0) then
#ifdef CUSOLVER

        ndim=nelecs
        mdim=mbasel
        
        if(isolver.eq.1) then
!$acc data copyin(w,ta) create(dev_info_d)
       
!$acc host_data use_device(ta)
          cusolver_status = cusolverDnDgesvd_bufferSize
     &         (cusolver_handle,ndim,ndim,worksize)
!$acc end host_data

          if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &         write(*,*) 'cusolverDnDgesvd_bufferSize failed'

!$acc end data   
        elseif(isolver.eq.2) then
          
          ndim=nelecs
          mdim=mbasel
          tol = 1.0d-07
          max_sweeps = 15
          worksize = 0
          
          jobz=CUSOLVER_EIG_MODE_VECTOR
        
!$acc data create(dev_info_d,u,w,ev,ta)
        
          cusolver_status = cudaStreamCreateWithFlags
     &         (stream,cudaStreamNonBlocking)
          
          cusolver_status = cusolverDnSetStream(cusolver_handle,stream)
          
          cusolver_status = cusolverDnCreateGesvdjInfo(gesvdj_params)
          
          cusolver_status = cusolverDnXgesvdjSetTolerance
     &         (gesvdj_params,tol)
          
          cusolver_status = cusolverDnXgesvdjSetMaxSweeps
     &         (gesvdj_params,max_sweeps)

!$acc host_data use_device(ta,ev,u,w)
          
          cusolver_status = cusolverDnDgesvdj_bufferSize
     &         (cusolver_handle,jobz,econ,
     &         ndim,ndim,ta,mdim,ev,u,ndim,w,ndim,worksize,
     &         gesvdj_params)
        
!$acc end host_data
          
          if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &         print *,"cusolverDnDgesvdj_bufferSize failed",
     &         cusolver_status

          worksize=10*worksize
!$acc end data
        endif

        if(jsolver.eq.1.or.jsolver.eq.2) then
!$acc data copyin(w,ta) create(dev_info_d)
       
!$acc host_data use_device(ta,w)
          cusolver_status = cusolverDnDsyevd_bufferSize
     &         (cusolver_handle,CUSOLVER_EIG_MODE_NOVECTOR,
     &         CUBLAS_FILL_MODE_LOWER,
     &         ndim,ta,ndim,w,worksize2)
!$acc end host_data

          if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &         print *,"cusolverDnDsyevd_bufferSize failed",
     &         cusolver_status

!$acc end data
        elseif(jsolver.eq.2) then

          jsolver = 1
! Jacobi EVD
          
        endif
      
        worksize=max(worksize,worksize2)
        
        allocate(workspace_d(worksize))
        allocate(rwork(nelecs))

#endif
       endif

      do i=1,8
       buffer(i)=0.0d0
      enddo

      if(idbg.gt.0) then
       call swatch(date,time)
       write(lfndbg,'(a,1x,a,1x,a,5i5)') date(1:8),time(1:8),
     & ' iamhead, numdev, master, mygroup =',
     & iamhead,numdev,master,mygroup
       call swatch(date,time)
       write(lfndbg,130) date(1:8),time(1:8),
     & ' thisgroup=',(thisgroup(i),i=1,mgr+1)
 130   format(a,1x,a,1x,a,t30,11i5,/,(t35,10i5))
       flush(lfndbg)
      endif

!     If head thread signal master thread to start sending tasks

      if(iamhead.eq.1) then
       call MPI_iSend(buffer,8,MPI_REAL8,master,1,
     & MPI_COMM_WORLD,ireq,ierr)
       if(idbg.gt.0) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &      ' Head signalled master'
        flush(lfndbg)
       endif
       if(idbg.gt.2) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &     me,' sent buffer to',master
        flush(lfndbg)
       endif
      endif

      ibase=-1

      do while(ibase.ne.0)

       call timer_start(39)

       if(iamhead.eq.1) then

!     Receive next task from master on head thread

        call MPI_Recv(ibuf,4,MPI_INTEGER8,master,2,
     &  MPI_COMM_WORLD,status,ierr)
        if(idbg.gt.2) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &      me,' received task from ',master
         flush(lfndbg)
        endif

!     Send task to other worker threads in the same group as current head thread

        if(mgr.gt.1) then

         do i=1,mgr-1

          call MPI_iSend(ibuf,4,MPI_INTEGER8,thisgroup(i+2),15,
     &                  MPI_COMM_WORLD,ireq,ierr)
         if(idbg.gt.2) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &       me,' sent task to group rank ',thisgroup(i+2)
          flush(lfndbg)
         endif

         enddo
        endif

       else

!      Receive task from head thread

        if(idbg.gt.2) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &      me,' waiting for task from head rank ',thisgroup(2)
         flush(lfndbg)
        endif
        call MPI_Recv(ibuf,4,MPI_INTEGER8,thisgroup(2),15,
     &  MPI_COMM_WORLD,status,ierr)
        if(idbg.gt.2) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &      me,' received task from head rank ',thisgroup(2)
         flush(lfndbg)
        endif

       endif

       call timer_stop(39)
       
       if(idbg.gt.0) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,i5,a,f12.6)') date(1:8),time(1:8),
     &        me,' Cumulative COMM1 Wait Time ',timer_wall_total(39)
         flush(lfndbg)
       endif
       
       ibase=ibuf(1)
       jbase=ibuf(2)
       idet=ibuf(3)
       jdet=ibuf(4)

       call timer_start(46)
       if(ibase.ne.0) then
        if(idbg.gt.2) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,i5,a,4i10)') date(1:8),time(1:8),
     &        me,' Entering gronor_calculate with ',
     &        ibase,jbase,idet,jdet,ibatch
         flush(lfndbg)
        endif
        call timer_start(47)
        call gronor_calculate(ibase,jbase,idet,jdet,ibatch)
        call timer_stop(47)
        buffer(3)=timer_wall(47)
        if(idbg.gt.2) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,i5,a)') date(1:8),time(1:8),
     &      me,' Returned from gronor_calculate '
         flush(lfndbg)
        endif
       call timer_start(48)
        if(iamhead.eq.1) then
         call MPI_iSend(buffer,8,MPI_REAL8,master,1,
     &   MPI_COMM_WORLD,ireq,ierr)
         if(idbg.gt.2) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &       me,' sent results to master ',master
          flush(lfndbg)
         endif
        endif
       call timer_stop(48)
       endif
       call timer_stop(46)

      enddo

      deallocate(cdiag,csdiag,bdiag,bsdiag,diag,sdiag,ev,w,u,a)
      deallocate(aat,tt,sm)
      deallocate(s12d,tb,aaa,taa,ta,vecb,veca,vnb,vna,vb,va)
      deallocate(s,st)
      deallocate(temp)

      if(iamacc.gt.0) then
#ifdef CUSOLVER
        cusolver_status = cusolverDnDestroy(cusolver_handle)
        if (cusolver_status /= CUSOLVER_STATUS_SUCCESS)
     &        write(*,*) 'cusolver_handle destruction failed'
#endif
      endif
      
      return
      end
