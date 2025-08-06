
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
  use gnome_solvers
  use omp_lib
  use gronor_calculate_mod, only: gronor_calculate

  implicit none

  external :: gronor_solver_init,gronor_solver_final
  external :: gronor_abort
  external :: swatch,timer_start,timer_stop

!  external :: MPI_Recv,MPI_iRecv,MPI_iSend

  real(kind=8), external :: timer_wall_total, timer_wall

  integer :: ibase,jbase,idet,jdet,nidet,njdet
  integer :: i,j,k,l2,n,iact
  integer (kind=4) :: ireq, ierr, ncount, mpitag, mpidest
  integer (kind=8) :: ibuf(4)
  integer (kind=4) :: status(MPI_STATUS_SIZE)
  real (kind=8) :: tbuf(18)
  integer :: thread_id, lfnmpi
  character(len=128) :: mpifile
  integer (kind=4) :: mpi_err_len, ierr2
  character(len=MPI_MAX_ERROR_STRING) :: mpi_err_str
  character(len=10) :: today, now

  real (kind=8), allocatable :: va(:,:),vb(:,:),tb(:,:),ta(:,:),a(:,:)
  real (kind=8), allocatable :: u(:,:),w(:,:),wt(:,:),ev(:)
  real (kind=8), allocatable :: sdiag(:),diag(:),bsdiag(:),bdiag(:)
  real (kind=8), allocatable :: csdiag(:),cdiag(:)
  real (kind=8), allocatable :: w1(:),w2(:,:)
  real (kind=8), allocatable :: taa(:,:),sm(:,:),aaa(:,:),aat(:,:),tt(:,:)

  logical (kind=4) :: flag

  ! Manager layer removed; master rank stored globally in mstr
  
  if(ntask.eq.0) return

  oterm = .false.
  otreq = .false.
  odupl = .false.
  itreq = 0
  irbuf = 0_8

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
  nveca_max=0
  n=0
  do ibase=1,nbase
    nveca_max=max(nveca_max,inactb(ibase)+nactb(ibase))
    n=2*inactb(ibase)
    do iact=1,nactb(ibase)
      n=n+iabs(int(iocc(1,ibase,iact),kind=kind(n)))
    enddo
    nelecs=max(nelecs,n)
  enddo
  nstdim=max(1,nelecs*nelecs,nbas*(nbas+1)/2)
  mbasel=max(nelecs,nbas)

  icur=0
  jcur=0
  ndeti=0
  ndetj=0
  nacti=0
  nactj=0
  inacti=0
  inactj=0
  lsvcpu=.false.
  levcpu=.false.
  lsvtrns=.false.

#ifdef _OPENMP
  call omp_set_num_threads(num_threads)

!$omp parallel private(thread_id,va,vb,tb,ta,a,u,w,wt,ev,w1,w2,taa,sm,aaa,aat,tt,sdiag,diag,bsdiag,bdiag,csdiag,cdiag, &
!$omp& ibase,jbase,idet,jdet,nidet,njdet,i,j,k,l2,n,iact,ibuf,status,tbuf,lfnmpi,mpifile,ireq,ierr,ncount,mpitag,mpidest, &
!$omp& mpi_err_len,ierr2,mpi_err_str,today,now,flag) &
!$omp& copyin(oterm,otreq,odupl,itreq,irbuf,icur,jcur,lsvcpu,levcpu,lsvtrns, &
!$omp&        ndeti,ndetj,nacti,nactj,inacti,inactj,nelecs,nveca_max,nstdim,mbasel, &
!$omp&        ntcl,ntop,nclose,nopen,nelec,nact,ninact)

  thread_id = omp_get_thread_num()

  allocate(a(nelecs,nelecs))
  allocate(ta(mbasel,max(mbasel,nveca_max)))
  allocate(tb(mbasel,nveca_max))
  allocate(va(nveca_max,mbasel))
  allocate(vb(nveca_max,mbasel))
  allocate(w1(max(nelecs,nbas,mbasel)))
  allocate(w2(max(nelecs,nbas,mbasel),max(nelecs,nbas,mbasel)))
  allocate(taa(mbasel,max(mbasel,nveca_max)))

  allocate(u(nelecs,nelecs))
  allocate(w(nelecs,nelecs))
  allocate(wt(nelecs,nelecs))
  allocate(ev(nelecs))
  allocate(diag(max(nelecs,nbas,mbasel)))
  allocate(bdiag(max(nelecs,nbas,mbasel)))
  allocate(cdiag(max(nelecs,nbas,mbasel)))
  allocate(bsdiag(max(nelecs,nbas,mbasel)))
  allocate(csdiag(max(nelecs,nbas,mbasel)))
  allocate(sdiag(max(nelecs,nbas,mbasel)))
  allocate(aaa(mbasel,max(mbasel,nveca_max)))
  allocate(tt(mbasel,max(mbasel,nveca_max)))
  allocate(aat(mbasel,max(mbasel,nveca_max)))
  allocate(sm(mbasel,max(mbasel,nveca_max)))
  allocate(ioccup(mnact,2))
  allocate(iocopen(mnact,2))  ! thread-local open-shell occupation
  allocate(vec(mvec,mbasel,2))
  allocate(vtemp(mvec,mbasel,2))
  allocate(ioccn(nsrep,2))
  allocate(veca(mbasel))
  allocate(vecb(mbasel))
  allocate(melist(memax,2))

#ifdef ACC
!$acc data create(va,vb,tb,ta,a,u,w,wt,ev,w1,w2,taa,sm,aaa,aat,tt,sdiag,diag,bsdiag,bdiag,csdiag,cdiag)
#endif

  if(idbg.gt.50 .and. thread_id==0) then
    call swatch(today,now)
    write(lfndbg,'(a,1x,a,a)') today(1:8),now(1:8)," Entering solver initialization"
    flush(lfndbg)
  endif

  call gronor_solver_init(nelecs, a, u, w, ev)

  if(idbg.gt.50 .and. thread_id==0) then
    call swatch(today,now)
    write(lfndbg,'(a,1x,a,a)') today(1:8),now(1:8)," Solver initialization completed"
    flush(lfndbg)
  endif

  if(idbg.gt.10 .and. thread_id==0) then
    call swatch(today,now)
    write(lfndbg,'(a,1x,a,a)') today(1:8),now(1:8), " Array dimensions check in gronor_worker_process:"
    ! 输出二维数组维度
    write(lfndbg,'(a,2i10)') " va:    ", size(va,1), size(va,2)
    write(lfndbg,'(a,2i10)') " vb:    ", size(vb,1), size(vb,2)
    write(lfndbg,'(a,2i10)') " tb:    ", size(tb,1), size(tb,2)
    write(lfndbg,'(a,2i10)') " ta:    ", size(ta,1), size(ta,2)
    write(lfndbg,'(a,2i10)') " a:     ", size(a,1), size(a,2)
    write(lfndbg,'(a,2i10)') " u:     ", size(u,1), size(u,2)
    write(lfndbg,'(a,2i10)') " w:     ", size(w,1), size(w,2)
    write(lfndbg,'(a,2i10)') " wt:    ", size(wt,1), size(wt,2)
    write(lfndbg,'(a,2i10)') " sm:    ", size(sm,1), size(sm,2)
    write(lfndbg,'(a,2i10)') " aaa:   ", size(aaa,1), size(aaa,2)
    write(lfndbg,'(a,2i10)') " aat:   ", size(aat,1), size(aat,2)
    write(lfndbg,'(a,2i10)') " tt:    ", size(tt,1), size(tt,2)
    ! 输出一维数组维度
    write(lfndbg,'(a,i10)')  " ev:    ", size(ev)
    write(lfndbg,'(a,i10)')  " w1:    ", size(w1)
    write(lfndbg,'(a,2i10)') " w2:    ", size(w2,1), size(w2,2)  ! 注意w2是二维
    write(lfndbg,'(a,i10)')  " sdiag: ", size(sdiag)
    write(lfndbg,'(a,i10)')  " diag:  ", size(diag)
    write(lfndbg,'(a,i10)')  " bsdiag:", size(bsdiag)
    write(lfndbg,'(a,i10)')  " bdiag: ", size(bdiag)
    write(lfndbg,'(a,i10)')  " csdiag:", size(csdiag)
    write(lfndbg,'(a,i10)')  " cdiag: ", size(cdiag)
    write(lfndbg,'(a,2i10)') " taa:   ", size(taa,1), size(taa,2)
    write(lfndbg,'(a,i10)')  " nelecs:", nelecs
    write(lfndbg,'(a,i10)')  " nveca_max: ", nveca_max
    write(lfndbg,'(a,i10)')  " mbasel:", mbasel
    write(lfndbg,'(a,i10)')  " nstdim:", nstdim
    flush(lfndbg)
  endif


  call gronor_worker_process(va,vb,tb,ta,a,u,w,wt,ev,w1,w2,taa,sm,aaa,aat,tt,sdiag,diag,bsdiag,bdiag,csdiag,cdiag)

  call gronor_solver_finalize()

  deallocate(a)
  deallocate(ta)
  deallocate(tb)
  deallocate(va)
  deallocate(vb)
  deallocate(w1)
  deallocate(w2)
  deallocate(taa)
  deallocate(u)
  deallocate(w)
  deallocate(wt)
  deallocate(ev)
  deallocate(diag)
  deallocate(bdiag)
  deallocate(cdiag)
  deallocate(bsdiag)
  deallocate(csdiag)
  deallocate(sdiag)
  deallocate(aaa)
  deallocate(tt)
  deallocate(aat)
  deallocate(sm)
  deallocate(ioccup)
  deallocate(iocopen)          ! thread-local open-shell occupation
  deallocate(vec)
  deallocate(vtemp)
  deallocate(veca)
  deallocate(vecb)
  deallocate(ioccn)
  if(allocated(melist)) deallocate(melist)

#ifdef ACC
!$acc end data
#endif

!$omp end parallel
#endif

  call gronor_update_device_info()

!  if(otreq) then
!    call MPI_Test(itreq,flag,status,ierr)
!    if(.not.flag) call MPI_Request_free(itreq,ierr)
!  endif

  return

contains

  subroutine gronor_worker_process(va,vb,tb,ta,a,u,w,wt,ev,w1,w2,taa,sm,aaa,aat,tt,sdiag,diag,bsdiag,bdiag,csdiag,cdiag)

  use mpi
  use cidef
  use cidist
  use gnome_integrals
  use gnome_data
  use gnome_parameters
  use gnome_solvers
  use omp_lib

  implicit none

  real (kind=8), intent(inout) :: va(:,:),vb(:,:),tb(:,:),ta(:,:),a(:,:)
  real (kind=8), intent(inout) :: u(:,:),w(:,:),wt(:,:),ev(:)
  real (kind=8), intent(inout) :: w1(:),w2(:,:),taa(:,:),sm(:,:),aaa(:,:),aat(:,:),tt(:,:)
  real (kind=8), intent(inout) :: sdiag(:),diag(:),bsdiag(:),bdiag(:),csdiag(:),cdiag(:)

!  external :: MPI_Recv,MPI_iRecv,MPI_iSend

  thread_id = 0
  lfnmpi    = 0
  mpifile   = ' '
  ireq      = 0
  ierr      = 0
  ncount    = 0
  mpitag    = 0
  mpidest   = 0
  ibuf      = 0_8
  status    = 0
  tbuf      = 0.0d0
  flag      = .false.

  thread_id = omp_get_thread_num()
  if(idbg.gt.0) then
    write(lfndbg,'(a,i0,a,i0,a,i0,a,i0,a,i0)') 'thread_id=',thread_id,&
         ' len_work_dbl=',len_work_dbl,' len_work_int=',len_work_int,&
         ' me=',me,' mstr=',mstr
    flush(lfndbg)
  endif

  if(thread_id.lt.0 .or. len_work_dbl.lt.0_8 .or. len_work_int.lt.0_8 .or.&
     me.lt.0 .or. mstr.lt.0) then
    write(*,'(a,5(1x,i0))') 'Error: invalid worker parameters', thread_id,&
         len_work_dbl,len_work_int,me,mstr
    call gronor_abort(910,'Invalid worker parameters')
  endif
  do i=1,18
    tbuf(i)=0.0d0
  enddo
  tbuf(16)=dble(len_work_dbl)
  tbuf(17)=dble(len_work_int)
  tbuf(18)=dble(thread_id)

  write(mpifile,'("mpi_log_rank",i0,"_thread",i0,".log")') me,thread_id
  open(newunit=lfnmpi,file=mpifile,status='replace',action='write',iostat=ierr)
  write(lfnmpi,'(a,i0,a,i0,a)') 'rank ',me,' thread ',thread_id,' starting'
  flush(lfnmpi)
  write(lfnmpi,'(a,2i8)') 'len_work_dbl len_work_int ',len_work_dbl,len_work_int
  flush(lfnmpi)
  write(lfnmpi,'("va=",i0,"x",i0," vb=",i0,"x",i0," tb=",i0,"x",i0,&
  " ta=",i0,"x",i0," a=",i0,"x",i0)') size(va,1),size(va,2),size(vb,1),size(vb,2), &
  size(tb,1),size(tb,2),size(ta,1),size(ta,2),size(a,1),size(a,2)
  flush(lfnmpi)
  
  if(idbg.gt.0) then
    call swatch(today,now)
    write(lfndbg,'(a,1x,a,1x,a,5i5)') today(1:8),now(1:8), &
        ' iamhead, numdev, master, mygroup =',iamhead,numdev,mstr,mygroup
    call swatch(today,now)
    write(lfndbg,130) today(1:8),now(1:8),' thisgroup=',(thisgroup(i),i=1,mgr+1)
130 format(a,1x,a,1x,a,t30,11i5,/,(t35,10i5))
    flush(lfndbg)
  endif

  !     Each OpenMP thread signals the master it is ready to receive tasks

  ncount=18
  mpitag=1
  call MPI_iSend(tbuf,ncount,MPI_REAL8,mstr,mpitag,MPI_COMM_WORLD,ireq,ierr)
  if(ierr .ne. MPI_SUCCESS) then
    call MPI_Error_string(ierr, mpi_err_str, mpi_err_len, ierr2)
    write(lfndbg,'(a)') 'MPI_iSend failed: '//mpi_err_str(1:mpi_err_len)
    flush(lfndbg)
    call MPI_Abort(MPI_COMM_WORLD, ierr, ierr2)
  endif
  call MPI_Request_free(ireq,ierr)
  if(ierr .ne. MPI_SUCCESS) then
    call MPI_Error_string(ierr, mpi_err_str, mpi_err_len, ierr2)
    write(lfndbg,'(a)') 'MPI_Request_free failed: '//mpi_err_str(1:mpi_err_len)
    flush(lfndbg)
    call MPI_Abort(MPI_COMM_WORLD, ierr, ierr2)
  endif
  write(lfnmpi,'("send ready len_work_dbl=",i0," len_work_int=",i0)') int(tbuf(16)),int(tbuf(17))
  flush(lfnmpi)
  if(idbg.gt.20) then
    call swatch(today,now)
    write(lfndbg,'(a,1x,a,1x,a)') today(1:8),now(1:8),' Head signalled master'
    flush(lfndbg)
  endif
  if(idbg.gt.10) then
    call swatch(today,now)
    write(lfndbg,'(a,1x,a,i5,a,4i7)') today(1:8),now(1:8),me,' sent buffer   ',mstr
    flush(lfndbg)
  endif

  ibase=1

  do while(ibase.gt.0)

    call timer_start(39)

    !     Receive next task directly from master
    ncount=4
    mpitag=100+thread_id
    call MPI_Recv(ibuf,ncount,MPI_INTEGER8,mstr,mpitag,MPI_COMM_WORLD,status,ierr)
    if(ierr .ne. MPI_SUCCESS) then
      call MPI_Error_string(ierr, mpi_err_str, mpi_err_len, ierr2)
      write(lfndbg,'(a)') 'MPI_Recv failed: '//mpi_err_str(1:mpi_err_len)
      flush(lfndbg)
      call MPI_Abort(MPI_COMM_WORLD, ierr, ierr2)
    endif
    write(lfnmpi,'("recv task ibuf=",4i12)') ibuf
    flush(lfnmpi)

    if(idbg.gt.10) then
      call swatch(today,now)
      write(lfndbg,'(a,1x,a,i5,a,7i7)') today(1:8),now(1:8), &
          me,' received task ',mstr,mpitag,(ibuf(i),i=1,4),ierr
      flush(lfndbg)
    endif

!     Generate the ME list for ibase=ibuf(1) and jbase=ibuf(2)

    if((icur.ne.iabs(ibuf(1)).or.jcur.ne.ibuf(2)).and.ibuf(2).gt.0) then
      if(.not.allocated(melist)) allocate(melist(memax,2))
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
              if(dabs(civb(i,icur)*civb(j,jcur)).lt.tau_CI_off) exit
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
!          if(otreq) call MPI_Cancel(itreq,ierr)
      call timer_stop(39)
!      if(otreq) then
!        call MPI_Test(itreq,flag,status,ierr)
!        if(.not.flag) call MPI_Request_free(itreq,ierr)
!      endif
      close(lfnmpi)
      return
    endif
    odupl=ibuf(1).lt.0
    ibuf(1)=iabs(ibuf(1))
    if(odupl.and.iint.ne.0) then
      ncount=4
      mpitag=99
      if(.not.otreq) then
        call MPI_iRecv(irbuf,ncount,MPI_INTEGER8,MPI_ANY_SOURCE, &
            mpitag,MPI_COMM_WORLD,itreq,ierr)
        if(ierr .ne. MPI_SUCCESS) then
          call MPI_Error_string(ierr, mpi_err_str, mpi_err_len, ierr2)
          write(lfndbg,'(a)') 'MPI_iRecv failed: '//mpi_err_str(1:mpi_err_len)
          flush(lfndbg)
          call MPI_Abort(MPI_COMM_WORLD, ierr, ierr2)
        endif
        !            call MPI_Request_free(itreq,ierr)
        if(idbg.gt.10) then
          call swatch(today,now)
          write(lfndbg,'(a,1x,a,a)') &
              today(1:8),now(1:8),' Terminate iRecv posted '
        endif
        otreq=.true.
      endif
      call MPI_Test(itreq,flag,status,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        call MPI_Error_string(ierr, mpi_err_str, mpi_err_len, ierr2)
        write(lfndbg,'(a)') 'MPI_Test failed: '//mpi_err_str(1:mpi_err_len)
        flush(lfndbg)
        call MPI_Abort(MPI_COMM_WORLD, ierr, ierr2)
      endif
      if(flag) then
!            call MPI_Cancel(itreq,ierr)
        if(idbg.gt.10) then
          call swatch(today,now)
          write(lfndbg,'(a,1x,a,a)') today(1:8),now(1:8), &
              ' Terminating in gronor_worker'
        endif
        call timer_stop(39)
!        call MPI_Request_free(itreq,ierr)
        oterm=.true.
        close(lfnmpi)
        return
      endif
    endif

    call timer_stop(39)
    
    if(idbg.gt.15) then
      call swatch(today,now)
      write(lfndbg,'(a,1x,a,a,f12.6)') today(1:8),now(1:8), &
          ' Cumulative COMM1 Wait Time ',timer_wall_total(39)
      flush(lfndbg)
    endif
    
    ibase=ibuf(1)
    jbase=ibuf(2)
    idet=ibuf(3)
    jdet=ibuf(4)
    
    call timer_start(46)
    if(ibase.ne.0.and..not.oterm) then
      if(idbg.gt.30) then
        call swatch(today,now)
        write(lfndbg,'(a,1x,a,i5,a,6i10)') today(1:8),now(1:8), &
            me,' Entering gronor_calculate with ',ibase,jbase,idet,jdet,ntask,nbatch
        flush(lfndbg)
      endif
      call timer_start(47)

      call gronor_calculate(ibase,jbase,idet,jdet,va,vb,tb,ta,a,u,w,wt,ev,w1,w2,taa,sm,aaa,aat,tt,sdiag,diag,bsdiag,bdiag,csdiag,cdiag)

      call timer_stop(47)

      if(oterm) then
!        if(otreq) then
!          call MPI_Test(itreq,flag,status,ierr)
!          if(.not.flag) call MPI_Request_free(itreq,ierr)
!        endif
        close(lfnmpi)
        return
      endif
      
      buffer(3)=timer_wall(47)
      if(idbg.gt.30) then
        call swatch(today,now)
        write(lfndbg,'(a,1x,a,i5,a)') today(1:8),now(1:8), &
            me,' Returned from gronor_calculate '
        flush(lfndbg)
      endif
      if(idbg.ge.12) then
        write(lfndbg,*)'Multipoles after multiplying the coeffs',(buffer(i),i=9,17)          
      endif
      call timer_start(48)
      !     Send results back to master
      do i=1,17
        tbuf(i)=buffer(i)
      enddo
      tbuf(18)=dble(thread_id)
      ncount=18
      mpitag=1
      call MPI_iSend(tbuf,ncount,MPI_REAL8,mstr,mpitag,MPI_COMM_WORLD,ireq,ierr)
      if(ierr .ne. MPI_SUCCESS) then
        call MPI_Error_string(ierr, mpi_err_str, mpi_err_len, ierr2)
        write(lfndbg,'(a)') 'MPI_iSend failed: '//mpi_err_str(1:mpi_err_len)
        flush(lfndbg)
        call MPI_Abort(MPI_COMM_WORLD, ierr, ierr2)
      endif
    call MPI_Request_free(ireq,ierr)
    if(ierr .ne. MPI_SUCCESS) then
      call MPI_Error_string(ierr, mpi_err_str, mpi_err_len, ierr2)
      write(lfndbg,'(a)') 'MPI_Request_free failed: '//mpi_err_str(1:mpi_err_len)
      flush(lfndbg)
      call MPI_Abort(MPI_COMM_WORLD, ierr, ierr2)
    endif
      write(lfnmpi,'("send result buffer=",17(1x,e16.8))') (buffer(i),i=1,17)
      flush(lfnmpi)
      if(idbg.gt.10) then
        call swatch(today,now)
        write(lfndbg,'(a,1x,a,i5,a,7i7)') today(1:8),now(1:8), &
            me,' sent results  ',mstr,(ibuf(i),i=1,4)
        flush(lfndbg)
      endif
      call timer_stop(48)
    endif
    call timer_stop(46)
    
  enddo

  close(lfnmpi)

  return
  end subroutine gronor_worker_process

end subroutine gronor_worker
