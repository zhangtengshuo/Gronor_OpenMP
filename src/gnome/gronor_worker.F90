
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

  implicit none

  external :: gronor_solver_init,gronor_solver_final
  external :: gronor_calculate
  external :: swatch,timer_start,timer_stop

  external :: MPI_Recv,MPI_iRecv,MPI_iSend

  real(kind=8), external :: timer_wall_total, timer_wall

  integer :: ibase,jbase,idet,jdet,nidet,njdet
  integer :: i,j,k,l2,n,iact
  integer (kind=4) :: ireq,ierr,ncount,mpitag,mpidest
  integer (kind=8) :: ibuf(4)
  integer (kind=4) :: status(MPI_STATUS_SIZE)
  real (kind=8) :: rbuf(17)

  logical (kind=4) :: flag

  flush(lfnout)

  if(managers.gt.0) then
    mstr=map2(me+1,9)
  endif
  
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

  call gronor_solver_init()

  do i=1,17
    rbuf(i)=0.0d0
  enddo
  
  if(idbg.gt.0) then
    call swatch(date,time)
    write(lfndbg,'(a,1x,a,1x,a,5i5)') date(1:8),time(1:8), &
        ' iamhead, numdev, master, mygroup =',iamhead,numdev,mstr,mygroup
    call swatch(date,time)
    write(lfndbg,130) date(1:8),time(1:8),' thisgroup=',(thisgroup(i),i=1,mgr+1)
130 format(a,1x,a,1x,a,t30,11i5,/,(t35,10i5))
    flush(lfndbg)
  endif

  !     If head thread signal master thread to start sending tasks

  if(iamhead.eq.1) then
    ncount=17
    mpitag=1
    call MPI_iSend(rbuf,ncount,MPI_REAL8,mstr,mpitag,MPI_COMM_WORLD,ireq,ierr)
    if(idbg.gt.20) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),' Head signalled master'
      flush(lfndbg)
    endif
    if(idbg.gt.10) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,i5,a,4i7)') date(1:8),time(1:8),me,' sent buffer   ',mstr
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
      call MPI_Recv(ibuf,ncount,MPI_INTEGER8,mstr,mpitag,MPI_COMM_WORLD,status,ierr)

      if(idbg.gt.10) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,i5,a,7i7)') date(1:8),time(1:8), &
            me,' received task ',mstr,mpitag,(ibuf(i),i=1,4),ierr
        flush(lfndbg)
      endif

      !     Send task to other worker threads in the same group as current head thread

      if(mgr.gt.1) then
        
        do i=1,mgr-1
          ncount=4
          mpidest=thisgroup(i+2)
          mpitag=15
          call MPI_iSend(ibuf,ncount,MPI_INTEGER8,mpidest,mpitag,MPI_COMM_WORLD,ireq,ierr)
          if(idbg.gt.10) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8), &
                me,' sent task to group rank ',thisgroup(i+2)
            flush(lfndbg)
          endif

        enddo
      endif

    else

      !     Receive task from head thread

      if(idbg.gt.30) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8), &
            me,' waiting for task from head rank ',thisgroup(2)
        flush(lfndbg)
      endif
      ncount=4
      mpidest=thisgroup(2)
      mpitag=15
      call MPI_Recv(ibuf,ncount,MPI_INTEGER8,mpidest,mpitag,MPI_COMM_WORLD,status,ierr)
      if(idbg.gt.10) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8), &
            me,' received task from head rank ',thisgroup(2)
        flush(lfndbg)
      endif

    endif

!     Generate the ME list for ibase=ibuf(1) and jbase=ibuf(2)

        if((icur.ne.iabs(ibuf(1)).or.jcur.ne.ibuf(2)).and.ibuf(2).gt.0) then
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
!          if(otreq) call MPI_Cancel(itreq,ierr)
          call timer_stop(39)
          return
        endif
        odupl=ibuf(1).lt.0
        ibuf(1)=iabs(ibuf(1))
        if(odupl) then
          ncount=4
          mpitag=99
          if(.not.otreq) then
            call MPI_iRecv(irbuf,ncount,MPI_INTEGER8,MPI_ANY_SOURCE, &
                mpitag,MPI_COMM_WORLD,itreq,ierr)
            if(idbg.gt.10) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,a)') &
                  date(1:8),time(1:8),' Terminate iRecv posted '
            endif
            otreq=.true.
          endif
          call MPI_Test(itreq,flag,status,ierr)
          if(flag) then
!            call MPI_Cancel(itreq,ierr)
            if(idbg.gt.10) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,a)') date(1:8),time(1:8), &
                  ' Terminating in gronor_worker'
            endif
            call timer_stop(39)
            oterm=.true.
            return
          endif
        endif

        call timer_stop(39)

        if(idbg.gt.15) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,a,f12.6)') date(1:8),time(1:8), &
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
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,i5,a,6i10)') date(1:8),time(1:8), &
                me,' Entering gronor_calculate with ',ibase,jbase,idet,jdet,ntask,nbatch
            flush(lfndbg)
          endif
          call timer_start(47)
          call gronor_calculate(ibase,jbase,idet,jdet)
          call timer_stop(47)

          if(oterm) return

          buffer(3)=timer_wall(47)
          if(idbg.gt.30) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,i5,a)') date(1:8),time(1:8), &
                me,' Returned from gronor_calculate '
            flush(lfndbg)
          endif
          if(idbg.ge.12) then
            write(lfndbg,*)'Multipoles after multiplying the coeffs',(buffer(i),i=9,17)          
          endif
          call timer_start(48)
          if(iamhead.eq.1) then
!     call MPI_Wait(ireq,status,ierr)
            do i=1,17
              rbuf(i)=buffer(i)
            enddo
            ncount=17
            mpitag=1
            call MPI_iSend(rbuf,ncount,MPI_REAL8,mstr,mpitag,MPI_COMM_WORLD,ireq,ierr)
            if(idbg.gt.10) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,i5,a,7i7)') date(1:8),time(1:8), &
                  me,' sent results  ',mstr,(ibuf(i),i=1,4)
              flush(lfndbg)
            endif
          endif
          call timer_stop(48)
        endif
        call timer_stop(46)

      enddo

      call gronor_solver_final()
      
      return
    end subroutine gronor_worker
