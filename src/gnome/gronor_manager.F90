!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ basis, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!> @brief
!! Manager rank assignments
!! @author  Tjerk P. Straatsma, ORNL
!! @author  Coen de Graaf, URV
!! @date    2016
!>

subroutine gronor_manager()

  use mpi
  use cidist
  use cidef
  use gnome_parameters

  implicit none

  integer (kind=4) :: ireq,ierr,ncount,mpitag,mpidest,iremote
  integer (kind=8) :: ibuf(4)
  integer (kind=4) :: status(MPI_STATUS_SIZE)
  real (kind=8) :: rbuf(17),tbuf(17),buffer(17)
  integer (kind=8) :: numsnd,numrcv,numdets,idet,jdet,ibase,jbase,numtsk

  integer (kind=8) :: i,j,k,m
  
  ! Signal the master to start sending tasks

  do i=1,17
    rbuf(i)=0.0d0
    tbuf(i)=0.0d0
  enddo

  numtsk=0
  do i=1,np
    if(map2(i,9).eq.me) numtsk=numtsk+1
  enddo
  
  ncount=17
  mpitag=1
  call MPI_iSend(rbuf,ncount,MPI_REAL8,mstr,mpitag,MPI_COMM_WORLD,ireq,ierr)
!  write(*,'(a,4f12.3)') "to mstr",(rbuf(k),k=1,4)
  if(idbg.gt.20) then
    call swatch(date,time)
    write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),' Head signaled master'
    flush(lfndbg)
  endif
  if(idbg.gt.10) then
    call swatch(date,time)
    write(lfndbg,'(a,1x,a,i5,a,4i7)') date(1:8),time(1:8),me,' sent buffer   ',mstr
    flush(lfndbg)
  endif

  ibase=1

  call timer_start(50)
  
  do while(ibase.gt.0)
    
    call timer_start(51)
    ! Receive task from the master
    ncount=4
    mpitag=2
    call MPI_Recv(ibuf,ncount,MPI_INTEGER8,mstr,mpitag,MPI_COMM_WORLD,status,ierr)
!    write(*,'(a,4i5)') "from mstr ",(ibuf(k),k=1,4)
    if(idbg.gt.10) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,i5,a,7i7)') date(1:8),time(1:8), &
          me,' received task ',mstr,mpitag,(ibuf(i),i=1,4),ierr
      flush(lfndbg)
    endif
    call timer_stop(51)

    if(ibuf(1).le.0) exit

    call timer_start(52)
    
    do i=1,17
      tbuf(i)=0.0d0
    enddo

    ibase=ibuf(1)
    jbase=ibuf(2)
    idet=ibuf(3)
    jdet=ibuf(4)

    numdets=jdet-idet+1

    m=numdets/numtsk
    do i=1,numtsk
      mgrbuf(i,5)=m
    enddo
    if(mod(numdets,numtsk).eq.0) then
      mgrbuf(i,5)=mgrbuf(i,5)+1
    endif
    numbuf=numtsk
    if(numdets.le.numtsk) then
      do i=1,numtsk
      mgrbuf(i,5)=0
      enddo
      do i=1,numdets
      mgrbuf(i,5)=1
      enddo
    numbuf=numdets
    endif
    
    ! Split the determinant list into numbuf pieces

!    write(*,'(i5,a,3i8)') me,' Size    ',idet,jdet,numdets
    mgrbuf(1,1)=idet
    do i=1,numbuf-1
      mgrbuf(i,2)=mgrbuf(i,1)+mgrbuf(i,5)-1
      mgrbuf(i+1,1)=mgrbuf(i,2)+1
    enddo
    mgrbuf(numbuf,2)=jdet
    do i=1,numbuf
      mgrbuf(i,3)=0
!      write(*,'(i5,a,3i8)') me,' Subtask ',i,mgrbuf(i,1),mgrbuf(i,2)
    enddo

    ! set number of tasks sent in numsnd, received in numrcv

    numsnd=0
    numrcv=0

    ! for i running over the number of sub-tasks (1..numbuf<maxbuf)
    ! mgrbuf(i,1) : first index for this part of the buffer
    ! mgrbuf(i,2) : last index for this part of the buffer
    ! mgrbuf(i,3) : 0=not sent; 1=sent; 2=received results
    ! mgrbuf(i,4) : worker rank this buffer was sent to

    ! for i running over the number of worker ranks of this manager (1..numwrk)
    ! mgrwrk(i,1) : rank of the worker
    ! mgrwrk(i,2) : 0=worker not ready; 1=worker ready to receive next sub-task

    ! first send sub tasks to workers ready to receive a new sub-task
    
    do i=1,numbuf
      if(mgrbuf(i,3).eq.0.and.mgrbuf(i,5).gt.0) then
        do j=1,numwrk
          if(mgrwrk(j,2).eq.1) then
            mipbuf(1,j)=ibase
            mipbuf(2,j)=jbase
            mipbuf(3,j)=mgrbuf(i,1)
            mipbuf(4,j)=mgrbuf(i,2)
            ncount=4
            mpitag=2
            iremote=mgrwrk(j,1)
            call MPI_iSend(mipbuf(1,j),ncount,MPI_INTEGER8, &
                iremote,mpitag,MPI_COMM_WORLD,ireq,ierr)
!            write(*,'(a,5i5)') "to wrkr1 ",iremote,(mipbuf(k,j),k=1,4)
            mgrbuf(i,3)=1
            mgrwrk(j,2)=0
            numsnd=numsnd+1
            exit
          endif
        enddo
      endif
    enddo

    call timer_stop(52)
    
    call timer_start(53)
    
    do while(numrcv.lt.numbuf)
      ncount=17
      mpitag=1
      call MPI_Recv(buffer,ncount,MPI_REAL8,MPI_ANY_SOURCE,mpitag,MPI_COMM_WORLD,status,ierr)
      iremote=status(MPI_SOURCE)
!      write(*,'(a,i5,4f12.3)') "from wrkr",iremote,(buffer(k),k=1,4)
      do j=1,numwrk
        if(mgrwrk(j,1).eq.iremote) then
          if(mgrwrk(j,2).eq.0) then
            numrcv=numrcv+1
            do i=1,17
              tbuf(i)=tbuf(i)+buffer(i)
            enddo
            mgrwrk(j,2)=1
          endif
          if(mgrwrk(j,2).lt.0) mgrwrk(j,2)=1
          do i=1,numbuf
            if(mgrbuf(i,3).eq.0.and.mgrwrk(j,2).eq.1) then
              mipbuf(1,j)=ibase
              mipbuf(2,j)=jbase
              mipbuf(3,j)=mgrbuf(i,1)
              mipbuf(4,j)=mgrbuf(i,2)
              ncount=4
              mpitag=2
              call MPI_iSend(mipbuf(1,j),ncount,MPI_INTEGER8, &
                  iremote,mpitag,MPI_COMM_WORLD,ireq,ierr)
!              write(*,'(a,5i5)') "to wrkr2 ",iremote,(mipbuf(k,j),k=1,4)
              mgrbuf(i,3)=1
              mgrwrk(j,2)=0
              numsnd=numsnd+1
            endif
          enddo
        endif
      enddo
    enddo
    
    call timer_stop(53)
    
    ! Send results buffer to master
    call timer_start(54)
    ncount=17
    mpitag=1
    call MPI_iSend(tbuf,ncount,MPI_REAL8,mstr,mpitag,MPI_COMM_WORLD,ireq,ierr)
!    write(*,'(a,4f12.3)') "to mstr",(tbuf(k),k=1,4)
    if(idbg.gt.10) then
      call swatch(date,time)
      write(lfndbg,'(a,1x,a,i5,a,7i7)') date(1:8),time(1:8), &
          me,' sent results  ',mstr,(ibuf(i),i=1,4)
      flush(lfndbg)
    endif
    do i=1,17
      rbuf(i)=0.0d0
      tbuf(i)=0.0d0
    enddo
    call timer_stop(54)
    
  enddo
  call timer_stop(50)

  call timer_start(55)
  ibuf(1)=0
  ibuf(2)=0
  ibuf(3)=0
  ibuf(4)=0
  ncount=4
  mpitag=2
  do j=1,numwrk
    iremote=mgrwrk(j,1)
    call MPI_iSend(ibuf,ncount,MPI_INTEGER8, &
        iremote,mpitag,MPI_COMM_WORLD,ireq,ierr)
!    write(*,'(a,4i5)') "to wrkr ",(ibuf(k),k=1,4)
  enddo
  call timer_stop(55)
  
  return
    
end subroutine gronor_manager
