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
!! Main driver routine for the master rank
!!
!! @author T. P. Straatsma, ORNL      
!! @author C. de Graaf, URV
!! @date 2016
!!
!! The master rank orchestrates the calculation by reponding to task requests received 
!! from worker ranks.
!! The master rank requests matrix element contributions from worker ranks as tasks
!! consisting of batches of contributions with user defined bacth size.
!! When all tasks have been dispatched by the master rank, a duplicate of still outstanding
!! tasks is sent to worker ranks requesting a new task. The master rank will ignore results
!! from a worker rank if these results were already received from a different rank.
!! This is the primary mechanism to achieve fault resiliency.
!! A consequence of this implementation is that when all results have been received
!! and the calculation has normally completed, the main program terminates with an MPI-Abort,
!! rather than an MPI_Finalize. This avoids the job from hanging if a harware fault occurred. 
      
      subroutine gronor_master()

      use mpi
      use cidef
      use cidist
      use gnome_data
      use gnome_parameters

      implicit none

      integer :: ibase,jbase,mdet
      integer :: l2,i,j,k,ndone,nleft,nrb
      integer (kind=4) :: ierr,ireq,iremote
      integer :: igrp, ltemp, last
      logical :: osame

      real (kind=8) :: tsum
      real (kind=8), allocatable :: fday(:,:),pnrb(:,:)

      integer :: ibuf(4),nsingt(5)
      integer, allocatable :: lgroup(:,:), lactive(:)
      integer (kind=4) :: status(MPI_STATUS_SIZE)

      real(kind=8), external :: timer_wall_total

      integer, allocatable :: ntasks(:,:),ndets(:,:,:)

      integer :: ioff
      
      allocate(lgroup(npg,5))
      allocate(lactive(npg))
      allocate(ntasks(nbase,nbase))
      allocate(ndets(nbase,nbase,2))
      allocate(pnrb(nbase,nbase))
      allocate(fday(nbase,nbase))

!     nsing(i,j,k) collects the number of singularities (0:k=1, 1:k=2, 2:k=3 3+:k=4)
!     for contributions to Hamiltonian matrix element H(i,j)
!     nsingt(k) sums the singularities over the entire Hamiltonain matrix
      
      do i=1,5
        nsingt(i)=0
      enddo
      tsum=0.0d0

!     Contributions to H(i,j) are sent as tasks to groups and tracked in lgroup:

!     lgroup(igrp,1) : current i in H(i,j) for current task sent to group igrp
!     lgroup(igrp,2) : current j in H(i,j) for current task sent to group igrp
!     lgroup(igrp,3) : index to first determinant pair for H(i,j) in the current task sent to group igrp
!     lgroup(igrp,4) : index to last determinant pair for H(i,j) in the current task sent to group igrp
!     lgroup(igrp,5) : 0:group is done, 1:group has task, 2:group has duplicate task 
      
      do j=1,5
        do i=1,npg
          lgroup(i,j)=0
        enddo
      enddo

!     Initialize the results matrices to zero

!     hbase(i,j)   collects the Hamiltonain matrix results for basestate combination i j
!     sbase(i,j)   collects the overlap matrix results for basestate combination i j
!     tbase(i,j)   collects the wall clock time spent basestate combination i j
!     ntasks(i,j)  contains the task status for base state combination i j
!     -2 : base pair calculation has completed
!     -1 : base pair calculation has not yet started
!     0 : base pair contributions have been calculated and accumulated
!     >0 : number of base pair contributions calculated
!     nsing(i,j,k) sums the number of singularities for base state combination i j

      do ibase=1,nbase
        do jbase=1,nbase
          if(.not.bpdone(ibase,jbase)) then
            hbase(ibase,jbase)=0.0d0
            sbase(ibase,jbase)=0.0d0
            tbase(ibase,jbase)=0.0d0
            ntasks(ibase,jbase)=-1
            do i=1,5
              nsing(ibase,jbase,i)=0
            enddo
          else
            ntasks(ibase,jbase)=-2
          endif
        enddo
      enddo

!     Calculation of the number of matrix element pairs for each base state combination

!     idetb(i)     has the number of determinants in base state i
!     ndets(i,j,1) has the number of determinant pairs for base state combination i j
!     ndets(i,j,2) has the number of completed determinant pairs for base state combination i j
!     mdet         has the total number of determinant pairs to be evaluated
      
      mdet=0
      do ibase=1,nbase
        do jbase=1,ibase
          ndeti=idetb(ibase)
          ndetj=idetb(jbase)
          if(ibase.eq.jbase) then
            mdet=mdet+idetb(ibase)*(idetb(jbase)+1)/2
          else
            mdet=mdet+idetb(ibase)*idetb(jbase)
          endif
          if(ibase.eq.jbase) then
            ndets(ibase,jbase,1)=idetb(ibase)*(idetb(ibase)+1)/2
          else
            ndets(ibase,jbase,1)=idetb(ibase)*idetb(jbase)
          endif
          ndets(ibase,jbase,2)=0
        enddo
      enddo

!     ndest(i,j,1) is replaced by offset in matrix element pair list for base state combination i j
!     ndets(i,j,2) has the number of completed determinant pairs for base state combination i j

      ioff=0
      do ibase=1,nbase
        do jbase=1,ibase
          ndets(ibase,jbase,1)=ndxdet(ibase,jbase)-ioff
          ioff=ndxdet(ibase,jbase)
          ndets(ibase,jbase,2)=0
        enddo
      enddo
      
      if(ipr.ge.40) write(lfnout,601) mdet,numdet
 601  format(' Number of determinant pairs',t50,i16,/,
     &     ' Number of determinants above CI threshold',t50,i16)
      
      if(ipr.gt.0) then
        write(*,602)
        write(lfnout,602)
      endif
 602  format(//,' Calculation of Hamiltonian matrix elements',//,
     &     '   Element       ------- Determinants -------',
     &     '          -------------- Singularities -------------',
     &     '        ----- Unnormalized -----   ',
     &     '         Normalized ',
     &     '  Wall Clock  Wall Clock',/
     &     '  idet jdet      left     right         pairs',
     &     '          none           1           2            3+',
     &     '        Hamiltonian      Overlap   ',
     &     '        Hamiltonian ',
     &     '        Cost        Time',/)

      write(lfnarx,521) nbase
 521  format('Data ',i10)

!     Loop over all base state combinations

      ioff=0
      do ibase=1,nbase
        do jbase=1,ibase

!     Determine offset and number of determinant pairs for this base pair

          osame=ibase.eq.jbase
          ndeti=idetb(ibase)
          ndetj=idetb(jbase)
          nacti=nactb(ibase)
          nactj=nactb(jbase)
          call iloop_init(ndeti,ndetj,osame)
          if(osame) then
            l2=ndeti*(ndeti+1)/2
            ijend=l2
          else
            l2=ndeti*ndetj
            ijend=l2
          endif
          ijend=ndxdet(ibase,jbase)-ioff
          ioff=ndxdet(ibase,jbase)

!     fday is the progress counter for messages to the dayfile
          
          fday(ibase,jbase)=dble(iday)
          
          if(ipr.ge.40) then
            write(lfnout,603) ibase,jbase
 603        format(/,' Base state pair ',2i5)
            write(lfnout,604) ndeti
 604        format(/,' Left wavefunction has ',i8,
     &           ' determinants:',/)
            do i=1,ndeti
              call unpack(ioccb(i,ibase),iocch,nacti)
              write(lfnout,605) i,civb(i,ibase),
     &             (iocch(k),k=1,nacti)
 605          format(i5,f25.15,35i3)
            enddo
            write(lfnout,606) ndetj
 606        format(//,' Right wavefunction has ',i8,
     &           ' determinants:',/)
            do i=1,ndetj
              call unpack(ioccb(i,jbase),iocch,nactj)
              write(lfnout,605) i,civb(i,jbase),
     &             (iocch(k),k=1,nactj)
            enddo
          endif
          flush(lfnout)
          
          call swatch(date,time)
          call timer_stop(99)
          if(ntasks(ibase,jbase).ne.-2) then
            write(lfnday,802) date(1:8),time(1:8),
     &           timer_wall_total(99),
     &           '  :  ',ibase,jbase,' started, dimensions ',
     &           ndeti,ndetj
          else
            write(lfnday,802) date(1:8),time(1:8),
     &           timer_wall_total(99),
     &           '  :  ',ibase,jbase,
     &           ' read from checkpoint restart'
          endif
 802      format(a8,2x,a8,f12.3,a,2i6,a,2i12)
          flush(lfnday)
          call timer_start(99)
          
          nrb=0
          ndone=0
          ibuf(1)=0
          ibuf(2)=0
          ibuf(3)=0
          ibuf(4)=0

          if(ntasks(ibase,jbase).eq.-2) ndone=ijend

!     Loop over the number of determinant pairs for the current base pair

          call timer_start(95)
          do while(ndone.lt.ijend)

!     Wait for a request for a task from any worker group
            
            call timer_start(94)
            call MPI_Recv(buffer,8,MPI_REAL8,MPI_ANY_SOURCE,1,
     &           MPI_COMM_WORLD,status,ierr)
            call timer_stop(94)

!     Determine rank and group from which the request came
            
            iremote=status(MPI_SOURCE)
            igrp=map2(iremote+1,3)
            
            if(idbg.gt.2) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,i5,a,7i7)')
     &             date(1:8),time(1:8),
     &             master,' received buffer from',iremote,igrp,
     &             (lgroup(igrp,i),i=1,5)
              flush(lfndbg)
            endif

!     If the request group completed a task accumulate into appropriate arrays
!     The very first request comes without results so this is skipped
            
            if(lgroup(igrp,5).ne.0) then

!     Determine the number of determinants completed for the base pair in ndets
!     Determine the percent complete for the base pair in pnrb

              ndets(lgroup(igrp,1),lgroup(igrp,2),2)=
     &             ndets(lgroup(igrp,1),lgroup(igrp,2),2)+
     &             lgroup(igrp,4)-lgroup(igrp,3)+1
              pnrb(lgroup(igrp,1),lgroup(igrp,2))=
     &             dble(100*ndets(lgroup(igrp,1),
     &             lgroup(igrp,2),2))/
     &             ndets(lgroup(igrp,1),lgroup(igrp,2),1)

!     Write entry to progress file

              if(ipro.eq.1.or.ipro.eq.3) then
                call swatch(date,time)
                if(ipro.eq.1) rewind(unit=lfnpro)
                write(lfnpro,680) date(1:8),time(1:8),
     &               iremote,igrp,(lgroup(igrp,j),j=1,4),
     &               iloop(lgroup(igrp,3),1),
     &               iloop(lgroup(igrp,3),2),
     &               lgroup(igrp,4)-lgroup(igrp,3)+1,nrb,
     &               pnrb(lgroup(igrp,1),lgroup(igrp,2)),
     &               buffer(1),buffer(2),(int(buffer(j)),j=4,8)
 680            format(a,1x,a,' Rcvd ',2i6,' : ',2i5,2i10,
     &               3i6,i12,f8.3,'% ',2e16.8,5i8)
                flush(lfnpro)
              elseif(ipro.eq.2.or.ipro.eq.4) then
                call swatch(date,time)
                if(ipro.eq.2) rewind(unit=lfnpro)
                write(lfnpro,681) date(1:8),time(1:8), 
     &               iremote,igrp,(lgroup(igrp,j),j=1,4),
     &               iloop(lgroup(igrp,3),1),
     &               iloop(lgroup(igrp,3),2),
     &               lgroup(igrp,4)-lgroup(igrp,3)+1,nrb,
     &               pnrb(lgroup(igrp,1),lgroup(igrp,2))
 681            format(a,1x,a,' Rcvd ',2i6,' : ',2i5,2i10,
     &               3i6,i12,f8.3,'% ')
                flush(lfnpro)
              endif
              if(me.eq.master.and.
     &             pnrb(lgroup(igrp,1),lgroup(igrp,2)).ge.
     &             fday(lgroup(igrp,1),lgroup(igrp,2))) then
                call timer_stop(99)
                call swatch(date,time)
                write(lfnday,702) date(1:8),time(1:8),
     &               timer_wall_total(99),
     &               '  :  ',lgroup(igrp,1),
     &               lgroup(igrp,2),' at ',
     &               pnrb(lgroup(igrp,1),lgroup(igrp,2)),'%'
 702            format(a8,2x,a8,f12.3,a,2i6,a,f5.1,a)
                flush(lfnday)
                call timer_start(99)
                fday(lgroup(igrp,1),lgroup(igrp,2))=
     &               fday(lgroup(igrp,1),
     &               lgroup(igrp,2))+dble(iday)
              endif
              
!     Count the number of tasks completed by the remote rank in numrecs
              
              numrecs(iremote+1)=numrecs(iremote+1)+1
              
!     Accumulate the returned contributions in buffer
              
              hbase(lgroup(igrp,1),lgroup(igrp,2))=
     &             hbase(lgroup(igrp,1),lgroup(igrp,2))+buffer(1)
              sbase(lgroup(igrp,1),lgroup(igrp,2))=
     &             sbase(lgroup(igrp,1),lgroup(igrp,2))+buffer(2)
              tbase(lgroup(igrp,1),lgroup(igrp,2))=
     &             tbase(lgroup(igrp,1),lgroup(igrp,2))+buffer(3)
              nsing(lgroup(igrp,1),lgroup(igrp,2),1)=
     &             nsing(lgroup(igrp,1),
     &             lgroup(igrp,2),1)+int(buffer(4))
              nsing(lgroup(igrp,1),lgroup(igrp,2),2)=
     &             nsing(lgroup(igrp,1),
     &             lgroup(igrp,2),2)+int(buffer(5))
              nsing(lgroup(igrp,1),lgroup(igrp,2),3)=
     &             nsing(lgroup(igrp,1),
     &             lgroup(igrp,2),3)+int(buffer(6))
              nsing(lgroup(igrp,1),lgroup(igrp,2),4)=
     &             nsing(lgroup(igrp,1),
     &             lgroup(igrp,2),4)+int(buffer(7))
              nsing(lgroup(igrp,1),lgroup(igrp,2),5)=
     &             nsing(lgroup(igrp,1),
     &             lgroup(igrp,2),5)+int(buffer(8))

!     Decrement the tasks still to be done for base pair
              
              ntasks(lgroup(igrp,1),lgroup(igrp,2))=
     &             ntasks(lgroup(igrp,1),lgroup(igrp,2))-1

              if(ntasks(lgroup(igrp,1),lgroup(igrp,2)).eq.0) then
                i=lgroup(igrp,1)
                j=lgroup(igrp,2)
                
                if(i.eq.j) then
                  ltemp=idetb(i)*(idetb(i)+1)/2
                else
                  ltemp=idetb(i)*idetb(j)
                endif
                
                write(lfnarx,522) i,j,idetb(i),idetb(j),
     &               ltemp,
     &               (nsing(i,j,k),k=1,4),
     &               hbase(i,j),sbase(i,j)
                
                call timer_stop(99)
                
                if(ipr.ge.1) then
                  if(i.ne.j) then
                    write(*,607) i,j,idetb(i),idetb(j),
     &                   ltemp,
     &                   (nsing(i,j,k),k=1,4),
     &                   hbase(i,j),sbase(i,j),
     &                   tbase(i,j),timer_wall_total(99)
                    write(lfnout,607) i,j,idetb(i),
     &                   idetb(j),ltemp,
     &                   (nsing(i,j,k),k=1,4),
     &                   hbase(i,j),sbase(i,j),
     &                   tbase(i,j),timer_wall_total(99)
                  else
                    write(*,608) i,j,idetb(i),idetb(j),
     &                   ltemp,
     &                   (nsing(i,j,k),k=1,4),
     &                   hbase(i,j),sbase(i,j),
     &                   hbase(i,j)/sbase(i,j),
     &                   tbase(i,j),timer_wall_total(99)
                    write(lfnout,608) i,j,idetb(i),
     &                   idetb(j),ltemp,
     &                   (nsing(i,j,k),k=1,4),
     &                   hbase(i,j),sbase(i,j),
     &                   hbase(i,j)/sbase(i,j),
     &                   tbase(i,j),timer_wall_total(99)
                  endif
                  flush(lfnout)
 522              format(1x,2i5,2i10,2i14,2i12,i14,2e20.13)
 608              format(1x,2i5,2i10,2i14,2i12,i14,f20.10,
     &                 f15.10,f20.10,2f12.1)
 607              format(1x,2i5,2i10,2i14,2i12,i14,f20.10,
     &                 f15.10,20x,2f12.1)
                endif

!     Write entry in the checkpoint restart file
                
                write(lfncpr) i,j,
     &               hbase(i,j),sbase(i,j),tbase(i,j),
     &               (nsing(i,j,k),k=1,4)
                flush(lfncpr)

!     Write completion message to dayfile

                call swatch(date,time)
                if(i.ne.j) then
                  write(lfnday,703) date(1:8),time(1:8),
     &                 timer_wall_total(99),
     &                 '  :  ',i,j,' completed, H,S   :',
     &                 hbase(i,j),sbase(i,j)
                else
                  write(lfnday,703) date(1:8),time(1:8),
     &                 timer_wall_total(99),
     &                 '  :  ',i,j,' completed, H,S,H*:',
     &                 hbase(i,j),sbase(i,j),
     &                 hbase(i,j)/sbase(i,j)
                endif
 703            format(a8,2x,a8,f12.3,a,2i6,a,3f20.10)
                flush(lfnday)
                call timer_start(99)

!     Setting ntasks to -2 identifies this base pair as completed and reported
                
                ntasks(i,j)=-2
              endif
              
!     Reset the task information for the group from which request came
              
              do k=1,5
                lgroup(igrp,k)=0
              enddo

!     Loop over all base pairs to check if all tasks have been completed
!     This will be the case if ntasks(i,j) is zero

            endif

!     Setup the next task to be sent to requesting rank in an integer buffer with length 4
!     buffer(1) : ibase of the base pair
!     buffer(2) : jbase of the base pair
!     buffer(3) : first index into the determinant pair list
!     buffer(4) : last index into the determinant pair list

!     ndone becomes the last index of the determinant pairs that have been put in a task

            ibuf(1)=ibase
            ibuf(2)=jbase
            ibuf(3)=min(ndone+1,ijend)
            ibuf(4)=min(ndone+ntask,ijend)
            ndone=ibuf(4)

!     Send buffer to requesting rank without waiting for completion of mpi

            call MPI_iSend(ibuf,4,MPI_INTEGER8,iremote,2,
     &           MPI_COMM_WORLD,ireq,ierr)

!     Set ntasks(ibase,jbase) to 0 if this is the first task for the base pair

            if(ntasks(ibase,jbase).eq.-1) ntasks(ibase,jbase)=0
            ntasks(ibase,jbase)=ntasks(ibase,jbase)+1

!     Debug message

            if(idbg.gt.2) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,i5,a,6i5)')
     &             date(1:8),time(1:8),master,' sent ibuf to',
     &             iremote,igrp,(ibuf(i),i=1,4)
              flush(lfndbg)
            endif

!     Progress entry

            if(ipr.eq.2.or.ipro.eq.4) then
              call swatch(date,time)
              write(lfnpro,682) date(1:8),time(1:8),
     &             iremote,igrp,(ibuf(j),j=1,4),
     &             iloop(ibuf(3),1),iloop(ibuf(3),2)
 682          format(a,1x,a,' Sent ',2i6,' : ',2i5,2i10,2i6)
              flush(lfnpro)
            endif

!     Store the buffer information sent to remote rank in lgroup

            do k=1,4
              lgroup(igrp,k)=ibuf(k)
            enddo
            lgroup(igrp,5)=1
            
          enddo

          call timer_stop(95)

!     At this point the while loop over the determinant pairs for base pair ibase jabase has completed
          
        enddo
      enddo

!     At this point the loopover base pairs has completed and all tasks have been dispatched

!     Now the still outstanding tasks need to be collected
!     Assuming that some may be unresponsive, duplicates of outstanding tasks are sent

      call swatch(date,time)
      call timer_stop(99)
      write(lfnday,8023) date(1:8),time(1:8),
     &     timer_wall_total(99),
     &     '  :  Switching to duplicate tasks'
 8023 format(a8,2x,a8,f12.3,a)
      flush(lfnday)
      call timer_start(99)

      last=0

!     Count the outstanding tasks in nleft
      
      nleft=0
      do i=1,npg
        if(lgroup(i,5).gt.0) nleft=nleft+1
        lactive(i)=0
      enddo

!     Loop until no tasks are left to be received and accumulated

      call timer_start(97)
      do while(nleft.gt.0)

!     Wait for request from any remote rank

        call timer_start(96)
        call MPI_Recv(buffer,8,MPI_REAL8,MPI_ANY_SOURCE,1,
     &       MPI_COMM_WORLD,status,ierr)
        call timer_stop(96)

!     Determine the rank, group and base pair of received request

        iremote=status(MPI_SOURCE)
        igrp=map2(iremote+1,3)
        ibase=lgroup(igrp,1)
        jbase=lgroup(igrp,2)

!     Flag this group as still active

        lactive(igrp)=1
        if(idbg.gt.2) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &         master,' received last buffer from',iremote
          flush(lfndbg)
        endif
        
!     Accumulate results into appropriate arrays

        if(lgroup(igrp,5).gt.0) then

!     Determine the number of determinants completed for the base pair in ndets
!     Determine the percent complete for the base pair in pnrb

          ndets(lgroup(igrp,1),lgroup(igrp,2),2)=
     &         ndets(lgroup(igrp,1),lgroup(igrp,2),2)+
     &         lgroup(igrp,4)-lgroup(igrp,3)+1
          pnrb=dble(100*ndets(lgroup(igrp,1),lgroup(igrp,2),2))/
     &         ndets(lgroup(igrp,1),lgroup(igrp,2),1)

!     Write entry to progress file
          
          if(ipro.eq.1.or.ipro.eq.3) then
            if(ipro.eq.1) rewind(unit=lfnpro)   
            write(lfnpro,680) date(1:8),time(1:8),
     &           iremote,igrp,(lgroup(igrp,j),j=1,4),
     &           iloop(lgroup(igrp,3),1),iloop(lgroup(igrp,3),2),
     &           lgroup(igrp,4)-lgroup(igrp,3)+1,nrb,
     &           pnrb(lgroup(igrp,1),lgroup(igrp,2)),
     &           buffer(1),buffer(2),(int(buffer(j)),j=4,8)
            flush(lfnpro)
          elseif(ipro.eq.2.or.ipro.eq.4) then
            if(ipro.eq.2) rewind(unit=lfnpro)
            write(lfnpro,681) date(1:8),time(1:8),
     &           iremote,igrp,(lgroup(igrp,j),j=1,4),
     &           iloop(lgroup(igrp,3),1),iloop(lgroup(igrp,3),2),
     &           lgroup(igrp,4)-lgroup(igrp,3)+1,nrb,
     &           pnrb(lgroup(igrp,1),lgroup(igrp,2))
            flush(lfnpro)
          endif
          
!     Write entry to dayfile
          
          if(me.eq.master.and.
     &         pnrb(lgroup(igrp,1),lgroup(igrp,2)).ge.
     &         fday(lgroup(igrp,1),lgroup(igrp,2))) then
            call timer_stop(99)
            call swatch(date,time)
            write(lfnday,702) date(1:8),time(1:8),
     &           timer_wall_total(99),
     &           '  :  ',lgroup(igrp,1),lgroup(igrp,2),' at ',
     &           pnrb(lgroup(igrp,1),lgroup(igrp,2)),'%'
            flush(lfnday)
            call timer_start(99)
            fday(lgroup(igrp,1),lgroup(igrp,2))=
     &           fday(lgroup(igrp,1),lgroup(igrp,2))+dble(iday)
          endif
        endif
        
!     Count the number of tasks completed by the remote rank in numrecs
        
        if(lgroup(igrp,5).gt.0) then
          numrecs(iremote+1)=numrecs(iremote+1)+1
          
!     Accumulate the returned contributions in buffer
          
          hbase(lgroup(igrp,1),lgroup(igrp,2))=
     &         hbase(lgroup(igrp,1),lgroup(igrp,2))+buffer(1)
          sbase(lgroup(igrp,1),lgroup(igrp,2))=
     &         sbase(lgroup(igrp,1),lgroup(igrp,2))+buffer(2)
          tbase(lgroup(igrp,1),lgroup(igrp,2))=
     &         tbase(lgroup(igrp,1),lgroup(igrp,2))+buffer(3)
          nsing(lgroup(igrp,1),lgroup(igrp,2),1)=
     &         nsing(lgroup(igrp,1),lgroup(igrp,2),1)+int(buffer(4))
          nsing(lgroup(igrp,1),lgroup(igrp,2),2)=
     &         nsing(lgroup(igrp,1),lgroup(igrp,2),2)+int(buffer(5))
          nsing(lgroup(igrp,1),lgroup(igrp,2),3)=
     &         nsing(lgroup(igrp,1),lgroup(igrp,2),3)+int(buffer(6))
          nsing(lgroup(igrp,1),lgroup(igrp,2),4)=
     &         nsing(lgroup(igrp,1),lgroup(igrp,2),4)+int(buffer(7))
          nsing(lgroup(igrp,1),lgroup(igrp,2),5)=
     &         nsing(lgroup(igrp,1),lgroup(igrp,2),5)+int(buffer(8))

!     Decrement the tasks still to be done for base pair
          
          ntasks(lgroup(igrp,1),lgroup(igrp,2))=
     &         ntasks(lgroup(igrp,1),lgroup(igrp,2))-1

!     Loop over all base pairs to check if all tasks have been completed
!     This will be the case if ntasks(i,j) is zero

          i=lgroup(igrp,1)
          j=lgroup(igrp,2)

          if(ntasks(i,j).eq.0) then
            call timer_stop(99)
            call swatch(date,time)
            if(i.eq.j) then
              ltemp=idetb(i)*(idetb(i)+1)/2
            else
              ltemp=idetb(i)*idetb(j)
            endif
            write(lfnarx,522) i,j,idetb(i),idetb(j),ltemp,
     &           (nsing(i,j,k),k=1,4),
     &           hbase(i,j),sbase(i,j)
            if(ipr.ge.1) then
              if(i.ne.j) then
                write(*,607) i,j,idetb(i),idetb(j),ltemp,
     &               (nsing(i,j,k),k=1,4),
     &               hbase(i,j),sbase(i,j),tbase(i,j),
     &               timer_wall_total(99)
                write(lfnout,607) i,j,idetb(i),idetb(j),
     &               ltemp,
     &               (nsing(i,j,k),k=1,4),
     &               hbase(i,j),sbase(i,j),tbase(i,j),
     &               timer_wall_total(99)
              else
                write(*,608) i,j,idetb(i),idetb(j),ltemp,
     &               (nsing(i,j,k),k=1,4),
     &               hbase(i,j),sbase(i,j),
     &               hbase(i,j)/sbase(i,j),
     &               tbase(i,j),timer_wall_total(99)
                write(lfnout,608) i,j,idetb(i),idetb(j),
     &               ltemp,
     &               (nsing(i,j,k),k=1,4),
     &               hbase(i,j),sbase(i,j),
     &               hbase(i,j)/sbase(i,j),
     &               tbase(i,j),timer_wall_total(99)
              endif
              flush(lfnout)
            endif
            
            call swatch(date,time)
            write(lfncpr) i,j,
     &           hbase(i,j),sbase(i,j),tbase(i,j),
     &           (nsing(i,j,k),k=1,4)
            flush(lfncpr)
            
            if(i.ne.j) then
              write(lfnday,703) date(1:8),time(1:8),
     &             timer_wall_total(99),
     &             '  :  ',i,j,' completed, H,S   :',
     &             hbase(i,j),sbase(i,j)
            else
              write(lfnday,703) date(1:8),time(1:8),
     &             timer_wall_total(99),
     &             '  :  ',i,j,' completed, H,S,H*:',
     &             hbase(i,j),sbase(i,j),
     &             hbase(i,j)/sbase(i,j)
            endif
            flush(lfnday)
            call timer_start(99)

!     Setting ntasks to -2 identifies this base pair as completed and reported
            
            ntasks(i,j)=-2
          endif
          
!     Remove duplicates that have already been sent
          
          do i=1,npg
            if(lgroup(i,1).eq.lgroup(igrp,1).and.
     &           lgroup(i,2).eq.lgroup(igrp,2).and.
     &           lgroup(i,3).eq.lgroup(igrp,3).and.
     &           lgroup(i,4).eq.lgroup(igrp,4).and.
     &           lgroup(i,5).gt.0.and.i.ne.igrp) lgroup(i,5)=0
          enddo

!     Reset the current entry

          do i=1,5
            lgroup(igrp,i)=0
          enddo

!     Send a duplicate of an outstanding entry

          j=0
          do i=1,npg
            if(j.eq.0) then
              last=last+1
              if(last.gt.npg) last=last-npg
              if(lgroup(last,5).eq.1) j=last
            endif
          enddo

          if(j.gt.0) then
            do i=1,4
              ibuf(i)=lgroup(j,i)
            enddo

!     Send the duplicate buffer to requesting rank without waiting for completion of mpi

            call MPI_iSend(ibuf,4,MPI_INTEGER8,iremote,2,
     &           MPI_COMM_WORLD,ireq,ierr)

!     Debug message

            if(idbg.gt.2) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &             master,' sent duplicate ibuf to',iremote
              flush(lfndbg)
            endif

!     Progress entry

            if(ipro.eq.2.or.ipro.eq.4) then
              call swatch(date,time)
              write(lfnpro,682) date(1:8),time(1:8),
     &             iremote,igrp,(ibuf(k),k=1,4),
     &             iloop(ibuf(3),1),iloop(ibuf(3),2)
!     682          format(a,1x,a,' Sent ',2i6,' : ',2i5,2i10,2i6)
              flush(lfnpro)
            endif

!     Store the buffer information sent to remote rank in lgroup

            do k=1,4
              lgroup(igrp,k)=ibuf(k)
            enddo
            lgroup(igrp,5)=2
            lgroup(j,5)=2

          endif

        else

!     Received data other than matrix elements, should currently be impossible

          if(ipro.eq.1.or.ipro.eq.3) then
            write(lfnpro,683) date(1:8),time(1:8), 
     &           iremote,igrp,(lgroup(igrp,j),j=1,4),
     &           iloop(lgroup(igrp,3),1),iloop(lgroup(igrp,3),2),
     &           lgroup(igrp,4)-lgroup(igrp,3)+1
 683        format(a,1x,a,' RCVD ',2i6,' : ',2i5,2i10,
     &           3i6)
          endif
        endif

!     Count the outstanding tasks in nleft
        
        nleft=0
        do i=1,npg
          if(lgroup(i,5).gt.0) nleft=nleft+1
        enddo
      enddo

!     At this point all outstanding tasks completed

      call timer_stop(97)

!     Count the number of groups that are still active in nalive

      nalive=0
      do i=1,npg
        if(lactive(i).eq.1) nalive=nalive+1
      enddo
      
!     Signal worker processes to return from gronor_master
      
      ibuf(1)=0
      ibuf(2)=0
      ibuf(3)=0
      ibuf(4)=0
      do i=1,npg
        
        iremote=allgroups(i,2)
        
        call MPI_iSend(ibuf,4,MPI_INTEGER8,iremote,2,
     &       MPI_COMM_WORLD,ireq,ierr)
        if(idbg.gt.2) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &         master,' sent return to',iremote
          flush(lfndbg)
        endif
      enddo

!     Fill lower half of the Hamiltonian and overlap matrices
      
      do i=1,nbase
        do j=i+1,nbase
          hbase(i,j)=hbase(j,i)
          sbase(i,j)=sbase(j,i)
        enddo
      enddo

!     Output of the un-normalized results
      
      if(ipr.ge.40) then
        write(lfnout,609)
 609    format(//,' Hamiltonian Matrix Unnormalized',/)
        if(nbase.le.7) then
          write(lfnout,610) (i,i=1,nbase)
 610      format(6x,7(6x,i8,6x))
        else
          write(lfnout,611) (i,i=1,nbase)
 611      format(6x,7(6x,i8,6x))
        endif
        do j=1,nbase
          if(nbase.le.7) then
            write(lfnout,612) j,(hbase(i,j),i=1,nbase)
 612        format(i5,1x,7f20.10)
          else
            write(lfnout,613) j,(hbase(i,j),i=1,nbase)
 613        format(i5,1x,7f20.10,/,(6x,7f20.10))
          endif
        enddo
        write(lfnout,614)
 614    format(//,' Overlap Matrix Unnormalized',/)
        write(lfnout,610) (i,i=1,nbase)
        do j=1,nbase
          if(nbase.le.7) then
            write(lfnout,610) j,(sbase(i,j),i=1,nbase)
          else
            write(lfnout,611) j,(sbase(i,j),i=1,nbase)
          endif
        enddo
      endif

!     Sum the different singularities
      
      do i=1,nbase
        do j=1,nbase
          do k=1,5
            nsingt(k)=nsingt(k)+nsing(i,j,k)
          enddo
        enddo
        do j=1,i
          tsum=tsum+tbase(i,j)
        enddo
      enddo

!     Output the total singularites and wall clock times

      call timer_stop(99)
      write(*,615) mdet,(nsingt(i),i=1,4),tsum,
     &     dble(mgr)*tsum/dble(np-1),timer_wall_total(99),
     &     100.0*dble(mgr)*tsum/dble(np-1)/timer_wall_total(99)
      if(ipr.gt.0) then
        write(lfnout,615) mdet,(nsingt(i),i=1,4),tsum,
     &       dble(mgr)*tsum/dble(np-1),timer_wall_total(99),
     &       100.0*dble(mgr)*tsum/dble(np-1)/timer_wall_total(99)
      endif
 615  format(/,' Total',25x,2i14,2i12,i14,55x,1f12.1,/,
     &     /,' Total averaged per rank',t153,2f12.1,/,
     &     ' Task execution efficiency',t153,f12.1,' %')
      call timer_start(99)
      
      deallocate(pnrb,fday,lgroup,lactive,ntasks,ndets)
      
      return
      end
      
