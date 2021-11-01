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
!! Collect and print timings to the output file
!!
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
!! The code is instrumented with timers to collect CPU and wall-clock times for the 
!! computationally demanding parts of the calculation. These timings are presented for each 
!! of the ranks used in the calculation.
!! Because of the fault resilient implementation, the master ranks only waits for a limited 
!! time for timning data from other ranks. In that case values of zero are reported for those
!! ranks. This is not an error!
!!
 
      subroutine gronor_timings(lfnout,lfnday,lfntim)

      use mpi
      use cidist
      use gnome_parameters

      implicit none

      integer (kind=8) :: lfnout,lfnday,lfntim,i,j,idum
      external :: timer_wall_total,timer_calls
      real (kind=8) :: timer_wall_total
      integer :: timer_calls

      real (kind=8), allocatable :: timings(:,:)
      real (kind=8), allocatable :: tdata(:),tdat(:),taver(:)
      integer (kind=4), allocatable :: isync(:),ireq(:)
      integer (kind=4) :: ierr,iremote,status(MPI_STATUS_SIZE)
      integer (kind=4) :: ncount,mpitag,mpidest,mpireq
      integer (kind=8) :: irtim,nreqso,nreqsc
      integer (kind=8) :: nrecav
      character (len=255) :: date,time
      logical (kind=4) :: flag
      logical :: openrcv

      if(itim.gt.0) then

       allocate(timings(np,68),tdata(68))

       if(me.eq.master) then

         nalive=nalive*mgr

         allocate(isync(nalive))
         allocate(ireq(nalive))
         allocate(tdat(68),taver(68))
         
         do i=1,nalive
           isync(i)=-1
         enddo
         
         do j=1,68
           do i=1,np
             timings(i,j)=0.0d0
           enddo
         enddo
         do j=1,60
           timings(me+1,j)=timer_wall_total(j)
         enddo
         timings(me+1,61)=dble(timer_calls(6))
         timings(me+1,62)=dble(timer_calls(31))
         timings(me+1,63)=dble(timer_calls(32))
         timings(me+1,64)=dble(timer_calls(34))
         timings(me+1,65)=dble(timer_calls(35))
         timings(me+1,66)=dble(nacc0)
         timings(me+1,67)=dble(nacc1)
         timings(me+1,68)=dble(memavail)

         nreqso=0
         nreqsc=0
         openrcv=.false.
         irtim=1
         do i=1,np

           ncount=1
           mpitag=20
           mpidest=i-1
           if(mpidest.ne.master)
     &          call MPI_iSend(irtim,ncount,MPI_INTEGER8,mpidest,
     &          mpitag,MPI_COMM_WORLD,mpireq,ierr)
           nreqso=nreqso+1
             
           if(.not.openrcv) then
             ncount=68
             mpitag=11
             call MPI_iRecv(tdat,ncount,MPI_REAL8,MPI_ANY_SOURCE,
     &            mpitag,MPI_COMM_WORLD,itreq,ierr)
             openrcv=.true.
           endif           

           do while(nreqso.ge.itim.or.(i.eq.np.and.nreqsc.lt.nalive))
             call MPI_Test(itreq,flag,status,ierr)
             if(flag) then
               iremote=status(MPI_SOURCE)
               isync(i)=iremote
               do j=1,68
                 timings(iremote+1,j)=tdat(j)
               enddo
               nreqso=nreqso-1
               nreqsc=nreqsc+1

               ncount=68
               mpitag=11
               call MPI_iRecv(tdat,ncount,MPI_REAL8,MPI_ANY_SOURCE,
     &              mpitag,MPI_COMM_WORLD,itreq,ierr)
               openrcv=.true.

             endif
           enddo

         enddo

       else
         
         do j=1,60
           tdata(j)=timer_wall_total(j)
         enddo
         tdata(61)=dble(timer_calls(6))
         tdata(62)=dble(timer_calls(31))
         tdata(63)=dble(timer_calls(32))
         tdata(64)=dble(timer_calls(34))
         tdata(65)=dble(timer_calls(35))
         tdata(66)=dble(nacc0)
         tdata(67)=dble(nacc1)
         tdata(68)=dble(memavail)

         ncount=1
         mpitag=20
         call MPI_Recv(irtim,ncount,MPI_INTEGER8,master,
     &        mpitag,MPI_COMM_WORLD,status,ierr)
         
         ncount=68
         mpitag=11
         call MPI_Send(tdata,ncount,MPI_REAL8,master,
     &        mpitag,MPI_COMM_WORLD,ierr)

       endif
       
       if(me.eq.master) then

         do i=1,68
           taver(i)=0.0d0
           do j=1,np
             taver(i)=taver(i)+timings(j,i)
           enddo
           taver(i)=taver(i)/dble(nalive)
         enddo
         
         write(lfntim,500) np,68
 500     format(2i10)
         write(lfntim,501) ((timings(j,i),i=1,68),j=1,np)
 501     format(10f12.3)
         write(lfntim,501) (taver(j),j=1,68)

         write(lfnout,600)
 600     format(/,' Wallclock Timing Analysis Main Program',/)
         write(lfnout,601)
 601     format('  Proc',
     &        '       Total',
     &        '       Input',
     &        '   Integrals',
     &        '    H Matrix',
     &        '    Elements',
     &        '       Gnome',
     &        '       Calls',
     &        '  Unused Dev Mem'/)

         flush(lfnout)
         do i=1,np
           if(i-1.eq.master) then
             write(lfnout,602) i-1,timings(i,1),timings(i,2),
     &            timings(i,9),timings(i,4)
           else
             write(lfnout,602) i-1,timings(i,1),timings(i,2),
     &            timings(i,9),timings(i,4),timings(i,5),timings(i,6),
     &            int(timings(i,61)),int(timings(i,68))
 602         format(1x,i5,6f12.3,i12,i16,4f12.3)
           endif
         enddo
         write(lfnout,622) taver(1),taver(2),taver(9),taver(4),
     &        taver(5),taver(6),int(taver(61)),int(taver(68))
 622     format(1x,105('-'),/,' Avrg ',6f12.3,i12,i16,4f12.3)
         
         flush(lfnout)
         
         write(lfnout,603)
 603     format(//,' Wallclock Timing Analysis Gnome',/)
         write(lfnout,604)
 604     format('  Proc',
     &        '     TransVc',
     &        '       Order',
     &        '     TranOut',
     &        '      MOOver',
     &        '      CoFac1',
     &        '      CorOrb',
     &        '      TraMat',
     &        '      Dipole',
     &        '       TrSym',
     &        '     TraMat2',
     &        '       GnOne',
     &        '       GnTwo',/)
         
         flush(lfnout)
         
         do i=1,np
           if(i-1.ne.master) then
             write(lfnout,605) i-1,(timings(i,j),j=11,22)
 605         format(1x,i5,12f12.3)
           endif
         enddo
         write(lfnout,625) (taver(j),j=11,22)
 625     format(1x,149('-'),/,' Avrg ',12f12.3)
         
         flush(lfnout)
         
         write(lfnout,606)
 606     format(//,' Wallclock Timing Analysis CoFac',/)
         write(lfnout,607)
 607     format('  Proc',
     &        '         SVD',
     &        '         Sum',
     &        '         EVD',
     &        '        Rest',/)
         
         flush(lfnout)
         
         do i=1,np
           if(i-1.ne.master) then
             write(lfnout,608) i-1,(timings(i,j),j=41,44)
 608         format(1x,i5,4f12.3)
           endif
         enddo
         write(lfnout,628) (taver(j),j=41,44)
 628     format(1x,53('-'),/,' Avrg ',12f12.3)
         
         flush(lfnout)
         
         write(lfnout,609)
 609     format(//,' Wallclock Timing Analysis GnTwo',/)
         
         write(lfnout,610)
 610     format('  Proc',
     &        '   Transpose',
     &        '      Two0-a',
     &        '      Two0-n',
     &        '      Two0-w',
     &        '      Two1-a',
     &        '      Two1-n',
     &        '      Two1-w',
     &        '       Comm1',
     &        '       Comm2',/)
         
         flush(lfnout)
         
         do i=1,np
           if(i-1.ne.master) then
             write(lfnout,611) i-1,(timings(i,j),j=30,36),
     &            timings(i,39),timings(i,37)
 611         format(1x,i5,9f12.3)
           endif
         enddo
         write(lfnout,631) (taver(j),j=30,36),
     &            taver(39),taver(37)
 631     format(1x,113('-'),/,' Avrg ',9f12.3)
         
         flush(lfnout)
         
         write(lfnout,612)
 612     format(//,'  Proc',
     &        '     NumRecv',
     &        '        Num0',
     &        '        Num1',
     &        '        Ndx0',
     &        '        Ndx1',/)
         
         do i=1,np
           if(i-1.ne.master) then
             write(lfnout,613) i-1,numrecs(i),
     &            int(timings(i,62))+int(timings(i,63)),
     &            int(timings(i,64))+int(timings(i,65)),
     &            int(timings(i,66)),int(timings(i,67))
 613         format(1x,i5,5i12)
           endif
         enddo
         nrecav=0
         do i=1,np
           if(i-1.ne.master) nrecav=nrecav+numrecs(i)
         enddo
         nrecav=int(dble(nrecav)/dble(nalive-1))
         write(lfnout,633) nrecav,
     &        int(taver(62))+int(taver(63)),
     &        int(taver(64))+int(taver(65)),
     &        int(taver(66)),int(taver(67))
 633     format(1x,65('-'),/,' Avrg ',5i12)
         
         flush(lfnout)
         
       endif
       
c       if(me.eq.master) then
c         idum=1
c         do i=1,nalive
c           ncount=1
c           mpitag=12
c           mpidest=isync(i)
c           if (isync(i).ne.-1)
c     &          call MPI_iSend(idum,ncount,MPI_INTEGER8,mpidest,
c     &          mpitag,MPI_COMM_WORLD,mpireq,ierr)
c           ireq(i)=mpireq
c         enddo
c       else
c         ncount=1
c         mpitag=12
c         call MPI_Recv(idum,ncount,MPI_INTEGER8,master,
c     &        mpitag,MPI_COMM_WORLD,status,ierr)
c       endif

       deallocate(tdata,timings)
       if(me.eq.master) deallocate(tdat,taver,ireq,isync)
       
      endif
      
      call timer_stop(99)

      if(me.eq.master) then
        if(ipr.ge.40) then
          write(lfnout,620) mgr*nalive+1,np
 620      format(/,' At the end of the run ',i8,' of ',i8,
     &         ' processes are still active')
        endif
        if(ipr.ge.20) then
          write(lfnout,621)
 621      format(/,' Timings summary')
          if(ipr.ge.2) then
            write(lfnout,617) numtasks,
     &           int(dble(numtasks)/timer_wall_total(98)),
     &           timer_wall_total(2),
     &           timer_wall_total(3),
     &           timer_wall_total(7),
     &           timer_wall_total(8),
     &           timer_wall_total(9), tmax
 617        format(/,
     &           ' Total number of tasks',t43,i20,/,
     &           ' Number of tasks per second',t43,i20,//,
     &           ' Reading and distribution input',t55,f12.3,/,
     &           ' Assignment of ranks',t55,f12.3,/,
     &           ' Generation base states',t55,f12.3,/,
     &           ' Generation matrix element list',t55,f12.3,/,
     &           ' Reading and distribution of integrals',t55,f12.3,/,
     &           ' Maximum task execution time',t55,f12.3)
            write(lfnout,615) 
     &           timer_wall_total(94),
     &           100.0d0*(timer_wall_total(94)/timer_wall_total(95)),
     &           timer_wall_total(95)-timer_wall_total(94),
     &           100.0d0*(timer_wall_total(95)-timer_wall_total(94))/
     &           timer_wall_total(95),
     &           timer_wall_total(96),
     &           100.0d0*(timer_wall_total(96)/timer_wall_total(97)),
     &           timer_wall_total(97)-timer_wall_total(96),
     &           100.0d0*(timer_wall_total(97)-timer_wall_total(96))/
     &           timer_wall_total(97)
 615        format(/,
     &           ' Request  time on master',t55,f12.3,' (',f7.3,'%)',/,
     &           ' Response time on master',t55,f12.3,' (',f7.3,'%)',/,
     &           ' Requestf time on master',t55,f12.3,' (',f7.3,'%)',/,
     &           ' Finalize time on master',t55,f12.3,' (',f7.3,'%)')
          endif
          write(lfnout,614) 
     &         timer_wall_total(99)-timer_wall_total(98),
     &         100d0*(timer_wall_total(99)-timer_wall_total(98))/
     &         timer_wall_total(99),
     &         timer_wall_total(98),
     &         100d0*timer_wall_total(98)/timer_wall_total(99)
 614      format(/,' Setup wall clock time on master',t55,f12.3,
     &         ' (',f7.3,'%)',/,
     &         ' Matrix element calculation wall clock time on master',
     &         t55,f12.3,' (',f7.3,'%)')
        endif
        
        if(ipr.ge.0) then
          write(lfnout,619) timer_wall_total(99)
 619      format(/,' Total execution wall clock time on master ',
     &         t55,f12.3)
        endif
        call swatch(date,time)
        write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &       '  :  Timings reported'
 702    format(a8,2x,a8,f12.3,a)
        call swatch(date,time)
        write(lfnday,703) date(1:8),time(1:8),timer_wall_total(99),
     &       nalive,npg
 703    format(a8,2x,a8,f12.3,'  :  Completion of job with',i8,' of ',
     &       i8,' groups still active')
        flush(lfnday)
      endif
      
      return
      end
      
