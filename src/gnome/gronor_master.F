      subroutine gronor_master()

      use mpi
      use cidef
      use cidist
      use gnome_data
      use gnome_parameters

      implicit none

      integer :: ibase,jbase,nidet,njdet,mdet
      integer :: l2,i,j,k,ij,n,ndone,nleft,nrb
      integer (kind=4) :: ierr,ireq,iremote
      integer :: igrp, ltemp, last
      logical :: osame

      real (kind=8) :: fac, tsum
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
      
      do ibase=1,nbase
        do jbase=1,ibase
          hbase(ibase,jbase)=0.0d0
          sbase(ibase,jbase)=0.0d0
          tbase(ibase,jbase)=0.0d0
        enddo
        do jbase=1,nbase
          ntasks(ibase,jbase)=-1
          do i=1,5
            nsing(ibase,jbase,i)=0
          enddo
        enddo
      enddo
      
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

      ioff=0
      do ibase=1,nbase
        do jbase=1,ibase
          ndets(ibase,jbase,1)=ndxdet(ibase,jbase)-ioff
          ioff=ndxdet(ibase,jbase)
          ndets(ibase,jbase,2)=0
        enddo
      enddo
      
      if(ipr.ge.0) write(lfnout,618) mdet,numdet
 618  format(' Number of determinant pairs',t50,i16,/,
     &     ' Number of determinants above CI threshold',t50,i16)
      
      if(ipr.gt.1) then
        write(*,621)
        write(lfnout,621)
      endif
 621      format(//,' Calculation of Hamiltonian matrix elements',//,
     &         '   Element       ------- Determinants -------',
     &         '          -------------- Singularities -------------',
     &         '        ----- Unnormalized -----   ',
     &         '  Wall Clock  Wall Clock',/
     &         '  idet jdet      left     right         pairs',
     &         '          none           1           2            3+',
     &         '        Hamiltonian      Overlap   ',
     &         '        Cost        Time',/)

          ioff=0
          do ibase=1,nbase
            do jbase=1,ibase
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
              
              fday(ibase,jbase)=dble(iday)
              
              if(ipr.ge.8) then
                write(lfnout,600) ibase,jbase
 600            format(/,' Base state pair ',2i5)
                write(lfnout,601) ndeti
 601            format(/,' Left wavefunction has ',i8,
     &               ' determinants:',/)
                do i=1,ndeti
                  call unpack(ioccb(i,ibase),iocch,nacti)
                  write(lfnout,602) i,civb(i,ibase),(iocch(k),k=1,nacti)
 602              format(i5,f15.5,10i3)
                enddo
                write(lfnout,603) ndetj
 603            format(//,' Right wavefunction has ',i8,
     &               ' determinants:',/)
                do i=1,ndetj
                  call unpack(ioccb(i,jbase),iocch,nactj)
                  write(lfnout,602) i,civb(i,jbase),(iocch(k),k=1,nactj)
                enddo
              endif
              flush(lfnout)
              
              call swatch(date,time)
              call timer_stop(99)
              write(lfnday,802) date(1:8),time(1:8),
     &             timer_wall_total(99),
     &             '  :  ',ibase,jbase,' started, dimensions ',
     &             ndeti,ndetj
 802          format(a8,2x,a8,f12.3,a,2i6,a,2i12)
              flush(lfnday)
              call timer_start(99)
              
              nrb=0
              ndone=0
              ibuf(1)=0
              ibuf(2)=0
              ibuf(3)=0
              ibuf(4)=0

          do while(ndone.lt.ijend)
            
            call timer_start(94)
            call MPI_Recv(buffer,8,MPI_REAL8,MPI_ANY_SOURCE,1,
     &           MPI_COMM_WORLD,status,ierr)
            call timer_stop(94)
            call timer_start(95)
            
            iremote=status(MPI_SOURCE)
            igrp=map2(iremote+1,3)
            
            if(idbg.gt.2) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,i5,a,7i7)') date(1:8),time(1:8),
     &             master,' received buffer from',iremote,igrp,
     &             (lgroup(igrp,i),i=1,5)
              flush(lfndbg)
            endif
            
            if(lgroup(igrp,5).ne.0) then
              ndets(lgroup(igrp,1),lgroup(igrp,2),2)=
     &             ndets(lgroup(igrp,1),lgroup(igrp,2),2)+
     &             lgroup(igrp,4)-lgroup(igrp,3)+1
              pnrb(lgroup(igrp,1),lgroup(igrp,2))=
     &             dble(100*ndets(lgroup(igrp,1),lgroup(igrp,2),2))/
     &             ndets(lgroup(igrp,1),lgroup(igrp,2),1)
              if(ipro.eq.1.or.ipro.eq.3) then
                call swatch(date,time)
                if(ipro.eq.1) rewind(unit=lfnpro)
                write(lfnpro,680) date(1:8),time(1:8),
     &               iremote,igrp,(lgroup(igrp,j),j=1,4),
     &               iloop(lgroup(igrp,3),1),iloop(lgroup(igrp,3),2),
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
     &               iloop(lgroup(igrp,3),1),iloop(lgroup(igrp,3),2),
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
     &               '  :  ',lgroup(igrp,1),lgroup(igrp,2),' at ',
     &               pnrb(lgroup(igrp,1),lgroup(igrp,2)),'%'
 702            format(a8,2x,a8,f12.3,a,2i6,a,f5.1,a)
                flush(lfnday)
                call timer_start(99)
                fday(lgroup(igrp,1),lgroup(igrp,2))=
     &               fday(lgroup(igrp,1),lgroup(igrp,2))+dble(iday)
              endif
            endif
            
            if(lgroup(igrp,5).ne.0) then
              numrecs(iremote+1)=numrecs(iremote+1)+1
              
!     Accumulate the returned contributions in buffer
              
              hbase(lgroup(igrp,1),lgroup(igrp,2))=
     &             hbase(lgroup(igrp,1),lgroup(igrp,2))+buffer(1)
              sbase(lgroup(igrp,1),lgroup(igrp,2))=
     &             sbase(lgroup(igrp,1),lgroup(igrp,2))+buffer(2)
              tbase(lgroup(igrp,1),lgroup(igrp,2))=
     &             tbase(lgroup(igrp,1),lgroup(igrp,2))+buffer(3)
              nsing(lgroup(igrp,1),lgroup(igrp,2),1)=
     &             nsing(lgroup(igrp,1),lgroup(igrp,2),1)+int(buffer(4))
              nsing(lgroup(igrp,1),lgroup(igrp,2),2)=
     &             nsing(lgroup(igrp,1),lgroup(igrp,2),2)+int(buffer(5))
              nsing(lgroup(igrp,1),lgroup(igrp,2),3)=
     &             nsing(lgroup(igrp,1),lgroup(igrp,2),3)+int(buffer(6))
              nsing(lgroup(igrp,1),lgroup(igrp,2),4)=
     &             nsing(lgroup(igrp,1),lgroup(igrp,2),4)+int(buffer(7))
              nsing(lgroup(igrp,1),lgroup(igrp,2),5)=
     &             nsing(lgroup(igrp,1),lgroup(igrp,2),5)+int(buffer(8))
              
              ntasks(lgroup(igrp,1),lgroup(igrp,2))=
     &             ntasks(lgroup(igrp,1),lgroup(igrp,2))-1
              
              do k=1,5
                lgroup(igrp,k)=0
              enddo
              
              do i=1,ibase
                do j=1,jbase
                  if((i.ne.ibase.or.j.ne.jbase).and.
     &                 ntasks(i,j).eq.0) then
                    call timer_stop(99)
                    call swatch(date,time)
                    
                    
                    if(ipr.ge.1) then
                      if(i.eq.j) then
                        ltemp=idetb(i)*(idetb(i)+1)/2
                      else
                        ltemp=idetb(i)*idetb(j)
                      endif
                      write(*,622) i,j,idetb(i),idetb(j),ltemp,
     &                     (nsing(i,j,k),k=1,4),
     &                     hbase(i,j),sbase(i,j),tbase(i,j),
     &                     timer_wall_total(99)
                      write(lfnout,622) i,j,idetb(i),idetb(j),ltemp,
     &                     (nsing(i,j,k),k=1,4),
     &                     hbase(i,j),sbase(i,j),tbase(i,j),
     &                     timer_wall_total(99)
                      flush(lfnout)
 622                  format(1x,2i5,2i10,2i14,2i12,i14,f20.10,
     &                     f15.10,2f12.1)
                    endif
                    call swatch(date,time)
                    write(lfnday,703) date(1:8),time(1:8),
     &                   timer_wall_total(99),
     &                   '  :  ',i,j,' completed, H,S:',
     &                   hbase(i,j),sbase(i,j)
 703                format(a8,2x,a8,f12.3,a,2i6,a,2f20.10)
                    flush(lfnday)
                    call timer_start(99)
                    ntasks(i,j)=-2
                  endif
                enddo
              enddo
            endif

            ibuf(1)=ibase
            ibuf(2)=jbase
            ibuf(3)=min(ndone+1,ijend)
            ibuf(4)=min(ndone+ibatch,ijend)
            ndone=ibuf(4)

            call timer_stop(95)
            
            call MPI_iSend(ibuf,4,MPI_INTEGER8,iremote,2,MPI_COMM_WORLD,
     &           ireq,ierr)
          
            if(ntasks(ibase,jbase).eq.-1) ntasks(ibase,jbase)=0
            ntasks(ibase,jbase)=ntasks(ibase,jbase)+1
            
            if(idbg.gt.2) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,i5,a,6i5)') date(1:8),time(1:8),
     &             master,' sent ibuf to',iremote,igrp,(ibuf(i),i=1,4)
              flush(lfndbg)
            endif
            if(ipr.eq.2.or.ipro.eq.4) then
              call swatch(date,time)
              write(lfnpro,682) date(1:8),time(1:8),
     &             iremote,igrp,(ibuf(j),j=1,4),
     &             iloop(ibuf(3),1),iloop(ibuf(3),2)
 682          format(a,1x,a,' Sent ',2i6,' : ',2i5,2i10,2i6)
              flush(lfnpro)
            endif
            
            do k=1,4
              lgroup(igrp,k)=ibuf(k)
            enddo
            lgroup(igrp,5)=1
            
          enddo
          
        enddo
      enddo
      
      call swatch(date,time)
      call timer_stop(99)
      write(lfnday,8023) date(1:8),time(1:8),
     &     timer_wall_total(99),
     &     '  :  Switching to duplicate tasks'
 8023 format(a8,2x,a8,f12.3,a)
      flush(lfnday)
      call timer_start(99)

      last=0
      
      nleft=0
      do i=1,npg
        if(lgroup(i,5).gt.0) nleft=nleft+1
        lactive(i)=0
      enddo
      
      call timer_start(97)
      do while(nleft.gt.0)
        call timer_start(96)
        call MPI_Recv(buffer,8,MPI_REAL8,MPI_ANY_SOURCE,1,
     &       MPI_COMM_WORLD,status,ierr)
        call timer_stop(96)
        
        iremote=status(MPI_SOURCE)
        igrp=map2(iremote+1,3)
        
        lactive(igrp)=1
        if(idbg.gt.2) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &         master,' received last buffer from',iremote
          flush(lfndbg)
        endif
        
        if(lgroup(igrp,5).gt.0) then
          ndets(lgroup(igrp,1),lgroup(igrp,2),2)=
     &         ndets(lgroup(igrp,1),lgroup(igrp,2),2)+
     &         lgroup(igrp,4)-lgroup(igrp,3)+1
          pnrb=dble(100*ndets(lgroup(igrp,1),lgroup(igrp,2),2))/
     &         ndets(lgroup(igrp,1),lgroup(igrp,2),1)
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
          if(me.eq.master.and.
     &         pnrb(lgroup(igrp,1),lgroup(igrp,2)).ge.
     &         fday(lgroup(igrp,1),lgroup(igrp,2))) then
            call timer_stop(99)
            call swatch(date,time)
            write(lfnday,702) date(1:8),time(1:8),timer_wall_total(99),
     &           '  :  ',lgroup(igrp,1),lgroup(igrp,2),' at ',
     &           pnrb(lgroup(igrp,1),lgroup(igrp,2)),'%'
            flush(lfnday)
            call timer_start(99)
            fday(lgroup(igrp,1),lgroup(igrp,2))=
     &           fday(lgroup(igrp,1),lgroup(igrp,2))+dble(iday)
          endif
        endif
        
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
          
          ntasks(lgroup(igrp,1),lgroup(igrp,2))=
     &         ntasks(lgroup(igrp,1),lgroup(igrp,2))-1
          
          do i=1,nbase
            do j=1,nbase
              if(i.ne.ibase.and.j.ne.jbase.and.ntasks(i,j).eq.0) then
                call timer_stop(99)
                call swatch(date,time)
                
                if(ipr.ge.1) then
                  if(i.eq.j) then
                    ltemp=idetb(i)*(idetb(i)+1)/2
                  else
                    ltemp=idetb(i)*idetb(j)
                  endif
                  write(*,622) i,j,idetb(i),idetb(j),ltemp,
     &                 (nsing(i,j,k),k=1,4),
     &                 hbase(i,j),sbase(i,j),tbase(i,j),
     &                 timer_wall_total(99)
                  write(lfnout,622) i,j,idetb(i),idetb(j),ltemp,
     &                 (nsing(i,j,k),k=1,4),
     &                 hbase(i,j),sbase(i,j),tbase(i,j),
     &                 timer_wall_total(99)
                  flush(lfnout)
                endif
                call swatch(date,time)
                write(lfnday,703) date(1:8),time(1:8),
     &               timer_wall_total(99),
     &               '  :  ',i,j,'         completed, H,S:',
     &               hbase(i,j),sbase(i,j)
                flush(lfnday)
                call timer_start(99)
                ntasks(i,j)=-2
              endif
            enddo
          enddo
          
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
!              if(lgroup(i,5).eq.1) j=i
            endif
          enddo
          if(j.gt.0) then
            do i=1,4
              ibuf(i)=lgroup(j,i)
            enddo
            
            call MPI_iSend(ibuf,4,MPI_INTEGER8,iremote,2,MPI_COMM_WORLD,
     &           ireq,ierr)
            if(idbg.gt.2) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,i5,a,4i5)') date(1:8),time(1:8),
     &             master,' sent duplicate ibuf to',iremote
              flush(lfndbg)
            endif
            if(ipr.eq.2.or.ipro.eq.4) then
              call swatch(date,time)
              write(lfnpro,682) date(1:8),time(1:8),
     &             iremote,igrp,(ibuf(j),j=1,4),
     &             iloop(ibuf(3),1),iloop(ibuf(3),2)
! 682          format(a,1x,a,' Sent ',2i6,' : ',2i5,2i10,2i6)
              flush(lfnpro)
            endif
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
        nleft=0
        do i=1,npg
          if(lgroup(i,5).gt.0) nleft=nleft+1
        enddo
      enddo
      
      call timer_stop(97)
      
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
      
      do i=1,nbase
        do j=i+1,nbase
          hbase(i,j)=hbase(j,i)
          sbase(i,j)=sbase(j,i)
        enddo
      enddo
      
      if(ipr.ge.3) then
        write(lfnout,670)
 670    format(//,' Hamiltonian Matrix Unnormalized',/)
        if(nbase.le.7) then
          write(lfnout,671) (i,i=1,nbase)
 671      format(6x,7(6x,i8,6x))
        else
          write(lfnout,1671) (i,i=1,nbase)
 1671   format(6x,7(6x,i8,6x))
        endif
        do j=1,nbase
          if(nbase.le.7) then
            write(lfnout,672) j,(hbase(i,j),i=1,nbase)
 672        format(i5,1x,7f20.10)
          else
            write(lfnout,1672) j,(hbase(i,j),i=1,nbase)
 1672       format(i5,1x,7f20.10,/,(6x,7f20.10))
          endif
        enddo
        write(lfnout,673)
 673    format(//,' Overlap Matrix Unnormalized',/)
        write(lfnout,671) (i,i=1,nbase)
        do j=1,nbase
          if(nbase.le.7) then
            write(lfnout,672) j,(sbase(i,j),i=1,nbase)
          else
            write(lfnout,1672) j,(sbase(i,j),i=1,nbase)
          endif
        enddo
      endif
      
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
      write(*,623) mdet,(nsingt(i),i=1,4),tsum,
     &     dble(mgr)*tsum/dble(np-1)
      write(lfnout,623) mdet,(nsingt(i),i=1,4),tsum,
     &     dble(mgr)*tsum/dble(np-1)
 623  format(/,' Total',25x,2i14,2i12,i14,35x,2f12.1)
      
      deallocate(pnrb,fday,lgroup,lactive,ntasks,ndets)
      
      return
      end
      
