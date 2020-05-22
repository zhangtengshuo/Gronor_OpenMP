      subroutine gronor_calculate(ib,jb,id1,id2,lbuf)

      use mpi
      use cidist
      use cidef
      use gnome_data
      use gnome_parameters

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ireq,ierr
      integer (kind=4) :: iremote,status(MPI_STATUS_SIZE)
      integer :: ib,jb,id1,id2,lbuf,ii,jj,ig0

      integer :: i,ibv,idet,iin,ij,j,jin,k,l2,m,n,nck,nop,nvc,indxx
      integer :: idown,iup,iv,ncleff,nopeff,mspin,iclose,nb,ncl,index
      logical :: osamij,osame
      integer :: ipopcnt
      integer :: ntvc,ivc,ibas

      real (kind=8) :: fac

      ndeti=idetb(ib)
      ndetj=idetb(jb)
      nacti=nactb(ib)
      nactj=nactb(jb)
      if(idbg.ge.10) then
       write(lfndbg,600) ib,jb,id1,id2
 600   format(/,20('*'),' ',i5,' -',i5,' : ',i5,' -',i5,' ',20('*'))
       write(lfndbg,601) ndeti
 601   format(/,' Left wavefunction has ',i8,' determinants:',/)
       do i=1,ndeti
        call unpack(ioccb(i,ib),iocch,nacti)
        write(lfndbg,602) i,civb(i,ib),(iocch(k),k=1,nacti)
 602    format(I5,F15.5,10I3)
       enddo
       write(lfndbg,603) ndetj
 603   format(//,' Right wavefunction has ',i8,' determinants:',/)
       do i=1,ndetj
        call unpack(ioccb(i,jb),iocch,nactj)
        write(lfndbg,602) i,civb(i,jb),(iocch(k),k=1,nactj)
       enddo
      endif
      osame=ib.eq.jb
      call iloop_init(ndeti,ndetj,osame)
      if(ib.eq.jb) then
       l2=ndeti*(ndeti+1)/2
       ijend=l2
      else
       l2=ndeti*ndetj
       ijend=l2
      endif

      do i=1,8
       buffer(i)=0.0d0
      enddo
      e2summ=0.0d0
      e2buff=0.0d0

      do ij=id1,id2
       call timer_start(5)

       if(idbg.ge.10)
     &  write(lfndbg,608) ij, iloop(ij,1),iloop(ij,2)
 608    format(//,
     &  25('='),' Matrix element',i6,' (',i4,'-',i4,')',25('='),//,
     &  '         Left determinant      Right determinant',/,
     &  ' Irrep   inactive   active     inactive   active',//)

       nelec(1)=0
       nelec(2)=0
       ntop(1)=0
       ntop(2)=0
       ntcl(1)=0
       ntcl(2)=0

       i=iloop(ij,1)
       j=iloop(ij,2)
       hh=0.0d0
       ss=0.0d0
       
       if(dabs(civb(i,ib)*civb(j,jb)).lt.thresh_CI) then
          ising=4
          goto 105
       endif

       inacti=inactb(ib)
       inactj=inactb(jb)
       nacti=nactb(ib)
       nactj=nactb(jb)

       nclose(1)=inacti
       nclose(2)=inactj
       nopen(1)=nacti
       nopen(2)=nactj
       nelec(1)=2*nclose(1)
       nelec(2)=2*nclose(2)
       ntop(1)=ntop(1)+nopen(1)
       ntop(2)=ntop(2)+nopen(2)
       ntcl(1)=ntcl(1)+nclose(1)
       ntcl(2)=ntcl(2)+nclose(2)

!     the following if statement should be removed once
!     we have more than one symmetry representation

       nelec(1)=nelec(1)+ipopcnt(ioccb(1,ib),nacti)
       nelec(2)=nelec(2)+ipopcnt(ioccb(1,jb),nactj)

       if(nelec(1).ne.nelec(2)) stop

       i=iloop(ij,1)
       j=iloop(ij,2)

       hh=0.0d0
       ss=0.0d0
       if(ib.eq.jb) then
        if(i.eq.j) then
         fac=1.0d0
        else
         fac=2.0d0
        endif
       else
        fac=1.0d0
       endif

!     copy the occupation of the active orbitals for determinants i and j

       ninact(1)=inacti
       ninact(2)=inactj
       nact(1)=nacti
       nact(2)=nactj

       call unpack(ioccb(i,ib),iocch,nact(1))
       do k=1,nact(1)
        ioccup(k,1)=iocch(k)
        iocopen(k,1)=iocch(k)
       enddo
       call unpack(ioccb(j,jb),iocch,nact(2))
       do k=1,nact(2)
        ioccup(k,2)=iocch(k)
        iocopen(k,2)=iocch(k)
       enddo

       if(idbg.ge.10) then
        write(lfndbg,610)
 610    format(//,' Active orbital occupation',/)
        write(lfndbg,611) (ioccup(k,1),k=1,nact(1))
        write(lfndbg,612) (ioccup(k,2),k=1,nact(2))
 611    format(' Left determinant :',(20i4))
 612    format(' Right determinant:',(20i4))
       endif

       ncl=0
       nop=0
       nvc=0
       indxx=ib
       do idet=1,2

        if(idet.eq.2) indxx=jb
        index=0

        iv=1
        ncl=nclose(idet)
        nop=nopen(idet)

        nvc=ncl+nop
        if(nvc.gt.0) then
         nb=nbas
         ncleff=ncl
         nopeff=nop
         mspin=0
         do k=1,nop
          if(abs(iocopen(k,idet)).eq.1) then
           mspin=mspin+iocopen(k,idet)
          else
           nopeff=nopeff-1
           if(iocopen(k,idet).eq.2) ncleff=ncleff+1
          endif
         enddo
         iclose=ncl
         iup=ncleff
         idown=ncleff+(nopeff+mspin)/2

         do k=1,nvc
          if(k.le.ncl) then
           m=k
          else
           m=0
           if(iocopen(k-ncl,idet).eq.1) then
            iup=iup+1
            m=iup
           else
            if(iocopen(k-ncl,idet).eq.-1) then
             idown=idown+1
             m=idown
            else
             if(iocopen(k-ncl,idet).eq.2) then
              iclose=iclose+1
              m=iclose
             endif
            endif
           endif
          endif
          if(m.eq.0) then
           iv=iv+1
          else
           do ibv=1,nb
            vec(index+m,ibv,idet)=vecsb(ibv,iv,indxx)
           enddo
           iv=iv+1
          endif
         enddo

         index=index+ncleff+nopeff
         nclose(idet)=ncleff
         nopen(idet)=nopeff
         do k=1,nopen(idet)
          iocopen(k,idet)=1
          if(k.gt.(nopeff+mspin)/2) iocopen(k,idet)=-1
         enddo
         endif
       enddo

       if(idbg.ge.11) then
        do idet=1,2
         ntvc=ntcl(idet)+ntop(idet)
         write(lfndbg,1603) ntvc
 1603     format(/,' Closed shell M.O.''s',i5,/)
         do ivc=1,ntvc
          if(ivc.eq.ntcl(idet)+1) write(lfndbg,1604)
 1604     format(/,' Open shell M.O.'' s:')
          write(lfndbg,1605)  ' (',ivc,')',
     &      (vec(ivc,ibas,idet),ibas=1,nbas)
 1605     format(a2,i3,a1,(t9,10f12.8))
         enddo
        enddo
       endif

       if(idbg.ge.15) write(lfndbg,618) ij,i,j
 618   format(//,' Entering GNOME :',i8,', for determinants',2i8)
 
       call timer_start(6)
       call gronor_gnome(lfndbg)
       call timer_stop(6)

 105   continue

       if((icalc.eq.2.or.icalc.eq.0).and.ising.le.2) then

        if(icalc.eq.0) then
         if(mgr.le.1) then
          etot=e1+e2
         else
          etot=e1
          e2buff=e2buff+fac*civb(i,ib)*civb(j,jb)*e2
         endif
         if(iamhead.eq.1) then
          hh=etot-bias*deta
          ss=deta
         endif
        endif
       else
       endif

!      Upon return from gronor_gnome the head thread has the accumulated
!      results from all threads in the group

!      hh contains the Hamiltonian matrix element
!      ss contains the Overlap matrix element

       if(idbg.ge.15) write(lfndbg,619) ij,i,j
 619   format(/,' Exiting GNOME :',i8,', for determinants',2i8,//)

       if(idbg.ge.11) then
        write(lfndbg,651) me,iamhead,hh,ss
 651    format(' Calculated values from gronor_gnome ',2i5,/,             &
     & ' hh    ',f20.10,/,                                              &
     & ' ss    ',f20.10)
       endif

       buffer(1)=buffer(1)+fac*civb(i,ib)*civb(j,jb)*hh
       buffer(2)=buffer(2)+fac*civb(i,ib)*civb(j,jb)*ss
       buffer(ising+4)=buffer(ising+4)+1.0d0

       call timer_stop(5)
      enddo

      call timer_start(37)
      if(mgr.gt.1) then
       if(iamhead.eq.1) then
        do ii=1,mgr-1
         call MPI_Recv(e2summ,1,MPI_REAL8,MPI_ANY_SOURCE,5,
     &                 MPI_COMM_WORLD,status,ierr)
         if(idbg.gt.3) then
         iremote=status(MPI_SOURCE)
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,i5,a,i5)') date(1:8),time(1:8),
     &       me,' received e2buf from',iremote
          flush(lfndbg)
         endif
         e2buff=e2buff+e2summ
        enddo
       else
        call MPI_iSend(e2buff,1,MPI_REAL8,thisgroup(2),5,
     &                MPI_COMM_WORLD,ireq,ierr)

        if(idbg.gt.3) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,i5,a,i5)') date(1:8),time(1:8),
     &      me,' sent e2buf to',thisgroup(2)
         flush(lfndbg)
        endif

       endif
       buffer(1)=buffer(1)+e2buff
      endif
      call timer_stop(37)

      return
      end
      
