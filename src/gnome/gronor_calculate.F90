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

!>    @brief
!!    Driver for calculation range of Hamiltonian matrix element contributions
!!    on worker ranks
!!    @details
!!    This routine calculates all Slater determinant contributions id1 to id2
!!    for <ib|H|jb>, and is called by all worker ranks.
!!    If the calling rank is part of a group of ranks assigned to calculating
!!    this contribution, the head rank of the group has the accumulated result
!!    upon exit of this routine.
!!    @author
!!    T. P. Straatsma (ORNL)
!!    @param[in] ib  : left basestate
!!    @param[in] jb  : right basestate
!!    @param[in] id1 : index of initial matrix contribution
!!    @param[in] id2 : index of final matrix contribution
!!    @date    2016

subroutine gronor_calculate(ib,jb,id1,id2)

  use mpi
  use cidist
  use cidef
  use gnome_data
  use gnome_parameters

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  external :: swatch,timer_start,timer_stop
!  external :: MPI_iSend,MPI_Recv
  external :: gronor_gnome,gronor_abort

  integer (kind=4) :: ireq,ierr
  integer (kind=4) :: iremote,status(MPI_STATUS_SIZE)
  integer (kind=4) :: ncount,mpitag,mpidest
  integer :: ib,jb,id1,id2,ii,ihc,nhc,iact

  integer :: i,ibv,idet,ij,j,k,l2,m,nop,nvc,indxx
  integer :: idown,iup,iv,ncleff,nopeff,mspin,iclose,nb,ncl,indexv
  integer :: ntvc,ivc,ibas

  integer :: ioff

  real (kind=8) :: btemp,btemp2

  logical (kind=4) :: flag

  ndeti=idetb(ib)
  ndetj=idetb(jb)
  nacti=nactb(ib)
  nactj=nactb(jb)
  if(idbg.ge.50) then
    write(lfndbg,600) ib,jb,id1,id2
600 format(/,20('*'),' ',i5,' -',i5,' : ',i5,' -',i5,' ',20('*'))
    write(lfndbg,601) ndeti
601 format(/,' Left wavefunction has ',i8,' determinants:',/)
    if(idbg.gt.90) then
      do i=1,ndeti
        write(lfndbg,602) i,civb(i,ib),(iocc(i,ib,k),k=1,nacti)
602     format(i5,f15.5,32i3)
      enddo
    else
      do i=1,min(ndeti,10)
        write(lfndbg,602) i,civb(i,ib),(iocc(i,ib,k),k=1,nacti)
      enddo
      do i=max(11,ndeti-10),ndeti
        write(lfndbg,602) i,civb(i,ib),(iocc(i,ib,k),k=1,nacti)
      enddo
    endif
    write(lfndbg,603) ndetj
603 format(//,' Right wavefunction has ',i8,' determinants:',/)
    do i=1,ndetj
      write(lfndbg,602) i,civb(i,jb),(iocc(i,jb,k),k=1,nactj)
    enddo
  endif

  ioff=ndxdet(ib,jb)

  if(ib.eq.jb) then
    l2=ndeti*(ndeti+1)/2
    ijend=l2
  else
    l2=ndeti*ndetj
    ijend=l2
  endif

  btemp=0.0d0
  btemp2=0.0d0
  do i=1,17
    buffer(i)=0.0d0
  enddo
  e2summ=0.0d0
  e2buff=0.0d0

  nhc=id2-id1+1
  e1tot=0.0d0
  e2tot=0.0d0
  etotb=0.0d0
  sstot=0.0d0

  do ij=id1,id2

    if(otreq.and.iint.ne.0) then
      call MPI_Test(itreq,flag,status,ierr)

      if(flag) then
        if(idbg.gt.10) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,a)') date(1:8),time(1:8),' Terminating in gronor_calculate'
        endif
        oterm=.true.
        return
      endif
    endif

    call timer_start(5)

    if(idbg.ge.50) write(lfndbg,608) ij
608 format(//,25('='),' Matrix element',i6,25('='),//, &
        '         Left determinant      Right determinant',/, &
        ' Irrep   inactive   active     inactive   active',//)

    nelec(1)=0
    nelec(2)=0
    ntop(1)=0
    ntop(2)=0
    ntcl(1)=0
    ntcl(2)=0

    i=melist(ij,1)
    j=melist(ij,2)

    ihc=ij-id1+1

    hh=0.0d0
    ss=0.0d0
    do k=1,9
      mpoles(k)=0.0d0
    enddo

    if(dabs(civb(i,ib)*civb(j,jb)).lt.tau_CI) then
      ising=4
    else

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

      do iact=1,nacti
        nelec(1)=nelec(1)+iabs(int(iocc(1,ib,iact),kind=kind(nelec(1))))
      enddo
      do iact=1,nactj
        nelec(2)=nelec(2)+iabs(int(iocc(1,jb,iact),kind=kind(nelec(2))))
      enddo

      if(nelec(1).ne.nelec(2)) call gronor_abort(300,"Inconsistent number electrons")

      i=melist(ij,1)
      j=melist(ij,2)

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
      fctr=fac*civb(i,ib)*civb(j,jb)

      !     copy the occupation of the active orbitals for determinants i and j

      ninact(1)=inacti
      ninact(2)=inactj
      nact(1)=nacti
      nact(2)=nactj

      do k=1,nact(1)
        ioccup(k,1)=iocc(i,ib,k)
        iocopen(k,1)=iocc(i,ib,k)
      enddo

      do k=1,nact(2)
        ioccup(k,2)=iocc(j,jb,k)
        iocopen(k,2)=iocc(j,jb,k)
      enddo

      if(idbg.ge.40) then
        write(lfndbg,610)
610     format(//,' Active orbital occupation',/)
        write(lfndbg,611) (ioccup(k,1),k=1,nact(1))
        write(lfndbg,612) (ioccup(k,2),k=1,nact(2))
611     format(' Left determinant :',(20i4))
612     format(' Right determinant:',(20i4))
      endif

      ncl=0
      nop=0
      nvc=0
      indxx=ib
      do idet=1,2
        ncl=nclose(idet)
        nop=nopen(idet)
        do iv=1,ncl+nop
          do k=1,nbas
            vec(iv,k,idet)=0.0d0
          enddo
        enddo
      enddo
      
      do idet=1,2
        if(idet.eq.2) indxx=jb
        indexv=0
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
                vec(indexv+m,ibv,idet)=vecsb(ibv,iv,indxx)
              enddo
              iv=iv+1
            endif
          enddo

          indexv=indexv+ncleff+nopeff
          nclose(idet)=ncleff
          nopen(idet)=nopeff
          do k=1,nopen(idet)
            iocopen(k,idet)=1
            if(k.gt.(nopeff+mspin)/2) iocopen(k,idet)=-1
          enddo
        endif
      enddo

      if(idbg.ge.60) then
        do idet=1,2
          ntvc=ntcl(idet)+ntop(idet)
          write(lfndbg,1603) ntvc
1603      format(/,' Closed shell M.O.''s',i5,/)
          do ivc=1,ntvc
            if(ivc.eq.ntcl(idet)+1) write(lfndbg,1604)
1604        format(/,' Open shell M.O.'' s:')
            write(lfndbg,1605)  ' (',ivc,')',(vec(ivc,ibas,idet),ibas=1,nbas)
1605        format(a2,i3,a1,(t9,10f12.8))
          enddo
        enddo
      endif

      if(idbg.ge.40) write(lfndbg,618) ij,i,j
618   format(//,' Entering GNOME :',i8,', for determinants',2i8)

      call timer_start(6)
      call gronor_gnome(lfndbg,ihc,nhc)
      call timer_stop(6)

      if(oterm) then
        call timer_stop(5)
        return
      endif

      if(idbg.ge.40) write(lfndbg,619) ij,i,j
619   format(/,' Exiting GNOME :',i8,', for determinants',2i8,//)

    endif

    !     hh contains the Hamiltonian matrix element
    !     ss contains the Overlap matrix element

    if(nbatch.eq.0) then
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
      endif
      if(idbg.ge.12) then
        call swatch(date,time)
        write(lfndbg,651) date(1:8),time(1:8),ij,i,j,hh,ss
651     format(a,1x,a,' Calculated values from gronor_gnome ',3i7,'  H:',f20.10,'  S:',f20.10)
      endif
      buffer(1)=buffer(1)+fac*civb(i,ib)*civb(j,jb)*hh
      buffer(2)=buffer(2)+fac*civb(i,ib)*civb(j,jb)*ss
      if(iamhead.eq.1)then
        do k=1,9
          buffer(k+8)=buffer(k+8)+fac*civb(i,ib)*civb(j,jb)*mpoles(k)
        enddo
        if(idbg.ge.12)then
          write(lfndbg,*)'Gronor calculate head rank contributions', &
              fac,(civb(i,ib)*civb(j,jb)*mpoles(k),k=1,9)
        endif
      endif
    else

    endif
    buffer(ising+4)=buffer(ising+4)+1.0d0

    if(iamhead.eq.1) e1tot=e1tot+e1*fctr

    call timer_stop(5)

  enddo

  if(nbatch.ne.0) then
    if(iamhead.eq.1) then
      buffer(1)=e1tot+etotb
      buffer(2)=sstot
      do k=1,9
        buffer(k+8)=buffer(k+8)+fac*civb(i,ib)*civb(j,jb)*mpoles(k)
      enddo
      e2buff=0.0d0
    else
      e2buff=e1tot+etotb
    endif

  endif

  !     Upon return from gronor_calculate the head thread has the
  !     accumulated results from all threads in the group

  call timer_start(37)
  if(mgr.gt.1) then
    if(iamhead.eq.1) then
      do ii=1,mgr-1
        ncount=1
        mpitag=5
        call MPI_Recv(e2summ,ncount,MPI_REAL8,MPI_ANY_SOURCE,mpitag,MPI_COMM_WORLD,status,ierr)
        iremote=status(MPI_SOURCE)
        if(idbg.gt.10) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,i5,a,i5)') date(1:8),time(1:8),me,' received e2buf from',iremote
          flush(lfndbg)
        endif
        e2buff=e2buff+e2summ
      enddo
    else
      ncount=1
      mpitag=5
      mpidest=thisgroup(2)
      call MPI_iSend(e2buff,ncount,MPI_REAL8,mpidest,mpitag,MPI_COMM_WORLD,ireq,ierr)
      call MPI_Request_free(ireq,ierr)
      if(idbg.gt.10) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,i5,a,i5)') date(1:8),time(1:8),me,' sent e2buf to',thisgroup(2)
        flush(lfndbg)
      endif

    endif
    buffer(1)=buffer(1)+e2buff
  endif
  call timer_stop(37)

  return
end subroutine gronor_calculate
