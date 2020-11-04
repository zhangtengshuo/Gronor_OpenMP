      subroutine gronor_read_integrals()

!      use mpi
      use inp
      use cidef
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none
#include 'mpif.h'

!tps#include 'mpif.h'

      integer (kind=4) :: ierr
      integer :: i,ie,im,is,ke,km,ks,le,lm,ne,ig0,j,ni,ng,nz
      integer :: nm,ns,kl,k,ll,l,ls,nl,n,nn,igr,ngi,ngir
      integer :: nmax,kbloc,mdum,itarget,noff
      integer, allocatable :: idum(:)
      integer :: ielem,jkoff
      real (kind=8) :: fik,fln,fikln
      real (kind=8) :: sdum
      real (kind=8) :: b(510)
      integer :: status(MPI_STATUS_SIZE)

      integer :: kint, kintbatch, mint,nint,jint,lint
      integer :: nacc, itoff,itnum
      integer :: icol,ncol,irow,nrow
      integer :: incint
      integer :: nsmall
      integer :: ilast,ii,jj
      integer(kind=2)::ibuf(1020)

      allocate(ig(mgr))
      allocate(ntarget(mgr))

      call MPI_Comm_Group(MPI_COMM_WORLD,group_world,ierr)

      allocate(new_comm(mgr))

      allocate(ranks_list(numgrp))

      do i=1,mgr
       do j=1,numgrp
         ranks_list(j)=allgroups(j,i+1)
       enddo

       call MPI_Group_Incl(group_world,int(numgrp,kind=4),
     & ranks_list,new_group,ierr)
       call MPI_Comm_Create(MPI_COMM_WORLD,new_group,
     &                      new_comm(i),ierr)
      enddo

      deallocate(ranks_list)

      allocate(ranks_list(mgr))
      allocate(grp_comm(numgrp))

      do i=1,numgrp
       do j=1,mgr
        ranks_list(j)=allgroups(i,j+1)
       enddo

       call MPI_Group_Incl(group_world,int(mgr,kind=4),
     &  ranks_list,grp_group,ierr)
       call MPI_Comm_Create(MPI_COMM_WORLD,grp_group,
     &                      grp_comm(i),ierr)
      enddo

      deallocate(ranks_list)

      if(me.eq.master) then

       intfil=0
       open(unit=lfnone,file=filone,form='unformatted',status='old',
     &  err=99)
       intfil=1
       goto 199
 99    continue

       open(unit=lfnint,file=filint,form='unformatted',status='old',
     &  err=994)
 199   continue
      endif
      call MPI_Bcast(intfil,1,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)

      if(intfil.eq.1) then
       call gronor_parallel_integral_input()
       flush(lfnout)
       return
      endif
      write(lfnout,*)'Do we ever get here?'
      flush(lfnout)

      if(me.eq.master) then
       rewind(lfnint)
       if(ione.eq.0) then
        read(lfnint,err=993) namint,intone,potnuc,nbas,mbuf,mdum
       else
        read(lfnint,err=993) namint,nbas,intone,potnuc
       endif
       if(ipr.ge.1) then
        write(lfnout,601) namint
        flush(lfnout)
  601   format(/,' Integral file title: ',t25,a,/,t25,a)
        write(lfnout,602) nbas,intone,potnuc
  602   format(/,' Number of basis functions',t50,i16,                    &
     &  /,' Number of one electron integrals',t50,i16,                    &
     &  //,' Nuclear potential energy',t57,f20.10)
       endif
       int1=(nbas*(nbas+1))/2
       if(ipr.ge.1) then
        write(lfnout,603) intone,int1
 603    format(/,' Number of one electron integrals',t50,i16,             &
     &  /,' Number of calculated one electron integrals ',t50,i16)
       endif
       allocate(s(nbas,nbas))
       allocate(t(int1))
       allocate(v(int1))
       if(ione.eq.0) then                                            !  Next simplification: restrict ione to 0 
        read(lfnint,err=993) ((s(i,j),i=1,j),j=1,nbas)
        do i=1,nbas
         do j=1,i
           s(i,j)=s(j,i)
         end do
        end do
        read(lfnint,err=993) (t(i),i=1,int1)
        read(lfnint,err=993) (v(i),i=1,int1)
       else
        call gronor_read1(lfnint,t,intone)
        call gronor_read1(lfnint,t,intone)
        call gronor_read1(lfnint,t,intone)
        call gronor_read1(lfnint,v,intone)
       endif
       int2=0
       read(lfnint,end=993) nt

       if (nt(2).eq.-1) read(lfnint) nt(2)

       if(ipr.ge.1) then
        write(lfnout,604) int1*(int1+1)/2,nt(2),
     &   dble(200*nt(2))/dble(int1*(int1+1))
 604    format(/,' Number of two electron integrals',t50,i16,             &
     &   /,' Number of non-zero two electron integrals',t50,i16,          &
     &     ' (',f6.2,'%)')
#ifdef SINGLEP
        write(lfnout,614)
 614    format(' Integrals are used in single precision')
#else
        write(lfnout,615)
 615    format(' Integrals are used in double precision')
#endif
       endif 

!     int1  : number of one electron integrals
!     int2  : number of two electron integrals
!     mint2 : dimension of the two electron integral array
!     mgr   : number of ranks in a group
!     ig(i) : number of integrals on rank i

       int2=nt(2)
       mint2=int2/mgr
       if(mod(mint2,mgr).ne.0) mint2=mint2+1
       do i=1,mgr
        ig(i)=mint2
       enddo
       ig(mgr)=int2-(mgr-1)*mint2
      endif

!     broadcast the dimensions to all ranks

          write(lfnout,*) 'neither here, are we?'
          flush(lfnout)

      if(np.gt.1) then
       allocate(idum(mgr+3))
       if(me.eq.master) then
        idum(1)=int1
        idum(2)=int2
        idum(3)=mint2
        do i=1,mgr
         idum(3+i)=ig(i)
        enddo
       endif
       call MPI_Bcast(idum,mgr+3,MPI_INTEGER8,master,
     &                MPI_COMM_WORLD,ierr)
       if(me.ne.master) then
        int1=idum(1)
        int2=idum(2)
        mint2=idum(3)
        do i=1,mgr
         ig(i)=idum(3+i)
        enddo
       endif

       deallocate(idum)

!      broadcast the nuclear potential energy

       call MPI_Bcast(potnuc,1,MPI_REAL8,master,
     & MPI_COMM_WORLD,ierr)

!      broadcast the overlap and one-electron integrals

       if(me.ne.master) then
        allocate(s(nbas,nbas))
        allocate(t(int1))
        allocate(v(int1))
       endif

       call MPI_Bcast(s,nbas*nbas,MPI_REAL8,master,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(t,int1,MPI_REAL8,master,MPI_COMM_WORLD,ierr)
       call MPI_Bcast(v,int1,MPI_REAL8,master,MPI_COMM_WORLD,ierr)

!      read and broadcast the two-electron integrals

       allocate(g(mint2))
       allocate(lab(4,mint2))

       if(idbg.gt.0) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &    ' Allocating g and lab with dimension ',mint2
       endif

!      for single rank groups all two-electron integrals can be read into
!      single array and broadcasted

       if(mgr.eq.1) then
        if(me.eq.master) then
         ilast=0
         ii=0
         do while(ilast.eq.0)
          write(lfnout,*) 'We aren-t here, are we?'
          flush(lfnout)
          read(lfnint,end=993) (b(n),n=1,255),(ibuf(nn),nn=1,1020),
     &                         jint,ilast
          jj=0
          do nn=1,jint
           ii=ii+1
           g(ii)=b(nn)
           do n=1,4
            jj=jj+1
            lab(n,ii)=ibuf(jj)
           enddo
          enddo
         enddo
        endif

        kint=mpibuf
        kintbatch=int2/kint
        if (mod(int2,kint).ne.0)kintbatch=kintbatch+1
        mint=int2
        nint=1
        do jint=1,kintbatch
         lint=min(mint,kint)
#ifdef SINGLEP
         call MPI_Bcast(g(nint),lint,MPI_REAL4,master,
     &                  MPI_COMM_WORLD,ierr)
#else
         call MPI_Bcast(g(nint),lint,MPI_REAL8,master,
     &                  MPI_COMM_WORLD,ierr)
#endif
         call MPI_Bcast(lab(1,nint),4*lint,MPI_INTEGER2,master,
     &                  MPI_COMM_WORLD,ierr)
         nint=nint+lint
         mint=mint-lint
        enddo
       endif

!      for multiple rank groups two electron integrals need to be read in
!      individual rank parts

       if(mgr.gt.1) then
        if(me.eq.master) then
         ilast=0
         ii=0
         igr=1
         do while(ilast.eq.0)
          read(lfnint,end=993) (b(n),n=1,255),(ibuf(nn),nn=1,1020),
     &                         jint,ilast
          jj=0
          do nn=1,jint
           ii=ii+1
           g(ii)=b(nn)
           do n=1,4
            jj=jj+1
            lab(n,ii)=ibuf(jj)
           enddo
           if(ii.eq.ig(igr)) then

!          send to appropriate rank in first group

#ifdef SINGLEP
            call MPI_Send(g,ii,MPI_REAL4,allgroups(1,igr+1),3,
     &       MPI_COMM_WORLD,ierr)
#else
            call MPI_Send(g,ii,MPI_REAL8,allgroups(1,igr+1),3,
     &       MPI_COMM_WORLD,ierr)
#endif
            call MPI_Send(lab,4*ii,MPI_INTEGER2,allgroups(1,igr+1),3,
     &                  MPI_COMM_WORLD,ierr)
           ii=0
           igr=igr+1
           endif
          enddo
         enddo
        else
         igr=0
         do i=1,mgr
          if(me.eq.allgroups(1,i+1)) igr=i
         enddo
         if(igr.gt.0) then
          ii=ig(igr)

#ifdef SINGLEP
          call MPI_Recv(g,ii,MPI_REAL4,master,3,
     &                  MPI_COMM_WORLD,status,ierr)
#else
          call MPI_Recv(g,ii,MPI_REAL8,master,3,
     &                  MPI_COMM_WORLD,status,ierr)
#endif
          call MPI_Recv(lab,4*ii,MPI_INTEGER2,master,3,
     &                  MPI_COMM_WORLD,status,ierr)
         endif
        endif

!      broadcast to corresponding ranks in other groups

        if(numgrp.gt.1) then

         ig0=-1
         do i=1,numgrp
          do j=1,mgr
           if(me.eq.allgroups(i,j+1)) ig0=j
          enddo
         enddo
         if(ig0.ge.0) then

          ngi=ig(ig0)
          kint=mpibuf
          kintbatch=ngi/kint
          if (mod(ngi,kint).ne.0)kintbatch=kintbatch+1
          mint=ngi
          nint=1

          do jint=1,kintbatch
           lint=min(mint,kint)

#ifdef SINGLEP
          call MPI_Bcast(g(nint),lint,MPI_REAL4,0,new_comm(ig0),ierr)
#else
          call MPI_Bcast(g(nint),lint,MPI_REAL8,0,new_comm(ig0),ierr)
#endif
          call MPI_Bcast(lab(1,nint),4*lint,MPI_INTEGER2,0,
     &                  new_comm(ig0),ierr)

           nint=nint+lint
           mint=mint-lint
          enddo

         endif

        endif

       endif

      endif

      close(unit=lfnint,status='keep')

      deallocate(new_comm)

!     Scaling of the one electron integrals

      ielem=0
      jkoff=0

      do j=jkoff+1,jkoff+nbas
       do k=jkoff+1,j
        ielem=ielem+1
       enddo
       t(ielem)=0.5d0*t(ielem)
       v(ielem)=0.5d0*v(ielem)
      enddo

      nacc=min(int(((myints(6)-myints(5))*0.05))+1,naccel)
      nacc0=nacc0+myints(5)
      nacc1=nacc1+myints(5)
      if (naccel.le.0) then
         nacc0=naccel
         nacc1=naccel
      endif

!     Calculate the integral index needed for openacc processing

      allocate(ndxtv(nbas))

      ndxtv(1)=0
      do i=2,nbas
      ndxtv(i)=ndxtv(i-1)+i-1
      enddo

      deallocate(ntarget)

      return
 993  write(lfnout,983) filint(1:inp_strlen(filint))
      call errquit()
 994  write(lfnout,984) filint(1:inp_strlen(filint))
      call errquit()
 983  format('Error reading integral file ',a)
 984  format('Unable to open integral file ',a)
      end
