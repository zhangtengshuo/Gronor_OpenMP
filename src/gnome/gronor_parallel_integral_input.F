
!> @brief
!! Reading integral fiels in parallel mode
!!
!! @author  T. P. Straatsma, ORNL
!! @date    2020
!!
!! Integrals are read in parallel by all ranks the first worker group. 
!! If worker groups have only one rank, one rank is doing the reading of integrals.
!! After the integrals are read, an MPI broadcast copies them to the corresponding ranks in other worker groups.
!!
!! The integrals may be read from multiple files, with or witjout integral labels. If used without labels, the
!! canonical order as defined by OpenMolcas is assumed.
      
      subroutine gronor_parallel_integral_input()

      use mpi
      use inp
      use cidef
      use cidist
      use gnome_data
      use gnome_integrals
      use gnome_parameters

      implicit none

      integer (kind=4) :: ierr
      integer (kind=4) :: nrecb4
      integer :: nrecb
      integer :: i,ie,im,is,ke,km,ks,le,lm,ne,ig0,j,ni,ng,nz
      integer :: nm,ns,kl,k,ll,l,ls,nl,n,nn,igr,ngi,ngir
      integer :: nmax,kbloc,mclab,itarget,noff
      integer :: numfiles,nrecfil,nrectot,irecs(4),myrecs,myrec,myfile
      integer :: nrfil,id,idfil,jdfil,nrec,ilast,ninttot
      integer (kind=8), allocatable :: igrrecs(:,:), ifil(:,:)
      integer :: nrecl
      real (kind=8) :: sdum
      real (kind=8),allocatable :: b(:)
      integer (kind=2),allocatable :: ibuf(:)
      integer :: kint, kintbatch, mint,nint,jint,lint
      integer :: ielem,jkoff,nacc

      integer (kind=8), allocatable :: indxr(:,:,:)
      integer (kind=8) :: numigr,idi,idf,idr,igp,jdi,ii,jj,kk

      if(me.eq.master) then
       rewind(lfnone)

       mclab=0
       if(ione.eq.0) then
         read(lfnone,err=993) namint,intone,potnuc,nbas,mbuf,mclab
       else
        read(lfnone,err=993) namint,nbas,intone,potnuc
      endif
      if(mclab.eq.1) labels=1

       if(ipr.ge.1) then
        write(lfnout,601) namint
        flush(lfnout)
  601   format(/,' Integral file title: ',t25,a,/,t25,a)
        write(lfnout,602) nbas,intone,potnuc
  602   format(/,' Number of basis functions',t50,i16,
     &  /,' Number of one electron integrals',t50,i16,
     &  //,' Nuclear potential energy',t57,f20.10)
       endif
       int1=int1+(nbas*(nbas+1))/2
       if(ipr.ge.1) then
        write(lfnout,603) intone,int1,mclab,labels
 603    format(/,' Number of one electron integrals',t50,i16,
     &  /,' Number of calculated one electron integrals ',t50,i16,
     &  /,' Two electron integral labels option on file ',t50,i16,
     &  /,' Two electron integral labels input option ',t50,i16)
       endif
       allocate(s(nbas,nbas))
       allocate(t(int1))
       allocate(v(int1))
       if(ione.eq.0) then
        read(lfnone,err=993) ((s(i,j),i=1,j),j=1,nbas)
        do i=1,nbas
         do j=1,i
           s(i,j)=s(j,i)
         end do
        end do
        read(lfnone,err=993) (t(i),i=1,int1)
        read(lfnone,err=993) (v(i),i=1,int1)
       else
        call gronor_read1(lfnone,t,intone)
        call gronor_read1(lfnone,t,intone)
        call gronor_read1(lfnone,t,intone)
        call gronor_read1(lfnone,v,intone)
       endif

       read(lfnone,end=993) nt
       if (nt(2).eq.-1) read(lfnone) nt(2)
       int2=nt(2)

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

       read(lfnone,err=993) numfiles

      endif

      call MPI_Bcast(int1,1,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(int2,1,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(root,255,MPI_CHAR,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(mclab,1,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(labels,1,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(numfiles,1,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(mbuf,1,MPI_INTEGER8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(potnuc,1,MPI_REAL8,master,
     & MPI_COMM_WORLD,ierr)

!     Broadcast the overlap and one-electron integrals

      if(me.ne.master) then
       allocate(s(nbas,nbas))
       allocate(t(int1))
       allocate(v(int1))
      endif

      call MPI_Bcast(s,nbas*nbas,MPI_REAL8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t,int1,MPI_REAL8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(v,int1,MPI_REAL8,master,MPI_COMM_WORLD,ierr)

      allocate(ifil(numfiles,2))

      allocate( b(mbuf) )
      allocate( ibuf(4*mbuf) )

      if(me.eq.master) then
       read(lfnone,err=993) (ifil(i,1),i=1,numfiles)
       read(lfnone,err=993) (ifil(i,2),i=1,numfiles)
      endif

      call MPI_Bcast(ifil,2*numfiles,MPI_INTEGER8,master,
     & MPI_COMM_WORLD,ierr)

      nrectot=0
      ninttot=0
      do i=1,numfiles
       nrectot=nrectot+ifil(i,1)
       ninttot=ninttot+ifil(i,2)
      enddo

      do i=1,mgr
       ig(i)=0
      enddo
      igr=0
      do i=1,mgr
       if(me.eq.thisgroup(i+1)) igr=i
      enddo

      if(me.eq.master.or.iamactive.eq.0) then
       allocate(g(1),lab(1,1),ndx(1))
      endif

!     Only the first group reads the integrals from file
      
      if(mygroup.eq.1) then
       myrecs=nrectot/mgr
       myrec=(igr-1)*myrecs+1
       myfile=1
       nrfil=1
       do i=1,numfiles-1
        nrfil=nrfil+ifil(i,1)
        if(myrec.ge.nrfil) myfile=i+1
       enddo
       idfil=myrec
       do i=1,myfile-1
        idfil=idfil-ifil(i,1)
       enddo

       mlab=nbas*(nbas+1)/2

       allocate(indxr(mgr,8,2))
       
       n=0
       nn=0
       idf=1
       idr=1
       idi=1
       igp=1
       numigr=int2/mgr
       indxr(igp,3,1)=1
       indxr(igp,4,1)=1
       do i=1,mlab
         do j=i,mlab
           n=n+1
         enddo
         if(n.ge.numigr.or.i.eq.mlab) then
           nn=nn+n
           n=0
           indxr(igp,3,2)=i
           indxr(igp,4,2)=nn
           igp=igp+1
           if(igp.le.mgr) then
           indxr(igp,3,1)=i+1
           indxr(igp,4,1)=nn+1
           endif
         endif
       enddo
       mint2=indxr(igr,4,2)
       if(igr.gt.1) mint2=indxr(igr,4,2)-indxr(igr-1,4,2)
       
       numigr=int2/mgr

       igp=1
       
       indxr(1,1,1)=1
       indxr(1,1,2)=1
       indxr(1,2,1)=1
       indxr(1,2,2)=ifil(1,1)
       indxr(1,5,1)=1

       n=0
       nn=0
       igp=1
       do i=1,numfiles
         do j=1,ifil(i,1)
           n=n+mbuf
           nn=nn+1
           if(n.gt.indxr(igp,4,2)) then
             indxr(igp,1,2)=i
             indxr(igp,2,2)=nn
             indxr(igp,5,2)=indxr(igp,4,2)-i*(j-1)*mbuf
             igp=igp+1
             if(igp.le.mgr) then
               indxr(igp,1,1)=i
               indxr(igp,2,1)=nn
               indxr(igp,5,1)=indxr(igp-1,5,2)+1
               if(indxr(igp-1,4,2)-i*j*mbuf.eq.0) then
                 indxr(igp,1,1)=i
                 indxr(igp,2,1)=nn+1
                 indxr(igp,5,1)=1
                 if(j.eq.ifil(i,1)) indxr(igp,1,1)=i+1
               endif
             endif
           endif
         enddo
         nn=0
       enddo

       if(labels.eq.0) then
         mint2 = myrecs * mbuf
         if(me.eq.thisgroup(mgr+1)) then
           mint2 = ninttot - (mgr-1) * myrecs * mbuf
           myrecs=nrectot-(mgr-1)*myrecs
         endif
       endif

       inquire(iolength=nrecb4) b(1)
       nrecb = nrecb4
       if(mclab.ne.1) then
         nrecl=(2*mbuf+2)*nrecb
       else
         nrecl=(mbuf+2)*nrecb
       endif
       nrec=0
       
       allocate(g(mint2))
       if(labels.eq.0) then
         allocate(lab(4,mint2),ndx(1))
       else
         mlab=nbas*(nbas+1)/2
         allocate(lab(2,mlab),ndx(mlab))
         kl=0
         do i=1,nbas
           do k=1,i
             kl=kl+1
             lab(1,kl)=k
             lab(2,kl)=i
           enddo
         enddo
       endif

       if(idbg.gt.0) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &    ' Allocating g and lab with dimension ',mint2
       endif
       
       if(idbg.gt.3) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,5i12,4i10)') date(1:8),time(1:8),
     & ' Allocated g and lab for mint2 (1)= ',
     &   mint2,ninttot,int2,
     &   ninttot-(myrec-1)*mbuf,myrecs,myfile,idfil,mbuf
        flush(lfndbg)
       endif

       if(labels.eq.0) then
         mint2=0
         do while(nrec.lt.myrecs)
           nrec=nrec+ifil(myfile,1)-idfil+1
           jdfil=ifil(myfile,1)
           if(nrec.gt.myrecs) jdfil=jdfil-nrec+myrecs
           write(filtwo,'(a,a,i3.3,a)')
     &          root(1:inp_strlen(root)),'_',myfile,'.two'
           open(unit=lfntwo,file=filtwo,access='direct',
     &          form='unformatted',recl=nrecl,status='old')
           if(idbg.gt.0) then
             call swatch(date,time)
             write(lfndbg,'(a,1x,a,1x,a,a)') date(1:8),time(1:8),
     &            ' Opened file ',filtwo(1:inp_strlen(filtwo))
             flush(lfndbg)
           endif
           if(idbg.gt.0) then
             call swatch(date,time)
             write(lfndbg,'(a,1x,a,1x,a,3i10)') date(1:8),time(1:8),
     &            ' Reading records ',idfil,jdfil,jdfil-idfil+1
             flush(lfndbg)
           endif
           flush(lfnout)
           do id=idfil,jdfil
             if(mclab.eq.0) then
               read(lfntwo,rec=id)
     &              (b(n),n=1,mbuf),(ibuf(nn),nn=1,4*mbuf),jint,ilast
             else
               read(lfntwo,rec=id)
     &              (b(n),n=1,mbuf),jint,ilast
               stop 'Error in integral options'
             endif
             j=0
             do i=1,jint
               g(mint2+i)=b(i)
               do n=1,4
                 j=j+1
                 lab(n,mint2+i)=ibuf(j)
               enddo
             enddo
             mint2=mint2+jint
             if(idbg.gt.0) then
               call swatch(date,time)
               write(lfndbg,'(a,1x,a,1x,a,3i10,i5)')
     &              date(1:8),time(1:8),' Read ',id,mint2,jint,ilast
               flush(lfndbg)
             endif
           enddo
           close(unit=lfntwo,status='keep')
           if(idbg.gt.0) then
             call swatch(date,time)
             write(lfndbg,'(a,1x,a,1x,a,a,i5,i10)')
     &            date(1:8),time(1:8),' Closed file ',
     &            filtwo(1:inp_strlen(filtwo)),ilast,jint
             flush(lfndbg)
           endif
           myfile=myfile+1
           idfil=1
         enddo
         
       else

         idi=indxr(igr,5,1)
         jdi=mbuf

         mint2=0

         do idf=indxr(igr,1,1),indxr(igr,1,2)
           
           write(filtwo,'(a,a,i3.3,a)')
     &          root(1:inp_strlen(root)),'_',idf,'.two'

           open(unit=lfntwo,file=filtwo,access='direct',
     &          form='unformatted',recl=nrecl,status='old')

           do idr=indxr(igr,2,1),indxr(igr,2,2)
             
             if(mclab.eq.0) then
               read(lfntwo,rec=idr)
     &              (b(n),n=1,mbuf),(ibuf(nn),nn=1,4*mbuf),jint,ilast
             else
               read(lfntwo,rec=idr)
     &              (b(n),n=1,mbuf),jint,ilast
             endif
             if(idf.eq.indxr(igr,1,2).and.idr.eq.indxr(igr,2,2))
     &            jdi=indxr(igr,5,2)
             do i=idi,jdi
               mint2=mint2+1
               g(mint2)=b(i)
             enddo
             idi=1
           enddo

           close(unit=lfntwo,status='keep')
           
         enddo
         
       endif

         intndx=indxr(igr,3,1)
         jntndx=indxr(igr,3,2)

         deallocate(indxr)
         
         ig(igr)=mint2

         if(idbg.gt.0) then
           call swatch(date,time)
           write(lfndbg,'(a,1x,a,1x,a,i12)') date(1:8),time(1:8),
     &          ' mint2= ',mint2
           flush(lfndbg)
         endif
         
       endif

!     Dimensions need to be broadcasted to all worker ranks
!     if there are more than one group

!     int1  : number of one electron integrals
!     int2  : number of two electron integrals
!     mint2 : number of two electron integrals on this rank

      if(numgrp.gt.1) then

       ig0=-1
       do i=1,numgrp
        do j=1,mgr
         if(me.eq.allgroups(i,j+1)) ig0=j
        enddo
       enddo

       if(ig0.ge.0) then
        call MPI_Bcast(intndx,1,MPI_INTEGER8,0,new_comm(ig0),ierr)
        call MPI_Bcast(jntndx,1,MPI_INTEGER8,0,new_comm(ig0),ierr)
        call MPI_Bcast(mint2,1,MPI_INTEGER8,0,new_comm(ig0),ierr)
        call MPI_Bcast(myrec,1,MPI_INTEGER8,0,new_comm(ig0),ierr)
       endif


       if(idbg.gt.0) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i12)') date(1:8),time(1:8),
     &  ' mint2= ',mint2
        flush(lfndbg)
       endif

       mlab=nbas*(nbas+1)/2

       
       if(mygroup.gt.1) then
         allocate(g(mint2))
         if(labels.eq.0) then
           allocate(lab(4,mint2),ndx(1))
         else
           allocate(lab(2,mlab),ndx(mlab))
           kl=0
           do i=1,nbas
             do k=1,i
               kl=kl+1
               lab(1,kl)=k
               lab(2,kl)=i
             enddo
           enddo
         endif

       if(idbg.gt.0) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &    ' Allocating g and lab with dimension ',mint2
       endif

        if(idbg.gt.3) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,1x,a,5i12,4i10)') date(1:8),time(1:8),
     &  ' Allocated g and lab for mint2 (2)= ',
     &    mint2,ninttot,int2,ninttot-(myrec-1)*mbuf,myrec,0,0,mbuf
         flush(lfndbg)
        endif
       endif

       if(idbg.gt.0) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,3i12)') date(1:8),time(1:8),
     &  ' ig0 (1)= ',ig0
        flush(lfndbg)
       endif

       if(ig0.ge.0) then
        ngi=mint2
        kint=mpibuf
        kintbatch=ngi/kint
        if(mod(ngi,kint).ne.0) kintbatch=kintbatch+1
        mint=ngi
        nint=1

        if(idbg.gt.0) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,1x,a,5i12)') date(1:8),time(1:8),
     &   ' ig0 (2)= ',ig0,ngi,mint,kint,kintbatch
         flush(lfndbg)
        endif

        do jint=1,kintbatch
         lint=min(mint,kint)

#ifdef SINGLEP
         call MPI_Bcast(g(nint),lint,MPI_REAL4,0,new_comm(ig0),ierr)
#else
         call MPI_Bcast(g(nint),lint,MPI_REAL8,0,new_comm(ig0),ierr)
#endif
         if(labels.eq.0) then
           call MPI_Bcast(lab(1,nint),4*lint,MPI_INTEGER2,0,
     &          new_comm(ig0),ierr)
         endif

         nint=nint+lint
         mint=mint-lint
        enddo

        if(labels.eq.1) then
           call MPI_Bcast(lab(1,1),2*mlab,MPI_INTEGER2,0,
     &          new_comm(ig0),ierr)
        endif
        
        if(idbg.gt.0) then
         call swatch(date,time)
         write(lfndbg,'(a,1x,a,1x,a,3i12)') date(1:8),time(1:8),
     &   ' Integrals and labels broadcasted ',mint2,nint,lint
         flush(lfndbg)
        endif

       endif

      endif

      close(unit=lfntwo,status='keep')

      deallocate(new_comm)

!     Scaling of the one electron integrals

      ielem=0
      do j=1,nbas
       do k=1,j
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

!      Need to determine the integral file for this rank based on
!      information read by master from the one electron integral file

!      What is needed:

!       intfil: id for the first integral file
!       intrec: index of the first record to be read from this file
!       maxrec: total number of records on an integral file
!       numrec: number of records to be read on the current rank

      intfil=me+1

!      write(filtwo,'(a,i5.5,a)') root(1:inp_strlen(root)),intfil,'.two'

      ig(igr)=mint2

      if(me.ne.master.and.iamactive.eq.1.and.labels.eq.1) then
        kk=1
        do ii=intndx,jntndx
          ndx(ii)=kk-ii
          do jj=ii,mlab
            kk=kk+1
          enddo
        enddo
      endif
      
      return
 993  write(lfnout,983) filone(1:inp_strlen(filone))
      call errquit()
 994  write(lfnout,984) filone(1:inp_strlen(filone))
      call errquit()
 983  format('Error reading one electron integral file ',a)
 984  format('Unable to open one electron integral file ',a)
 995  write(lfnout,985) filtwo(1:inp_strlen(filtwo))
      call errquit()
 996  write(lfnout,986) filtwo(1:inp_strlen(filtwo))
      call errquit()
 985  format('Error reading two electron integral file ',a)
 986  format('Unable to open two electron integral file ',a)
      end
