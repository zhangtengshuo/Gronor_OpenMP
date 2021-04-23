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
!! Reading integral files in parallel mode
!!
!! @author  T. P. Straatsma, ORNL
!! @date    2020
!!
!! Integrals are read in parallel by all ranks the first worker group. 
!! If worker groups have only one rank, one rank is doing the reading of integrals.
!! After the integrals are read, an MPI broadcast copies them to the corresponding
!! ranks in other worker groups.
!!
!! The integrals may be read from multiple files, with or witjout integral labels.
!! If used without labels, the canonical order as defined by OpenMolcas is assumed.
      
      subroutine gronor_parallel_integral_input()

      use mpi
      use inp
      use cidef
      use cidist
      use gnome_data
      use gnome_integrals
      use gnome_parameters

      implicit none

      integer (kind=4) :: ierr,nrecb4
      integer (kind=8) :: nrecb,i,j,kl,k,n,nn,igr,ngi
      integer (kind=8) :: nrectot,ilast,ninttot,nrecl,ielem
      integer (kind=8), allocatable :: ifil(:,:)
      real (kind=8), allocatable :: b(:)
      integer (kind=2), allocatable :: ibuf(:)
      integer (kind=8), allocatable :: ijbuf(:,:)
      integer :: kint,kintbatch,mint,nint,jint,lint,jbuf(6)

      integer (kind=8) :: numigr,idf,idr,ii,jj,kk,inode

      integer (kind=8) :: idri,idrf,idg,idgi,idgf
      integer (kind=8), allocatable :: ndxf(:,:)
      integer (kind=8), allocatable :: nag(:),nig(:),nng(:)

      inode=map2(me+1,6)

!     Read the one-electron integral file on the master

      if(me.eq.master) then
        open(unit=lfnone,file=filone,form='unformatted',status='old',
     &       err=994)
c     rewind(lfnone)
        mclab=0
        read(lfnone,err=993) namint,intone,potnuc,nbas,mbuf,mclab
        if(mclab.eq.1) labels=1

        if(ipr.gt.0) then
          write(lfnout,601) namint
          flush(lfnout)
 601      format(/,' Integral file title: ',t25,a,/,t25,a)
          write(lfnout,602) potnuc
 602      format(/,' Nuclear potential energy',t57,f20.10)
        endif
        int1=(nbas*(nbas+1))/2
        if(ipr.ge.30) then
          write(lfnout,603) mclab,labels
 603      format(/,
     &         ' Two electron integral labels option on file ',t50,i16,
     &         /,' Two electron integral labels input option ',t50,i16)
        endif
Coen: maybe here a check to be sure that mclab and labels are coherent
        
        allocate(s(nbas,nbas))
        allocate(t(int1))
        allocate(v(int1))
        
        read(lfnone,err=993) ((s(i,j),i=1,j),j=1,nbas)
        do i=1,nbas
          do j=1,i
            s(i,j)=s(j,i)
          end do
        end do
        read(lfnone,err=993) (t(i),i=1,int1)
        read(lfnone,err=993) (v(i),i=1,int1)

Coen: this part is not valid anymore.
C     although we'll loose backward compatibility, we might 
C     want to replace it with just "read(lfnone) int2"
C     The number of small integrals can be written by rdtraint
C     and retrieved here, if considered useful info
        
        read(lfnone,end=993) nt
        if (nt(2).eq.-1) read(lfnone) nt(2)
        int2=nt(2)

#ifdef SINGLEP
        if(ipr.ge.30) write(lfnout,605)
 605    format(' Integrals are used in single precision')
#else
        if(ipr.ge.30) write(lfnout,606)
 606    format(' Integrals are used in double precision')
#endif
        
        read(lfnone,err=993) numfiles
      endif

!     Allocate memory for the overlap and one-electron integrals

      if(me.ne.master) then
        allocate(s(nbas,nbas))
        allocate(t(int1))
        allocate(v(int1))
      endif

!     Broadcast the overlap and one-electron integrals

      call MPI_Bcast(s,nbas*nbas,MPI_REAL8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(t,int1,MPI_REAL8,master,MPI_COMM_WORLD,ierr)
      call MPI_Bcast(v,int1,MPI_REAL8,master,MPI_COMM_WORLD,ierr)

      allocate(ifil(numfiles,2))
      
      if(me.eq.master) then
        read(lfnone,err=993) (ifil(i,1),i=1,numfiles)
        read(lfnone,err=993) (ifil(i,2),i=1,numfiles)
        nrectot=0
        ninttot=0
        do i=1,numfiles
          nrectot=nrectot+ifil(i,1)
          ninttot=ninttot+ifil(i,2)
        enddo
        if(ipr.ge.2) then
          write(lfnout,607)
          do i=1,numfiles
            write(filtwo,'(a,a,i3.3,a)')
     &           root(1:inp_strlen(root)),'_',i,'.two'
            write(lfnout,608) i,ifil(i,1),ifil(i,2),trim(filtwo)
          enddo
          write(lfnout,609) nrectot,ninttot
 607      format(/,' Dimensions on two-electron integral files',//,
     &         ' File     Number Records    Number Integrals'
     &         '     Filename(s)',/)
 608      format(i4,5x,i12,i20,8x,a)
 609      format(11x,'----------',10x,'----------',/,
     &         ' Total   ',i12,i20)
        endif
      endif

      call MPI_Bcast(ifil,2*numfiles,MPI_INTEGER8,master,
     &     MPI_COMM_WORLD,ierr)

      
!     Determine id of rank within its group

      igr=0
      do i=1,mgr
        if(me.eq.thisgroup(i+1)) igr=i
      enddo

!     Determine dimension of super-index to integrals
      
      mlab=nbas*(nbas+1)/2

!     For each rank i within a group determine
!     Initial and final file in              : ndxf(i,1) and ndxf(i,2)
!     Initial and final record in            : ndxf(i,3) and ndxf(i,4)
!     Initial and final index into record in : ndxf(i,5) and ndxf(i,6)

!     Target number of integrals per rank is numigr
      
      allocate(ndxf(mgr,6),nag(mgr),nig(mgr),nng(mgr))
      if(me.eq.master) then
        
        numigr=int2/mgr
        
        k=0
        j=1
        do i=1,mlab
          k=k+mlab-i+1
          if(k.le.j*numigr) then
            nag(j)=k
            nig(j)=i
          else
            j=j+1
          endif
        enddo
        nag(mgr)=int2
        nig(mgr)=mlab
        nng(1)=nag(1)
        do j=2,mgr
          nng(j)=nag(j)-nag(j-1)
        enddo
        if(ipr.ge.20.and.mgr.gt.1) then
          write(lfnout,610)
 610      format(/,' Integrals divided by ranks in a group',/)
          write(lfnout,611) (nig(j),j=1,mgr)
 611      format(' Super index integrals ',10i16)
          write(lfnout,612) (nng(j),j=1,mgr)
 612      format(' Number of integrals   ',10i16)
          write(lfnout,613) (nag(j),j=1,mgr)
 613      format(' Accumulated integrals ',10i16)
        endif

        ndxf(1,1)=1
        ndxf(1,3)=1
        ndxf(1,5)=1
        k=0
        j=1
        do idf=1,numfiles
          do idr=1,ifil(idf,1)
            k=k+1
            if(k*mbuf.ge.nag(j)) then
              ndxf(j,2)=idf
              ndxf(j,4)=idr
              ndxf(j,6)=nag(j)-(k-1)*mbuf
              if(j.lt.mgr) then
                if(k.eq.nag(j)) then
                  ndxf(j+1,5)=1
                  if(idr.eq.ifil(idf,1)) then
                    ndxf(j+1,1)=ndxf(j,2)+1
                    ndxf(j+1,3)=1
                    ndxf(j+1,5)=1
                  else
                    ndxf(j+1,1)=ndxf(j,2)
                    ndxf(j+1,3)=ndxf(j,4)+1
                    ndxf(j+1,5)=1
                  endif
                else
                  ndxf(j+1,1)=ndxf(j,2)
                  ndxf(j+1,3)=ndxf(j,4)
                  ndxf(j+1,5)=ndxf(j,6)+1
                endif
              endif
              j=j+1
            endif
          enddo
        enddo
        if(ipr.ge.40) then
          write(lfnout,614)
 614      format(/,' Index into file buffers by ranks in a group',//,
     &         t12,' ------- Initial -------',
     &         t50,' -------- Final --------',/,
     &         t12,' File   Record     Index',
     &         t50,' File   Record     Index',/)
          do j=1,mgr
            write(lfnout,615) j,ndxf(j,1),ndxf(j,3),ndxf(j,5),
     &           ndxf(j,2),ndxf(j,4),ndxf(j,6)
 615        format(' Rank',i5,t12,i4,2i10,t50,i4,2i10)
          enddo
        endif

      endif

      allocate(ijbuf(mgr,8))

      do j=1,6
        do i=1,mgr
          ijbuf(i,j)=ndxf(i,j)
        enddo
      enddo
      do i=1,mgr
        ijbuf(i,7)=nig(i)
        ijbuf(i,8)=nng(i)
      enddo
     
      call MPI_Bcast(ijbuf,8*mgr,MPI_INTEGER8,master,
     &     MPI_COMM_WORLD,ierr)
      
      do j=1,6
        do i=1,mgr
          ndxf(i,j)=ijbuf(i,j)
        enddo
      enddo
      do i=1,mgr
        nig(i)=ijbuf(i,7)
        nng(i)=ijbuf(i,8)
      enddo
      
      deallocate(ijbuf)

!     Determine the correct record length of the integral files
      
      inquire(iolength=nrecb4) b(1)
      nrecb = nrecb4
      if(mclab.ne.1) then
        nrecl=(2*mbuf+2)*nrecb
      else
        nrecl=(mbuf+2)*nrecb
      endif
      

      if(igr.gt.0) mint2=nng(igr)
      
!     Allocate the arrays to hold integrals and labels

      if(me.eq.master.or.iamactive.eq.0) then
        allocate(g(1),lab(1,1),ndx(1))
      else
        allocate(g(mint2))
        if(labels.eq.0) then
          allocate(lab(4,mint2),ndx(1))
        else
          mlab=nbas*(nbas+1)/2
          allocate(lab(2,mlab),ndx(mlab))
        endif
      endif
      
!     Only the first group reads the integrals from file

      if(mygroup.eq.1) then
!     Allocate the buffers needed to read the two-electron integral files
        allocate(b(mbuf),ibuf(4*mbuf))
        k=0
        do idf=ndxf(igr,1),ndxf(igr,2)
          idri=ndxf(igr,3)
          idrf=ifil(idf,1)
          if(idf.ne.ndxf(igr,1)) idri=1
          if(idf.eq.ndxf(igr,2)) idrf=ndxf(igr,4)
          write(filtwo,'(a,a,i3.3,a)')
     &         root(1:inp_strlen(root)),'_',idf,'.two'
          open(unit=lfntwo,file=filtwo,access='direct',
     &         form='unformatted',recl=nrecl,status='old')
          do idr=idri,idrf
            if(mclab.eq.0) then
              read(lfntwo,rec=idr)
     &             (b(n),n=1,mbuf),(ibuf(nn),nn=1,4*mbuf),jint,ilast
            else
              read(lfntwo,rec=idr)
     &             (b(n),n=1,mbuf),jint,ilast
            endif
            idgi=ndxf(igr,5)
            idgf=mbuf
            if(idf.ne.ndxf(igr,1).or.idr.ne.idri) idgi=1
            if(idf.eq.ndxf(igr,2).and.idr.eq.idrf) idgf=ndxf(igr,6)
            do idg=idgi,idgf
              k=k+1
              g(k)=b(idg)
              if(labels.eq.0) then
                do n=1,4
                  lab(n,k)=ibuf((k-1)*4+n)
                enddo
              endif
            enddo

          enddo  
          close(unit=lfntwo,status='keep')
        enddo
        deallocate(b,ibuf)
        if(labels.ne.0) then
          kl=0
          do i=1,nbas
            do k=1,i
              kl=kl+1
              lab(1,kl)=k
              lab(2,kl)=i
            enddo
          enddo
        endif
      endif

      intndx=1
      jntndx=mlab
      if(igr.gt.1) intndx=nig(igr-1)+1
      if(igr.lt.mgr) jntndx=nig(igr)
      
!     Distribute integrals and labels if more than one group

      if(numgrp.gt.1) then
        if(igr.ge.1) then
          ngi=mint2
          kint=mpibuf
          kintbatch=ngi/kint
          if(mod(ngi,kint).ne.0) kintbatch=kintbatch+1
          mint=ngi
          nint=1
          do jint=1,kintbatch
            lint=min(mint,kint)
            if(itest.eq.0) then
#ifdef SINGLEP
              call MPI_Bcast(g(nint),lint,MPI_REAL4,0,
     &             new_comm(igr),ierr)
#else
              call MPI_Bcast(g(nint),lint,MPI_REAL8,0,
     &             new_comm(igr),ierr)
#endif
            else
#ifdef SINGLEP
              if(lcomm1.and.nnodes.gt.1)
     &             call MPI_Bcast(g(nint),lint,MPI_REAL4,0,
     &             new_comm1(igr),ierr)
              call MPI_Bcast(g(nint),lint,MPI_REAL4,0,
     &             new_comm2(igr,inode),ierr)
              
#else
              if(lcomm1.and.nnodes.gt.1)
     &             call MPI_Bcast(g(nint),lint,MPI_REAL8,0,
     &             new_comm1(igr),ierr)
              call MPI_Bcast(g(nint),lint,MPI_REAL8,0,
     &             new_comm2(igr,inode),ierr)
#endif
            endif
            if(labels.eq.0) then
              call MPI_Bcast(lab(1,nint),4*lint,MPI_INTEGER2,0,
     &             new_comm(igr),ierr)
            endif

            if(idbg.gt.10) then
              call swatch(date,time)
              write(lfndbg,130) date(1:8),time(1:8),jint
 130          format(a,1x,a,1x,' Integrals broadcasted for batch ',i5)
              flush(lfndbg)
            endif

            nint=nint+lint
            mint=mint-lint
          enddo

          if(labels.eq.1) then
            call MPI_Bcast(lab(1,1),2*mlab,MPI_INTEGER2,0,
     &           new_comm(igr),ierr)
          endif

        endif
      endif

!     Calculate the integral index needed for openacc processing

      allocate(ndxtv(nbas))

      ndxtv(1)=0
      do i=2,nbas
        ndxtv(i)=ndxtv(i-1)+i-1
      enddo

      if(me.ne.master.and.iamactive.eq.1.and.labels.eq.1) then
        kk=1
        do ii=intndx,jntndx
          ndx(ii)=kk-ii
          do jj=ii,mlab
            kk=kk+1
          enddo
        enddo
      endif

!     Scaling of the one electron integrals

      ielem=0
      do j=1,nbas
        do k=1,j
          ielem=ielem+1
        enddo
        t(ielem)=0.5d0*t(ielem)
        v(ielem)=0.5d0*v(ielem)       
      enddo
      
      flush(lfnout)
      
      return
      
 993  write(lfnout,983) filone(1:inp_strlen(filone))
      call errquit(me,4,"Error in reading integrals")
 994  write(lfnout,984) filone(1:inp_strlen(filone))
      call errquit(me,4,"Error in reading integrals")
 983  format('Error reading one electron integral file ',a)
 984  format('Unable to open one electron integral file ',a)
 995  write(lfnout,985) filtwo(1:inp_strlen(filtwo))
      call errquit(me,4,"Error in reading integrals")
 996  write(lfnout,986) filtwo(1:inp_strlen(filtwo))
      call errquit(me,4,"Error in reading integrals")
 985  format('Error reading two electron integral file ',a)
 986  format('Unable to open two electron integral file ',a)
      end
