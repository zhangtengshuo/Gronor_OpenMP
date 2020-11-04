      subroutine productf
     &     (civm,ioccm,nactm,idetm,nmol,ncomb,civb,ioccb,
     &     maxci,mstates,maxnact,mactb,maxcib,ndetout,thresh_CI,lfnout)
      implicit none
      integer :: nmol,maxci,mstates,maxnact,mactb,maxcib
      real (kind=8) :: civm(maxci,mstates)
      integer :: ioccm(maxnact,maxci,mstates)
      integer :: nactm(mstates)
      integer :: idetm(mstates)
      real (kind=8) :: civb(maxcib),ci2,thresh_CI
      integer :: ioccb(maxcib),ncomb(nmol)
      integer, allocatable :: kloop(:),jdet(:)
      integer :: ndets,i,j,jact,k,isetsign4,iocch(32),ndetout
      integer :: lfnout
      integer, external :: setsign_mebf

      ndets=1
      allocate(kloop(nmol),jdet(nmol))

      do i=1,nmol
       ndets=ndets*idetm(ncomb(i))
       jdet(i)=idetm(ncomb(i))
      enddo
      if (ndets.gt.maxcib) then
         print *,'error ndets',ndets,maxcib
         stop 'ndets error'
      endif

      ndetout=0
      do i=1,ndets
       call calciloop(i,kloop,nmol,jdet)
       ci2=1.0d0
       do j=1,nmol
        ci2=ci2*civm(kloop(j),ncomb(j))
       enddo
       if (dabs(ci2).gt.thresh_CI) then
         ndetout = ndetout + 1
         civb(ndetout)=ci2
         jact=1
         do j=1,nmol
            do k=1,nactm(ncomb(j))
             iocch(jact)=ioccm(k,kloop(j),ncomb(j))
             jact=jact+1
            enddo
         enddo
         jact=jact-1
         if(nmol.gt.1)then
            civb(ndetout)=
     &         setsign_mebf(iocch,nactm(ncomb(1)),nactm(ncomb(2)),jact)*
     &         civb(ndetout)
         endif
         call pack(ioccb(ndetout),iocch,jact)
       endif
      enddo

      deallocate(kloop,jdet)

      return
      end

      integer function isetsign4(jocc,nact,lfnout)
      implicit none
      integer :: lfnout
      integer :: nact
      integer :: jocc(nact),iocca(32),ioccb(32)
      integer :: iclosed(32)
      integer, allocatable :: idet(:,:)
      integer :: ia,ib,ibita,ibitb,is,kk,jj,j,i,k,jperm
      integer :: nelec,icl,nclosed,index,nperm,ih1,ih2,na,nb,ialp,ih
      integer, external :: ibitset,ibitc

      isetsign4=0
      ia=0
      ib=0
      ibita=0
      ibitb=0
      iclosed=0
      is=0

      do i=1,nact
       iocca(i)=0
       ioccb(i)=0
      enddo

      do i=1,nact
       iclosed(i)=0
      enddo

      do j=1,nact
       if(jocc(j).eq.2) then
        ia=ia+1
        iocca(ia)=j
        ib=ib+1
        ioccb(ib)=j
       elseif(jocc(j).eq.-1) then
        ib=ib+1
        ioccb(ib)=j
       elseif(jocc(j).eq.1) then
        ia=ia+1
        iocca(ia)=j
       endif
      enddo

      ibita=ibitset(iocca,ia)
      ibitb=ibitset(ioccb,ib)
      nclosed=ibitc(iand(ibita,ibitb))
      call seticlosed(iclosed,iand(ibita,ibitb),nclosed)
      nperm=0
      icl=1
      do i=1,nclosed
       do j=1,ia
        if(iocca(j).eq.iclosed(i)) then
         if(icl.ne.j) then
          ih=iocca(j)
          iocca(j)=iocca(icl)
          iocca(icl)=ih
          nperm=nperm+1
         endif
         icl=icl+1
         goto 101
        endif
       enddo
 101   continue
      enddo

      icl=1

      do i=1,nclosed
       do j=1,ib
        if(ioccb(j).eq.iclosed(i)) then
         if(icl.ne.j) then
          ih=ioccb(j)
          ioccb(j)=ioccb(icl)
          ioccb(icl)=ih
          nperm=nperm+1
         endif
         icl=icl+1
         goto 102
        endif
       enddo
 102   continue
      enddo

      na=ia-nclosed
      nb=ib-nclosed
      call bubble_sort(iocca(nclosed+1),na,jperm)
      nperm=nperm+jperm
      call bubble_sort(ioccb(nclosed+1),nb,jperm)
      nperm=nperm+jperm
      
      if(mod(nperm,2).eq.1)then
       isetsign4=-1
      else
       isetsign4=1
      endif
      return
      end

      subroutine calciloop(ind,kloop,nmol,jdet)
      implicit none
      integer :: nmol, kloop(nmol), jdet(nmol)
      integer :: ind, index, ihelp, i, k

      index=ind

      do i=1,nmol
       ihelp=1
       do k=i+1,nmol
        ihelp=ihelp*jdet(k)
       enddo
       kloop(i)=(index-1)/ihelp+1
       index=index-(kloop(i)-1)*ihelp
      enddo

      return
      end

      integer function ifind(nnn,iocc,nact)
      implicit none
      integer :: iocc(nact)
      integer :: i,nact,nnn
      ifind=0
      do i=1,nact
       if(iocc(i).eq.nnn) then
        ifind=i
        return
       endif
      enddo
      stop 'I should not be here!'
      return
      end

      subroutine seticlosed(iclosed,in,nclosed)
      implicit none
      integer :: nclosed,in,iclosed(nclosed)
      integer :: mask,jj,i,mask1
      mask=1
      jj=0
      do i=0,31
       mask1=ishft(mask,i)
       if(iand(in,mask1).ne.0) then
        jj=jj+1
        iclosed(jj)=i+1
       endif
      enddo
      if(jj.ne.nclosed)stop 'you have to learn how to count'
      return
      end
      
      subroutine old_productf
     &     (civm,ioccm,nactm,idetm,nmol,ncomb,civb,ioccb,
     &     maxci,mstates,maxnact,mactb,maxcib,ndetout,thresh_CI,
     &     lfnout)
      implicit none
      integer :: nmol,maxci,mstates,maxnact,mactb,maxcib
      real (kind=8) :: civm(maxci,mstates)
      integer :: ioccm(maxnact,maxci,mstates)
      integer :: nactm(mstates)
      integer :: idetm(mstates)
      real (kind=8) :: civb(maxcib),ci2,thresh_CI
      integer :: ioccb(maxcib),ncomb(nmol)
      integer, allocatable :: kloop(:),jdet(:)
      integer :: ndets,i,j,jact,k,isetsign4,iocch(32),ndetout
      integer :: lfnout
      external :: isetsign4

      ndets=1
      allocate(kloop(nmol),jdet(nmol))

      do i=1,nmol
       ndets=ndets*idetm(ncomb(i))
       jdet(i)=idetm(ncomb(i))
      enddo
      if (ndets.gt.maxcib) then
         print *,'error ndets',ndets,maxcib
         stop 'ndets error'
      endif

      ndetout=0
      do i=1,ndets
       call calciloop(i,kloop,nmol,jdet)
       ci2=1.0d0
       do j=1,nmol
        ci2=ci2*civm(kloop(j),ncomb(j))
       enddo
       if (dabs(ci2).gt.thresh_CI) then
         ndetout = ndetout + 1
         civb(ndetout)=ci2
         jact=1
         do j=1,nmol
            do k=1,nactm(ncomb(j))
             iocch(jact)=ioccm(k,kloop(j),ncomb(j))
             jact=jact+1
            enddo
         enddo
         jact=jact-1
         write(lfnout,*)'ndet',ndetout
         civb(ndetout)=isetsign4(iocch,jact,lfnout)*civb(ndetout)
         call pack(ioccb(ndetout),iocch,jact)
       endif
      enddo
      deallocate(kloop,jdet)
      return
      end

      subroutine productf2
     &     (civm,ioccm,nactm,idetm,nmol,ncomb,civb,ioccb,
     &     maxci,mstates,maxnact,mactb,maxcib,ndetout,thresh_CI,lfnout)
      implicit none
      integer :: nmol,maxci,mstates,maxnact,mactb,maxcib
      real (kind=8) :: civm(maxci,mstates)
      integer :: ioccm(maxnact,maxci,mstates)
      integer :: nactm(mstates)
      integer :: idetm(mstates)
      real (kind=8) :: civb(maxcib),ci2,thresh_CI
      integer :: ioccb(maxcib),ncomb(nmol)
      integer, allocatable :: kloop(:),jdet(:)
      integer :: ndets,i,j,jact,k,iocch(32),ndetout
      integer :: lfnout
      integer, external :: perm_ab
      integer, external :: isetsign4
      ndets=1
      allocate(kloop(nmol),jdet(nmol))
      do i=1,nmol
         ndets=ndets*idetm(ncomb(i))         
         jdet(i)=idetm(ncomb(i))
      enddo
      if (ndets.gt.maxcib) then
         write(lfnout,*)'error ndets',ndets,maxcib
         flush(lfnout)
         stop 'ndets error'
      endif
      ndetout=0
      do i=1,ndets
         call calciloop(i,kloop,nmol,jdet)
         ci2=1.0d0
         do j=1,nmol
            ci2=ci2*civm(kloop(j),ncomb(j))
         enddo
         if (dabs(ci2).gt.thresh_CI) then
            ndetout = ndetout + 1
            civb(ndetout)=ci2
            jact=1
            do j=1,nmol
               
               do k=1,nactm(ncomb(j))
                  iocch(jact)=ioccm(k,kloop(j),ncomb(j))
                  
                  jact=jact+1
               enddo
            enddo            
            jact=jact-1                        
            civb(ndetout)=civb(ndetout)*perm_ab(iocch,jact,lfnout)  
            civb(ndetout)=isetsign4(iocch,jact,lfnout)*civb(ndetout)       
            call pack(ioccb(ndetout),iocch,jact)
         endif
      enddo      
         
      deallocate(kloop,jdet)
      return
      end

