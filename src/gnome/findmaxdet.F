
      subroutine findmaxdet(ioccm,maxnact,maxci,mstates,nactm,ifac)

      implicit none

      integer :: maxnact,maxci,mstates,ifac
      integer :: ioccm(maxnact,maxci,mstates),nactm(mstates)

      integer :: i,j,na,nb,ni,ndif

      ndif=0
      do i=1,mstates
       na=0
       nb=0
       do j=1,nactm(i)
        if(ioccm(j,1,i).eq.2) then
         na=na+1
         nb=nb+1
        elseif(ioccm(j,1,i).eq.1) then
         na=na+1
        elseif(ioccm(j,1,i).eq.-1) then
         nb=nb+1
        endif
       enddo
       if(na.ne.nb) then
        ni=nactm(i)-na
        if(ni.lt.nb) then
         nb=nb-ni
        else
         nb=0
        endif
        ndif=max(ndif,na-nb)
       endif
      enddo

      if(ndif.eq.0) then
       ifac=1
      else
       ifac=ndif
      endif

      return
      end

      integer function ibitc(ibit)
      implicit none
      integer :: ibit,i,mask
      ibitc=0
      do i=0,31
       mask=1
       mask=lshift(mask,i)
       if(iand(mask,ibit).ne.0) ibitc=ibitc+1
      enddo
      return
      end

      integer function ibitset(iocc,n)
      implicit none
      integer :: n,iocc(n)
      integer :: mask,i
      ibitset=0
      do i=1,n
       mask=1
       mask=lshift(mask,iocc(i)-1)
       ibitset=ior(ibitset,mask)
      enddo
      return
      end
