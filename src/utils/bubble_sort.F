
      subroutine bubble_sort(idata,icount,nperm)
      implicit none
      integer :: idata(icount),icount,nperm
      integer :: ipass,isorted,i,itemp
      ipass=1
      nperm=0
 1    continue
      isorted=1
      do i=1,icount-ipass
       if(idata(i) .gt. idata(i+1)) then
        itemp=idata(i)
        idata(i)=idata(i+1)
        idata(i+1)=itemp
        isorted=0
        nperm=nperm+1
       endif
      enddo
      ipass=ipass +1
      if(isorted.eq.0) goto 1
      return
      end
