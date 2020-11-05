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
