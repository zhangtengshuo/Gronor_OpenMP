
      subroutine pack(ipack,iocc,n)
      implicit none
      integer :: ipack,n,iocc(n)
      integer :: i,ior
      integer :: mask(-1:2),m2
      data mask/1,0,2,3/
      data m2/2/

      ipack=0
      do i=1,n
       ipack=lshift(ipack,m2)
       ipack=ior(ipack,mask(iocc(i)))
      enddo
      return
      end

      subroutine unpack(ipackin,iocc,n)
      implicit none
      integer :: ipack,n,iocc(n),ipackin
      integer :: i,iand,rshift
      integer :: mask,m2,umask(0:3)
      data mask/3/
      data umask/0,-1,1,2/
      data m2/2/
      ipack=ipackin
      do i=n,1,-1
       iocc(i)=umask(iand(ipack,mask))
       ipack=rshift(ipack,m2)
      enddo
      return
      end
