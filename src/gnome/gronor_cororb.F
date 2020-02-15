
      subroutine gronor_cororb(lfndbg)
      use gnome_parameters
      use gnome_data
      implicit none
      integer :: lfndbg
      integer :: i,j,k

      do i=1,nelecs
       if(ev(i).lt.0.99d0) then
        do j=1,mbasel
         veca(j)=0.0d0
         vecb(j)=0.0d0
        enddo
        do j=1,nveca
         do k=1,mbasel
          veca(k)=veca(k)+u(j,i)*va(j,k)
          vecb(k)=vecb(k)+w(j,i)*vb(j,k)
         enddo
        enddo
       endif
      enddo

      return
      end
