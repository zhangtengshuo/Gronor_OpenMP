

      subroutine Check_Integrals()

      use gnome_integrals
      use cidef

      implicit none

      integer :: i
      real (kind=8) :: rg(10)

!$acc kernels present(g) copyout(rg)

      do i=1,10
       rg(i)=g(i)
      enddo
!$acc end kernels

      do i=1,10
       write(lfnout,300) i,g(i),rg(i)
 300   format('Check Integrals',i5,2f12.6)
      enddo

      return
      end
