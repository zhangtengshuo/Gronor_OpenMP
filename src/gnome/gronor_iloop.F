
!> @brief
!! iloop functions
!!     
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  

      subroutine iloop_init(ni,nj,os)
      implicit none
      integer :: ni,nj
      logical os
      integer :: ndeti,ndetj
      logical osame
      common/iloopvar/ndeti,ndetj,osame
      ndeti=ni
      ndetj=nj
      osame=os
      return
      end

      integer function iloop(ij,mode)
      implicit none
      integer ij,mode,ni,nj
      integer :: ndeti,ndetj
      logical osame
      common/iloopvar/ndeti,ndetj,osame
      real(kind = 8) :: p,q
      integer :: l,i,j
      if (osame) then
        p=-0.5+dsqrt(8.0d0*ij+1)/2
        i=p
        l=i*(i+1)/2
        q=p-i
        if (q.gt.1.0d-8) i=i+1
         j=ij-l
        if (j.eq.0) j=i
        ni=i
        nj=j
      else
        ni=ij/ndetj+1
        if (mod(ij,ndetj).eq.0) ni=ni-1
        nj=ij-(ni-1)*ndetj
      endif
      if (mode.eq.1) then
         iloop=ni
      elseif(mode.eq.2) then
         iloop=nj
      else
         stop 'wrong mode for iloop'
      endif
      return
      end
