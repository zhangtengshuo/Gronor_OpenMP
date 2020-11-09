
!> @brief
!! Transform vectors
!!     
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  
      
      subroutine gronor_transvc(lfndbg,idet)
      use cidist
      use gnome_data
      use gnome_parameters
      implicit none
      integer :: idet, m, n, n1, iop, ib, iv, lfndbg
      integer :: k, ke, ks, i, is, ie, mstand, m1, nn, ns, ne, l
      integer :: im, m2, ls,nd
      if(idbg.ge.17) write(lfndbg,601)
 601  format(/,' Determinant Irrep  Nclose  Nopen OccActive',/)

      ntop(idet)=0
      ntcl(idet)=0

      itemp(1)=0
      m=0

      itemp(2)=itemp(1)+nclose(idet)+nopen(idet)
      n1=0
      do iop=1,nopen(idet)
       n1=n1+1
       if(abs(iocopen(n1,idet)).eq.1) then
        m=m+1
        ioccn(m,idet)=iocopen(n1,idet)
       endif
      enddo
      if(idbg.ge.13) write(lfndbg,602) idet,
     & nclose(idet),nopen(idet),
     & (ioccn(i,idet),i=1,nopen(idet))
 602  format(3i8,20i4)

       ntop(idet)=ntop(idet)+nopen(idet)
       ntcl(idet)=ntcl(idet)+nclose(idet)

!     Initialize the output M.O.'s

      do iv=1,mvec
       do ib=1,nbas
        vtemp(iv,ib,idet)=0.0d0
       enddo
      enddo

!     Transform the M.O.'s   (no so much to transform if there is no symmetry)

      do im = 1, itemp(2)           ! number of non-empty orbitals
        do k = 1, nbas              ! number of basis functions
          vtemp(im,k,idet)=vec(im,k,idet)
        end do
      end do
      if(idbg.ge.13) write(lfndbg,604)
 604  format(/,' M.O.''s are transformed')

!     Put transformed M.O's in correct order: closed, alpha, beta
!     First order closed shell M.O.'s of all subspecies

      l=0
      m=0

      do k=1,nclose(idet)
       m=m+1
       l=l+1
       do ib=1,nbas
        vec(m,ib,idet)=vtemp(l,ib,idet)
       enddo
      enddo
      l=l+nopen(idet)

      if(idbg.ge.13) write(lfndbg,606)
 606  format(/,' Closed shells ordered')

!     Then order open alpha shell M.O.'s of all subspecies

      ls=0

      l=ls+nclose(idet)
      do k=1,nopen(idet)
       l=l+1
       if(iocopen(k,idet).eq.1) then
        m=m+1
        do ib=1,nbas
         vec(m,ib,idet)=vtemp(l,ib,idet)
        enddo
       endif
      enddo
      ls=l

      if(idbg.ge.13) write(lfndbg,608)
 608  format(/,' Open alpha shells ordered')

!     Then order open beta shell M.O.'s of all subspecies

      ls=0

      l=ls+nclose(idet)
      do k=1,nopen(idet)
       l=l+1
       if(iocopen(k,idet).eq.-1) then
        m=m+1
        do ib=1,nbas
         vec(m,ib,idet)=vtemp(l,ib,idet)
        enddo
       endif
      enddo
      ls=l

      if(idbg.ge.13) write(lfndbg,610)
 610  format(/,' Open beta shells ordered')

      return
      end
