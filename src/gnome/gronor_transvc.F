
      subroutine gronor_transvc(lfndbg,idet)
      use cidist
      use gnome_data
      use gnome_parameters
      implicit none
      integer :: idet, m, n, n1, iop, ib, iv, lfndbg
      integer :: k, ke, ks, i, is, ie, mstand, m1, nn, ns, ne, l
      integer :: im, m2, ls,nd
      integer :: nsymrp
      if(idbg.ge.17) write(lfndbg,601)
 601  format(/,' Determinant Irrep  Nclose  Nopen OccActive',/)

      ntop(idet)=0
      ntcl(idet)=0

      nsymrp = 1     ! Emergency solution, should be changed

      itemp(1)=0
      m=0
      do n=1,nsymrp
       itemp(n+1)=itemp(n)+nclose(n,idet)+nopen(n,idet)
       n1=0
       do iop=1,nopen(n,idet)
        n1=n1+1
        if(abs(iocopen(n,n1,idet)).eq.1) then
         m=m+1
         ioccn(m,idet)=iocopen(n,n1,idet)
        endif
       enddo
       if(idbg.ge.13) write(lfndbg,602) idet,n,
     & nclose(n,idet),nopen(n,idet),
     & (ioccn(i,idet),i=1,nopen(n,idet))
 602   format(4i8,20i4)

       ntop(idet)=ntop(idet)+nopen(n,idet)
       ntcl(idet)=ntcl(idet)+nclose(n,idet)

      enddo

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
      do n=1,nsymrp
       do k=1,nclose(n,idet)
        m=m+1
        l=l+1
        do ib=1,nbas
         vec(m,ib,idet)=vtemp(l,ib,idet)
        enddo
       enddo
       l=l+nopen(n,idet)
      enddo
      if(idbg.ge.13) write(lfndbg,606)
 606  format(/,' Closed shells ordered')

!     Then order open alpha shell M.O.'s of all subspecies

      ls=0
      do n=1,nsymrp
       l=ls+nclose(n,idet)
       do k=1,nopen(n,idet)
        l=l+1
        if(iocopen(n,k,idet).eq.1) then
         m=m+1
         do ib=1,nbas
          vec(m,ib,idet)=vtemp(l,ib,idet)
         enddo
        endif
       enddo
       ls=l
      enddo
      if(idbg.ge.13) write(lfndbg,608)
 608  format(/,' Open alpha shells ordered')

!     Then order open beta shell M.O.'s of all subspecies

      ls=0
      do n=1,nsymrp
       l=ls+nclose(n,idet)
       do k=1,nopen(n,idet)
        l=l+1
        if(iocopen(n,k,idet).eq.-1) then
         m=m+1
         do ib=1,nbas
          vec(m,ib,idet)=vtemp(l,ib,idet)
         enddo
        endif
       enddo
       ls=l
      enddo

      if(idbg.ge.13) write(lfndbg,610)
 610  format(/,' Open beta shells ordered')

      return
      end
