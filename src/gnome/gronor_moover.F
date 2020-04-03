      subroutine gronor_moover(lfndbg)

      use mpi
      use cidist
      use gnome_integrals
      use gnome_parameters
      use gnome_data

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg
      integer :: i, nopala, nopalb, nalfab, i1, i2
      integer :: ict,j,ib,kb,iv,ie,ke,le,kk,ii,k,l,m1,k1
      real (kind=8) :: sint, si, sum

      if(idbg.ge.13) write(lfndbg,601)
 601  format(/,' Calculation of the overlap matrix')

!     calculate number of open shells with alpha spin

      nopala=0
      nopalb=0
      do i=1,ntopa
       if(ioccn(i,1).eq.1) nopala=nopala+1
      enddo
      do i=1,ntopb
       if(ioccn(i,2).eq.1) nopalb=nopalb+1
      enddo

      nalfa=ntcla+nopala
      nalfab=ntclb+nopalb
      if(nalfa.ne.nalfab) stop'i don''t know spin-flip operators'

      if(idbg.ge.14) write(lfndbg,602)
 602  format('  itypen icentn jtype jcent   icount      indbas     '
     &        ,' indbas         overlap',//)

!$acc data copyin(s,va,vb) copyout(ta,tb)

!$acc kernels present(ta,tb,s,vb)

      do iv=1,nvecb
       do ib=1,nbas
        tb(ib,iv)=0.0d0
       enddo
      enddo

       do iv=1,nvecb
        do ib=1,nbas
         sum=0.0d0
         do kb=1,nbas
           sum=sum+s(ib,kb)*vb(iv,kb)
         enddo
         tb(ib,iv)=tb(ib,iv)+sum
        enddo
       enddo

!     calculation of a **  a(i,k)=sigma(l)[ va(i,l)*b(l,k) ]
!     calculation for alpha spin in open- and closed shell-m.o.'s

      do le=1,nelecs
       do ie=1,nelecs
        ta(ie,le)=0.0d0
       enddo
      enddo

!$acc end kernels

      if(nalfa.ne.0) then
        
!$acc kernels present(va,tb)
        
       do ke=1,nalfa
        do ie=1,nalfa
         sum=0.0d0
         do ib=1,nbas
          sum=sum+va(ie,ib)*tb(ib,ke)
         enddo
         ta(ie,ke)=sum
        enddo
       enddo

!$acc end kernels

      endif
      
      
!     generate overlapmatrix-elements for beta spin in closed-shell m.o.

      if(ntcla.ne.0 .and. ntclb.ne.0) then
       
!$acc parallel loop present(ta)
        do i=1,ntcla
          ii=i+nalfa
!$acc loop vector    
        do k=1,ntclb
         kk=k+nalfa
         ta(ii,kk)=ta(i,k)
        enddo
       enddo
      
      endif

      
      m1=nalfa+1

!     calculation for beta spin in open-shell-m.o.'s

      if(nveca.ne.nalfa) then

        if(ntclb.gt.0) then
!$acc kernels present(va,ta,tb)
          do k=1,ntclb
            kk=k+nalfa
            do i=m1,nveca
              sum=0.0d0
              do l=1,nbas
                sum=sum+va(i,l)*tb(l,k)
              enddo
              ta(i+ntcla,kk)=sum
            enddo
          enddo
!$acc end kernels
        endif

        if(nvecb.gt.nalfa) then
!$acc kernels present(va,ta,tb)
          do k=nalfa+1,nvecb
            kk=k+ntclb
            do i=m1,nveca
              sum=0.0d0
              do l=1,nbas
                sum=sum+va(i,l)*tb(l,k)
              enddo
              ta(i+ntcla,kk)=sum
            enddo
          enddo
!$acc end kernels
        endif

      endif

      
      if(nvecb.ne.nalfa .and. ntcla.ne.0) then
!$acc kernels present(va,ta,tb)
       do i=1,ntcla
        do k=m1,nvecb
         sum=0.0d0
         do l=1,nbas
          sum=sum+va(i,l)*tb(l,k)
         enddo
         ta(i+nalfa,k+ntclb)=sum
        enddo
       enddo
!$acc end kernels
      endif

!     There are problems with svd if many off diagonal elements occur
!     of size 1e-13. therefore all elements < 1e-10 are set to zero

!$acc kernels present(ta)
      do i1=1,nelecs
       do i2=1,nelecs
        if(abs(ta(i1,i2)).lt.1.0d-10) ta(i1,i2)=0.0d0
       enddo
      enddo
!$acc end kernels

!$acc end data
      
      return
      end
