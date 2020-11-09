
!> @brief
!! Transformation of MOs
!!     
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  
      subroutine gronor_tramat2(lfndbg)

!      Transformation of the  m.o.'s
!      the new  m.o.'s are adapted to the basis set of the two electon
!      integrals

      use cidist
      use gnome_parameters
      use gnome_data
      implicit none
      integer :: lfndbg

      integer :: i,j,k,kk,m1
      real (kind=8) :: sum

      if(idbg.ge.13) write(lfndbg,600) 
 600  format(/,' Cofactor matrix transform to symmetry functions')

!$acc kernels present(va,vb,diag,bdiag,cdiag,bsdiag,csdiag,ta,aaa,w1,w2)

      if(ising .ne. 0.) then


!     transformation of diag and sdiag

       if(nalfa.ne.0) then
*        do j=1,nsymf
        do j=1,nbas 
         sum=0.0d0
!$acc loop reduction(+:sum)
         do k=1,nalfa
          sum=sum+cdiag(k)*va(k,j)
         enddo
         diag(j)=sum
        enddo
       endif

       if(ntcla.ne.0) then
*        do j=1,nsymf
        do j=1,nbas 
         sum=0.0d0
!$acc loop reduction(+:sum)
         do k=1,ntcla
          kk=k+nalfa
          sum=sum+cdiag(kk)*va(k,j)
         enddo
         bdiag(j)=sum
        enddo
       endif

       if(nalfa.ne.nveca) then
        m1=nalfa+1
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
!$acc loop reduction(+:sum)
         do k=m1,nveca
          kk=k+ntcla
          sum=sum+cdiag(kk)*va(k,j)
         enddo
         bdiag(j)=bdiag(j)+sum
        enddo
       endif

       if(nalfa.ne.0) then
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
!$acc loop reduction(+:sum)
         do k=1,nalfa
          sum=sum+csdiag(k)*vb(k,j)
         enddo
         w1(j)=sum
        enddo
       endif

       if(ntclb.ne.0) then
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
!$acc loop reduction(+:sum)
         do k=1,ntclb
          kk=k+nalfa
          sum=sum+csdiag(kk)*vb(k,j)
         enddo
         bsdiag(j)=sum
        enddo
       endif

       if(nalfa.ne.nvecb) then
        m1=nalfa+1
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
!$acc loop reduction(+:sum)
         do k=m1,nvecb
          kk=k+ntclb
          sum=sum+csdiag(kk)*vb(k,j)
         enddo
         bsdiag(j)=bsdiag(j)+sum
        enddo
       endif

*        do j=1,nsymf
        do j=1,nbas
        csdiag(j)=w1(j)
       enddo
      endif

!     transformation of the input-matrix a
!                  first:
!     calculation of the intermediate matrices

      if(nalfa.ne.0) then
       do i=1,nalfa
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
!$acc loop reduction(+:sum)
         do k=1,nalfa
          sum=sum+ta(i,k)*vb(k,j)
         enddo
         w2(i,j)=sum
        enddo
       enddo
*        do j=1,nsymf
        do j=1,nbas
        do i=1,nalfa
         ta(i,j)=w2(i,j)
        enddo
       enddo
      endif

      m1=nalfa+1
      if(nalfa.ne.nelecs) then
       do i=m1,nelecs
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
         if(ntclb.ne.0) then
!$acc loop reduction(+:sum)
          do k=1,ntclb
           kk=k+nalfa
           sum=sum+ta(i,kk)*vb(k,j)
          enddo
         endif
         if(nalfa.ne.nvecb) then
!$acc loop reduction(+:sum)
          do k=m1,nvecb
           kk=k+ntclb
           sum=sum+ta(i,kk)*vb(k,j)
          enddo
         endif
         w2(i,j)=sum
        enddo
       enddo
*        do j=1,nsymf
        do j=1,nbas
        do i=m1,nelecs
         ta(i,j)=w2(i,j)
        enddo
       enddo
      endif

!     calculation of the two final matrices

      if(nalfa.ne.0) then
*       do j=1,nsymf
*        do i=1,nsymf
       do j=1,nbas
        do i=1,nbas
         sum=0.0
!$acc loop reduction(+:sum)
         do k=1,nalfa
          sum=sum+ta(k,j)*va(k,i)
         enddo
         aaa(i,j)=sum
        enddo
       enddo
      endif

      if(nalfa.ne.nelecs) then
*       do j=1,nsymf
*        do i=1,nsymf
       do j=1,nbas
        do i=1,nbas
         sum=0.0d0
         if(ntcla.ne.0) then
!$acc loop reduction(+:sum)
          do k=1,ntcla
           kk=k+nalfa
           sum=sum+ta(kk,j)*va(k,i)
          enddo
         endif
         if(nalfa.ne.nveca) then
!$acc loop reduction(+:sum)
          do k=m1,nveca
           kk=k+ntcla
           sum=sum+ta(kk,j)*va(k,i)
          enddo
         endif
         w2(i,j)=sum
        enddo
       enddo
*       do j=1,nsymf
*        do i=1,nsymf
       do j=1,nbas
        do i=1,nbas
         ta(i,j)=w2(i,j)
        enddo
       enddo
      else
*       do j=1,nsymf
*        do i=1,nsymf
       do j=1,nbas
        do i=1,nbas
         ta(i,j)=0.0d0
        enddo
       enddo
      endif

!$acc end kernels

c!$acc update host(ta,w2)

c        write(*,'(a,i5)') 'TA tramat2 '
c        do k=1,min(5,nbas)
c          write(*,'(10f12.6)') (ta(i,k),i=1,nbas)
c        enddo
c        write(*,'(a,i5)') 'W2 tramat'
c        do k=1,min(5,nbas)
c          write(*,'(10f12.6)') (w2(i,k),i=1,nbas)
c        enddo

      return
      end

      subroutine gronor_tramat2_omp(lfndbg)

!      Transformation of the  m.o.'s
!      the new  m.o.'s are adapted to the basis set of the two electon
!      integrals

      use cidist
      use gnome_parameters
      use gnome_data
      implicit none
      integer :: lfndbg

      integer :: i,j,k,kk,m1
      real (kind=8) :: sum

      if(idbg.ge.13) write(lfndbg,600) 
 600  format(/,' Cofactor matrix transform to symmetry functions')

cc!$acc kernels present(va,vb,diag,bdiag,cdiag,bsdiag,csdiag,ta,aaa,w1,w2)

      if(ising .ne. 0.) then


!     transformation of diag and sdiag

       if(nalfa.ne.0) then
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
cc!$acc loop reduction(+:sum)
         do k=1,nalfa
          sum=sum+cdiag(k)*va(k,j)
         enddo
         diag(j)=sum
        enddo
       endif

       if(ntcla.ne.0) then
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
cc!$acc loop reduction(+:sum)
         do k=1,ntcla
          kk=k+nalfa
          sum=sum+cdiag(kk)*va(k,j)
         enddo
         bdiag(j)=sum
        enddo
       endif

       if(nalfa.ne.nveca) then
        m1=nalfa+1
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
cc!$acc loop reduction(+:sum)
         do k=m1,nveca
          kk=k+ntcla
          sum=sum+cdiag(kk)*va(k,j)
         enddo
         bdiag(j)=bdiag(j)+sum
        enddo
       endif

       if(nalfa.ne.0) then
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
!$acc loop reduction(+:sum)
         do k=1,nalfa
          sum=sum+csdiag(k)*vb(k,j)
         enddo
         w1(j)=sum
        enddo
       endif

       if(ntclb.ne.0) then
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
cc!$acc loop reduction(+:sum)
         do k=1,ntclb
          kk=k+nalfa
          sum=sum+csdiag(kk)*vb(k,j)
         enddo
         bsdiag(j)=sum
        enddo
       endif

       if(nalfa.ne.nvecb) then
        m1=nalfa+1
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
cc!$acc loop reduction(+:sum)
         do k=m1,nvecb
          kk=k+ntclb
          sum=sum+csdiag(kk)*vb(k,j)
         enddo
         bsdiag(j)=bsdiag(j)+sum
        enddo
       endif

*        do j=1,nsymf
        do j=1,nbas
        csdiag(j)=w1(j)
       enddo
      endif

!     transformation of the input-matrix a
!                  first:
!     calculation of the intermediate matrices

      if(nalfa.ne.0) then
       do i=1,nalfa
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
cc!$acc loop reduction(+:sum)
         do k=1,nalfa
          sum=sum+ta(i,k)*vb(k,j)
         enddo
         w2(i,j)=sum
        enddo
       enddo
*        do j=1,nsymf
        do j=1,nbas
        do i=1,nalfa
         ta(i,j)=w2(i,j)
        enddo
       enddo
      endif

      m1=nalfa+1
      if(nalfa.ne.nelecs) then
       do i=m1,nelecs
*        do j=1,nsymf
        do j=1,nbas
         sum=0.0d0
         if(ntclb.ne.0) then
cc!$acc loop reduction(+:sum)
          do k=1,ntclb
           kk=k+nalfa
           sum=sum+ta(i,kk)*vb(k,j)
          enddo
         endif
         if(nalfa.ne.nvecb) then
cc!$acc loop reduction(+:sum)
          do k=m1,nvecb
           kk=k+ntclb
           sum=sum+ta(i,kk)*vb(k,j)
          enddo
         endif
         w2(i,j)=sum
        enddo
       enddo
*        do j=1,nsymf
        do j=1,nbas
        do i=m1,nelecs
         ta(i,j)=w2(i,j)
        enddo
       enddo
      endif

!     calculation of the two final matrices

      if(nalfa.ne.0) then
*       do j=1,nsymf
*        do i=1,nsymf
       do j=1,nbas 
        do i=1,nbas 
         sum=0.0
cc!$acc loop reduction(+:sum)
         do k=1,nalfa
          sum=sum+ta(k,j)*va(k,i)
         enddo
         aaa(i,j)=sum
        enddo
       enddo
      endif

      if(nalfa.ne.nelecs) then
*       do j=1,nsymf
*        do i=1,nsymf
       do j=1,nbas 
        do i=1,nbas
         sum=0.0d0
         if(ntcla.ne.0) then
cc!$acc loop reduction(+:sum)
          do k=1,ntcla
           kk=k+nalfa
           sum=sum+ta(kk,j)*va(k,i)
          enddo
         endif
         if(nalfa.ne.nveca) then
cc!$acc loop reduction(+:sum)
          do k=m1,nveca
           kk=k+ntcla
           sum=sum+ta(kk,j)*va(k,i)
          enddo
         endif
         w2(i,j)=sum
        enddo
       enddo
*       do j=1,nsymf
*        do i=1,nsymf
       do j=1,nbas 
        do i=1,nbas
         ta(i,j)=w2(i,j)
        enddo
       enddo
      else
*       do j=1,nsymf
*        do i=1,nsymf
       do j=1,nbas 
        do i=1,nbas
         ta(i,j)=0.0d0
        enddo
       enddo
      endif

cc!$acc end kernels

      return
      end
