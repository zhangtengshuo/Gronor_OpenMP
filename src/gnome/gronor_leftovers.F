!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ bases, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!> @brief
!! Several routines used in generating basestates
!!
!! @author  T. P. Straatsma, ORNL
!! @author  R. Havenith, RuG
!! @author  C. de Graaf, URV
!! @date    2016, 2020
!!
!! @detail  Last leftovers from the old code. Might be worth
!!          the effort to rewrite in a cleaner way.

      integer function isetsign4(jocc,nact)
      implicit none

      integer :: jocc(nact),iocca(32),ioccb(32),iclosed(32)
      integer :: nact,ia,ib,ibita,ibitb,is,j,i,jperm
      integer :: icl,nclosed,nperm,na,nb,ih
      integer, external :: ibitset,ibitc

      isetsign4=0
      ia=0
      ib=0
      ibita=0
      ibitb=0
      iclosed=0
      is=0

      do i=1,nact
       iocca(i)=0
       ioccb(i)=0
      enddo

      do i=1,nact
       iclosed(i)=0
      enddo

      do j=1,nact
       if(jocc(j).eq.2) then
        ia=ia+1
        iocca(ia)=j
        ib=ib+1
        ioccb(ib)=j
       elseif(jocc(j).eq.-1) then
        ib=ib+1
        ioccb(ib)=j
       elseif(jocc(j).eq.1) then
        ia=ia+1
        iocca(ia)=j
       endif
      enddo

      ibita=ibitset(iocca,ia)
      ibitb=ibitset(ioccb,ib)
      nclosed=ibitc(iand(ibita,ibitb))
      call gronor_seticlosed(iclosed,iand(ibita,ibitb),nclosed)
      nperm=0
      icl=1
      do i=1,nclosed
       do j=1,ia
        if(iocca(j).eq.iclosed(i)) then
         if(icl.ne.j) then
          ih=iocca(j)
          iocca(j)=iocca(icl)
          iocca(icl)=ih
          nperm=nperm+1
         endif
         icl=icl+1
         goto 101
        endif
       enddo
 101   continue
      enddo

      icl=1

      do i=1,nclosed
       do j=1,ib
        if(ioccb(j).eq.iclosed(i)) then
         if(icl.ne.j) then
          ih=ioccb(j)
          ioccb(j)=ioccb(icl)
          ioccb(icl)=ih
          nperm=nperm+1
         endif
         icl=icl+1
         goto 102
        endif
       enddo
 102   continue
      enddo

      na=ia-nclosed
      nb=ib-nclosed
      call bubble_sort(iocca(nclosed+1),na,jperm)
      nperm=nperm+jperm
      call bubble_sort(ioccb(nclosed+1),nb,jperm)
      nperm=nperm+jperm

      if(mod(nperm,2).eq.1)then
       isetsign4=-1
      else
       isetsign4=1
      endif

      return
      end

      subroutine gronor_seticlosed(iclosed,in,nclosed)
      use cidist
      implicit none
      integer :: nclosed,in,iclosed(nclosed)
      integer :: mask,jj,i,mask1
      mask=1
      jj=0
      do i=0,31
       mask1=ishft(mask,i)
       if(iand(in,mask1).ne.0) then
        jj=jj+1
        iclosed(jj)=i+1
       endif
      enddo
      if(jj.ne.nclosed) then
        call gronor_abort(290,"Problem in seticlosed")
      end if
      return
      end
