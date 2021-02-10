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
      
      subroutine check_duplicate(list,det,korb,nact,ndetlist,dup,i_dup,
     &     lfnout)

      implicit none
      
      integer :: nact,ndetlist,i,j,korb
      integer :: list(nact,ndetlist),det(nact),newdet(nact)
      integer :: i_dup,lfnout
      logical :: dup
      do i=1,nact
         newdet(i)=det(i)
      enddo
      newdet(korb)=-1
      dup=.false.
      do i=1,ndetlist         
         if(all(newdet.eq.list(:,i)))then
            dup=.true.
            i_dup=i
            goto 99
         else
            dup=.false.
         endif
      enddo
 99   continue
      return
      end
