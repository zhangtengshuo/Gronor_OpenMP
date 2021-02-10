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

      subroutine remove_duplicates_fildet
     &     (ndetin,ndetout,nact,cicoef,occ,thresh_CI)
      
      implicit none
      
      integer      :: ndetin,ndetout,nact
      integer      :: ndetaux
      integer      :: i,j
      real(kind=8) :: cicoef(ndetin)
      real(kind=8) :: ciaux(ndetin)
      real(kind=8) :: thresh_CI
      character(len=255)    :: occ(ndetin)
      character(len=255)    :: occaux(ndetin)      
      
      do i=1,ndetin
         ciaux(i)=0.0d0
         occaux(i)=''
      enddo
      ndetout=1
      ciaux(1)=cicoef(1)
      occaux(1)=occ(1)
      outer : do i=2,ndetin
         do j=1,ndetout
            if(occaux(j).eq.occ(i))then
               ciaux(j)=ciaux(j)+cicoef(i)               
               cycle outer
            endif
         enddo
         ndetout=ndetout+1
         ciaux(ndetout)=cicoef(i)
         occaux(ndetout)=occ(i)
      enddo outer
 
      do i=1,ndetout   
         cicoef(i)=0.0d0
         occ(i)=''
      enddo
      ndetaux=ndetout
      ndetout=0      
      do i=1,ndetaux
         if(dabs(ciaux(i)).ge.thresh_CI)then
            ndetout=ndetout+1
            cicoef(ndetout)=ciaux(i)
            occ(ndetout)=occaux(i)
         endif
      enddo            
      return
      end
