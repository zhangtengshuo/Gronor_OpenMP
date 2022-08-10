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
      
!> @brief   Application of the S- operator to generate all the
!!          determinants with M_S smaller than M_S,max
!! @author  Aitor Sanchez, URV
!! @date    October 2020

      subroutine sminop(cicoef_in,cicoef_out,occ_in,occ_out,
     &     ndetin,ndetout,spin,nact,nci)      
      use gnome_parameters
      implicit none
      
      integer          :: nact,spin,ndetin,nci
      integer          :: ndetout(spin+1),occ_in(nact,ndetin)      
      integer          :: occ_out(nact,nci,spin+1)  
      integer          :: i,j,k,l,m,idup
      integer,allocatable :: checklist(:,:)
      real (kind=8)    :: cicoef_in(ndetin),dnorm
      real (kind=8)    :: cicoef_out(nci,spin+1)
      logical          :: duplicated
      
      ndetout=0
      ndetout(1)=ndetin
      do i=1,ndetout(1)
         cicoef_out(i,1)=cicoef_in(i)
         do j=1,nact
            occ_out(j,i,1)=occ_in(j,i)         
         enddo
      enddo      
*     generate ms states with the s minus operator
      do i=2,spin+1
         do j=1,ndetout(i-1)            
            do k=1,nact
               if(occ_out(k,j,i-1).eq.1)then
                  allocate(checklist(nact,ndetout(i)))
                  do l=1,ndetout(i)
                     do m=1,nact
                        checklist(m,l)=occ_out(m,l,i)
                     enddo
                  enddo
                  call check_duplicate(checklist,occ_out(:,j,i-1),
     &                 k,nact,ndetout(i),duplicated,idup)
                  deallocate(checklist)
                  if(.not.duplicated)then
                     ndetout(i)=ndetout(i)+1                     
                     do l=1,nact
                        occ_out(l,ndetout(i),i)=occ_out(l,j,i-1)
                     enddo
                     occ_out(k,ndetout(i),i)=-1
                     cicoef_out(ndetout(i),i)=cicoef_out(j,i-1)
                  else                       
                     cicoef_out(idup,i)=cicoef_out(idup,i)+
     &                    cicoef_out(j,i-1)
                     if( abs(cicoef_out(idup,i)).le.tau_CI) then
                        do l=idup,ndetout(i)
                           cicoef_out(l,i)=cicoef_out(l+1,i)
                           do m=1,nact
                              occ_out(m,l,i)=occ_out(m,l+1,i)
                           enddo
                        enddo
                        ndetout(i)=ndetout(i)-1                        
                     endif                     
                  endif
               endif               
            enddo            
         enddo
         dnorm=0.0d0
         do j=1,ndetout(i)
            dnorm=dnorm+cicoef_out(j,i)**2                        
         enddo
         do j=1,ndetout(i)
            cicoef_out(j,i)=cicoef_out(j,i)*(1.0d0/dsqrt(dnorm))
         enddo
      enddo
* Normalize the MS(max) component
      dnorm = 0.0
      do j = 1, ndetout(1)
        dnorm = dnorm + cicoef_out(j,1)**2
      end do
      dnorm = 1.0/dsqrt(dnorm)
      do j = 1, ndetout(1)
        cicoef_out(j,1) = cicoef_out(j,1)*dnorm
      end do

      return
      end
