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
      
!> @brief   Counts the number of permutations to bring a determinant in
!!          the alpha-beta-alpha-beta-... order.
!!
!! @details  When M_S larger than zero, the excess alphas are moved to
!!          the end.
!! @details  Doubly occupied orbitals are first expanded and empty
!!          orbitals eliminated: a2b0 --> aabb
!! @details  Adaptation of the Python script "civec.py"
!! @author  Aitor Sanchez, URV
!! @author  Coen de Graaf, URV
!! @date    October 2020
!! @param   nact number of active orbitals
!! @param   nperm number of permutations
!! @param   occ input occupation string
!! @param   occ2 expanded occupation string

      integer function perm_ab(occ,nact)
*     permutes the determinant in the alpha-beta order
      implicit none

      integer :: nact,nelec
      integer :: nperm,i,j
      integer :: occ(nact)
      integer,allocatable  :: occ2(:)
      integer :: ms,pos1,pos2,spin,oppspin,aux
      logical :: permutation
      
      nelec = 0
      do i = 1, nact
        nelec = nelec + abs(occ(i))
      end do
      allocate( occ2(nelec) )

      j=0
      ms=0
      nperm=0
      occ2=0
      do i=1,nact
         if(occ(i).eq.2)then
            j=j+1
            occ2(j)=1
            ms=ms+1
            j=j+1            
            occ2(j)=-1
            ms=ms-1
         else if(occ(i).eq.1)then
            j=j+1
            occ2(j)=1
            ms=ms+1
         else if(occ(i).eq.-1)then
            j=j+1
            occ2(j)=-1
            ms=ms-1
         endif
      enddo
      if(ms.gt.0)then
         call alphastoend(occ2,nelec,nperm,ms)
      endif      
      do i=1,nact-ms-1
         pos1=mod(i-1,2)
         if(occ2(i).eq.1)then
            spin=0
         else
            spin=1
         endif
         if(pos1.ne.spin)then
            oppspin=abs(spin-1)
            permutation=.false.
            j=i
            do while(.not.permutation)
               j=j+1
               if(occ2(j).ne.occ2(i))then
                  pos2=mod(j-1,2)
                  if(pos2.ne.oppspin)then
                     aux=occ2(i)
                     occ2(i)=occ2(j)
                     occ2(j)=aux
                     permutation=.true.
                     nperm=nperm+j-i
                  endif
               endif
            enddo
         endif
      enddo
      if(mod(nperm,2).eq.1)then
         perm_ab=-1
      else
         perm_ab=1
      endif
      deallocate (occ2)
      return     
      end

      subroutine alphastoend(occ,nact,nperm,ms)
      
      implicit none

      integer,intent(in)  :: nact,ms
      integer             :: index,i,nperm
      integer             :: occ(nact)
      logical             :: permutation
      
      do i=nact-ms+1,nact
         if(occ(i).eq.-1)then
            index=1
            permutation=.false.
            do while(.not.permutation)
               if(occ(index).eq.1)then
                  occ(index)=-1
                  occ(i)=1
                  nperm=nperm+i-index
                  permutation=.true.
               endif
               index=index+1
            enddo
         endif
      enddo
      return      
      end
