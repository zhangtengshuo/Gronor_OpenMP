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
!! GNOME input of sym file
!!     
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  

      subroutine gronor_gnome_molcas_input
      use gnome_parameters
      use gnome_data
      use cidef
      implicit none

      integer               :: i,l,iAtom
*      integer, allocatable  :: nCntr(:,:)
      character (len = 19)  :: lAO

      read(lfnsym,*) text
!     first line of $Project.sym is a dummy
      
      read(lfnsym,'(a72)') name(1)
      read(lfnsym,'(a72)') name(2)
      if(ipr.gt.0) write(lfnout,601) name(1),name(2)
 601  format(/,' OpenMolcas header from Seward (integral module)',/,
     &     a,/,a)
      flush(lfnout) 
      read(lfnsym,'(a17,i4)') text,nNucl

      allocate( centn(nNucl)   )
      allocate( xCord(nNucl)   )
      allocate( yCord(nNucl)   )
      allocate( zCord(nNucl)   )
      allocate(  zNuc(nNucl)   )
      allocate( nCntr(nNucl,7) )

      read(lfnsym,*) text
      lAO = 's  p  d  f  g  h  i'
      
      if(ipr.gt.0) write(lfnout,602) lAO

      write(lfnarx,603) nNucl
 603  format('Atoms',i10)
      write(lfnarx,602) lAO
 602  format(/,6x,'Atom',10x,'x',13x,'y',13x,'z',10x,
     &     'Znuc',5x,a,/)

      zNucTot=0.0d0
      do iAtom = 1, nNucl
        read(lfnsym,*) i,centn(iAtom),xCord(iAtom),yCord(iAtom),
     &       zCord(iAtom),zNuc(iAtom),(nCntr(iAtom,l),l=1,7)
        if(ipr.gt.0)
     &       write(lfnout,604)i,centn(iAtom),xCord(iAtom),yCord(iAtom),
     &       zCord(iAtom),zNuc(iAtom),(nCntr(iAtom,l),l=1,7)
        write(lfnarx,604)i,centn(iAtom),xCord(iAtom),yCord(iAtom),
     &       zCord(iAtom),zNuc(iAtom),(nCntr(iAtom,l),l=1,7)
        zNucTot=zNucTot+zNuc(iAtom)
      enddo
      
      flush(lfnout)
 604  format (i4,2x,a4,2x,3f14.8,2x,f7.2,3x,7i3)
*      deallocate( nCntr )
      
      return
      end
