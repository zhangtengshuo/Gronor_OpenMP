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
!! Nuclear contribution to multipole moments
!! Center for dipole moment is the origin, for quadrupole it is the charge center
!! @author  A. Sanchez-Mansilla, URV
!! @date    2022

      subroutine gronor_multipoles_nuclear()
      use cidef
      use gnome_data
      
      implicit none
            
      integer      :: i
      real(kind=8) :: xorig,yorig,zorig

      xorig=0.0d0
      yorig=0.0d0
      zorig=0.0d0
      do i=1,9
        mnuc(i)=0.0d0
      enddo      
      do i=1,nnucl
        xorig=xorig+znuc(i)*xcord(i)/zNucTot
        yorig=yorig+znuc(i)*xcord(i)/zNucTot
        zorig=zorig+znuc(i)*xcord(i)/zNucTot
      enddo      
      do i=1,nnucl
        mnuc(1)=mnuc(1)+znuc(i)*(xcord(i)-xorig)
        mnuc(2)=mnuc(2)+znuc(i)*(ycord(i)-yorig)
        mnuc(3)=mnuc(3)+znuc(i)*(zcord(i)-zorig)
        mnuc(4)=mnuc(4)+
     &       znuc(i)*(xcord(i)-xorig)*(xcord(i)-xorig)
        mnuc(5)=mnuc(5)+
     &       znuc(i)*(xcord(i)-xorig)*(ycord(i)-yorig)
        mnuc(6)=mnuc(6)+
     &       znuc(i)*(xcord(i)-xorig)*(zcord(i)-zorig)
        mnuc(7)=mnuc(7)+
     &       znuc(i)*(ycord(i)-yorig)*(ycord(i)-yorig)
        mnuc(8)=mnuc(8)+
     &       znuc(i)*(ycord(i)-yorig)*(zcord(i)-zorig)
        mnuc(9)=mnuc(9)+
     &       znuc(i)*(zcord(i)-zorig)*(zcord(i)-zorig)
      enddo
      
      return
      end
