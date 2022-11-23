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
!! Center for dipole moment intagrals is the origin, for quadrupole integrals is the charge center
!! @author  A. Sanchez-Mansilla, URV
!! @date    2022

      subroutine gronor_multipoles_nuclear()
      use cidef
      use gnome_data
      use gnome_parameters
      implicit none
            
      integer      :: i
      real(kind=8) :: xorig,yorig,zorig

      xorig=0.0d0
      yorig=0.0d0
      zorig=0.0d0
      do i=1,9
        mnuc(i)=0.0d0
      enddo
      print*,"Origin for dipole: x y z"
      print*,xorig,yorig,zorig      
      do i=1,nnucl
        mnuc(1)=mnuc(1)+znuc(i)*(xcord(i)-xorig)
        mnuc(2)=mnuc(2)+znuc(i)*(ycord(i)-yorig)
        mnuc(3)=mnuc(3)+znuc(i)*(zcord(i)-zorig)
      enddo
      xorig=CoM(1)
      yorig=CoM(2)
      zorig=CoM(3)
      print*,"Origin for quadrupole: x y z"
      print*,xorig,yorig,zorig
      do i=1,nnucl
        mnuc(4)=mnuc(4)+                                                &
     &       znuc(i)*(xcord(i)-xorig)*(xcord(i)-xorig)                  
        mnuc(5)=mnuc(5)+                                                &
     &       znuc(i)*(xcord(i)-xorig)*(ycord(i)-yorig)
        mnuc(6)=mnuc(6)+                                                &
     &       znuc(i)*(xcord(i)-xorig)*(zcord(i)-zorig)
        mnuc(7)=mnuc(7)+                                                &
     &       znuc(i)*(ycord(i)-yorig)*(ycord(i)-yorig)
        mnuc(8)=mnuc(8)+                                                &
     &       znuc(i)*(ycord(i)-yorig)*(zcord(i)-zorig)
        mnuc(9)=mnuc(9)+                                                &
     &       znuc(i)*(zcord(i)-zorig)*(zcord(i)-zorig)
      enddo      

      if (idbg.ge.12) then
        write(lfndbg,120)
        write(lfndbg,130)mnuc(1),mnuc(2),mnuc(3)
        write(lfndbg,140)
        write(lfndbg,150)mnuc(4),mnuc(5),mnuc(6),mnuc(7),mnuc(8),mnuc(9)
 120    format(///,1x,'The nuclear contributions to dipole are :',//)
 130    format(15x,'X component              : ',f20.12,//,                  &
     &       15x,'Y component              : ',f20.12,//,                    &
     &       15x,'Z component              : ',f20.12,//)
 140   format(///,1x,'The nuclear contributions to quadrupole are :',//)
 150    format(15x,'XX component              : ',f20.12,//,                 &
     &      15x,'XY component              : ',f20.12,//,                    &
     &      15x,'XZ component              : ',f20.12,//,                    &
     &      15x,'YY component              : ',f20.12,//,                    &
     &      15x,'YZ component              : ',f20.12,//,                    &
     &      15x,'ZZ component              : ',f20.12,//)
      endif
      return
      end
