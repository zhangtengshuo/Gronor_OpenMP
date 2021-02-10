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
!! Determine number of integrals
!!     
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  
      subroutine gronor_number_integrals(numone,numtwo)

      use inp
      use cidef
      use gnome_data
      use gnome_integrals
      use gnome_parameters

      implicit none

      integer :: mdum, numone, numtwo
      real (kind=8) :: sdum

      open(unit=lfnint,file=filone,form='unformatted',status='old',
     & err=1994)
      goto 1995
 1994 continue
      open(unit=lfnint,file=filint,form='unformatted',status='old',
     & err=994)
 1995 continue
      rewind(lfnint)
      if(ione.eq.0) then
       read(lfnint,err=993) namint,intone,potnuc,nbas,mbuf,mdum
      else
       read(lfnint,err=993) namint,nbas,intone,potnuc
      endif
      if(ione.eq.0) then
       read(lfnint,err=993) sdum
       read(lfnint,err=993) sdum
       read(lfnint,err=993) sdum
      else
       call gronor_read0(lfnint,intone)
       call gronor_read0(lfnint,intone)
       call gronor_read0(lfnint,intone)
       call gronor_read0(lfnint,intone)
      endif
      int2=0
      read(lfnint,end=993) nt

      numone=nt(1)
      numtwo=nt(2)

      rewind(unit=lfnint)
      close(unit=lfnint)

      return
 993  write(lfnout,983) filint(1:inp_strlen(filint))
      call errquit()
 994  write(lfnout,984) filint(1:inp_strlen(filint))
      call errquit()
 983  format('Error reading integral file ',a)
 984  format('Unable to open integral file ',a)
      end
