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
      use cidist
      use gnome_data
      use gnome_integrals
      use gnome_parameters

      implicit none

      external :: gronor_abort

      integer :: numone, numtwo
      real (kind=8) :: sdum

      open(unit=lfnint,file=filone,form='unformatted',status='old',             &
     & err=1994)
      goto 1995
 1994 continue
      open(unit=lfnint,file=filint,form='unformatted',status='old',             &
     & err=994)
 1995 continue
      rewind(lfnint)
      read(lfnint,err=993) namint,intone,potnuc,nbas,mbuf,mclab
      if(mclab.eq.1) labels=1

      read(lfnint,err=993) sdum
      read(lfnint,err=993) sdum
      read(lfnint,err=993) sdum

      int2=0
      read(lfnint,end=993) nt
      if(nt(2).eq.-1) read(lfnint) nt(2)
      read(lfnint,err=993) numfiles

      numone=nt(1)
      numtwo=nt(2)

      int1=(nbas*(nbas+1))/2
      int2=numtwo

      rewind(unit=lfnint)
      close(unit=lfnint)

      return
 993  write(lfnout,983) trim(filint)
      call gronor_abort(230,trim(filint))
 994  write(lfnout,984) trim(filint)
      call gronor_abort(231,trim(filint))
 983  format('Error reading integral file ',a)
 984  format('Unable to open integral file ',a)
      
      end subroutine gronor_number_integrals
