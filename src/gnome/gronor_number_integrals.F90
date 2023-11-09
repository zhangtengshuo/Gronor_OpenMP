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

  open(unit=lfnone,file=filone,form='unformatted',status='old',err=994)

  rewind(lfnone)
  read(lfnone,err=993) namint,intone,potnuc,nbas,mbuf,mclab

  read(lfnone,err=993) sdum
  read(lfnone,err=993) sdum
  read(lfnone,err=993) sdum

  int2=0
  read(lfnone,end=993) nt
  if(nt(2).eq.-1) read(lfnone) nt(2)
  read(lfnone,err=993) numfiles

  numone=nt(1)
  numtwo=nt(2)

  int1=(nbas*(nbas+1))/2
  int2=numtwo

  rewind(unit=lfnone)
  close(unit=lfnone)

  return
993 call gronor_abort(230,trim(filone))
994 call gronor_abort(231,trim(filone))
983 format('Error reading integral file ',a)
984 format('Unable to open integral file ',a)

end subroutine gronor_number_integrals
