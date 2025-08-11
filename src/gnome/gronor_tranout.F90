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
!! Write out transformed MOs
!!
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!

subroutine gronor_tranout(lfndbg,idet)
  use cidist
  use gnome_parameters
  use gnome_data
#ifdef DEBUG_HDF5
  use debug_hdf5
#endif
  implicit none
  integer :: idet, lfndbg
  integer :: ivc,ntvc,ibas,i

  ntvc=ntcl(idet)+ntop(idet)

  
#ifdef DEBUG_HDF5
  if(idbg.ge.25) then
    call dbg_write_int_scalar('tranout','nclose',ntcl(idet))
    call dbg_write_int_scalar('tranout','nopen',ntop(idet))
    call dbg_log_msg('tranout','M.O.s transformed and ordered')
  endif
#endif

  
#ifdef DEBUG_HDF5
  if(idbg.gt.23) then
    do ivc=1,ntvc
      if(ivc.eq.ntcl(idet)+1) call dbg_log_msg('tranout','open shell M.O.s:')
      if(idbg.gt.90.or.ivc.lt.11.or.ivc.gt.ntvc-10) then
        call dbg_write_array('tranout','vec',vec(ivc,1:nbas,idet),step=ivc)
      endif
    enddo
  endif
#endif

  return
end subroutine gronor_tranout