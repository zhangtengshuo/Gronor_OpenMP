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
!! General Non-Orthogonal Matrix Element calculation using GNOME structure
!!
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!


subroutine gronor_gnome(lfndbg,ihc,nhc)

  use mpi
  use cidist
  use gnome_integrals
  use gnome_parameters
  use gnome_data
#ifdef DEBUG_HDF5
  use debug_hdf5
#endif
  use iso_fortran_env, only: int32
  !      use nvtx

  implicit none

  external :: timer_start,timer_stop
  external :: gronor_dipole
  external :: gronor_cororb
  external :: gronor_gntwo
  external :: gronor_gntwo_canonical
!  external :: gronor_gntwo_batch_indexed
  external :: gronor_gnone
  external :: gronor_tramat2
  external :: gronor_cofac1
  external :: gronor_moover
  external :: gronor_abort
  external :: gronor_tranout
  external :: gronor_transvc

  integer :: lfndbg,idet,k,iv,ib,ihc,nhc,ntvc,ivc,ibas

  logical (kind=4) :: flag
  integer (kind=4) :: ierr,status(MPI_STATUS_SIZE)


  e1=0.0d0
  e2=0.0d0
  e2c=0.0d0
  etot=0.0d0
  ttest=0

  !     If duplicate check for terminate signal

  if(odupl.and.iint.ne.0) then
    call MPI_Test(itreq,flag,status,ierr)
    oterm=flag
    if(oterm) return
  endif

  
#ifdef DEBUG_HDF5
  if(idbg.ge.20) then
    call dbg_write_int_scalar('gnome','nbasis',int(nbasis,int32))
    call dbg_write_int_scalar('gnome','nelec',int(nelec(1),int32))
  endif
#endif

  do idet=1,2

    
#ifdef DEBUG_HDF5
    if(idbg.ge.20) then
      call dbg_write_int_scalar('gnome','idet',int(idet,int32),step=idet)
      call dbg_write_int_array('gnome','ioccup',int(ioccup(1:nact(idet),idet),int32),step=idet)
    endif
#endif

    call timer_start(11)
    call gronor_transvc(lfndbg,idet)
    call timer_stop(11)

    call timer_start(12)
    !       call gronor_order(lfndbg,idet)
    call timer_stop(12)

    call timer_start(13)
    call gronor_tranout(lfndbg,idet)
    call timer_stop(13)

#ifdef DEBUG_HDF5
    if(idbg.ge.20) call dbg_log_msg('gnome','Construction of M.O.set completed',step=idet)
#endif

  enddo

  ntcla=ntcl(1)
  ntclb=ntcl(2)
  ntopa=ntop(1)
  ntopb=ntop(2)

  nveca=ntcla+ntopa
  nvecb=ntclb+ntopb
  ntesta=nveca+ntcla
  ntestb=nvecb+ntclb

#ifdef DEBUG_HDF5
  if(idbg.ge.20) then
    call dbg_write_int_scalar('gnome','ntcla',int(ntcla,int32))
    call dbg_write_int_scalar('gnome','ntclb',int(ntclb,int32))
    call dbg_write_int_scalar('gnome','nveca',int(nveca,int32))
    call dbg_write_int_scalar('gnome','nvecb',int(nvecb,int32))
  endif
#endif

  if(ntesta.ne.ntestb) call gronor_abort(305,"Number of electrons is inconsistent")

  nelecs=ntesta
  n1bas=nbas*(nbas+1)/2
  nstdim=max(1,nelecs*nelecs,n1bas)
  mbasel=max(nelecs,nbas)

  if(nveca.ne.ntcl(1)+ntop(1)) call gronor_abort(306,"Incompatible nveca")
  if(nvecb.ne.ntcl(2)+ntop(2)) call gronor_abort(307,"Incompatible nvecb")

  
#ifdef DEBUG_HDF5
  if(idbg.gt.40) then
    do idet=1,2
      ntvc=ntcl(idet)+ntop(idet)
      call dbg_log_msg('gnome','Closed shell M.O.s',step=idet)
      call dbg_write_int_scalar('gnome','ntvc',int(ntvc,int32),step=idet)
      do ivc=1,ntvc
        if(ivc.eq.ntcl(idet)+1) call dbg_log_msg('gnome','Open shell M.O.s:',step=idet)
        if(idbg.gt.90.or.ivc.lt.11.or.ivc.gt.ntvc-10) then
          call dbg_write_array('gnome','vec',vec(ivc,1:nbas,idet),step=ivc)
        endif
      enddo
    enddo
  endif
#endif

  do ib=1,nbas
    do iv=1,nveca
      va(iv,ib)=vec(iv,ib,1)
    enddo
  enddo
  do ib=1,nbas
    do iv=1,nvecb
      vb(iv,ib)=vec(iv,ib,2)
    enddo
  enddo

#ifdef DEBUG_HDF5
  if(idbg.ge.10) call dbg_write_int_scalar('gnome','iamacc',int(iamacc,int32))
#endif

  if(iamacc.gt.0) then

!$acc data copyin(va,vb)
    
    !  Calculations of the overlap matrices

    call timer_start(14)
    call gronor_moover(lfndbg)
    call timer_stop(14)

    !  Calculation of the cofactor matrices and arrays corresponding to the total overlap

    call timer_start(15)
    call gronor_cofac1(lfndbg)
    call timer_stop(15)
    
#ifdef DEBUG_HDF5
    if(idbg.ge.20) then
      select case(ising)
      case(0)
        call dbg_log_msg('gnome','A has no singularities')
      case(1)
        call dbg_log_msg('gnome','A has a single singularity')
      case(2)
        call dbg_log_msg('gnome','A has two singularities: one electron matrix elements are zero')
      case(3)
        call dbg_log_msg('gnome','A has more than two singularities: one and two electron matrix elements are zero')
      end select
    endif
#endif

    if(ising.lt.3) then

      if(corres) then
        call timer_start(16)
        call gronor_cororb()
        call timer_stop(16)
      endif

      call timer_start(17)
      !          call gronor_tramat()
      call timer_stop(17)

      if(idipole.ne.0) then
        call timer_start(18)
        call gronor_dipole(lfndbg)
        call timer_stop(18)
      endif

      
#ifdef DEBUG_HDF5
      if(idbg.ge.30) then
        call dbg_log_msg('gnome','No calculation of two-electron matrix elements')
      endif
#endif

      !     Transformation of the  m.o.'s into the bassisset of the two
      !     electron integrals

      call timer_start(19)
      !         call gronor_trsym(lfndbg)
      call timer_stop(19)

      !     transformation of the first order cofactor matrix
      !     (x-matrix to f-matrix in terms of the basis set of the 2-el.integr

      call timer_start(20)
      call gronor_tramat2(lfndbg)
      call timer_stop(20)

      !     Calculation of the one electron Hamiltonian matrix elements

      call timer_start(21)
      if(icalc.le.1.and.ising.le.1) call gronor_gnone()
      call timer_stop(21)
    endif

    !     Calculation of the two-electron matrix elements

    if((icalc.eq.2.or.icalc.eq.0)) then
      call timer_start(22)
      !           call nvtxStartRange("gntwo")
      ! if(nbatch.gt.1) then
      !   call gronor_gntwo_batch_indexed(lfndbg,ihc,nhc)
      ! else
      ! Thanks! but we do not use batch 
        if(idevel.eq.0.or.mgr.gt.1) then
          if(ising.le.2) call gronor_gntwo()
        else
          if(ising.le.2) call gronor_gntwo_canonical()
        endif
      ! endif
      !         call nvtxEndRange
      call timer_stop(22)
    endif

!$acc end data

  else
! error! not in acc
  endif

#ifdef DEBUG_HDF5
  if(idbg.ge.20) then
    call dbg_write_scalar('gnome','e1',e1)
    call dbg_write_scalar('gnome','e2',e2)
    call dbg_write_scalar('gnome','e2c',e2c)
    call dbg_write_scalar('gnome','etot',etot)
  endif
#endif

  return
end subroutine gronor_gnome