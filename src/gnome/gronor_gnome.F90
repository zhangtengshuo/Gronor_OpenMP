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


module gronor_gnome_mod
  use mpi
  use cidist
  use gnome_integrals
  use gnome_parameters
  use gnome_data
#ifdef _OPENMP
  use omp_lib
#endif
  use gronor_moover_mod,    only: gronor_moover
  use gronor_cofac1_mod,    only: gronor_cofac1
  use gronor_cororb_mod,    only: gronor_cororb
  use gronor_gntwo_mod,     only: gronor_gntwo, gronor_gntwo_canonical
  use gronor_gnone_mod,     only: gronor_gnone
  use gronor_tramat2_mod,   only: gronor_tramat2
  use gronor_dipole_mod,    only: gronor_dipole
  implicit none
contains
subroutine gronor_gnome(lfndbg,ihc,nhc,va,vb,tb,ta,a,u,w,wt,ev,w1,w2,taa,sm,aaa,aat,tt,sdiag,diag,bsdiag,bdiag,csdiag,cdiag)

  implicit none

  real (kind=8), intent(inout) :: va(:,:),vb(:,:),tb(:,:),ta(:,:),a(:,:)
  real (kind=8), intent(inout) :: u(:,:),w(:,:),wt(:,:),ev(:)
  real (kind=8), intent(inout) :: w1(:),w2(:,:),taa(:,:),sm(:,:),aaa(:,:),aat(:,:),tt(:,:)
  real (kind=8), intent(inout) :: sdiag(:),diag(:),bsdiag(:),bdiag(:),csdiag(:),cdiag(:)

  external :: timer_start,timer_stop
  external :: gronor_abort
  external :: gronor_tranout
  external :: gronor_transvc
  external :: swatch

  integer :: lfndbg,ihc,nhc
  integer :: idet=0,k=0,iv=0,ib=0,ntvc=0,ivc=0,ibas=0

  logical (kind=4) :: flag=.false.
  integer (kind=4) :: ierr=0,status(MPI_STATUS_SIZE)=0

  integer :: thread_id
  character(len=10) :: today, now

#ifdef _OPENMP
  thread_id = omp_get_thread_num()
#else
  thread_id = 0
#endif

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

  if(idbg.ge.20) then
    write(lfndbg,600) nbasis
600 format(/,' Number of basis functions is',t45,i8)
    write(lfndbg,603) nelec(1)
603 format(' Number of electrons is',t45,i8)
  endif

  do idet=1,2

    if(idbg.ge.20) then
      write(lfndbg,604) idet
604   format(/,' Transformation of MO set',i8,/)
      write(lfndbg,605) (ioccup(k,idet),k=1,nact(idet))
605   format(' Active orbital occupation : ',32i3)
    endif

    call timer_start(11)
    call gronor_transvc(lfndbg, idet)
    call timer_stop(11)

    call timer_start(12)
    !       call gronor_order(lfndbg,idet)
    call timer_stop(12)

    call timer_start(13)
    call gronor_tranout(lfndbg,idet)
    call timer_stop(13)

    if(idbg.ge.20) write(lfndbg,606) idet
606 format(/,' Construction of M.O.set',i2,' completed')

  enddo

  ntcla=ntcl(1)
  ntclb=ntcl(2)
  ntopa=ntop(1)
  ntopb=ntop(2)

  nveca=ntcla+ntopa
  nvecb=ntclb+ntopb
  ntesta=nveca+ntcla
  ntestb=nvecb+ntclb

  if(ntesta.ne.ntestb) call gronor_abort(305,"Number of electrons is inconsistent")

  nelecs=ntesta
  n1bas=nbas*(nbas+1)/2
  nstdim=max(1,nelecs*nelecs,n1bas)
  mbasel=max(nelecs,nbas)

  if(idbg.gt.10 .and. thread_id==0) then
    call swatch(today,now)
    write(lfndbg,'(a,1x,a,a)') today(1:8),now(1:8), " Array dimensions check in gronor_gnome:"
    write(lfndbg,'(a,2i10)') " va:    ", size(va,1), size(va,2)
    write(lfndbg,'(a,2i10)') " vb:    ", size(vb,1), size(vb,2)
    write(lfndbg,'(a,2i10)') " tb:    ", size(tb,1), size(tb,2)
    write(lfndbg,'(a,2i10)') " ta:    ", size(ta,1), size(ta,2)
    write(lfndbg,'(a,2i10)') " a:     ", size(a,1), size(a,2)
    write(lfndbg,'(a,2i10)') " u:     ", size(u,1), size(u,2)
    write(lfndbg,'(a,2i10)') " w:     ", size(w,1), size(w,2)
    write(lfndbg,'(a,2i10)') " wt:    ", size(wt,1), size(wt,2)
    write(lfndbg,'(a,2i10)') " sm:    ", size(sm,1), size(sm,2)
    write(lfndbg,'(a,2i10)') " aaa:   ", size(aaa,1), size(aaa,2)
    write(lfndbg,'(a,2i10)') " aat:   ", size(aat,1), size(aat,2)
    write(lfndbg,'(a,2i10)') " tt:    ", size(tt,1), size(tt,2)
    write(lfndbg,'(a,i10)')  " ev:    ", size(ev)
    write(lfndbg,'(a,i10)')  " w1:    ", size(w1)
    write(lfndbg,'(a,2i10)') " w2:    ", size(w2,1), size(w2,2)
    write(lfndbg,'(a,i10)')  " sdiag: ", size(sdiag)
    write(lfndbg,'(a,i10)')  " diag:  ", size(diag)
    write(lfndbg,'(a,i10)')  " bsdiag:", size(bsdiag)
    write(lfndbg,'(a,i10)')  " bdiag: ", size(bdiag)
    write(lfndbg,'(a,i10)')  " csdiag:", size(csdiag)
    write(lfndbg,'(a,i10)')  " cdiag: ", size(cdiag)
    write(lfndbg,'(a,2i10)') " taa:   ", size(taa,1), size(taa,2)
    write(lfndbg,'(a,i10)')  " ihc:   ", ihc
    write(lfndbg,'(a,i10)')  " nhc:   ", nhc
    write(lfndbg,'(a,i10)')  " nveca: ", nveca
    write(lfndbg,'(a,i10)')  " nvecb: ", nvecb
    write(lfndbg,'(a,i10)')  " nelecs:", nelecs
    write(lfndbg,'(a,i10)')  " mbasel:", mbasel
    flush(lfndbg)
  endif

  if(nveca.ne.ntcl(1)+ntop(1)) call gronor_abort(306,"Incompatible nveca")
  if(nvecb.ne.ntcl(2)+ntop(2)) call gronor_abort(307,"Incompatible nvecb")

  if(idbg.gt.40) then
    do idet=1,2
      ntvc=ntcl(idet)+ntop(idet)
      write(lfndbg,1603) ntvc
1603  format(/,' Closed shell M.O.''s',i5,/)
      do ivc=1,ntvc
        if(ivc.eq.ntcl(idet)+1) write(lfndbg,1604)
1604    format(/,' Open shell M.O.'' s:')
        if(idbg.gt.90.or.ivc.lt.11.or.ivc.gt.ntvc-10) then
          write(lfndbg,1605)  ' (',ivc,')',(vec(ivc,ibas,idet),ibas=1,nbas)
1605      format(a2,i3,a1,(t9,10f12.8))
        endif
      enddo
    enddo
  endif

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

  if(idbg.ge.10) then
    write(lfndbg,3612) iamacc
3612 format(" GNOME with iamacc ",i5)
    flush(lfndbg)
  endif

  if(iamacc.gt.0) then

!$acc data present(va,vb,tb,ta,a,u,w,wt,ev,w1,w2,sm,aaa,aat,tt,sdiag,diag,bsdiag,bdiag,csdiag,cdiag)
    
    !  Calculations of the overlap matrices

    call timer_start(14)
    call gronor_moover(lfndbg,va,vb,tb,ta,a)
    call timer_stop(14)

    !  Calculation of the cofactor matrices and arrays corresponding to the total overlap

    call timer_start(15)
    call gronor_cofac1(lfndbg,a,u,w,wt,ev,ta,diag,sdiag,cdiag,csdiag)
    call timer_stop(15)

    if(idbg.ge.20) then
      if(ising.eq.0) write(lfndbg,609)
      if(ising.eq.1) write(lfndbg,610)
      if(ising.eq.2) write(lfndbg,611)
      if(ising.eq.3) write(lfndbg,612)
609   format(/,' A has no singularities')
610   format(/,' A has a single singularity')
611   format(/,' A has two singularities: one electron matrix elements are zero')
612   format(/,' A has more than two singularities: one and two electron matrix elements are zero')
    endif

    if(ising.lt.3) then

      if(corres) then
        call timer_start(16)
        call gronor_cororb(u,w,va,vb,ev)
        call timer_stop(16)
      endif

      call timer_start(17)
      !          call gronor_tramat()
      call timer_stop(17)

      if(idipole.ne.0) then
        call timer_start(18)
        call gronor_dipole(lfndbg,ta,diag,sdiag)
        call timer_stop(18)
      endif

      if(idbg.ge.30) then
        if(icalc.eq.1.or.icalc.eq.3) write(lfndbg,613) icalc
613     format(' No calculation of two-electron matrix elements',i4)
      endif

      !     Transformation of the  m.o.'s into the bassisset of the two
      !     electron integrals

      call timer_start(19)
      !         call gronor_trsym(lfndbg)
      call timer_stop(19)

      !     transformation of the first order cofactor matrix
      !     (x-matrix to f-matrix in terms of the basis set of the 2-el.integr

      call timer_start(20)
      call gronor_tramat2(lfndbg,va,vb,ta,aaa,w1,w2,diag,bdiag,bsdiag,cdiag,csdiag,sdiag)
      call timer_stop(20)

      !     Calculation of the one electron Hamiltonian matrix elements

      call timer_start(21)
      if(icalc.le.1.and.ising.le.1) call gronor_gnone(lfndbg,diag,bdiag,bsdiag,csdiag,ta,aaa)
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
          if(ising.le.2) call gronor_gntwo(lfndbg,aat,aaa,tt,ta,sm,diag,bdiag,bsdiag,csdiag)
        else
          if(ising.le.2) call gronor_gntwo_canonical(lfndbg,aat,aaa,tt,ta,sm,diag,bdiag,bsdiag,csdiag)
        endif
      ! endif
      !         call nvtxEndRange
      call timer_stop(22)
    endif

!$acc end data

  else
! error! not in acc
  endif

  return
end subroutine gronor_gnome

end module gronor_gnome_mod
