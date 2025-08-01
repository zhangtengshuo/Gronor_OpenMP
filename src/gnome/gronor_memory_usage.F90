!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ basis, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!> @brief
!! Collect and print timings to the output file
!!
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
!! The code is instrumented with timers to collect CPU and wall-clock times for the
!! computationally demanding parts of the calculation. These timings are presented for each
!! of the ranks used in the calculation.
!! Because of the fault resilient implementation, the master ranks only waits for a limited
!! time for timning data from other ranks. In that case values of zero are reported for those
!! ranks. This is not an error!
!!

subroutine gronor_memory_usage()

  use mpi
  use cidist
  use cidef
  use gnome_parameters
  use gnome_data
  use gnome_integrals

  implicit none

!  external :: MPI_iSend,MPI_Recv

  integer (kind=8) :: ni,nr,membuf(9)
  real (kind=8) ::gb
  integer (kind=4) :: ireq,ierr,istat(MPI_STATUS_SIZE)
  integer (kind=4) :: ncount,mpitag,mpisrc

  if(me.eq.0) then
    membuf(1)=mbasel
    membuf(2)=mvec
    membuf(3)=nveca
    membuf(4)=nvecb
    membuf(5)=nstdim
    membuf(6)=nelecs
    membuf(7)=nbas
    membuf(8)=ntask
    membuf(9)=nbatch
    ncount=9
    mpitag=12
    call MPI_iSend(membuf,ncount,MPI_INTEGER8,mstr,mpitag,MPI_COMM_WORLD,ireq,ierr)
    call MPI_Request_free(ireq,ierr)
  endif


  if(me.eq.mstr) then
    ncount=9
    mpitag=12
    mpisrc=0
    call MPI_Recv(membuf,ncount,MPI_INTEGER8,MPI_ANY_SOURCE,mpitag,MPI_COMM_WORLD,istat,ierr)
    mbasel=membuf(1)
    mvec=membuf(2)
    nveca=membuf(3)
    nvecb=membuf(4)
    nstdim=membuf(5)
    nelecs=membuf(6)
    nbas=membuf(7)
    ntask=membuf(8)
    nbatch=membuf(9)

    write(lfnout,600)
600 format(/,' Memory Usage Summary per Rank (Approximate)',//,' Array',t25,'    Size in GB',/)
601 format(1x,a,t25,f14.6)

    nr=nbas*nbas+2*int1
    gb=real(8*nr)/real(1073741824)
    write(lfnout,601) "One-electron integrals",gb

    ni=3*nbas*(nbas+1)/2
    gb=real(8*ni)/real(1073741824)
    write(lfnout,601) "Two-electron label index length",gb

    nr=int2
    gb=real(8*nr)/real(1073741824*mgr)
    write(lfnout,601) "Two-electron integrals",gb

    ni=nbase*(4+nmol+maxcib)+4*mstates+3*nmol
    nr=mstates*(maxci+maxvec*maxvec+maxnact*maxci)
    nr=nr+nbase*(1+nbasis*nbasis+maxcib+12*nbase)+maxcib*2
    gb=real(8*ni+8*nr)/real(1073741824)
    write(lfnout,601) "Fragment Lists",gb

    ni=memax*2
    gb=real(4*ni)/real(1073741824)
    write(lfnout,601) "ME List",gb

    ni=mbasel*(2+4*mvec+nveca+2*nvecb+7*max(mbasel,nveca)+mbasel)
    ni=ni+nstdim+nelecs*(1+3*nelecs+6*mbasel)
    if(nbatch.lt.0) then
      ni=ni+ntask*(1+4*nbas+6*nbas*nbas)+16
    elseif(nbatch.eq.0) then
      ni=ni+27
    else
      ni=ni+nbatch*(2+4*nbas+10*nbas*nbas)+13
    endif
    gb=real(8*ni)/real(1073741824)
    write(lfnout,601) "Cofac 1 Arrays",gb


  endif

  if(idbg.gt.50) then
    call swatch(date,time)
    write(lfndbg,'(a,1x,a,a)') date(1:8),time(1:8)," Memory usage analysis completed"
    flush(lfndbg)
  endif
      
  return
end subroutine gronor_memory_usage

