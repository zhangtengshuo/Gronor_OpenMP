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
!! Transform vectors
!!
!! This routine orders the vectors for determinant idet (value 1 or 2)
!!
!! Upon exit the array vec has the following order:
!!   1: closed shell orbitals
!!   2: alpha spin electron orbitals
!!   3: beta spin electron orbitals
!!
!! Upon exit:
!!   vec(idet)     : orbital coefficients ordered: closed, alpha, beta
!!   ntcl(idet)    : number of closed shell oirbitals
!!   ntop(idet)    : number of open shell orbitals
!!   ioccn(*,idet) : occupation of the open shells: 1=alpha, -1=beta
!!
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!

subroutine gronor_transvc(lfndbg,idet,ntcl,ntop,nclose,nopen,iocopen)
  
  use cidist
  use gnome_data
  use gnome_parameters
  
  implicit none

  integer, intent(in) :: lfndbg,idet
  integer, intent(inout) :: ntcl(:),ntop(:),nclose(:),nopen(:)
  integer, intent(inout) :: iocopen(:,:)
  integer :: m,n1,iop,ib,iv,k,i,l,im,ls,norbs
  if(idbg.ge.17) write(lfndbg,601)
601 format(/,' Determinant Irrep  Nclose  Nopen OccActive',/)

  ntop(idet)=0
  ntcl(idet)=0

  norbs=nclose(idet)+nopen(idet)
  n1=0
  m=0
  do iop=1,nopen(idet)
    n1=n1+1
    if(abs(iocopen(n1,idet)).eq.1) then
      m=m+1
      ioccn(m,idet)=iocopen(n1,idet)
    endif
  enddo
  if(idbg.ge.13) write(lfndbg,602) idet,nclose(idet),nopen(idet),(ioccn(i,idet),i=1,nopen(idet))
602 format(3i8,20i4)

  ntop(idet)=ntop(idet)+nopen(idet)
  ntcl(idet)=ntcl(idet)+nclose(idet)

  !     Initialize the output M.O.'s

  do iv=1,mvec
    do ib=1,nbas
      vtemp(iv,ib,idet)=0.0d0
    enddo
  enddo

  !     Transform the M.O.'s   (no so much to transform if there is no symmetry)

  do im=1,norbs           ! number of non-empty orbitals
    do k=1,nbas           ! number of basis functions
      vtemp(im,k,idet)=vec(im,k,idet)
    enddo
  enddo
  if(idbg.ge.13) write(lfndbg,604)
604 format(/,' M.O.''s are transformed')

  !     Put transformed M.O's in correct order: closed, alpha, beta
  !     First order closed shell M.O.'s of all subspecies

  l=0
  m=0
  do k=1,nclose(idet)
    m=m+1
    l=l+1
    do ib=1,nbas
      vec(m,ib,idet)=vtemp(l,ib,idet)
    enddo
  enddo
  l=l+nopen(idet)

  if(idbg.ge.13) write(lfndbg,606)
606 format(/,' Closed shells ordered')

  !     Then order open alpha shell M.O.'s of all subspecies

  ls=0

  l=ls+nclose(idet)
  do k=1,nopen(idet)
    l=l+1
    if(iocopen(k,idet).eq.1) then
      m=m+1
      do ib=1,nbas
        vec(m,ib,idet)=vtemp(l,ib,idet)
      enddo
    endif
  enddo
  ls=l

  if(idbg.ge.13) write(lfndbg,608)
608 format(/,' Open alpha shells ordered')

  !     Then order open beta shell M.O.'s of all subspecies

  ls=0

  l=ls+nclose(idet)
  do k=1,nopen(idet)
    l=l+1
    if(iocopen(k,idet).eq.-1) then
      m=m+1
      do ib=1,nbas
        vec(m,ib,idet)=vtemp(l,ib,idet)
      enddo
    endif
  enddo
  ls=l

  if(idbg.ge.13) write(lfndbg,610)
610 format(/,' Open beta shells ordered')

  return
end subroutine gronor_transvc
