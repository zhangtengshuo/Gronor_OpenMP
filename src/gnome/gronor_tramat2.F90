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
!! Transformation of MOs
!!
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!

subroutine gronor_tramat2(lfndbg)

  !      Transformation of the  m.o.'s
  !      the new  m.o.'s are adapted to the basis set of the two electon
  !      integrals

  use cidist
  use gnome_parameters
  use gnome_data
  implicit none
  integer :: lfndbg

  integer :: i,j,k,kk,m1
  real (kind=8) :: sum
  real (kind=8) :: aaamax,tamax
  real (kind=8) :: diagmax,bdiagmax,bsdiagmax,csdiagmax,sdiagmax

  if(idbg.ge.13) write(lfndbg,600)
600 format(/,' Cofactor matrix transform to symmetry functions')

#ifdef ACC
!$acc kernels present(va,vb,diag,bdiag,cdiag,bsdiag,csdiag,ta,aaa,w1,w2)
#endif
  do j=1,nbas
    diag(j)=0.0d0
    bdiag(j)=0.0d0
    bsdiag(j)=0.0d0
  enddo

  if(ising .ne. 0.) then


    !     transformation of diag and sdiag

    if(nalfa.ne.0) then

#ifdef ACC
!$acc loop private(sum)
#endif
      do j=1,nbas
        sum=0.0d0
#ifdef ACC
!$acc loop reduction(+:sum)
#endif
        do k=1,nalfa
          sum=sum+cdiag(k)*va(k,j)
        enddo
        diag(j)=sum
      enddo

    endif

    if(ntcla.ne.0) then

#ifdef ACC
!$acc loop private(sum)
#endif
      do j=1,nbas
        sum=0.0d0
#ifdef ACC
!$acc loop reduction(+:sum)
#endif
        do k=1,ntcla
          kk=k+nalfa
          sum=sum+cdiag(kk)*va(k,j)
        enddo
        bdiag(j)=sum
      enddo

    endif

    if(nalfa.ne.nveca) then
      m1=nalfa+1

#ifdef ACC
!$acc loop private(sum)
#endif
      do j=1,nbas
        sum=0.0d0
#ifdef ACC
!$acc loop reduction(+:sum)
#endif
        do k=m1,nveca
          kk=k+ntcla
          sum=sum+cdiag(kk)*va(k,j)
        enddo
        bdiag(j)=bdiag(j)+sum
      enddo

    endif

    if(nalfa.ne.0) then

#ifdef ACC
!$acc loop private(sum)
#endif
      do j=1,nbas
        sum=0.0d0
#ifdef ACC
!$acc loop reduction(+:sum)
#endif
        do k=1,nalfa
          sum=sum+csdiag(k)*vb(k,j)
        enddo
        w1(j)=sum
      enddo

    endif

    if(ntclb.ne.0) then

#ifdef ACC
!$acc loop private(sum)
#endif
      do j=1,nbas
        sum=0.0d0
#ifdef ACC
!$acc loop reduction(+:sum) private(kk)
#endif
        do k=1,ntclb
          kk=k+nalfa
          sum=sum+csdiag(kk)*vb(k,j)
        enddo
        bsdiag(j)=sum
      enddo

    endif

    if(nalfa.ne.nvecb) then
      m1=nalfa+1

#ifdef ACC
!$acc loop private(sum)
#endif
      do j=1,nbas
        sum=0.0d0
#ifdef ACC
!$acc loop reduction(+:sum) private(kk)
#endif
        do k=m1,nvecb
          kk=k+ntclb
          sum=sum+csdiag(kk)*vb(k,j)
        enddo
        bsdiag(j)=bsdiag(j)+sum
      enddo

    endif

#ifdef ACC
!$acc loop
#endif
    do j=1,nbas
      csdiag(j)=w1(j)
    enddo

  endif

!     transformation of the input-matrix a
!                  first:
!     calculation of the intermediate matrices

  if(nalfa.ne.0) then

#ifdef ACC
!$acc loop collapse(2) private(sum)
#endif
    do i=1,nalfa
      do j=1,nbas
        sum=0.0d0
#ifdef ACC
!$acc loop reduction(+:sum)
#endif
        do k=1,nalfa
          sum=sum+ta(i,k)*vb(k,j)
        enddo
        w2(i,j)=sum
      enddo
    enddo

#ifdef ACC
!$acc loop collapse(2)
#endif
    do j=1,nbas
      do i=1,nalfa
        ta(i,j)=w2(i,j)
      enddo
    enddo

  endif

  m1=nalfa+1
  if(nalfa.ne.nelecs) then

    do i=m1,nelecs
      do j=1,nbas
        sum=0.0d0
        if(ntclb.ne.0) then
#ifdef ACC
!$acc loop reduction(+:sum) private(kk)
#endif
          do k=1,ntclb
            kk=k+nalfa
            sum=sum+ta(i,kk)*vb(k,j)
          enddo
        endif
        if(nalfa.ne.nvecb) then
#ifdef ACC
!$acc loop reduction(+:sum) private(kk)
#endif
          do k=m1,nvecb
            kk=k+ntclb
            sum=sum+ta(i,kk)*vb(k,j)
          enddo
        endif
        w2(i,j)=sum
      enddo
    enddo
    do j=1,nbas
      do i=m1,nelecs
        ta(i,j)=w2(i,j)
      enddo
    enddo

  endif

  !     calculation of the two final matrices

  if(nalfa.ne.0) then

#ifdef ACC
!$acc loop collapse(2) private(sum)
#endif
    do j=1,nbas
      do i=1,nbas
        sum=0.0
#ifdef ACC
!$acc loop seq reduction(+:sum)
#endif
        do k=1,nalfa
          sum=sum+ta(k,j)*va(k,i)
        enddo
        aaa(i,j)=sum
      enddo
    enddo

  endif

  if(nalfa.ne.nelecs) then

#ifdef ACC
!$acc loop collapse(2) private(sum)
#endif
    do j=1,nbas
      do i=1,nbas
        sum=0.0d0
        if(ntcla.ne.0) then
#ifdef ACC
!$acc loop seq reduction(+:sum) private(kk)
#endif
          do k=1,ntcla
            kk=k+nalfa
            sum=sum+ta(kk,j)*va(k,i)
          enddo
        endif
        if(nalfa.ne.nveca) then
#ifdef ACC
!$acc loop seq reduction(+:sum) private(kk)
#endif
          do k=m1,nveca
            kk=k+ntcla
            sum=sum+ta(kk,j)*va(k,i)
          enddo
        endif
        w2(i,j)=sum
      enddo
    enddo

#ifdef ACC
!$acc loop collapse(2)
#endif
    do j=1,nbas
      do i=1,nbas
        ta(i,j)=w2(i,j)
      enddo
    enddo

  else

#ifdef ACC
!$acc loop collapse(2)
#endif
    do j=1,nbas
      do i=1,nbas
        ta(i,j)=0.0d0
      enddo
    enddo

  endif

#ifdef ACC
!$acc end kernels
#endif
  if(idbg.gt.90) then
#ifdef ACC
!$acc update host(aaa,ta,diag,bdiag,bsdiag,csdiag,sdiag)
#endif

    write(lfndbg,1601) nbas,nelecs
1601 format(//,' diag:',3i6,/)
    do i=1,nbas
      write(lfndbg,1602) diag(i),bdiag(i),bsdiag(i),csdiag(i),sdiag(i)
    enddo
1602 format((3x,6e20.12))
    tamax=0.0d0
    aaamax=0.0d0
    diagmax=0.0d0
    bdiagmax=0.0d0
    bsdiagmax=0.0d0
    csdiagmax=0.0d0
    sdiagmax=0.0d0
    do j=1,nbas
      diagmax=max(diagmax,diag(j))
      bdiagmax=max(bdiagmax,bdiag(j))
      bsdiagmax=max(bsdiagmax,bsdiag(j))
      csdiagmax=max(csdiagmax,csdiag(j))
      sdiagmax=max(sdiagmax,sdiag(j))
      do i=1,nbas
        tamax=max(tamax,ta(i,j))
        aaamax=max(aaamax,aaa(i,j))
      enddo
    enddo
    write(lfndbg,'(a,7f16.5)') 'Max=',tamax,aaamax,diagmax,bdiagmax,bsdiagmax,csdiagmax,sdiagmax
    flush(lfndbg)
  endif

  return
end subroutine gronor_tramat2
