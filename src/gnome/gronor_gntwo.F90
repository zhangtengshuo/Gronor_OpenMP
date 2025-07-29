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
!! Processing the two-electron integrals
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
!! These routines come in eight different forms
!!
!! <table>
!! <caption id="multi_row">gntwo routines</caption>
!! <tr><th> Routine name <th> OpenACC <th> OpenMP <th> Labeled <th> Batched
!! <tr><td> gronor_gntwo <td> \htmlonly &#x2714;\endhtmlonly <td>  <td>  <td>
!! <tr><td> gronor_gntwo_batch <td> \htmlonly &#x2714;\endhtmlonly <td>  <td>  <td>\htmlonly &#x2714;\endhtmlonly
!! <tr><td> gronor_gntwo_omp <td>  <td> \htmlonly &#x2714;\endhtmlonly  <td>  <td>
!! <tr><td> gronor_gntwo_omp_batch <td>  <td> \htmlonly &#x2714;\endhtmlonly <td>  <td>\htmlonly &#x2714;\endhtmlonly
!! </table>
!!

subroutine gronor_gntwo(lfndbg)

  use mpi
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals
  use iso_c_binding, only : c_loc, c_ptr

  use openacc
  use cuda_functions

  implicit none

  external :: timer_start,timer_stop

  integer :: lfndbg,i,ii,jj,k,l,n,kl,intg

  real(kind=8) :: e2n,tsn,sum2,ts,fourdet

#ifdef SINGLEP
  real (kind=4) :: e2t,tst
  real (kind=4) :: aai,abi,bai,bbi,aak,abk,bak,bbk
  real (kind=4) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#else
  real (kind=8) :: e2t,tst
  real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
  real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#endif

  real(kind=8), external :: timer_wall

  real (kind=8) :: aaamax,aatmax
  real (kind=8) :: diagmax,bdiagmax,bsdiagmax,csdiagmax,sdiagmax

  logical :: ldiag,lbdiag

  integer (kind=4) :: istat
  type(c_ptr) :: cpfre, cptot

  if(ising.ge.3) return

  call timer_start(30)

  if(idbg.ge.13) write(lfndbg,600)
600 format(' Calculating two electron matrix element',/)

  !      numint=ig(mod(me,mgr)+1)

  numint=mint2

  if(ising.gt.1) e1=0.0d0
  e2=0.0d0
  ts=0.0d0
  e2n=0.0d0
  tsn=0.0d0
  fourdet=4.0d0*deta

!$acc kernels present(aat,aaa,tt,ta,sm)
  do i=1,nbas
    do k=1,nbas
      tt(i,k)=ta(k,i)
      aat(i,k)=aaa(k,i)
    enddo
  enddo
!$acc end kernels

!$acc kernels present(aat,aaa,tt,ta,sm)
  do i=1,nbas
    do k=1,nbas
      sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+ta(k,i)
    enddo
  enddo
!$acc end kernels

  call timer_stop(30)

  call timer_start(38)

  if(ising.le.0) then

    call timer_start(31)

    tst=ts
    kl=nbas*(nbas+1)/2

!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) copyin(kl,intndx,jntndx)
    do ii=intndx,jntndx
      do jj=ii,kl
        intg=ndx(ii)+jj
        i=lab(1,ii)
        k=lab(2,ii)
        l=lab(1,jj)
        n=lab(2,jj)
        tst=tst+g(intg)*(sm(i,k)*sm(l,n) -aaa(i,n)*aaa(l,k)-ta(i,n)*ta(l,k) &
            -aat(i,n)*aat(l,k)-tt(i,n)*tt(l,k)-aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k) &
            -aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
      enddo
    enddo
!$acc end kernels

!!! 可选改动。上面的代码实际上外层循环就已经占满所有 sm。
! !$acc parallel loop gang vector collapse(2) &
! !$acc   reduction(+:tst)                    &
! !$acc   present(aat,aaa,tt,ta,sm,g,lab,ndx) &
! !$acc   copyin(kl,intndx,jntndx)
!     do ii=intndx,jntndx
!       do jj=intndx,kl
!         if (jj < ii) cycle
!         intg=ndx(ii)+jj
!         i=lab(1,ii)
!         k=lab(2,ii)
!         l=lab(1,jj)
!         n=lab(2,jj)
!         tst=tst+g(intg)*(sm(i,k)*sm(l,n) -aaa(i,n)*aaa(l,k)-ta(i,n)*ta(l,k) &
!             -aat(i,n)*aat(l,k)-tt(i,n)*tt(l,k)-aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k) &
!             -aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
!       enddo
!     enddo
! !$acc end parallel

    ts=tst

    call timer_stop(31)

  else

    call timer_start(34)
    
    diagmax=0.0d0
    bdiagmax=0.0d0
    
!$acc kernels present(diag,bdiag,bsdiag,csdiag)
    do k=1,nbas
      diagmax=max(diagmax,abs(diag(k)))
      bdiagmax=max(bdiagmax,abs(bdiag(k)))
    enddo
!$acc end kernels

    ldiag=diagmax.gt.1.0e-16
    lbdiag=bdiagmax.gt.1.0e-16
    
    e2t=e2
    kl=nbas*(nbas+1)/2
    
    if(ldiag.and.lbdiag) then

!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) &
!$acc& present(diag,bdiag,bsdiag,csdiag) copyin(kl,intndx,jntndx)
    do ii=intndx,jntndx
      do jj=ii,kl
        intg=ndx(ii)+jj
        i=lab(1,ii)
        k=lab(2,ii)
        l=lab(1,jj)
        n=lab(2,jj)
        aai=diag(i)
        abi=bdiag(i)
        bai=csdiag(i)
        bbi=bsdiag(i)
        aak=diag(k)
        abk=bdiag(k)
        bak=csdiag(k)
        bbk=bsdiag(k)
        sum2=aai*bak+aak*bai+abi*bbk+abk*bbi
        aaj=diag(l)
        abj=bdiag(l)
        baj=csdiag(l)
        bbj=bsdiag(l)
        aal=diag(n)
        abl=bdiag(n)
        bal=csdiag(n)
        bbl=bsdiag(n)
        e2t=e2t+g(intg)*(sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l) &
            -bai*(aat(n,k)*aaj+aat(l,k)*aal)-bbi*(tt(n,k)*abj+tt(l,k)*abl) &
            -baj*(aaa(n,i)*aak+aaa(n,k)*aai)-bbj*(ta(n,i)*abk+ta(n,k)*abi) &
            -bak*(aat(n,i)*aaj+aat(l,i)*aal)-bbk*(tt(n,i)*abj+tt(l,i)*abl) &
            -bal*(aaa(l,i)*aak+aaa(l,k)*aai)-bbl*(ta(l,i)*abk+ta(l,k)*abi))
      enddo
    enddo
!$acc end kernels

    elseif(.not.ldiag .and. lbdiag) then

!$acc update host(aat,aaa,tt,ta,sm,g,lab,ndx,diag,bdiag,bsdiag,csdiag)
    do ii=intndx,jntndx
      do jj=ii,kl
        intg=ndx(ii)+jj
        i=lab(1,ii)
        k=lab(2,ii)
        l=lab(1,jj)
        n=lab(2,jj)
        abi=bdiag(i)
        bai=csdiag(i)
        bbi=bsdiag(i)
        abk=bdiag(k)
        bak=csdiag(k)
        bbk=bsdiag(k)
        sum2=abi*bbk+abk*bbi
        abj=bdiag(l)
        baj=csdiag(l)
        bbj=bsdiag(l)
        abl=bdiag(n)
        bal=csdiag(n)
        bbl=bsdiag(n)
        e2t=e2t+g(intg)*(sm(k,i)*(abj*bbl+abl*bbj)+sum2*sm(n,l) &
            -bbi*(tt(n,k)*abj+tt(l,k)*abl)-bbj*(ta(n,i)*abk+ta(n,k)*abi) &
            -bbk*(tt(n,i)*abj+tt(l,i)*abl)-bbl*(ta(l,i)*abk+ta(l,k)*abi)) 
      enddo
    enddo

    elseif(ldiag .and. .not. lbdiag) then

!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) &
!$acc& present(diag,bdiag,bsdiag,csdiag) copyin(kl,intndx,jntndx)
    do ii=intndx,jntndx
      do jj=ii,kl
        intg=ndx(ii)+jj
        i=lab(1,ii)
        k=lab(2,ii)
        l=lab(1,jj)
        n=lab(2,jj)
        aai=diag(i)
        bai=csdiag(i)
        bbi=bsdiag(i)
        aak=diag(k)
        bak=csdiag(k)
        bbk=bsdiag(k)
        sum2=aai*bak+aak*bai
        aaj=diag(l)
        baj=csdiag(l)
        bbj=bsdiag(l)
        aal=diag(n)
        bal=csdiag(n)
        bbl=bsdiag(n)
        e2t=e2t+g(intg)*(sm(k,i)*(aaj*bal+aal*baj)+sum2*sm(n,l) &
            -bai*(aat(n,k)*aaj+aat(l,k)*aal)-baj*(aaa(n,i)*aak+aaa(n,k)*aai) &
            -bak*(aat(n,i)*aaj+aat(l,i)*aal)-bal*(aaa(l,i)*aak+aaa(l,k)*aai))
      enddo
    enddo
!$acc end kernels

    endif
    
    e2=e2t

    call timer_stop(34)

  endif

  if(numdev.gt.1) then
    cpfre=c_loc(memfre)
    cptot=c_loc(memtot)
    !        istat=cudaMemGetInfo(cpfre,cptot)
    memavail=min(memfre,memavail)
  endif

  ts=ts*fourdet
  e2=e2+ts

  call timer_stop(38)

  return
end subroutine gronor_gntwo

subroutine gronor_gntwo_canonical(lfndbg)

  use mpi
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals
  use iso_c_binding, only : c_loc, c_ptr

#ifdef _OPENACC
  use openacc
#ifdef CUDA
  use cuda_functions
#endif
#endif

  implicit none

  external :: timer_start,timer_stop

  integer :: lfndbg,i,k,l,n,kl,intg,intl,ls
  !     integer :: ii,jj

  real(kind=8) :: sum1,sum2,ts,fourdet

#ifdef SINGLEP
  real (kind=4) :: e2t,e2l,e2n,tst,tsl,tsn
  real (kind=4) :: aai,abi,bai,bbi,aak,abk,bak,bbk
  real (kind=4) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#else
  real (kind=8) :: e2t,e2l,e2n,tst,tsl,tsn
  real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
  real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#endif

  real(kind=8), external :: timer_wall

  real (kind=8) :: aaamax,aatmax
  real (kind=8) :: diagmax,bdiagmax,bsdiagmax,csdiagmax,sdiagmax

  logical :: ldiag,lbdiag

#ifdef _OPENACC
  integer (kind=4) :: istat
  type(c_ptr) :: cpfre, cptot
#endif

  if(ising.ge.3) return

  call timer_start(30)

  if(idbg.ge.13) write(lfndbg,600)
600 format(' Calculating two electron matrix element',/)

  !      numint=ig(mod(me,mgr)+1)

  numint=mint2

  if(ising.gt.1) e1=0.0d0
  e2=0.0d0
  ts=0.0d0
  e2n=0.0d0
  tsn=0.0d0
  fourdet=4.0d0*deta

!$acc kernels present(aat,aaa,tt,ta,sm)
  do i=1,nbas
    do k=1,nbas
      tt(i,k)=ta(k,i)
      aat(i,k)=aaa(k,i)
    enddo
  enddo
!$acc end kernels

!$acc kernels present(aat,aaa,tt,ta,sm)
  do i=1,nbas
    do k=1,nbas
      sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+ta(k,i)
    enddo
  enddo
!$acc end kernels

  call timer_stop(30)

  call timer_start(38)

  if(ising.le.0) then

    call timer_start(31)

    tst=ts
    kl=nbas*(nbas+1)/2
    intg=0

!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) copyin(intg,nbas)
    do k=1,nbas
      do i=1,k
        tsn=0.0d0
        do n=k,nbas
          ls=1
          if(n.eq.k) ls=i
          tsl=0.0d0
          do l=ls,n
            intg=intg+1
            tsl=tsl+g(intg)*(sm(i,k)*sm(l,n)-aaa(i,n)*aaa(l,k)-ta(i,n)*ta(l,k) &
                -aat(i,n)*aat(l,k)-tt(i,n)*tt(l,k)-aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k) &
                -aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
          enddo
          tsn=tsn+tsl
        enddo
        tst=tst+tsn
      enddo
    enddo
!$acc end kernels

    ts=tst

    call timer_stop(31)

  else

    call timer_start(34)
    
    diagmax=0.0d0
    bdiagmax=0.0d0
    
!$acc kernels present(diag,bdiag,bsdiag,csdiag)
    do k=1,nbas
      diagmax=max(diagmax,abs(diag(k)))
      bdiagmax=max(bdiagmax,abs(bdiag(k)))
    enddo
!$acc end kernels

    ldiag=diagmax.gt.1.0e-16
    lbdiag=bdiagmax.gt.1.0e-16

    e2t=e2
    kl=nbas*(nbas+1)/2
    intg=0

    if(ldiag.and.lbdiag) then

!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) present(diag,bdiag,bsdiag,csdiag) &
!$acc& copyin(intg,kl,intndx,jntndx)
    do k=1,nbas
      aak=diag(k)
      abk=bdiag(k)
      bak=csdiag(k)
      bbk=bsdiag(k)
      do i=1,k
        aai=diag(i)
        abi=bdiag(i)
        bai=csdiag(i)
        bbi=bsdiag(i)
        sum1=sm(k,i)
        sum2=aai*bak+aak*bai+abi*bbk+abk*bbi
        e2n=0.0d0
        do n=k,nbas
          aal=diag(n)
          abl=bdiag(n)
          bal=csdiag(n)
          bbl=bsdiag(n)
          ls=1
          if(n.eq.k) ls=i
          e2l=0.0d0
          intl=intg
          do l=ls,n
            intl=intl+1
            aaj=diag(l)
            abj=bdiag(l)
            baj=csdiag(l)
            bbj=bsdiag(l)
            e2l=e2l+g(intl)*(sum1*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l) &
                -bai*(aat(n,k)*aaj+aat(l,k)*aal)-bbi*(tt(n,k)*abj+tt(l,k)*abl) &
                -baj*(aaa(n,i)*aak+aaa(n,k)*aai)-bbj*(ta(n,i)*abk+ta(n,k)*abi) &
                -bak*(aat(n,i)*aaj+aat(l,i)*aal)-bbk*(tt(n,i)*abj+tt(l,i)*abl) &
                -bal*(aaa(l,i)*aak+aaa(l,k)*aai)-bbl*(ta(l,i)*abk+ta(l,k)*abi))
          enddo
          intg=intg+n+1-ls
          e2n=e2n+e2l
        enddo
        e2t=e2t+e2n
      enddo
    enddo
!$acc end kernels

  elseif(.not.ldiag .and. lbdiag) then
    
!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) present(diag,bdiag,bsdiag,csdiag) &
!$acc& copyin(intg,kl,intndx,jntndx)
    do k=1,nbas
      abk=bdiag(k)
      bak=csdiag(k)
      bbk=bsdiag(k)
      do i=1,k
        abi=bdiag(i)
        bai=csdiag(i)
        bbi=bsdiag(i)
        sum1=sm(k,i)
        sum2=abi*bbk+abk*bbi
        e2n=0.0d0
        do n=k,nbas
          abl=bdiag(n)
          bal=csdiag(n)
          bbl=bsdiag(n)
          ls=1
          if(n.eq.k) ls=i
          e2l=0.0d0
          intl=intg
          do l=ls,n
            intl=intl+1
            abj=bdiag(l)
            baj=csdiag(l)
            bbj=bsdiag(l)
            e2l=e2l+g(intl)*(sum1*(abj*bbl+abl*bbj)+sum2*sm(n,l) &
                -bbi*(tt(n,k)*abj+tt(l,k)*abl)-bbj*(ta(n,i)*abk+ta(n,k)*abi) &
                -bbk*(tt(n,i)*abj+tt(l,i)*abl)-bbl*(ta(l,i)*abk+ta(l,k)*abi))
          enddo
          intg=intg+n+1-ls
          e2n=e2n+e2l
        enddo
        e2t=e2t+e2n
      enddo
    enddo
!$acc end kernels

  elseif(ldiag .and. .not. lbdiag) then
    
!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) present(diag,bdiag,bsdiag,csdiag) &
!$acc& copyin(intg,kl,intndx,jntndx)
    do k=1,nbas
      aak=diag(k)
      bak=csdiag(k)
      bbk=bsdiag(k)
      do i=1,k
        aai=diag(i)
        bai=csdiag(i)
        bbi=bsdiag(i)
        sum1=sm(k,i)
        sum2=aai*bak+aak*bai
        e2n=0.0d0
        do n=k,nbas
          aal=diag(n)
          bal=csdiag(n)
          bbl=bsdiag(n)
          ls=1
          if(n.eq.k) ls=i
          e2l=0.0d0
          intl=intg
          do l=ls,n
            intl=intl+1
            aaj=diag(l)
            baj=csdiag(l)
            bbj=bsdiag(l)
            e2l=e2l+g(intl)*(sum1*(aaj*bal+aal*baj)+sum2*sm(n,l) &
                -bai*(aat(n,k)*aaj+aat(l,k)*aal)-baj*(aaa(n,i)*aak+aaa(n,k)*aai) &
                -bak*(aat(n,i)*aaj+aat(l,i)*aal)-bal*(aaa(l,i)*aak+aaa(l,k)*aai))
          enddo
          intg=intg+n+1-ls
          e2n=e2n+e2l
        enddo
        e2t=e2t+e2n
      enddo
    enddo
!$acc end kernels

    endif

    e2=e2t

    call timer_stop(34)

  endif

#ifdef _OPENACC
  if(numdev.gt.1) then
    cpfre=c_loc(memfre)
    cptot=c_loc(memtot)
    !        istat=cudaMemGetInfo(cpfre,cptot)
    memavail=min(memfre,memavail)
  endif
#endif

  ts=ts*fourdet
  e2=e2+ts

  call timer_stop(38)

  return
end subroutine gronor_gntwo_canonical
