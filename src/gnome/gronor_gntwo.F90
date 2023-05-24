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

#ifdef _OPENACC
  use openacc
#ifdef CUDA
  use cuda_functions
#endif
#endif

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

#ifdef ACC
!$acc kernels present(aat,aaa,tt,ta,sm)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop collapse(2) private(i,k)
#else
!$omp target teams distribute parallel do collapse(2) private(i,k)
#endif
#endif
  do i=1,nbas
    do k=1,nbas
      tt(i,k)=ta(k,i)
      aat(i,k)=aaa(k,i)
    enddo
  enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
#ifdef ACC
!$acc kernels present(aat,aaa,tt,ta,sm)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop collapse(2)
#else
!$omp target teams distribute parallel do collapse(2)
#endif
#endif
  do i=1,nbas
    do k=1,nbas
      sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+ta(k,i)
    enddo
  enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
  call timer_stop(30)

  call timer_start(38)

  if(ising.le.0) then

    call timer_start(31)

    tst=ts
    kl=nbas*(nbas+1)/2

#ifdef ACC
!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) copyin(kl,intndx,jntndx)
#endif
#ifdef OMPTGT
!$omp taskwait
#ifdef OMP5
!$omp target teams loop reduction(+:tst)
#else
!$omp target teams distribute parallel do reduction(+:tst)
#endif
#endif
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
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
    ts=tst

    call timer_stop(31)

  else

    call timer_start(34)

    e2t=e2
    kl=nbas*(nbas+1)/2

#ifdef ACC
!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) &
!$acc& present(diag,bdiag,bsdiag,csdiag) copyin(kl,intndx,jntndx)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop reduction(+:e2t) &
!$omp& private(intg,i,k,l,n,aai,abi,bai,bbi,aak,abk,bak,bbk,sum2) &
!$omp& private(aaj,abj,baj,bbj,aal,abl,bal,bbl)
#else
!$omp target teams distribute parallel do reduction(+:e2t) &
!$omp& private(intg,i,k,l,n,aai,abi,bai,bbi,aak,abk,bak,bbk,sum2) &
!$omp& private(aaj,abj,baj,bbj,aal,abl,bal,bbl)
#endif
#endif
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
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
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

#ifdef ACC
!$acc kernels present(aat,aaa,tt,ta,sm)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop collapse(2) private(i,k)
#else
!$omp target teams distribute parallel do collapse(2) private(i,k)
#endif
#endif
  do i=1,nbas
    do k=1,nbas
      tt(i,k)=ta(k,i)
      aat(i,k)=aaa(k,i)
    enddo
  enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
#ifdef ACC
!$acc kernels present(aat,aaa,tt,ta,sm)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop collapse(2)
#else
!$omp target teams distribute parallel do collapse(2)
#endif
#endif
  do i=1,nbas
    do k=1,nbas
      sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+ta(k,i)
    enddo
  enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
  call timer_stop(30)

  call timer_start(38)

  if(ising.le.0) then

    call timer_start(31)

    tst=ts
    kl=nbas*(nbas+1)/2
    intg=0

#ifdef ACC
!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) copyin(intg,nbas)
#endif
#ifdef OMPTGT
!$omp taskwait
#ifdef OMP5
!$omp target teams loop reduction(+:tst) private(tsn,tsl)
#else
!$omp target teams distribute parallel do reduction(+:tst) private(tsn,tsl)
#endif
#endif
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
    !        do ii=intndx,jntndx
    !          do jj=ii,kl
    !            intg=ndx(ii)+jj
    !            i=lab(1,ii)
    !            k=lab(2,ii)
    !            l=lab(1,jj)
    !            n=lab(2,jj)
    !            tst=tst+g(intg)*(sm(i,k)*sm(l,n)
    !     &           -aaa(i,n)*aaa(l,k)-ta(i,n)*ta(l,k)
    !     &           -aat(i,n)*aat(l,k)-tt(i,n)*tt(l,k)
    !     &           -aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k)
    !     &           -aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
    !          enddo
    !        enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
    ts=tst

    call timer_stop(31)

  else

    call timer_start(34)

    e2t=e2
    kl=nbas*(nbas+1)/2
    intg=0

#ifdef ACC
!$acc kernels present(aat,aaa,tt,ta,sm,g,lab,ndx) present(diag,bdiag,bsdiag,csdiag) &
!$acc& copyin(intg,kl,intndx,jntndx)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop reduction(+:e2t) &
!$omp& private(aak,abk,bak,bbk,aai,abi,bai,bbi,sum1,sum2,e2n) &
!$omp& private(aal,abl,bal,bbl,ls,e2l,intl,aaj,abj,baj,bbj)
#else
!$omp target teams distribute parallel do reduction(+:e2t) &
!$omp& private(aak,abk,bak,bbk,aai,abi,bai,bbi,sum1,sum2,e2n) &
!$omp& private(aal,abl,bal,bbl,ls,e2l,intl,aaj,abj,baj,bbj)
#endif
#endif
    !        do ii=intndx,jntndx
    !          do jj=ii,kl
    !            intg=ndx(ii)+jj
    !            i=lab(1,ii)
    !            k=lab(2,ii)
    !            l=lab(1,jj)
    !     n=lab(2,jj)

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
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
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

subroutine gronor_gntwo_batch_indexed(lfndbg,ihc,nhc)

  use mpi
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  external :: timer_start,timer_stop

  integer :: lfndbg,i,ii,jj,k,l,n,kl,intg,ihc,nhc,ibl

  real (kind=8) :: ts,fourdet

  real(kind=8), external :: timer_wall

  !     Copying arrays into batch arrays

  !     Initialize nb0 and nb1 when first matrix element in a batch

  !     nb0 : number of singularity 0 elements
  !     nb1 : number of singularity 1 or 2 elements

  if(ihc.eq.1) then
    nb0=0
    nb1=0
    etotb=0.0d0
  endif

  if(ising.lt.3) then

    if(iamhead.eq.1) sstot=sstot+fctr*deta

    if(idbg.ge.50) write(lfndbg,600)
600 format(' Calculating two electron matrix element',/)

    numint=mint2

    if(ising.gt.1) e1=0.0d0
    e2=0.0d0
    ts=0.0d0
    fourdet=4.0d0*deta

    call timer_start(30)

    if(ising.eq.0) then
      nb0=nb0+1
      prefac0(nb0)=4.0d0*deta*fctr

      !     ASSEMBLE ARRAYS ON GPU

#ifdef ACC
!$acc kernels present(aaa0,aat0,tt0,ta0,tt,aat,aaa,ta,sm0) copyin(nb0)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop collapse(2) private(k,i)
#else
!$omp target teams distribute parallel do collapse(2) private(k,i)
#endif
#endif
      do k=1,nbas
        do i=1,nbas
          sm0(i,k,nb0)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa0(i,k,nb0)=aaa(i,k)
          aat0(i,k,nb0)=aaa(k,i)
          ta0(i,k,nb0)=ta(i,k)
          tt0(i,k,nb0)=ta(k,i)
        enddo
      enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
      !     END ASSEMBLY ON GPU

      !!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt

      !!_ACCTGT_($acc update self(ta,aaa,aat,tt))
      !!_OMPTGT_($omp target update from(ta,aaa,aat,tt))
      !      do k=1,nbas
      !        do i=1,nbas
      !          sm0(i,k,nb0)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
      !          aaa0(i,k,nb0)=aaa(i,k)
      !          aat0(i,k,nb0)=aaa(k,i)
      !          tt0(i,k,nb0)=ta(k,i)
      !        enddo
      !      enddo

      !!     END ASSEMBLY ON CPU

    else
      nb1=nb1+1
      prefac1(nb1)=fctr

      !     ASSEMBLE ARRAYS ON GPU

#ifdef ACC
!$acc kernels present(aaa1,aat1,tt1,ta1,tt,aat,aaa,ta,sm1) copyin(nb1)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop collapse(2) private(k,i)
#else
!$omp target teams distribute parallel do collapse(2) private(k,i)
#endif
#endif
      do k=1,nbas
        do i=1,nbas
          sm1(i,k,nb1)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa1(i,k,nb1)=aaa(i,k)
          aat1(i,k,nb1)=aaa(k,i)
          ta1(i,k,nb1)=ta(i,k)
          tt1(i,k,nb1)=ta(k,i)
        enddo
      enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
      !     END ASSEMBLY ON GPU

      !!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt

      !!_ACCTGT_($acc update self(ta,aaa,aat,tt))
      !!_OMPTGT_($omp target update from(ta,aaa,aat,tt))
      !      do k=1,nbas
      !        do i=1,nbas
      !          sm1(i,k,nb1)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
      !          aaa1(i,k,nb1)=aaa(i,k)
      !          aat1(i,k,nb1)=aaa(k,i)
      !          tt1(i,k,nb1)=ta(k,i)
      !        enddo
      !      enddo

      !!     END ASSEMBLY ON CPU

#ifdef ACC
!$acc kernels present(diag,bdiag,bsdiag,csdiag) present(diag1,bdiag1,bsdiag1,csdiag1) copyin(nb1)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop
#else
!$omp target teams distribute parallel do
#endif
#endif
      do i=1,nbas
        diag1(nb1,i)=diag(i)
        bdiag1(nb1,i)=bdiag(i)
        bsdiag1(nb1,i)=bsdiag(i)
        csdiag1(nb1,i)=csdiag(i)
      enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
    endif

    call timer_stop(30)

  endif

  call timer_start(38)

  !     After the last element of a task the list of energies is generated

  if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.nbatch) then

    call timer_start(31)

    !!_ACCTGT_($acc update device(ta0,aaa0,aat0,tt0,sm0))
    !!_OMPTGT_($omp target update to(ta0,aaa0,aat0,tt0,sm0))

#ifdef ACC
!$acc data copyin(prefac0)
#endif
#ifdef OMPTGT
!$omp target data map(to:prefac0)
#endif
    kl=nbas*(nbas+1)/2

#ifdef ACC
!$acc kernels present(lab,ndx,g,sm0,aat0,aaa0,tt0,ta0,prefac0)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop reduction(+:etotb) private(i,k,l,n,intg)
#else
!$omp target teams distribute parallel do reduction(+:etotb) private(i,k,l,n,intg)
#endif
#endif
    do ibl=1,nb0
      do ii=intndx,jntndx
        do jj=ii,kl
          intg=ndx(ii)+jj
          i=lab(1,ii)
          k=lab(2,ii)
          l=lab(1,jj)
          n=lab(2,jj)
          etotb=etotb+g(intg)*prefac0(ibl)*(sm0(i,k,ibl)*sm0(l,n,ibl) &
              -aaa0(i,n,ibl)*aaa0(l,k,ibl)-ta0(i,n,ibl)*ta0(l,k,ibl) &
              -aat0(i,n,ibl)*aat0(l,k,ibl)-tt0(i,n,ibl)*tt0(l,k,ibl) &
              -aat0(l,i,ibl)*aaa0(n,k,ibl)-tt0(l,i,ibl)*ta0(n,k,ibl) &
              -aaa0(l,i,ibl)*aat0(n,k,ibl)-ta0(l,i,ibl)*tt0(n,k,ibl))
        enddo
      enddo
    enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
    !     Reset index in buffer to zero
    nb0=0

#ifdef ACC
!$acc end data
#endif
#ifdef OMPTGT
!$omp end target data
#endif
    call timer_stop(31)

  endif

  if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.nbatch) then

    call timer_start(34)

    !!_ACCTGT_($acc update device(ta1,aaa1,aat1,tt1,sm1))
    !!_OMPTGT_($omp target update to(ta1,aaa1,aat1,tt1,sm1))

#ifdef ACC
!$acc data copyin(prefac1)
#endif
#ifdef OMPTGT
!$omp target data map(to:prefac1)
#endif
    kl=nbas*(nbas+1)/2

#ifdef ACC
!$acc kernels present(lab,ndx,g,sm1,aat1,aaa1,tt1,ta1,prefac1) present(diag1,bdiag1,bsdiag1,csdiag1)
#endif
#ifdef OMPTGT
#ifdef OMP5
!$omp target teams loop reduction(+:etotb) private(i,k,l,n,intg)
#else
!$omp target teams distribute parallel do reduction(+:etotb) private(i,k,l,n,intg)
#endif
#endif
    do ii=intndx,jntndx
      do jj=ii,kl
        intg=ndx(ii)+jj
        i=lab(1,ii)
        k=lab(2,ii)
        l=lab(1,jj)
        n=lab(2,jj)
#ifdef ACC
!$acc loop vector(32)
#endif
        do ibl=1,nb1
          etotb=etotb+g(intg)*prefac1(ibl)*(sm1(k,i,ibl)*(diag1(l,ibl)*csdiag1(n,ibl) &
              +diag1(n,ibl)*csdiag1(l,ibl)+bdiag1(l,ibl)*bsdiag1(n,ibl) &
              +bdiag1(n,ibl)*bsdiag1(l,ibl))+(diag1(i,ibl)*csdiag1(k,ibl) &
              +diag1(k,ibl)*csdiag1(i,ibl)+bdiag1(i,ibl)*bsdiag1(k,ibl) &
              +bdiag1(k,ibl)*bsdiag1(i,ibl))*sm1(n,l,ibl) &
              -csdiag1(i,ibl)*(aat1(n,k,ibl)*diag1(l,ibl) &
              +aat1(l,k,ibl)*diag1(n,ibl))-bsdiag1(i,ibl)*(tt1(n,k,ibl)*bdiag1(l,ibl) &
              +tt1(l,k,ibl)*bdiag1(n,ibl))-csdiag1(l,ibl)*(aaa1(n,i,ibl)*diag1(k,ibl) &
              +aaa1(n,k,ibl)*diag1(i,ibl))-bsdiag1(l,ibl)*(ta1(n,i,ibl)*bdiag1(k,ibl) &
              +ta1(n,k,ibl)*bdiag1(i,ibl))-csdiag1(k,ibl)*(aat1(n,i,ibl)*diag1(l,ibl) &
              +aat1(l,i,ibl)*diag1(n,ibl))-bsdiag1(k,ibl)*(tt1(n,i,ibl)*bdiag1(l,ibl) &
              +tt1(l,i,ibl)*bdiag1(n,ibl))-csdiag1(n,ibl)*(aaa1(l,i,ibl)*diag1(k,ibl) &
              +aaa1(l,k,ibl)*diag1(i,ibl))-bsdiag1(n,ibl)*(ta1(l,i,ibl)*bdiag1(k,ibl) &
              +ta1(l,k,ibl)*bdiag1(i,ibl)))
        enddo
      enddo
    enddo
#ifdef OMPTGT
#ifdef OMP5
!$omp end target teams loop
#else
!$omp end target teams distribute parallel do
#endif
#endif
#ifdef ACC
!$acc end kernels
#endif
    !     Reset index into buffer
    nb1=0

    call timer_stop(34)

#ifdef ACC
!$acc end data
#endif
#ifdef OMPTGT
!$omp end target data
#endif
  endif

  call timer_stop(38)

  return
end subroutine gronor_gntwo_batch_indexed

subroutine gronor_gntwo_omp(lfndbg)

  use mpi
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  external :: timer_start,timer_stop

  integer :: lfndbg

  integer :: i,ii,jj,k,l
  integer :: n,kl
  integer :: intg

  real (kind=8) :: sum2,ts,fourdet
  real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
  real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

  real(kind=8), external :: timer_wall

  if(ising.ge.3) return

  if(idbg.ge.50) write(lfndbg,600)
600 format(' Calculating two electron matrix element',/)

  numint=mint2

  if(ising.gt.1) e1=0.0d0
  e2=0.0d0
  ts=0.0d0
  fourdet=4.0d0*deta

  call timer_start(30)

  do i=1,nbas
    do k=1,nbas
      tt(i,k)=ta(k,i)
      aat(i,k)=aaa(k,i)
    enddo
  enddo

  do i=1,nbas
    do k=1,nbas
      sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+ta(k,i)
    enddo
  enddo

  call timer_stop(30)

  call timer_start(38)

  if(ising.le.0) then

    call timer_start(32)

    kl=nbas*(nbas+1)/2
    intg=0
#ifdef OMP
!$omp parallel do reduction(+:ts) shared(sm,aat,aaa,tt,ta,g,lab,numint) private(intg,i,k,l,n,ii,jj)
#endif
    do ii=intndx,jntndx
      do jj=ii,kl
        intg=ndx(ii)+jj
        i=lab(1,ii)
        k=lab(2,ii)
        l=lab(1,jj)
        n=lab(2,jj)
        ts=ts+g(intg)*(sm(k,i)*sm(n,l) &
            -aat(n,i)*aaa(l,k)-tt(n,i)*ta(l,k)-aaa(n,i)*aat(l,k)-ta(n,i)*tt(l,k) &
            -aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k)-aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
      enddo
    enddo
#ifdef OMP
!$omp end parallel do
#endif
    call timer_stop(32)

  else

    call timer_start(35)

    kl=nbas*(nbas+1)/2
    intg=0
#ifdef OMP
!$omp parallel do reduction(+:e2) shared(diag,bdiag,csdiag,bsdiag,sm,aat,tt,aaa,ta,g,lab,numint) &
!$omp& private(intg,i,k,l,n,aai,abi,bai,bbi,aak,abk,bak,bbk) &
!$omp& private(aaj,abj,baj,bbj,sum2,aal,abl,bal,bbl,ii,jj)
#endif
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
        e2=e2+g(intg)*(sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l) &
            -bai*(aat(n,k)*aaj+aat(l,k)*aal)-bbi*(tt(n,k)*abj+tt(l,k)*abl) &
            -baj*(aaa(n,i)*aak+aaa(n,k)*aai)-bbj*(ta(n,i)*abk+ta(n,k)*abi) &
            -bak*(aat(n,i)*aaj+aat(l,i)*aal)-bbk*(tt(n,i)*abj+tt(l,i)*abl) &
            -bal*(aaa(l,i)*aak+aaa(l,k)*aai)-bbl*(ta(l,i)*abk+ta(l,k)*abi))
      enddo
    enddo
#ifdef OMP
!$omp end parallel do
#endif
    call timer_stop(35)
  endif

  ts=ts*fourdet
  e2=e2+ts

  call timer_stop(38)

  return
end subroutine gronor_gntwo_omp

subroutine gronor_gntwo_omp_batch_indexed(lfndbg,ihc,nhc)

  use mpi
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  external :: timer_start,timer_stop

  integer :: lfndbg

  integer :: i,ii,k,l,n,kl,intg,ihc,nhc,ibl,jj,igg

  real (kind=8) :: etemp

  real (kind=8) :: ts,fourdet

  real(kind=8), external :: timer_wall

  !     Copying arrays into batch arrays

  !     Initialize nb0 and nb1 when first matrix element in a batch

  !     nb0 : number of singularity 0 elements
  !     nb1 : number of singularity 1 or 2 elements

  if(ihc.eq.1) then
    nb0=0
    nb1=0
    etotb=0.0d0
  endif

  if(ising.lt.3) then

    call timer_start(30)

    if(iamhead.eq.1) sstot=sstot+fctr*deta

    if(idbg.ge.50) write(lfndbg,600)
600 format(' Calculating two electron matrix element',/)

    numint=mint2

    if(ising.gt.1) e1=0.0d0
    e2=0.0d0
    ts=0.0d0
    fourdet=4.0d0*deta

    if(ising.eq.0) then
      nb0=nb0+1
      prefac0(nb0)=4.0d0*deta*fctr
#ifdef OMP
!$omp parallel do shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas) private(i,k) collapse(2)
#endif
      do k=1,nbas
        do i=1,nbas
          sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa0(nb0,i,k)=aaa(i,k)
          !              aat0(nb0,i,k)=aaa(k,i)
          tt0(nb0,i,k)=ta(k,i)
        enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
    else
      nb1=nb1+1
      prefac1(nb1)=fctr
#ifdef OMP
!$omp parallel do shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas) private(i,k) collapse(2)
#endif
      do k=1,nbas
        do i=1,nbas
          sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa1(nb1,i,k)=aaa(i,k)
          !              aat1(nb1,i,k)=aaa(k,i)
          tt1(nb1,i,k)=ta(k,i)
        enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
#ifdef OMP
!$omp parallel do shared(diagl,bdiagl,csdiagl,bsdiagl,nbas) private(i)
#endif
      do i=1,nbas
        diag1(nb1,i)=diag(i)
        bdiag1(nb1,i)=bdiag(i)
        bsdiag1(nb1,i)=bsdiag(i)
        csdiag1(nb1,i)=csdiag(i)
      enddo
#ifdef OMP
!$omp end parallel do
#endif
    endif

    call timer_stop(30)

  endif

  call timer_start(38)

  !     After the last element of a task the list or when the buffer is full
  !     the energies are evaluated

  if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.nbatch) then

    call timer_start(32)

    kl=nbas*(nbas+1)/2
    do ii=intndx,jntndx
      igg=ndx(ii)
      i=lab(1,ii)
      k=lab(2,ii)
#ifdef OMP
!$omp parallel shared(igg,i,k,ii,sm0,aaa0,tt0,g,lab,nb0,prefac)
!$omp do reduction(+:etotb,etemp) private(intg,jj,l,n,ibl)
#endif
      do jj=ii,kl
        intg=igg+jj
        l=lab(1,jj)
        n=lab(2,jj)
        etemp=0.0d0
        do ibl=1,nb0
          etemp=etemp+prefac0(ibl)*(sm0(ibl,k,i)*sm0(ibl,n,l) &
              -aaa0(ibl,i,n)*aaa0(ibl,l,k)-aaa0(ibl,n,i)*aaa0(ibl,k,l) &
              -aaa0(ibl,i,l)*aaa0(ibl,n,k)-aaa0(ibl,l,i)*aaa0(ibl,k,n) &
              -tt0(ibl,i,n)*tt0(ibl,l,k)-tt0(ibl,n,i)*tt0(ibl,k,l) &
              -tt0(ibl,i,l)*tt0(ibl,n,k)-tt0(ibl,l,i)*tt0(ibl,k,n))
        enddo
        etotb=etotb+g(intg)*etemp
      enddo
#ifdef OMP
!$omp end do
#endif
#ifdef OMP
!$omp end parallel
#endif
    enddo

    !     Reset index in buffer to zero
    nb0=0

    call timer_stop(32)

  endif

  if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.nbatch) then

    call timer_start(35)

    kl=nbas*(nbas+1)/2
    do ii=intndx,jntndx
      igg=ndx(ii)
      i=lab(1,ii)
      k=lab(2,ii)
#ifdef OMP
!$omp parallel do reduction(+:etotb,etemp) shared(igg,i,k,diagl,bdiagl,csdiagl,bsdiagl) &
!$omp& shared(sm1,aaa1,tt1,g,lab,nb1,prefac) private(intg,l,n,ibl)
#endif
      do jj=ii,kl
        intg=igg+jj
        l=lab(1,jj)
        n=lab(2,jj)
        etemp=0.0d0
        do ibl=1,nb1
          etemp=etemp+prefac1(ibl)* &
              (sm1(ibl,k,i)*(diag1(ibl,l)*csdiag1(ibl,n)+diag1(ibl,n)*csdiag1(ibl,l) &
              +bdiag1(ibl,l)*bsdiag1(ibl,n)+bdiag1(ibl,n)*bsdiag1(ibl,l)) &
              +(diag1(ibl,i)*csdiag1(ibl,k)+diag1(ibl,k)*csdiag1(ibl,i) &
              +bdiag1(ibl,i)*bsdiag1(ibl,k)+bdiag1(ibl,k)*bsdiag1(ibl,i))*sm1(ibl,n,l) &
              -csdiag1(ibl,i)*(aaa1(ibl,k,n)*diag1(ibl,l)+aaa1(ibl,k,l)*diag1(ibl,n)) &
              -bsdiag1(ibl,i)*(tt1(ibl,n,k)*bdiag1(ibl,l)+tt1(ibl,l,k)*bdiag1(ibl,n)) &
              -csdiag1(ibl,l)*(aaa1(ibl,n,i)*diag1(ibl,k)+aaa1(ibl,n,k)*diag1(ibl,i)) &
              -bsdiag1(ibl,l)*(tt1(ibl,i,n)*bdiag1(ibl,k)+tt1(ibl,k,n)*bdiag1(ibl,i)) &
              -csdiag1(ibl,k)*(aaa1(ibl,i,n)*diag1(ibl,l)+aaa1(ibl,i,l)*diag1(ibl,n)) &
              -bsdiag1(ibl,k)*(tt1(ibl,n,i)*bdiag1(ibl,l)+tt1(ibl,l,i)*bdiag1(ibl,n)) &
              -csdiag1(ibl,n)*(aaa1(ibl,l,i)*diag1(ibl,k)+aaa1(ibl,l,k)*diag1(ibl,i)) &
              -bsdiag1(ibl,n)*(tt1(ibl,i,l)*bdiag1(ibl,k)+tt1(ibl,k,l)*bdiag1(ibl,i)))
        enddo
        etotb=etotb+g(intg)*etemp
      enddo
#ifdef OMP
!$omp end parallel do
#endif
    enddo

    !     Reset index into buffer
    nb1=0

    call timer_stop(35)

  endif

  return
end subroutine gronor_gntwo_omp_batch_indexed

subroutine gronor_gntwo_omp_batch_canonical(lfndbg,ihc,nhc)

  use mpi
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals

#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  external :: timer_start,timer_stop

  integer :: lfndbg

  integer :: i,ii,k,l,n,kl,intg,ihc,nhc,ibl,jj,igg,ls,noff

  real (kind=8) :: etemp

  real (kind=8) :: ts,fourdet

  real(kind=8), external :: timer_wall

  !     Copying arrays into batch arrays

  !     Initialize nb0 and nb1 when first matrix element in a batch

  !     nb0 : number of singularity 0 elements
  !     nb1 : number of singularity 1 or 2 elements

  if(ihc.eq.1) then
    nb0=0
    nb1=0
    etotb=0.0d0
  endif

  if(ising.lt.3) then

    call timer_start(30)

    if(iamhead.eq.1) sstot=sstot+fctr*deta

    if(idbg.ge.50) write(lfndbg,600)
600 format(' Calculating two electron matrix element',/)

    numint=mint2

    if(ising.gt.1) e1=0.0d0
    e2=0.0d0
    ts=0.0d0
    fourdet=4.0d0*deta

    if(ising.eq.0) then
      nb0=nb0+1
      prefac0(nb0)=4.0d0*deta*fctr
#ifdef OMP
!$omp parallel do shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas) private(i,k) collapse(2)
#endif
      do k=1,nbas
        do i=1,nbas
          sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa0(nb0,i,k)=aaa(i,k)
          aat0(nb0,i,k)=aaa(k,i)
          ta0(nb0,i,k)=ta(i,k)
          tt0(nb0,i,k)=ta(k,i)
        enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
    else
      nb1=nb1+1
      prefac1(nb1)=fctr
#ifdef OMP
!$omp parallel do shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas) private(i,k) collapse(2)
#endif
      do k=1,nbas
        do i=1,nbas
          sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa1(nb1,i,k)=aaa(i,k)
          !              aat1(nb1,i,k)=aaa(k,i)
          tt1(nb1,i,k)=ta(k,i)
        enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
#ifdef OMP
!$omp parallel do shared(diagl,bdiagl,csdiagl,bsdiagl,nbas) private(i)
#endif
      do i=1,nbas
        diag1(nb1,i)=diag(i)
        bdiag1(nb1,i)=bdiag(i)
        bsdiag1(nb1,i)=bsdiag(i)
        csdiag1(nb1,i)=csdiag(i)
      enddo
#ifdef OMP
!$omp end parallel do
#endif
    endif

    call timer_stop(30)

  endif

  call timer_start(38)

  !     After the last element of a task the list or when the buffer is full
  !     the energies are evaluated

  if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.nbatch) then

    call timer_start(32)

    intg=0
    
#ifdef OMP
!$omp parallel shared(sm0,aaa0,aat0,tt0,ta0,g,nb0,prefac,ndxk)
!$omp do reduction(+:etotb) private(intg,ibl,ls,noff) schedule(dynamic)
#endif
    do k=1,nbas
      intg=ndxk(k)
      do i=1,k
        do n=k,nbas
          ls=1
          if(n.eq.k) ls=i
          noff=intg+1-ls
          do l=ls,n
            do ibl=1,nb0
              etotb=etotb+g(noff+l)*prefac0(ibl)*(sm0(ibl,k,i)*sm0(ibl,n,l) &
                  -aat0(ibl,n,i)*aaa0(ibl,l,k)-aaa0(ibl,n,i)*aat0(ibl,l,k) &
                  -aat0(ibl,l,i)*aaa0(ibl,n,k)-aaa0(ibl,l,i)*aat0(ibl,n,k) &
                  -ta0(ibl,n,i)*tt0(ibl,l,k)-tt0(ibl,n,i)*ta0(ibl,l,k) &
                  -ta0(ibl,l,i)*tt0(ibl,n,k)-tt0(ibl,l,i)*ta0(ibl,n,k))
            enddo
          enddo
          intg=noff+n
        enddo
      enddo
    enddo
#ifdef OMP
!$omp end do
!$omp end parallel
#endif
    
    kl=nbas*(nbas+1)/2

    !     Reset index in buffer to zero

    nb0=0

    call timer_stop(32)

  endif

  if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.nbatch) then

    call timer_start(35)

    kl=nbas*(nbas+1)/2
    do ii=intndx,jntndx
      igg=ndx(ii)
      i=lab(1,ii)
      k=lab(2,ii)
#ifdef OMP
!$omp parallel do reduction(+:etotb,etemp) shared(igg,i,k,diagl,bdiagl,csdiagl,bsdiagl) &
!$omp& shared(sm1,aaa1,tt1,g,lab,nb1,prefac) private(intg,l,n,ibl)
#endif
      do jj=ii,kl
        intg=igg+jj
        l=lab(1,jj)
        n=lab(2,jj)
        etemp=0.0d0
        do ibl=1,nb1
          etemp=etemp+prefac1(ibl)* &
              (sm1(ibl,k,i)*(diag1(ibl,l)*csdiag1(ibl,n)+diag1(ibl,n)*csdiag1(ibl,l) &
              +bdiag1(ibl,l)*bsdiag1(ibl,n)+bdiag1(ibl,n)*bsdiag1(ibl,l)) &
              +(diag1(ibl,i)*csdiag1(ibl,k)+diag1(ibl,k)*csdiag1(ibl,i) &
              +bdiag1(ibl,i)*bsdiag1(ibl,k)+bdiag1(ibl,k)*bsdiag1(ibl,i))*sm1(ibl,n,l) &
              -csdiag1(ibl,i)*(aaa1(ibl,k,n)*diag1(ibl,l)+aaa1(ibl,k,l)*diag1(ibl,n)) &
              -bsdiag1(ibl,i)*(tt1(ibl,n,k)*bdiag1(ibl,l)+tt1(ibl,l,k)*bdiag1(ibl,n)) &
              -csdiag1(ibl,l)*(aaa1(ibl,n,i)*diag1(ibl,k)+aaa1(ibl,n,k)*diag1(ibl,i)) &
              -bsdiag1(ibl,l)*(tt1(ibl,i,n)*bdiag1(ibl,k)+tt1(ibl,k,n)*bdiag1(ibl,i)) &
              -csdiag1(ibl,k)*(aaa1(ibl,i,n)*diag1(ibl,l)+aaa1(ibl,i,l)*diag1(ibl,n)) &
              -bsdiag1(ibl,k)*(tt1(ibl,n,i)*bdiag1(ibl,l)+tt1(ibl,l,i)*bdiag1(ibl,n)) &
              -csdiag1(ibl,n)*(aaa1(ibl,l,i)*diag1(ibl,k)+aaa1(ibl,l,k)*diag1(ibl,i)) &
              -bsdiag1(ibl,n)*(tt1(ibl,i,l)*bdiag1(ibl,k)+tt1(ibl,k,l)*bdiag1(ibl,i)))
        enddo
        etotb=etotb+g(intg)*etemp
      enddo
#ifdef OMP
!$omp end parallel do
#endif
    enddo

    !     Reset index into buffer
    nb1=0

    call timer_stop(35)

  endif

  return
end subroutine gronor_gntwo_omp_batch_canonical
