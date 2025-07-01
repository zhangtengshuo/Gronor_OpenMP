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

subroutine gronor_gntwo_batch_indexed(lfndbg,ihc,nhc)

  use mpi
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals

  implicit none

  external :: timer_start,timer_stop

  integer :: lfndbg,i,ii,jj,k,l,n,kl,intg,ihc,nhc,ibl

  real (kind=8) :: ts,fourdet

  real(kind=8), external :: timer_wall

  real (kind=8) :: aaamax,aatmax
  real (kind=8) :: diagmax,bdiagmax,bsdiagmax,csdiagmax,sdiagmax

  logical :: ldiag,lbdiag

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
!$acc kernels present(aaa0,aat0,tt0,ta0,tt,aat,aaa,ta,sm0) copyin(nb0)
      do k=1,nbas
        do i=1,nbas
          sm0(i,k,nb0)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa0(i,k,nb0)=aaa(i,k)
          aat0(i,k,nb0)=aaa(k,i)
          ta0(i,k,nb0)=ta(i,k)
          tt0(i,k,nb0)=ta(k,i)
        enddo
      enddo
!$acc end kernels
!     END ASSEMBLY ON GPU

    else
      nb1=nb1+1
      prefac1(nb1)=fctr

!     ASSEMBLE ARRAYS ON GPU
!$acc kernels present(aaa1,aat1,tt1,ta1,tt,aat,aaa,ta,sm1) copyin(nb1)
      do k=1,nbas
        do i=1,nbas
          sm1(i,k,nb1)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa1(i,k,nb1)=aaa(i,k)
          aat1(i,k,nb1)=aaa(k,i)
          ta1(i,k,nb1)=ta(i,k)
          tt1(i,k,nb1)=ta(k,i)
        enddo
      enddo
!$acc end kernels
      !     END ASSEMBLY ON GPU

!$acc kernels present(diag,bdiag,bsdiag,csdiag) present(diag1,bdiag1,bsdiag1,csdiag1) copyin(nb1)
      do i=1,nbas
        diag1(nb1,i)=diag(i)
        bdiag1(nb1,i)=bdiag(i)
        bsdiag1(nb1,i)=bsdiag(i)
        csdiag1(nb1,i)=csdiag(i)
      enddo
!$acc end kernels

    endif

    call timer_stop(30)

  endif

  call timer_start(38)

  !     After the last element of a task the list of energies is generated

  if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.nbatch) then

    call timer_start(31)

    !!_ACCTGT_($acc update device(ta0,aaa0,aat0,tt0,sm0))
    !!_OMPTGT_($omp target update to(ta0,aaa0,aat0,tt0,sm0))

!$acc data copyin(prefac0)
    kl=nbas*(nbas+1)/2
!$acc kernels present(lab,ndx,g,sm0,aat0,aaa0,tt0,ta0,prefac0)
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
!$acc end kernels
    !     Reset index in buffer to zero
    nb0=0

!$acc end data

    call timer_stop(31)

  endif

  if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.nbatch) then

    call timer_start(34)

    !!_ACCTGT_($acc update device(ta1,aaa1,aat1,tt1,sm1))
    !!_OMPTGT_($omp target update to(ta1,aaa1,aat1,tt1,sm1))

!$acc data copyin(prefac1)

    kl=nbas*(nbas+1)/2

!$acc kernels present(lab,ndx,g,sm1,aat1,aaa1,tt1,ta1,prefac1) present(diag1,bdiag1,bsdiag1,csdiag1)
    do ii=intndx,jntndx
      do jj=ii,kl
        intg=ndx(ii)+jj
        i=lab(1,ii)
        k=lab(2,ii)
        l=lab(1,jj)
        n=lab(2,jj)
!$acc loop vector(32)
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
!$acc end kernels
    !     Reset index into buffer
    nb1=0

    call timer_stop(34)

!$acc end data

  endif

  call timer_stop(38)

  return
end subroutine gronor_gntwo_batch_indexed

