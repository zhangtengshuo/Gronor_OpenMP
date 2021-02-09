
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
!! <tr><td> gronor_gntwo <td> \htmlonly &#x2714;\endhtmlonly <td>  <td> \htmlonly &#x2714;\endhtmlonly <td>
!! <tr><td> gronor_gntwo_mlab <td> \htmlonly &#x2714;\endhtmlonly <td>  <td>  <td>
!! <tr><td> gronor_gntwo_batch <td> \htmlonly &#x2714;\endhtmlonly <td>  <td> \htmlonly &#x2714;\endhtmlonly <td>\htmlonly &#x2714;\endhtmlonly
!! <tr><td> gronor_gntwo_batch_mlab <td> \htmlonly &#x2714;\endhtmlonly <td>  <td>  <td>\htmlonly &#x2714;\endhtmlonly
!! <tr><td> gronor_gntwo_omp <td> <td> \htmlonly &#x2714;\endhtmlonly <td>  \htmlonly &#x2714;\endhtmlonly <td>
!! <tr><td> gronor_gntwo_omp_mlab <td>  <td> \htmlonly &#x2714;\endhtmlonly  <td>  <td>
!! <tr><td> gronor_gntwo_omp_batch <td> <td> \htmlonly &#x2714;\endhtmlonly  <td> \htmlonly &#x2714;\endhtmlonly <td>\htmlonly &#x2714;\endhtmlonly
!! <tr><td> gronor_gntwo_omp_batch_mlab <td>  <td> \htmlonly &#x2714;\endhtmlonly <td>  <td>\htmlonly &#x2714;\endhtmlonly
!! </table>
!!

#include "tgt_guard.fh"
      
      subroutine gronor_gntwo(lfndbg)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENACC
      use openacc
      use cuda_functions
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg
      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg

      integer :: iaf0,ial0,inf0,inl0,iaf1,ial1,inf1,inl1

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet

#ifdef SINGLEP
      real (kind=4) :: e2t,e2nt,tst,tsnt
      real (kind=4) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=4) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#else
      real (kind=8) :: e2t,e2nt,tst,tsnt
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
 600  format(' Calculating two electron matrix element',/)

!      numint=ig(mod(me,mgr)+1)

      numint=mint2

      if(ising.gt.1) e1=0.0d0
      e2=0.0d0
      ts=0.0d0
      e2n=0.0d0
      tsn=0.0d0
      fourdet=4.0d0*deta

!_ACCTGT_($acc kernels present(aat,aaa,tt,ta,tat,sm))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
      do i=1,nbas 
        do k=1,nbas
          tat(k,i)=ta(k,i)
          tt(i,k)=tat(k,i)
          aat(i,k)=aaa(k,i)
        enddo
      enddo
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
      do i=1,nbas 
        do k=1,nbas
          sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+tat(k,i)
        enddo
      enddo
!_ACCTGT_($acc end kernels)
      call timer_stop(30)

      call timer_start(38)

      if(ising.le.0) then
        
        call timer_start(31)
        
        tst=ts
!_ACCTGT_($acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab))
!_ACCTGT_($acc loop reduction (+:tst))
!_OMPTGT_($omp taskwait)
!_OMPTGT_($omp target teams distribute parallel do reduction(+:tst))
        do intg=1,numint
          i=lab(1,intg)
          k=lab(2,intg)
          l=lab(3,intg)
          n=lab(4,intg)
          tst=tst+g(intg)*(sm(k,i)*sm(n,l)
     &         -aat(n,i)*aaa(l,k)-tt(n,i)*tat(l,k)
     &         -aaa(n,i)*aat(l,k)-tat(n,i)*tt(l,k)
     &         -aat(l,i)*aaa(n,k)-tt(l,i)*tat(n,k)
     &         -aaa(l,i)*aat(n,k)-tat(l,i)*tt(n,k))
        enddo
!_ACCTGT_($acc end kernels)
        ts=tst

        call timer_stop(31)

      else
        
        call timer_start(34)

        e2t=e2
!_ACCTGT_($acc kernels present(aat,aaa,tt,ta,tat,sm,g))
!_ACCTGT_($acc& present(lab,diag,bdiag,bsdiag,csdiag))
!_ACCTGT_($acc loop reduction (+:e2t))
!_OMPTGT_($omp target teams distribute parallel do reduction(+:e2t))
        do intg=1,numint
          i=lab(1,intg) 
          k=lab(2,intg) 
          l=lab(3,intg) 
          n=lab(4,intg)
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
          e2t=e2t+g(intg)*
     &         (sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l)
     &         -bai*(aat(n,k)*aaj+aat(l,k)*aal)
     &         -bbi*(tt(n,k)*abj+tt(l,k)*abl)
     &         -baj*(aaa(n,i)*aak+aaa(n,k)*aai)
     &         -bbj*(tat(n,i)*abk+tat(n,k)*abi)
     &         -bak*(aat(n,i)*aaj+aat(l,i)*aal)
     &         -bbk*(tt(n,i)*abj+tt(l,i)*abl)
     &         -bal*(aaa(l,i)*aak+aaa(l,k)*aai)
     &         -bbl*(tat(l,i)*abk+tat(l,k)*abi))
        enddo
!_ACCTGT_($acc end kernels)
        e2=e2t

        call timer_stop(34)

      endif

#ifdef _OPENACC
      if(numdev.gt.1) then
        cpfre=c_loc(memfre)
        cptot=c_loc(memtot)
        istat=cudaMemGetInfo(cpfre,cptot)
        memavail=min(memfre,memavail)
      endif
#endif

      ts=ts*fourdet
      e2=e2+ts
      
      call timer_stop(38)
      
      return
      end
      
      subroutine gronor_gntwo_mlab(lfndbg)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENACC
      use openacc
      use cuda_functions
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg
      integer :: i,ii,jj,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg

      integer :: iaf0,ial0,inf0,inl0,iaf1,ial1,inf1,inl1

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet

#ifdef SINGLEP
      real (kind=4) :: e2t,e2nt,tst,tsnt
      real (kind=4) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=4) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#else
      real (kind=8) :: e2t,e2nt,tst,tsnt
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
 600  format(' Calculating two electron matrix element',/)

!      numint=ig(mod(me,mgr)+1)

      numint=mint2

      if(ising.gt.1) e1=0.0d0
      e2=0.0d0
      ts=0.0d0
      e2n=0.0d0
      tsn=0.0d0
      fourdet=4.0d0*deta

!_ACCTGT_($acc kernels present(aat,aaa,tt,ta,tat,sm))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
      do i=1,nbas 
        do k=1,nbas
          tat(k,i)=ta(k,i)
          tt(i,k)=tat(k,i)
          aat(i,k)=aaa(k,i)
        enddo
      enddo
!_OMPTGT_($omp target teams distribute parallel do collapse(2))      
      do i=1,nbas 
        do k=1,nbas
          sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+tat(k,i)
        enddo
      enddo
!_ACCTGT_($acc end kernels)

      call timer_stop(30)
      
      call timer_start(38)

      if(ising.le.0) then
        
        call timer_start(31)

        tst=ts
        kl=nbas*(nbas+1)/2

!_ACCTGT_($acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab,ndx))
!_ACCTGT_($acc& copyin(kl,intndx,jntndx))
!_OMPTGT_($omp taskwait)
!_OMPTGT_($omp target teams distribute parallel do reduction(+:tst))
        do ii=intndx,jntndx
          do jj=ii,kl
            intg=ndx(ii)+jj
            i=lab(1,ii)
            k=lab(2,ii)
            l=lab(1,jj)
            n=lab(2,jj)
            tst=tst+g(intg)*(sm(k,i)*sm(n,l)
     &           -aat(n,i)*aaa(l,k)-tt(n,i)*tat(l,k)
     &           -aaa(n,i)*aat(l,k)-tat(n,i)*tt(l,k)
     &           -aat(l,i)*aaa(n,k)-tt(l,i)*tat(n,k)
     &           -aaa(l,i)*aat(n,k)-tat(l,i)*tt(n,k))
          enddo
        enddo
!_ACCTGT_($acc end kernels)
        ts=tst

        call timer_stop(31)

      else
        
        call timer_start(34)

        e2t=e2
        kl=nbas*(nbas+1)/2

!_ACCTGT_($acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab,ndx))
!_ACCTGT_($acc& present(diag,bdiag,bsdiag,csdiag) copyin(kl,intndx,jntndx))
!_OMPTGT_($omp target teams distribute parallel do reduction(+:e2t))
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
            e2t=e2t+g(intg)*
     &           (sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l)
     &           -bai*(aat(n,k)*aaj+aat(l,k)*aal)
     &           -bbi*(tt(n,k)*abj+tt(l,k)*abl)
     &           -baj*(aaa(n,i)*aak+aaa(n,k)*aai)
     &           -bbj*(tat(n,i)*abk+tat(n,k)*abi)
     &           -bak*(aat(n,i)*aaj+aat(l,i)*aal)
     &           -bbk*(tt(n,i)*abj+tt(l,i)*abl)
     &           -bal*(aaa(l,i)*aak+aaa(l,k)*aai)
     &           -bbl*(tat(l,i)*abk+tat(l,k)*abi))
          enddo
        enddo
!_ACCTGT_($acc end kernels)
        e2=e2t

        call timer_stop(34)

      endif

#ifdef _OPENACC
      if(numdev.gt.1) then
        cpfre=c_loc(memfre)
        cptot=c_loc(memtot)
        istat=cudaMemGetInfo(cpfre,cptot)
        memavail=min(memfre,memavail)
      endif
#endif

      ts=ts*fourdet
      e2=e2+ts
      
      call timer_stop(38)
      
      return
      end
      
      subroutine gronor_gntwo_batch(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg,ihc,nhc,ibl

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

      etotb=0.0d0
      
!     Copying arrays into batch arrays

!     Initialize nb0 and nb1 when first matrix element in a batch

!     nb0 : number of singularity 0 elements
!     nb1 : number of singularity 1 or 2 elements
      
      if(ihc.eq.1) then
        nb0=0
        nb1=0
      endif
      
      if(ising.lt.3) then
        
        call timer_start(30)
        
        if(ising.eq.0) then
          nb0=nb0+1
          ibl=nb0        
          prefac(ibl)=4.0d0*deta*fctr
        else
          ibl=nhc-nb1
          nb1=nb1+1
          prefac(ibl)=fctr
        endif
        
        if(iamhead.eq.1) sstot=sstot+fctr*deta
        
        if(idbg.ge.13) write(lfndbg,600)
 600    format(' Calculating two electron matrix element',/)
        
        numint=mint2
        
        if(ising.gt.1) e1=0.0d0
        e2=0.0d0
        ts=0.0d0
        fourdet=4.0d0*deta
        
        if(ising.gt.0) then

!_ACCTGT_($acc kernels present(diag,bdiag,bsdiag,csdiag))
!_ACCTGT_($acc& present(diagl,bdiagl,bsdiagl,csdiagl) copyin(ibl))
!_OMPTGT_($omp target teams distribute parallel do)
          do i=1,nbas
            diagl(ibl,i)=diag(i)
            bdiagl(ibl,i)=bdiag(i)
            bsdiagl(ibl,i)=bsdiag(i)
            csdiagl(ibl,i)=csdiag(i)
          enddo
!_ACCTGT_($acc end kernels)

        endif
!_ACCTGT_($acc kernels present(aat,aaa,tt,ta))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
        do i=1,nbas 
          do k=1,nbas 
            tt(i,k)=ta(k,i)
            aat(i,k)=aaa(k,i)
          enddo
        enddo
!_ACCTGT_($acc end kernels)

!     ASSEMBLE ARRAYS ON GPU
      
!_ACCTGT_($acc kernels present(aaal,aatl,ttl,tatl,tt,aat,aaa,ta,sml) copyin(ibl))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
        do k=1,nbas
          do i=1,nbas
            sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
            aaal(ibl,i,k)=aaa(i,k)
            aatl(ibl,i,k)=aat(i,k)
            ttl(ibl,i,k)=tt(i,k)
            tatl(ibl,i,k)=ta(i,k)
          enddo
        enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt)) 

c      do k=1,nbas
c        do i=1,nbas
c          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
c          aaal(ibl,i,k)=aaa(i,k)
c          aatl(ibl,i,k)=aat(i,k)
c          ttl(ibl,i,k)=tt(i,k)
c          tatl(ibl,i,k)=ta(i,k)
c        enddo
c      enddo

c!     END ASSEMBLY ON CPU

        call timer_stop(30)

      endif

      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if(ihc.eq.nhc) then

c!_ACCTGT_($acc update device(tatl,aaal,aatl,ttl,sml))
c!_OMPTGT_($omp target update to(tatl,aaal,aatl,ttl,sml))

        etotb=0.0d0
        
!_ACCTGT_($acc data copyin(prefac))
!_OMPTGT_($omp target data map(to:prefac))        

        if(nb0.gt.0) then
          call timer_start(31)

!_ACCTGT_($acc kernels present(lab,g,sml,aatl,aaal,ttl,tatl,prefac))
!_OMPTGT_($omp target teams distribute parallel do)
          do intg=1,numint
            i=lab(1,intg)
            k=lab(2,intg)
            l=lab(3,intg)
            n=lab(4,intg)
!_ACCTGT_($acc loop private(ibl) vector(32))
            do ibl=1,nb0
              etotb=etotb+g(intg)*prefac(ibl)*(sml(ibl,k,i)*sml(ibl,n,l)
     &             -aatl(ibl,n,i)*aaal(ibl,l,k)
     &             -ttl(ibl,n,i)*tatl(ibl,l,k)
     &             -aaal(ibl,n,i)*aatl(ibl,l,k)
     &             -tatl(ibl,n,i)*ttl(ibl,l,k)
     &             -aatl(ibl,l,i)*aaal(ibl,n,k)
     &             -ttl(ibl,l,i)*tatl(ibl,n,k)
     &             -aaal(ibl,l,i)*aatl(ibl,n,k)
     &             -tatl(ibl,l,i)*ttl(ibl,n,k))
            enddo
          enddo
!_ACCTGT_($acc end kernels)

          call timer_stop(31)
          
        endif

        if(nb1.gt.0) then
          
          call timer_start(34)

!_ACCTGT_($acc kernels present(lab,g,sml,aatl,aaal,ttl,tatl,prefac))
!_ACCTGT_($acc& present(diagl,bdiagl,bsdiagl,csdiagl))
!_OMPTGT_($omp target teams distribute parallel do)
          do intg=1,numint
            i=lab(1,intg) 
            k=lab(2,intg) 
            l=lab(3,intg) 
            n=lab(4,intg)
!_ACCTGT_($acc loop private(ibl) vector(32))
            do ibl=nhc-nb1+1,nhc
              etotb=etotb+g(intg)*prefac(ibl)*
     &             (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
     &             +diagl(ibl,n)*csdiagl(ibl,l)
     &             +bdiagl(ibl,l)*bsdiagl(ibl,n)
     &             +bdiagl(ibl,n)*bsdiagl(ibl,l))
     &             +(diagl(ibl,i)*csdiagl(ibl,k)
     &             +diagl(ibl,k)*csdiagl(ibl,i)
     &             +bdiagl(ibl,i)*bsdiagl(ibl,k)
     &             +bdiagl(ibl,k)*bsdiagl(ibl,i))*sml(ibl,n,l)
     &             -csdiagl(ibl,i)*(aatl(ibl,n,k)*diagl(ibl,l)
     &             +aatl(ibl,l,k)*diagl(ibl,n))
     &             -bsdiagl(ibl,i)*(ttl(ibl,n,k)*bdiagl(ibl,l)
     &             +ttl(ibl,l,k)*bdiagl(ibl,n))
     &             -csdiagl(ibl,l)*(aaal(ibl,n,i)*diagl(ibl,k)
     &             +aaal(ibl,n,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,l)*(tatl(ibl,n,i)*bdiagl(ibl,k)
     &             +tatl(ibl,n,k)*bdiagl(ibl,i))
     &             -csdiagl(ibl,k)*(aatl(ibl,n,i)*diagl(ibl,l)
     &             +aatl(ibl,l,i)*diagl(ibl,n))
     &             -bsdiagl(ibl,k)*(ttl(ibl,n,i)*bdiagl(ibl,l)
     &             +ttl(ibl,l,i)*bdiagl(ibl,n))
     &             -csdiagl(ibl,n)*(aaal(ibl,l,i)*diagl(ibl,k)
     &             +aaal(ibl,l,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,n)*(tatl(ibl,l,i)*bdiagl(ibl,k)
     &             +tatl(ibl,l,k)*bdiagl(ibl,i)))
            enddo
          enddo
!_ACCTGT_($acc end kernels)
    
          call timer_stop(34)
          
        endif
        
!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)
      endif

      call timer_stop(38)
      
      return
      end

      subroutine gronor_gntwo_batch_mlab(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,jj,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg,ihc,nhc,ibl

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

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
        
        if(idbg.ge.13) write(lfndbg,600)
 600    format(' Calculating two electron matrix element',/)

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
      
!_ACCTGT_($acc kernels present(aaal,aatl,ttl,tatl,tt,aat,aaa,ta,sml) copyin(ibl))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
          do k=1,nbas
            do i=1,nbas
              sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa0(nb0,i,k)=aaa(i,k)
              aat0(nb0,i,k)=aaa(k,i)
              tt0(nb0,i,k)=ta(k,i)
              tat0(nb0,i,k)=ta(i,k)
            enddo
          enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt))     
c      do k=1,nbas
c        do i=1,nbas
c          sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa0(nb0,i,k)=aaa(i,k)
c          aat0(nb0,i,k)=aaa(k,i)
c          tt0(nb0,i,k)=ta(k,i)
c          tat0(nb0,i,k)=ta(i,k)
c        enddo
c      enddo

c!     END ASSEMBLY ON CPU
        else
          nb1=nb1+1
          prefac1(nb1)=fctr
          
!     ASSEMBLE ARRAYS ON GPU
      
!_ACCTGT_($acc kernels present(aaal,aatl,ttl,tatl,tt,aat,aaa,ta,sml) copyin(ibl))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
          do k=1,nbas
            do i=1,nbas
              sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa1(nb1,i,k)=aaa(i,k)
              aat1(nb1,i,k)=aaa(k,i)
              tt1(nb1,i,k)=ta(k,i)
              tat1(nb1,i,k)=ta(i,k)
            enddo
          enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt))     
c      do k=1,nbas
c        do i=1,nbas
c          sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa1(nb1,i,k)=aaa(i,k)
c          aat1(nb1,i,k)=aaa(k,i)
c          tt1(nb1,i,k)=ta(k,i)
c          tat1(nb1,i,k)=ta(i,k)
c        enddo
c      enddo

c!     END ASSEMBLY ON CPU
      
!_ACCTGT_($acc kernels present(diag,bdiag,bsdiag,csdiag))
!_ACCTGT_($acc& present(diagl,bdiagl,bsdiagl,csdiagl) copyin(ibl))
!_OMPTGT_($omp target teams distribute parallel do)
          do i=1,nbas
            diag1(nb1,i)=diag(i)
            bdiag1(nb1,i)=bdiag(i)
            bsdiag1(nb1,i)=bsdiag(i)
            csdiag1(nb1,i)=csdiag(i)
          enddo
!_ACCTGT_($acc end kernels)

        endif

        call timer_stop(30)

      endif

      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.lbatch) then

        call timer_start(31)

c!_ACCTGT_($acc update device(tatl,aaal,aatl,ttl,sml))
c!_OMPTGT_($omp target update to(tatl,aaal,aatl,ttl,sml)) 

!_ACCTGT_($acc data copyin(prefac0))
!_OMPTGT_($omp target data map(to:prefac0))

        kl=nbas*(nbas+1)/2

!_ACCTGT_($acc kernels present(lab,ndx,g,sml,aatl,aaal,ttl,tatl,prefac))
!_OMPTGT_($omp target teams distribute parallel do)
        do ii=intndx,jntndx
          do jj=ii,kl
            intg=ndx(ii)+jj
            i=lab(1,ii)
            k=lab(2,ii)
            l=lab(1,jj)
            n=lab(2,jj)
!_ACCTGT_($acc loop private(ibl) vector(32))
            do ibl=1,nb0
              etotb=etotb+g(intg)*
     &             prefac(ibl)*(sm0(ibl,k,i)*sm0(ibl,n,l)
     &             -aat0(ibl,n,i)*aaa0(ibl,l,k)
     &             -tt0(ibl,n,i)*tat0(ibl,l,k)
     &             -aaa0(ibl,n,i)*aat0(ibl,l,k)
     &             -tat0(ibl,n,i)*tt0(ibl,l,k)
     &             -aat0(ibl,l,i)*aaa0(ibl,n,k)
     &             -tt0(ibl,l,i)*tat0(ibl,n,k)
     &             -aaa0(ibl,l,i)*aat0(ibl,n,k)
     &             -tat0(ibl,l,i)*tt0(ibl,n,k))
            enddo
          enddo
        enddo
!_ACCTGT_($acc end kernels)

!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)

!     Reset index in buffer to zero
        nb0=0
        
        call timer_stop(31)
        
      endif

      if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.lbatch) then
        
        call timer_start(34)
        
!_ACCTGT_($acc data copyin(prefac0))
!_OMPTGT_($omp target data map(to:prefac0))

        kl=nbas*(nbas+1)/2

!_ACCTGT_($acc kernels present(lab,ndx,g,sml,aatl,aaal,ttl,tatl,prefac))
!_ACCTGT_($acc& present(diagl,bdiagl,bsdiagl,csdiagl))
!_OMPTGT_($omp target teams distribute parallel do)
        do ii=intndx,jntndx
          do jj=ii,kl
            intg=ndx(ii)+jj
            i=lab(1,ii)
            k=lab(2,ii)
            l=lab(1,jj)
            n=lab(2,jj)
!_ACCTGT_($acc loop private(ibl) vector(32))
            do ibl=1,nb1
              etotb=etotb+g(intg)*prefac(ibl)*
     &             (sm1(ibl,k,i)*(diag1(ibl,l)*csdiag1(ibl,n)
     &             +diag1(ibl,n)*csdiag1(ibl,l)
     &             +bdiag1(ibl,l)*bsdiag1(ibl,n)
     &             +bdiag1(ibl,n)*bsdiag1(ibl,l))
     &             +(diag1(ibl,i)*csdiag1(ibl,k)
     &             +diag1(ibl,k)*csdiag1(ibl,i)
     &             +bdiag1(ibl,i)*bsdiag1(ibl,k)
     &             +bdiag1(ibl,k)*bsdiag1(ibl,i))*sm1(ibl,n,l)
     &             -csdiag1(ibl,i)*(aat1(ibl,n,k)*diag1(ibl,l)
     &             +aat1(ibl,l,k)*diag1(ibl,n))
     &             -bsdiag1(ibl,i)*(tt1(ibl,n,k)*bdiag1(ibl,l)
     &             +tt1(ibl,l,k)*bdiag1(ibl,n))
     &             -csdiag1(ibl,l)*(aaa1(ibl,n,i)*diag1(ibl,k)
     &             +aaa1(ibl,n,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,l)*(tat1(ibl,n,i)*bdiag1(ibl,k)
     &             +tat1(ibl,n,k)*bdiag1(ibl,i))
     &             -csdiag1(ibl,k)*(aat1(ibl,n,i)*diag1(ibl,l)
     &             +aat1(ibl,l,i)*diag1(ibl,n))
     &             -bsdiag1(ibl,k)*(tt1(ibl,n,i)*bdiag1(ibl,l)
     &             +tt1(ibl,l,i)*bdiag1(ibl,n))
     &             -csdiag1(ibl,n)*(aaa1(ibl,l,i)*diag1(ibl,k)
     &             +aaa1(ibl,l,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,n)*(tat1(ibl,l,i)*bdiag1(ibl,k)
     &             +tat1(ibl,l,k)*bdiag1(ibl,i)))
            enddo
          enddo
        enddo
!_ACCTGT_($acc end kernels)
    
!     Reset index into buffer
        nb1=0

!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)

        call timer_stop(34)
        
      endif

      call timer_stop(38)
      
      return
      end

      subroutine gronor_gntwo_batch2(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg,ihc,nhc,ibl

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

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
        
        if(idbg.ge.13) write(lfndbg,600)
 600    format(' Calculating two electron matrix element',/)
        
        numint=mint2
        
        if(ising.gt.1) e1=0.0d0
        e2=0.0d0
        ts=0.0d0
        fourdet=4.0d0*deta
        
        if(ising.eq.0) then
          nb0=nb0+1
          prefac0(nb0)=4.0d0*deta*fctr
          
!     ASSEMBLE ARRAYS ON GPU
      
!_ACCTGT_($acc kernels present(aaal,aatl,ttl,tatl,tt,aat,aaa,ta,sml) copyin(ibl))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
          do k=1,nbas
            do i=1,nbas
              sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa0(nb0,i,k)=aaa(i,k)
              aat0(nb0,i,k)=aaa(k,i)
              tt0(nb0,i,k)=ta(k,i)
              tat0(nb0,i,k)=ta(i,k)
            enddo
          enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt)) 

c      do k=1,nbas
c        do i=1,nbas
c          sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa0(nb0,i,k)=aaa(i,k)
c          aat0(nb0,i,k)=aaa(k,i)
c          tt0(nb0,i,k)=ta(k,i)
c          tat0(nb0,i,k)=ta(i,k)
c        enddo
c      enddo

c!     END ASSEMBLY ON CPU
        else
          nb1=nb1+1
          prefac1(nb1)=fctr
!     ASSEMBLE ARRAYS ON GPU
      
!_ACCTGT_($acc kernels present(aaal,aatl,ttl,tatl,tt,aat,aaa,ta,sml) copyin(ibl))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
          do k=1,nbas
            do i=1,nbas
              sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa1(nb1,i,k)=aaa(i,k)
              aat1(nb1,i,k)=aaa(k,i)
              tt1(nb1,i,k)=ta(k,i)
              tat1(nb1,i,k)=ta(i,k)
            enddo
          enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt)) 

c      do k=1,nbas
c        do i=1,nbas
c          sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa1(nb1,i,k)=aaa(i,k)
c          aat1(nb1,i,k)=aaa(k,i)
c          tt1(nb1,i,k)=ta(k,i)
c          tat1(nb1,i,k)=ta(i,k)
c        enddo
c      enddo

c!     END ASSEMBLY ON CPU

!_ACCTGT_($acc kernels present(diag,bdiag,bsdiag,csdiag))
!_ACCTGT_($acc& present(diagl,bdiagl,bsdiagl,csdiagl) copyin(ibl))
!_OMPTGT_($omp target teams distribute parallel do)
          do i=1,nbas
            diag1(nb1,i)=diag(i)
            bdiag1(nb1,i)=bdiag(i)
            bsdiag1(nb1,i)=bsdiag(i)
            csdiag1(nb1,i)=csdiag(i)
          enddo
!_ACCTGT_($acc end kernels)

        endif

        call timer_stop(30)

      endif
      
      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.lbatch) then
        
        call timer_start(31)

c!_ACCTGT_($acc update device(tat0,aaa0,aat0,tt0,sm0))
c!_OMPTGT_($omp target update to(tat0,aaa0,aat0,tt0,sm0))
        
!_ACCTGT_($acc data copyin(prefac0))
!_OMPTGT_($omp target data map(to:prefac0))

!_ACCTGT_($acc kernels present(lab,g,sm0,aat0,aaa0,tt0,tat0,prefac0))
!_OMPTGT_($omp target teams distribute parallel do)
        do intg=1,numint
          i=lab(1,intg)
          k=lab(2,intg)
          l=lab(3,intg)
          n=lab(4,intg)
!_ACCTGT_($acc loop private(ibl) vector(32))
          do ibl=1,nb0
            etotb=etotb+g(intg)*prefac0(ibl)*(sm0(ibl,k,i)*sm0(ibl,n,l)
     &           -aat0(ibl,n,i)*aaa0(ibl,l,k)
     &           -tt0(ibl,n,i)*tat0(ibl,l,k)
     &           -aaa0(ibl,n,i)*aat0(ibl,l,k)
     &           -tat0(ibl,n,i)*tt0(ibl,l,k)
     &           -aat0(ibl,l,i)*aaa0(ibl,n,k)
     &           -tt0(ibl,l,i)*tat0(ibl,n,k)
     &           -aaa0(ibl,l,i)*aat0(ibl,n,k)
     &           -tat0(ibl,l,i)*tt0(ibl,n,k))
          enddo
        enddo
!_ACCTGT_($acc end kernels)
        
!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)

!     Reset index in buffer to zero
        nb0=0
        
        call timer_stop(31)
        
      endif
        
      if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.lbatch) then
        
        call timer_start(34)

c!_ACCTGT_($acc update device(tat1,aaa1,aat1,tt1,sm1))
c!_OMPTGT_($omp target update to(tat1,aaa1,aat1,tt1,sm1))
        
!_ACCTGT_($acc data copyin(prefac1))
!_OMPTGT_($omp target data map(to:prefac1))

!_ACCTGT_($acc kernels present(lab,g,sm1,aat1,aaa1,tt1,tat1,prefac1))
!_ACCTGT_($acc& present(diagl,bdiagl,bsdiagl,csdiagl))
!_OMPTGT_($omp target teams distribute parallel do)
        do intg=1,numint
          i=lab(1,intg) 
          k=lab(2,intg) 
          l=lab(3,intg) 
          n=lab(4,intg)
!_ACCTGT_($acc loop private(ibl) vector(32))
          do ibl=1,nb1
            etotb=etotb+g(intg)*prefac1(ibl)*
     &           (sm1(ibl,k,i)*(diag1(ibl,l)*csdiag1(ibl,n)
     &           +diag1(ibl,n)*csdiag1(ibl,l)
     &           +bdiag1(ibl,l)*bsdiag1(ibl,n)
     &           +bdiag1(ibl,n)*bsdiag1(ibl,l))
     &           +(diag1(ibl,i)*csdiag1(ibl,k)
     &           +diag1(ibl,k)*csdiag1(ibl,i)
     &           +bdiag1(ibl,i)*bsdiag1(ibl,k)
     &           +bdiag1(ibl,k)*bsdiag1(ibl,i))*sm1(ibl,n,l)
     &           -csdiag1(ibl,i)*(aat1(ibl,n,k)*diag1(ibl,l)
     &           +aat1(ibl,l,k)*diag1(ibl,n))
     &           -bsdiag1(ibl,i)*(tt1(ibl,n,k)*bdiag1(ibl,l)
     &           +tt1(ibl,l,k)*bdiag1(ibl,n))
     &           -csdiag1(ibl,l)*(aaa1(ibl,n,i)*diag1(ibl,k)
     &           +aaa1(ibl,n,k)*diag1(ibl,i))
     &           -bsdiag1(ibl,l)*(tat1(ibl,n,i)*bdiag1(ibl,k)
     &           +tat1(ibl,n,k)*bdiag1(ibl,i))
     &           -csdiag1(ibl,k)*(aat1(ibl,n,i)*diag1(ibl,l)
     &           +aat1(ibl,l,i)*diag1(ibl,n))
     &           -bsdiag1(ibl,k)*(tt1(ibl,n,i)*bdiag1(ibl,l)
     &           +tt1(ibl,l,i)*bdiag1(ibl,n))
     &           -csdiag1(ibl,n)*(aaa1(ibl,l,i)*diag1(ibl,k)
     &           +aaa1(ibl,l,k)*diag1(ibl,i))
     &           -bsdiag1(ibl,n)*(tat1(ibl,l,i)*bdiag1(ibl,k)
     &           +tat1(ibl,l,k)*bdiag1(ibl,i)))
          enddo
        enddo
!_ACCTGT_($acc end kernels)
    
!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)

!     Reset index into buffer
        nb1=0
          
        call timer_stop(34)
      endif

      call timer_stop(38)
      
      return
      end

      subroutine gronor_gntwo_batch2_mlab(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,jj,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg,ihc,nhc,ibl

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

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
        
        if(idbg.ge.13) write(lfndbg,600)
 600    format(' Calculating two electron matrix element',/)
        
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
      
!_ACCTGT_($acc kernels present(aaa0,aat0,tt0,tat0,tt,aat,aaa,ta,sm0) copyin(nb0))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
      do k=1,nbas
        do i=1,nbas
          sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa0(nb0,i,k)=aaa(i,k)
          aat0(nb0,i,k)=aaa(k,i)
          tt0(nb0,i,k)=ta(k,i)
          tat0(nb0,i,k)=ta(i,k)
        enddo
      enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt))     
c      do k=1,nbas
c        do i=1,nbas
c          sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa0(nb0,i,k)=aaa(i,k)
c          aat0(nb0,i,k)=aaa(k,i)
c          tt0(nb0,i,k)=ta(k,i)
c          tat0(nb0,i,k)=ta(i,k)
c        enddo
c      enddo

c!     END ASSEMBLY ON CPU

      else
        nb1=nb1+1
        prefac1(nb1)=fctr
      
!     ASSEMBLE ARRAYS ON GPU
      
!_ACCTGT_($acc kernels present(aaal,aatl,ttl,tatl,tt,aat,aaa,ta,sml) copyin(nb1))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
      do k=1,nbas
        do i=1,nbas
          sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa1(nb1,i,k)=aaa(i,k)
          aat1(nb1,i,k)=aaa(k,i)
          tt1(nb1,i,k)=ta(k,i)
          tat1(nb1,i,k)=ta(i,k)
        enddo
      enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt))     
c      do k=1,nbas
c        do i=1,nbas
c          sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa1(nb1,i,k)=aaa(i,k)
c          aat1(nb1,i,k)=aaa(k,i)
c          tt1(nb1,i,k)=ta(k,i)
c          tat1(nb1,i,k)=ta(i,k)
c        enddo
c      enddo

c!     END ASSEMBLY ON CPU

!_ACCTGT_($acc kernels present(diag,bdiag,bsdiag,csdiag))
!_ACCTGT_($acc& present(diagl,bdiagl,bsdiagl,csdiagl) copyin(nb1))
!_OMPTGT_($omp target teams distribute parallel do)
        do i=1,nbas
          diag1(nb1,i)=diag(i)
          bdiag1(nb1,i)=bdiag(i)
          bsdiag1(nb1,i)=bsdiag(i)
          csdiag1(nb1,i)=csdiag(i)
        enddo
!_ACCTGT_($acc end kernels)

      endif

      call timer_stop(30)

      endif

      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.lbatch) then
        
        call timer_start(31)

c!_ACCTGT_($acc update device(tatl,aaal,aatl,ttl,sml))
c!_OMPTGT_($omp target update to(tatl,aaal,aatl,ttl,sml)) 

!_ACCTGT_($acc data copyin(prefac0))
!_OMPTGT_($omp target data map(to:prefac0))
        
        kl=nbas*(nbas+1)/2

!_ACCTGT_($acc kernels present(lab,ndx,g,sm0,aat0,aaa0,tt0,tat0,prefac0))
!_OMPTGT_($omp target teams distribute parallel do)
        do ii=intndx,jntndx
          do jj=ii,kl
            intg=ndx(ii)+jj
            i=lab(1,ii)
            k=lab(2,ii)
            l=lab(1,jj)
            n=lab(2,jj)
!_ACCTGT_($acc loop private(ibl) vector(32))
            do ibl=1,nb0
              etotb=etotb+g(intg)*prefac0(ibl)*
     &             (sm0(ibl,k,i)*sm0(ibl,n,l)
     &             -aat0(ibl,n,i)*aaa0(ibl,l,k)
     &             -tt0(ibl,n,i)*tat0(ibl,l,k)
     &             -aaa0(ibl,n,i)*aat0(ibl,l,k)
     &             -tat0(ibl,n,i)*tt0(ibl,l,k)
     &             -aat0(ibl,l,i)*aaa0(ibl,n,k)
     &             -tt0(ibl,l,i)*tat0(ibl,n,k)
     &             -aaa0(ibl,l,i)*aat0(ibl,n,k)
     &             -tat0(ibl,l,i)*tt0(ibl,n,k))
            enddo
          enddo
        enddo
!_ACCTGT_($acc end kernels)

!     Reset index in buffer to zero
        nb0=0
        
!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)

        call timer_stop(31)

      endif

      if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.lbatch) then
        
        call timer_start(34)
        
c!_ACCTGT_($acc update device(tatl,aaal,aatl,ttl,sml))
c!_OMPTGT_($omp target update to(tatl,aaal,aatl,ttl,sml)) 

!_ACCTGT_($acc data copyin(prefac1))
!_OMPTGT_($omp target data map(to:prefac1))

        kl=nbas*(nbas+1)/2

!_ACCTGT_($acc kernels present(lab,ndx,g,sm1,aat1,aaa1,tt1,tat1,prefac1))
!_ACCTGT_($acc& present(diag1,bdiag1,bsdiag1,csdiag1))
!_OMPTGT_($omp target teams distribute parallel do)
        do ii=intndx,jntndx
          do jj=ii,kl
            intg=ndx(ii)+jj
            i=lab(1,ii)
            k=lab(2,ii)
            l=lab(1,jj)
            n=lab(2,jj)
!_ACCTGT_($acc loop private(ibl) vector(32))
            do ibl=1,nb1
              etotb=etotb+g(intg)*prefac1(ibl)*
     &             (sm1(ibl,k,i)*(diag1(ibl,l)*csdiag1(ibl,n)
     &             +diag1(ibl,n)*csdiag1(ibl,l)
     &             +bdiag1(ibl,l)*bsdiag1(ibl,n)
     &             +bdiag1(ibl,n)*bsdiag1(ibl,l))
     &             +(diag1(ibl,i)*csdiag1(ibl,k)
     &             +diag1(ibl,k)*csdiag1(ibl,i)
     &             +bdiag1(ibl,i)*bsdiag1(ibl,k)
     &             +bdiag1(ibl,k)*bsdiag1(ibl,i))*sm1(ibl,n,l)
     &             -csdiag1(ibl,i)*(aat1(ibl,n,k)*diag1(ibl,l)
     &             +aat1(ibl,l,k)*diag1(ibl,n))
     &             -bsdiag1(ibl,i)*(tt1(ibl,n,k)*bdiag1(ibl,l)
     &             +tt1(ibl,l,k)*bdiag1(ibl,n))
     &             -csdiag1(ibl,l)*(aaa1(ibl,n,i)*diag1(ibl,k)
     &             +aaa1(ibl,n,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,l)*(tat1(ibl,n,i)*bdiag1(ibl,k)
     &             +tat1(ibl,n,k)*bdiag1(ibl,i))
     &             -csdiag1(ibl,k)*(aat1(ibl,n,i)*diag1(ibl,l)
     &             +aat1(ibl,l,i)*diag1(ibl,n))
     &             -bsdiag1(ibl,k)*(tt1(ibl,n,i)*bdiag1(ibl,l)
     &             +tt1(ibl,l,i)*bdiag1(ibl,n))
     &             -csdiag1(ibl,n)*(aaa1(ibl,l,i)*diag1(ibl,k)
     &             +aaa1(ibl,l,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,n)*(tat1(ibl,l,i)*bdiag1(ibl,k)
     &             +tat1(ibl,l,k)*bdiag1(ibl,i)))
            enddo
          enddo
        enddo
!_ACCTGT_($acc end kernels)
    
!     Reset index into buffer
        nb1=0
          
        call timer_stop(34)
          
!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)

      endif

      call timer_stop(38)
      
      return
      end

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

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

      if(ising.ge.3) return
      
      if(idbg.ge.13) write(lfndbg,600)
 600  format(' Calculating two electron matrix element',/)

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

!$omp parallel do reduction(+:ts)
!$omp& shared(sm,aat,aaa,tt,ta,g,lab,numint)
!$omp& private(intg,i,k,l,n)
       do intg=1,numint
        i=lab(1,intg)
        k=lab(2,intg)
        l=lab(3,intg)
        n=lab(4,intg)
        ts=ts+g(intg)*(sm(k,i)*sm(n,l)
     &       -aat(n,i)*aaa(l,k)-tt(n,i)*ta(l,k)
     &       -aaa(n,i)*aat(l,k)-ta(n,i)*tt(l,k)
     &       -aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k)
     &       -aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
       enddo
!$omp end parallel do

       call timer_stop(32)

      else

       call timer_start(35)

!$omp parallel do reduction(+:e2)
!$omp& shared(diag,bdiag,csdiag,bsdiag,sm,aat,tt,aaa,ta,g,lab,numint)
!$omp& private(intg,i,k,l,n,aai,abi,bai,bbi,aak,abk,bak,bbk)
!$omp& private(aaj,abj,baj,bbj,sum2,aal,abl,bal,bbl)
       do intg=1,numint
        i=lab(1,intg) 
        k=lab(2,intg) 
        l=lab(3,intg) 
        n=lab(4,intg)
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
        e2=e2+g(intg)*
     &    (sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l)
     &          -bai*(aat(n,k)*aaj+aat(l,k)*aal)
     &          -bbi*(tt(n,k)*abj+tt(l,k)*abl)
     &          -baj*(aaa(n,i)*aak+aaa(n,k)*aai)
     &          -bbj*(ta(n,i)*abk+ta(n,k)*abi)
     &          -bak*(aat(n,i)*aaj+aat(l,i)*aal)
     &          -bbk*(tt(n,i)*abj+tt(l,i)*abl)
     &          -bal*(aaa(l,i)*aak+aaa(l,k)*aai)
     &          -bbl*(ta(l,i)*abk+ta(l,k)*abi))
       enddo
!$omp end parallel do

       call timer_stop(35)
      endif

      ts=ts*fourdet
      e2=e2+ts

      call timer_stop(38)

      return
      end

      subroutine gronor_gntwo_omp_mlab(lfndbg)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,jj,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

      if(ising.ge.3) return
      
      if(idbg.ge.13) write(lfndbg,600)
 600  format(' Calculating two electron matrix element',/)

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
!$omp parallel do reduction(+:ts)
!$omp& shared(sm,aat,aaa,tt,ta,g,lab,numint)
!$omp& private(intg,i,k,l,n,ii,jj)
       do ii=intndx,jntndx
         do jj=ii,kl
           intg=ndx(ii)+jj
           i=lab(1,ii)
           k=lab(2,ii)
           l=lab(1,jj)
           n=lab(2,jj)
           ts=ts+g(intg)*(sm(k,i)*sm(n,l)
     &          -aat(n,i)*aaa(l,k)-tt(n,i)*ta(l,k)
     &          -aaa(n,i)*aat(l,k)-ta(n,i)*tt(l,k)
     &          -aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k)
     &          -aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
         enddo
       enddo
!$omp end parallel do

       call timer_stop(32)

      else

       call timer_start(35)

       kl=nbas*(nbas+1)/2
       intg=0
!$omp parallel do reduction(+:e2)
!$omp& shared(diag,bdiag,csdiag,bsdiag,sm,aat,tt,aaa,ta,g,lab,numint)
!$omp& private(intg,i,k,l,n,aai,abi,bai,bbi,aak,abk,bak,bbk)
!$omp& private(aaj,abj,baj,bbj,sum2,aal,abl,bal,bbl,ii,jj)
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
           e2=e2+g(intg)*
     &          (sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l)
     &          -bai*(aat(n,k)*aaj+aat(l,k)*aal)
     &          -bbi*(tt(n,k)*abj+tt(l,k)*abl)
     &          -baj*(aaa(n,i)*aak+aaa(n,k)*aai)
     &          -bbj*(ta(n,i)*abk+ta(n,k)*abi)
     &          -bak*(aat(n,i)*aaj+aat(l,i)*aal)
     &          -bbk*(tt(n,i)*abj+tt(l,i)*abl)
     &          -bal*(aaa(l,i)*aak+aaa(l,k)*aai)
     &          -bbl*(ta(l,i)*abk+ta(l,k)*abi))
         enddo
       enddo
!$omp end parallel do

       call timer_stop(35)
      endif

      ts=ts*fourdet
      e2=e2+ts

      call timer_stop(38)

      return
      end

      subroutine gronor_gntwo_omp_batch(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg,ihc,nhc,ibl

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

      etotb=0.0d0
      
!     Copying arrays into batch arrays

!     Initialize nb0 and nb1 when first matrix element in a batch

!     nb0 : number of singularity 0 elements
!     nb1 : number of singularity 1 or 2 elements
      
      if(ihc.eq.1) then
        nb0=0
        nb1=0
      endif
      
      if(ising.lt.3) then

      if(ising.eq.0) then
        nb0=nb0+1
        ibl=nb0        
        prefac(ibl)=4.0d0*deta*fctr
      else
        ibl=nhc-nb1
        nb1=nb1+1
        prefac(ibl)=fctr
      endif

      if(iamhead.eq.1) sstot=sstot+fctr*deta

      if(idbg.ge.13) write(lfndbg,600)
 600  format(' Calculating two electron matrix element',/)

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

!$omp parallel do
!$omp& shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas)
!$omp& private(i,k)
      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo
!$omp end parallel do
          
      if(ising.gt.0) then

!$omp parallel do
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl,nbas)
!$omp& private(i)
        do i=1,nbas
          diagl(ibl,i)=diag(i)
          bdiagl(ibl,i)=bdiag(i)
          bsdiagl(ibl,i)=bsdiag(i)
          csdiagl(ibl,i)=csdiag(i)
        enddo
!$omp end parallel do

      endif
      
      call timer_stop(30)

      endif
      
      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if(ihc.eq.nhc) then

        call timer_start(38)
        
        if(nb0.gt.0.and.nb1.le.0) then
          call timer_start(32)

!$omp parallel do reduction(+:etotb)
!$omp& shared(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,prefac)
!$omp& private(intg,i,k,l,n,ibl)
          do intg=1,numint
            i=lab(1,intg)
            k=lab(2,intg)
            l=lab(3,intg)
            n=lab(4,intg)
            do ibl=1,nb0
              etotb=etotb+g(intg)*prefac(ibl)*(sml(ibl,k,i)*sml(ibl,n,l)
     &             -aatl(ibl,n,i)*aaal(ibl,l,k)
     &             -ttl(ibl,n,i)*tatl(ibl,l,k)
     &             -aaal(ibl,n,i)*aatl(ibl,l,k)
     &             -tatl(ibl,n,i)*ttl(ibl,l,k)
     &             -aatl(ibl,l,i)*aaal(ibl,n,k)
     &             -ttl(ibl,l,i)*tatl(ibl,n,k)
     &             -aaal(ibl,l,i)*aatl(ibl,n,k)
     &             -tatl(ibl,l,i)*ttl(ibl,n,k))
            enddo
          enddo
!$omp end parallel do

          call timer_stop(32)
        endif

        if(nb0.gt.0.and.nb1.gt.0) then
          call timer_start(36)
!$omp parallel do reduction(+:etotb)
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl)
!$omp& shared(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,nhc,nb1,prefac)
!$omp& private(intg,i,k,l,n,ibl)
          do intg=1,numint
            i=lab(1,intg)
            k=lab(2,intg)
            l=lab(3,intg)
            n=lab(4,intg)
            do ibl=1,nb0
              etotb=etotb+g(intg)*prefac(ibl)*(sml(ibl,k,i)*sml(ibl,n,l)
     &             -aatl(ibl,n,i)*aaal(ibl,l,k)
     &             -ttl(ibl,n,i)*tatl(ibl,l,k)
     &             -aaal(ibl,n,i)*aatl(ibl,l,k)
     &             -tatl(ibl,n,i)*ttl(ibl,l,k)
     &             -aatl(ibl,l,i)*aaal(ibl,n,k)
     &             -ttl(ibl,l,i)*tatl(ibl,n,k)
     &             -aaal(ibl,l,i)*aatl(ibl,n,k)
     &             -tatl(ibl,l,i)*ttl(ibl,n,k))
            enddo
            do ibl=nhc-nb1+1,nhc
              etotb=etotb+g(intg)*prefac(ibl)*
     &             (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
     &             +diagl(ibl,n)*csdiagl(ibl,l)
     &             +bdiagl(ibl,l)*bsdiagl(ibl,n)
     &             +bdiagl(ibl,n)*bsdiagl(ibl,l))
     &             +(diagl(ibl,i)*csdiagl(ibl,k)
     &             +diagl(ibl,k)*csdiagl(ibl,i)
     &             +bdiagl(ibl,i)*bsdiagl(ibl,k)
     &             +bdiagl(ibl,k)*bsdiagl(ibl,i))*sml(ibl,n,l)
     &             -csdiagl(ibl,i)*(aatl(ibl,n,k)*diagl(ibl,l)
     &             +aatl(ibl,l,k)*diagl(ibl,n))
     &             -bsdiagl(ibl,i)*(ttl(ibl,n,k)*bdiagl(ibl,l)
     &             +ttl(ibl,l,k)*bdiagl(ibl,n))
     &             -csdiagl(ibl,l)*(aaal(ibl,n,i)*diagl(ibl,k)
     &             +aaal(ibl,n,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,l)*(tatl(ibl,n,i)*bdiagl(ibl,k)
     &             +tatl(ibl,n,k)*bdiagl(ibl,i))
     &             -csdiagl(ibl,k)*(aatl(ibl,n,i)*diagl(ibl,l)
     &             +aatl(ibl,l,i)*diagl(ibl,n))
     &             -bsdiagl(ibl,k)*(ttl(ibl,n,i)*bdiagl(ibl,l)
     &             +ttl(ibl,l,i)*bdiagl(ibl,n))
     &             -csdiagl(ibl,n)*(aaal(ibl,l,i)*diagl(ibl,k)
     &             +aaal(ibl,l,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,n)*(tatl(ibl,l,i)*bdiagl(ibl,k)
     &             +tatl(ibl,l,k)*bdiagl(ibl,i)))
            enddo
          enddo
!$omp end parallel do
          call timer_stop(36)
        endif
        
        if(nb0.le.0.and.nb1.gt.0) then
          call timer_start(35)
!$omp parallel do reduction(+:etotb)
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl)
!$omp& shared(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,nhc,nb1,prefac)
!$omp& private(intg,i,k,l,n,ibl)
          do intg=1,numint
            i=lab(1,intg) 
            k=lab(2,intg) 
            l=lab(3,intg) 
            n=lab(4,intg)
            do ibl=nhc-nb1+1,nhc
              etotb=etotb+g(intg)*prefac(ibl)*
     &             (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
     &             +diagl(ibl,n)*csdiagl(ibl,l)
     &             +bdiagl(ibl,l)*bsdiagl(ibl,n)
     &             +bdiagl(ibl,n)*bsdiagl(ibl,l))
     &             +(diagl(ibl,i)*csdiagl(ibl,k)
     &             +diagl(ibl,k)*csdiagl(ibl,i)
     &             +bdiagl(ibl,i)*bsdiagl(ibl,k)
     &             +bdiagl(ibl,k)*bsdiagl(ibl,i))*sml(ibl,n,l)
     &             -csdiagl(ibl,i)*(aatl(ibl,n,k)*diagl(ibl,l)
     &             +aatl(ibl,l,k)*diagl(ibl,n))
     &             -bsdiagl(ibl,i)*(ttl(ibl,n,k)*bdiagl(ibl,l)
     &             +ttl(ibl,l,k)*bdiagl(ibl,n))
     &             -csdiagl(ibl,l)*(aaal(ibl,n,i)*diagl(ibl,k)
     &             +aaal(ibl,n,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,l)*(tatl(ibl,n,i)*bdiagl(ibl,k)
     &             +tatl(ibl,n,k)*bdiagl(ibl,i))
     &             -csdiagl(ibl,k)*(aatl(ibl,n,i)*diagl(ibl,l)
     &             +aatl(ibl,l,i)*diagl(ibl,n))
     &             -bsdiagl(ibl,k)*(ttl(ibl,n,i)*bdiagl(ibl,l)
     &             +ttl(ibl,l,i)*bdiagl(ibl,n))
     &             -csdiagl(ibl,n)*(aaal(ibl,l,i)*diagl(ibl,k)
     &             +aaal(ibl,l,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,n)*(tatl(ibl,l,i)*bdiagl(ibl,k)
     &             +tatl(ibl,l,k)*bdiagl(ibl,i)))
            enddo
          enddo
!$omp end parallel do
          call timer_stop(35)
        endif
        
      endif

      return
      end

      subroutine gronor_gntwo_omp_batch_mlab(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg,ihc,nhc,ibl,j,jj

      real (kind=8) :: e2n,tsn,gtemp,et, etemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

      character (len=1) :: ch

      etotb=0.0d0
      
!     Copying arrays into batch arrays

!     Initialize nb0 and nb1 when first matrix element in a batch

!     nb0 : number of singularity 0 elements
!     nb1 : number of singularity 1 or 2 elements
      
      if(ihc.eq.1) then
        nb0=0
        nb1=0
      endif
      
      if(ising.lt.3) then

      if(ising.eq.0) then
        nb0=nb0+1
        ibl=nb0        
        prefac(ibl)=4.0d0*deta*fctr
      else
        ibl=nhc-nb1
        nb1=nb1+1
        prefac(ibl)=fctr
      endif

      if(iamhead.eq.1) sstot=sstot+fctr*deta

      if(idbg.ge.13) write(lfndbg,600)
 600  format(' Calculating two electron matrix element',/)

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

!$omp parallel do
!$omp& shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas)
!$omp& private(i,k)
      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo
!$omp end parallel do
          
      if(ising.gt.0) then

!$omp parallel do
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl,nbas)
!$omp& private(i)
        do i=1,nbas
          diagl(ibl,i)=diag(i)
          bdiagl(ibl,i)=bdiag(i)
          bsdiagl(ibl,i)=bsdiag(i)
          csdiagl(ibl,i)=csdiag(i)
        enddo
!$omp end parallel do

      endif
      
      call timer_stop(30)

      endif
      
      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if(ihc.eq.nhc) then

        call timer_start(38)
        
        if(nb0.gt.0.and.nb1.le.0) then
          call timer_start(32)

          kl=nbas*(nbas+1)/2
          do ii=intndx,jntndx
!$omp parallel
!$omp& shared(ii,sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,prefac)
!$omp do reduction(+:etotb)
!$omp& private(intg,jj,i,k,l,n,ibl)
            do jj=ii,kl
              intg=ndx(ii)+jj
              i=lab(1,ii)
              k=lab(2,ii)
              l=lab(1,jj)
              n=lab(2,jj)
              etemp=0.0d0
              do ibl=1,nb0
                etemp=etemp+prefac(ibl)*
     &               (sml(ibl,k,i)*sml(ibl,n,l)
     &               -aatl(ibl,n,i)*aaal(ibl,l,k)
     &               -ttl(ibl,n,i)*tatl(ibl,l,k)
     &               -aaal(ibl,n,i)*aatl(ibl,l,k)
     &               -tatl(ibl,n,i)*ttl(ibl,l,k)
     &               -aatl(ibl,l,i)*aaal(ibl,n,k)
     &               -ttl(ibl,l,i)*tatl(ibl,n,k)
     &               -aaal(ibl,l,i)*aatl(ibl,n,k)
     &               -tatl(ibl,l,i)*ttl(ibl,n,k))
              enddo
              etotb=etotb+g(intg)*etemp
            enddo
!$omp end do
!$omp end parallel
          enddo

          call timer_stop(32)
        endif

        if(nb0.gt.0.and.nb1.gt.0) then
          call timer_start(36)
          kl=nbas*(nbas+1)/2
!$omp parallel do reduction(+:etotb)
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl)
!$omp& shared(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,nhc,nb1,prefac)
!$omp& private(intg,i,k,l,n,ibl)
          do ii=intndx,jntndx
            do jj=ii,kl
              intg=ndx(ii)+jj
              i=lab(1,ii)
              k=lab(2,ii)
              l=lab(1,jj)
              n=lab(2,jj)
              etemp=0.0d0
              do ibl=1,nb0
                etemp=etemp+prefac(ibl)*
     &               (sml(ibl,k,i)*sml(ibl,n,l)
     &               -aatl(ibl,n,i)*aaal(ibl,l,k)
     &               -ttl(ibl,n,i)*tatl(ibl,l,k)
     &               -aaal(ibl,n,i)*aatl(ibl,l,k)
     &               -tatl(ibl,n,i)*ttl(ibl,l,k)
     &               -aatl(ibl,l,i)*aaal(ibl,n,k)
     &               -ttl(ibl,l,i)*tatl(ibl,n,k)
     &               -aaal(ibl,l,i)*aatl(ibl,n,k)
     &               -tatl(ibl,l,i)*ttl(ibl,n,k))
              enddo
              do ibl=nhc-nb1+1,nhc
                etemp=etemp+prefac(ibl)*
     &               (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
     &               +diagl(ibl,n)*csdiagl(ibl,l)
     &               +bdiagl(ibl,l)*bsdiagl(ibl,n)
     &               +bdiagl(ibl,n)*bsdiagl(ibl,l))
     &               +(diagl(ibl,i)*csdiagl(ibl,k)
     &               +diagl(ibl,k)*csdiagl(ibl,i)
     &               +bdiagl(ibl,i)*bsdiagl(ibl,k)
     &               +bdiagl(ibl,k)*bsdiagl(ibl,i))*sml(ibl,n,l)
     &               -csdiagl(ibl,i)*(aatl(ibl,n,k)*diagl(ibl,l)
     &               +aatl(ibl,l,k)*diagl(ibl,n))
     &               -bsdiagl(ibl,i)*(ttl(ibl,n,k)*bdiagl(ibl,l)
     &               +ttl(ibl,l,k)*bdiagl(ibl,n))
     &               -csdiagl(ibl,l)*(aaal(ibl,n,i)*diagl(ibl,k)
     &               +aaal(ibl,n,k)*diagl(ibl,i))
     &               -bsdiagl(ibl,l)*(tatl(ibl,n,i)*bdiagl(ibl,k)
     &               +tatl(ibl,n,k)*bdiagl(ibl,i))
     &               -csdiagl(ibl,k)*(aatl(ibl,n,i)*diagl(ibl,l)
     &               +aatl(ibl,l,i)*diagl(ibl,n))
     &               -bsdiagl(ibl,k)*(ttl(ibl,n,i)*bdiagl(ibl,l)
     &               +ttl(ibl,l,i)*bdiagl(ibl,n))
     &               -csdiagl(ibl,n)*(aaal(ibl,l,i)*diagl(ibl,k)
     &               +aaal(ibl,l,k)*diagl(ibl,i))
     &               -bsdiagl(ibl,n)*(tatl(ibl,l,i)*bdiagl(ibl,k)
     &               +tatl(ibl,l,k)*bdiagl(ibl,i)))
              enddo
              etotb=etotb+g(intg)*etemp
            enddo
          enddo
!$omp end parallel do
          call timer_stop(36)
        endif
        
        if(nb0.le.0.and.nb1.gt.0) then
          call timer_start(35)
          kl=nbas*(nbas+1)/2
!$omp parallel do reduction(+:etotb)
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl)
!$omp& shared(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,nhc,nb1,prefac)
!$omp& private(intg,i,k,l,n,ibl)
          do ii=intndx,jntndx
            do jj=ii,kl
              intg=ndx(ii)+jj
              i=lab(1,ii)
              k=lab(2,ii)
              l=lab(1,jj)
              n=lab(2,jj)
              etemp=0.0d0
              do ibl=nhc-nb1+1,nhc
                etemp=etemp+prefac(ibl)*
     &               (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
     &               +diagl(ibl,n)*csdiagl(ibl,l)
     &               +bdiagl(ibl,l)*bsdiagl(ibl,n)
     &               +bdiagl(ibl,n)*bsdiagl(ibl,l))
     &               +(diagl(ibl,i)*csdiagl(ibl,k)
     &               +diagl(ibl,k)*csdiagl(ibl,i)
     &               +bdiagl(ibl,i)*bsdiagl(ibl,k)
     &               +bdiagl(ibl,k)*bsdiagl(ibl,i))*sml(ibl,n,l)
     &               -csdiagl(ibl,i)*(aatl(ibl,n,k)*diagl(ibl,l)
     &               +aatl(ibl,l,k)*diagl(ibl,n))
     &               -bsdiagl(ibl,i)*(ttl(ibl,n,k)*bdiagl(ibl,l)
     &               +ttl(ibl,l,k)*bdiagl(ibl,n))
     &               -csdiagl(ibl,l)*(aaal(ibl,n,i)*diagl(ibl,k)
     &               +aaal(ibl,n,k)*diagl(ibl,i))
     &               -bsdiagl(ibl,l)*(tatl(ibl,n,i)*bdiagl(ibl,k)
     &               +tatl(ibl,n,k)*bdiagl(ibl,i))
     &               -csdiagl(ibl,k)*(aatl(ibl,n,i)*diagl(ibl,l)
     &               +aatl(ibl,l,i)*diagl(ibl,n))
     &               -bsdiagl(ibl,k)*(ttl(ibl,n,i)*bdiagl(ibl,l)
     &               +ttl(ibl,l,i)*bdiagl(ibl,n))
     &               -csdiagl(ibl,n)*(aaal(ibl,l,i)*diagl(ibl,k)
     &               +aaal(ibl,l,k)*diagl(ibl,i))
     &               -bsdiagl(ibl,n)*(tatl(ibl,l,i)*bdiagl(ibl,k)
     &               +tatl(ibl,l,k)*bdiagl(ibl,i)))
              enddo
              etotb=etotb+g(intg)*etemp
            enddo
          enddo
!$omp end parallel do
          call timer_stop(35)
        endif
        
      endif

      return
      end

      subroutine gronor_gntwo_omp_batch2(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg,ihc,nhc,ibl

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

      etotb=0.0d0
      
!     Copying arrays into batch arrays

!     Initialize nb0 and nb1 when first matrix element in a batch

!     nb0 : number of singularity 0 elements
!     nb1 : number of singularity 1 or 2 elements
      
      if(ihc.eq.1) then
        nb0=0
        nb1=0
      endif
      
      if(ising.lt.3) then

      if(ising.eq.0) then
        nb0=nb0+1
        ibl=nb0        
        prefac(ibl)=4.0d0*deta*fctr
      else
        ibl=nhc-nb1
        nb1=nb1+1
        prefac(ibl)=fctr
      endif

      if(iamhead.eq.1) sstot=sstot+fctr*deta

      if(idbg.ge.13) write(lfndbg,600)
 600  format(' Calculating two electron matrix element',/)

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

!$omp parallel do
!$omp& shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas)
!$omp& private(i,k)
      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo
!$omp end parallel do
          
      if(ising.gt.0) then

!$omp parallel do
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl,nbas)
!$omp& private(i)
        do i=1,nbas
          diagl(ibl,i)=diag(i)
          bdiagl(ibl,i)=bdiag(i)
          bsdiagl(ibl,i)=bsdiag(i)
          csdiagl(ibl,i)=csdiag(i)
        enddo
!$omp end parallel do

      endif
      
      call timer_stop(30)

      endif
      
      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if(ihc.eq.nhc) then

        call timer_start(38)
        
        if(nb0.gt.0.and.nb1.le.0) then
          call timer_start(32)

!$omp parallel do reduction(+:etotb)
!$omp& shared(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,prefac)
!$omp& private(intg,i,k,l,n,ibl)
          do intg=1,numint
            i=lab(1,intg)
            k=lab(2,intg)
            l=lab(3,intg)
            n=lab(4,intg)
            do ibl=1,nb0
              etotb=etotb+g(intg)*prefac(ibl)*(sml(ibl,k,i)*sml(ibl,n,l)
     &             -aatl(ibl,n,i)*aaal(ibl,l,k)
     &             -ttl(ibl,n,i)*tatl(ibl,l,k)
     &             -aaal(ibl,n,i)*aatl(ibl,l,k)
     &             -tatl(ibl,n,i)*ttl(ibl,l,k)
     &             -aatl(ibl,l,i)*aaal(ibl,n,k)
     &             -ttl(ibl,l,i)*tatl(ibl,n,k)
     &             -aaal(ibl,l,i)*aatl(ibl,n,k)
     &             -tatl(ibl,l,i)*ttl(ibl,n,k))
            enddo
          enddo
!$omp end parallel do

          call timer_stop(32)
        endif

        if(nb0.gt.0.and.nb1.gt.0) then
          call timer_start(36)
!$omp parallel do reduction(+:etotb)
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl)
!$omp& shared(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,nhc,nb1,prefac)
!$omp& private(intg,i,k,l,n,ibl)
          do intg=1,numint
            i=lab(1,intg)
            k=lab(2,intg)
            l=lab(3,intg)
            n=lab(4,intg)
            do ibl=1,nb0
              etotb=etotb+g(intg)*prefac(ibl)*(sml(ibl,k,i)*sml(ibl,n,l)
     &             -aatl(ibl,n,i)*aaal(ibl,l,k)
     &             -ttl(ibl,n,i)*tatl(ibl,l,k)
     &             -aaal(ibl,n,i)*aatl(ibl,l,k)
     &             -tatl(ibl,n,i)*ttl(ibl,l,k)
     &             -aatl(ibl,l,i)*aaal(ibl,n,k)
     &             -ttl(ibl,l,i)*tatl(ibl,n,k)
     &             -aaal(ibl,l,i)*aatl(ibl,n,k)
     &             -tatl(ibl,l,i)*ttl(ibl,n,k))
            enddo
            do ibl=nhc-nb1+1,nhc
              etotb=etotb+g(intg)*prefac(ibl)*
     &             (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
     &             +diagl(ibl,n)*csdiagl(ibl,l)
     &             +bdiagl(ibl,l)*bsdiagl(ibl,n)
     &             +bdiagl(ibl,n)*bsdiagl(ibl,l))
     &             +(diagl(ibl,i)*csdiagl(ibl,k)
     &             +diagl(ibl,k)*csdiagl(ibl,i)
     &             +bdiagl(ibl,i)*bsdiagl(ibl,k)
     &             +bdiagl(ibl,k)*bsdiagl(ibl,i))*sml(ibl,n,l)
     &             -csdiagl(ibl,i)*(aatl(ibl,n,k)*diagl(ibl,l)
     &             +aatl(ibl,l,k)*diagl(ibl,n))
     &             -bsdiagl(ibl,i)*(ttl(ibl,n,k)*bdiagl(ibl,l)
     &             +ttl(ibl,l,k)*bdiagl(ibl,n))
     &             -csdiagl(ibl,l)*(aaal(ibl,n,i)*diagl(ibl,k)
     &             +aaal(ibl,n,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,l)*(tatl(ibl,n,i)*bdiagl(ibl,k)
     &             +tatl(ibl,n,k)*bdiagl(ibl,i))
     &             -csdiagl(ibl,k)*(aatl(ibl,n,i)*diagl(ibl,l)
     &             +aatl(ibl,l,i)*diagl(ibl,n))
     &             -bsdiagl(ibl,k)*(ttl(ibl,n,i)*bdiagl(ibl,l)
     &             +ttl(ibl,l,i)*bdiagl(ibl,n))
     &             -csdiagl(ibl,n)*(aaal(ibl,l,i)*diagl(ibl,k)
     &             +aaal(ibl,l,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,n)*(tatl(ibl,l,i)*bdiagl(ibl,k)
     &             +tatl(ibl,l,k)*bdiagl(ibl,i)))
            enddo
          enddo
!$omp end parallel do
          call timer_stop(36)
        endif
        
        if(nb0.le.0.and.nb1.gt.0) then
          call timer_start(35)
!$omp parallel do reduction(+:etotb)
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl)
!$omp& shared(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,nhc,nb1,prefac)
!$omp& private(intg,i,k,l,n,ibl)
          do intg=1,numint
            i=lab(1,intg) 
            k=lab(2,intg) 
            l=lab(3,intg) 
            n=lab(4,intg)
            do ibl=nhc-nb1+1,nhc
              etotb=etotb+g(intg)*prefac(ibl)*
     &             (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
     &             +diagl(ibl,n)*csdiagl(ibl,l)
     &             +bdiagl(ibl,l)*bsdiagl(ibl,n)
     &             +bdiagl(ibl,n)*bsdiagl(ibl,l))
     &             +(diagl(ibl,i)*csdiagl(ibl,k)
     &             +diagl(ibl,k)*csdiagl(ibl,i)
     &             +bdiagl(ibl,i)*bsdiagl(ibl,k)
     &             +bdiagl(ibl,k)*bsdiagl(ibl,i))*sml(ibl,n,l)
     &             -csdiagl(ibl,i)*(aatl(ibl,n,k)*diagl(ibl,l)
     &             +aatl(ibl,l,k)*diagl(ibl,n))
     &             -bsdiagl(ibl,i)*(ttl(ibl,n,k)*bdiagl(ibl,l)
     &             +ttl(ibl,l,k)*bdiagl(ibl,n))
     &             -csdiagl(ibl,l)*(aaal(ibl,n,i)*diagl(ibl,k)
     &             +aaal(ibl,n,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,l)*(tatl(ibl,n,i)*bdiagl(ibl,k)
     &             +tatl(ibl,n,k)*bdiagl(ibl,i))
     &             -csdiagl(ibl,k)*(aatl(ibl,n,i)*diagl(ibl,l)
     &             +aatl(ibl,l,i)*diagl(ibl,n))
     &             -bsdiagl(ibl,k)*(ttl(ibl,n,i)*bdiagl(ibl,l)
     &             +ttl(ibl,l,i)*bdiagl(ibl,n))
     &             -csdiagl(ibl,n)*(aaal(ibl,l,i)*diagl(ibl,k)
     &             +aaal(ibl,l,k)*diagl(ibl,i))
     &             -bsdiagl(ibl,n)*(tatl(ibl,l,i)*bdiagl(ibl,k)
     &             +tatl(ibl,l,k)*bdiagl(ibl,i)))
            enddo
          enddo
!$omp end parallel do
          call timer_stop(35)
        endif
        
      endif

      return
      end

      subroutine gronor_gntwo_omp_batch2_mlab(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,intg,ihc,nhc,ibl,j,jj

      real (kind=8) :: e2n,tsn,gtemp,et, etemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

      character (len=1) :: ch

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

        if(idbg.ge.13) write(lfndbg,600)
 600    format(' Calculating two electron matrix element',/)

        numint=mint2

        if(ising.gt.1) e1=0.0d0
        e2=0.0d0
        ts=0.0d0
        fourdet=4.0d0*deta

        if(ising.eq.0) then
          nb0=nb0+1      
          prefac0(nb0)=4.0d0*deta*fctr
!$omp parallel do
!$omp& shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas)
!$omp& private(i,k)
          do k=1,nbas
            do i=1,nbas
              sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa0(nb0,i,k)=aaa(i,k)
              aat0(nb0,i,k)=aaa(k,i)
              tt0(nb0,i,k)=ta(k,i)
              tat0(nb0,i,k)=ta(i,k)
            enddo
          enddo
!$omp end parallel do
        else
          nb1=nb1+1
          prefac1(nb1)=fctr
!$omp parallel do
!$omp& shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas)
!$omp& private(i,k)
          do k=1,nbas
            do i=1,nbas
              sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa1(nb1,i,k)=aaa(i,k)
              aat1(nb1,i,k)=aaa(k,i)
              tt1(nb1,i,k)=ta(k,i)
              tat1(nb1,i,k)=ta(i,k)
            enddo
          enddo
!$omp end parallel do

!$omp parallel do
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl,nbas)
!$omp& private(i)
          do i=1,nbas
            diag1(nb1,i)=diag(i)
            bdiag1(nb1,i)=bdiag(i)
            bsdiag1(nb1,i)=bsdiag(i)
            csdiag1(nb1,i)=csdiag(i)
          enddo
!$omp end parallel do
        endif
        
        call timer_stop(30)

      endif
      
      call timer_start(38)

!     After the last element of a task the list or when the buffer is full
!     the energies are evaluated

      if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.lbatch) then

        call timer_start(32)

        kl=nbas*(nbas+1)/2
        do ii=intndx,jntndx
!$omp parallel
!$omp& shared(ii,sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,prefac)
!$omp do reduction(+:etotb)
!$omp& private(intg,jj,i,k,l,n,ibl)
          do jj=ii,kl
            intg=ndx(ii)+jj
            i=lab(1,ii)
            k=lab(2,ii)
            l=lab(1,jj)
            n=lab(2,jj)
            etemp=0.0d0
            do ibl=1,nb0
              etemp=etemp+prefac0(ibl)*
     &             (sm0(ibl,k,i)*sm0(ibl,n,l)
     &             -aat0(ibl,n,i)*aaa0(ibl,l,k)
     &             -tt0(ibl,n,i)*tat0(ibl,l,k)
     &             -aaa0(ibl,n,i)*aat0(ibl,l,k)
     &             -tat0(ibl,n,i)*tt0(ibl,l,k)
     &             -aat0(ibl,l,i)*aaa0(ibl,n,k)
     &             -tt0(ibl,l,i)*tat0(ibl,n,k)
     &             -aaa0(ibl,l,i)*aat0(ibl,n,k)
     &             -tat0(ibl,l,i)*tt0(ibl,n,k))
            enddo
            etotb=etotb+g(intg)*etemp
          enddo
!$omp end do
!$omp end parallel
        enddo

!     Reset index in buffer to zero
        nb0=0

        call timer_stop(32)
        
      endif

      if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.lbatch) then
        
        call timer_start(35)
        
        kl=nbas*(nbas+1)/2
        do ii=intndx,jntndx
!$omp parallel do reduction(+:etotb)
!$omp& shared(diagl,bdiagl,csdiagl,bsdiagl)
!$omp& shared(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,nhc,nb1,prefac)
!$omp& private(intg,i,k,l,n,ibl)
          do jj=ii,kl
            intg=ndx(ii)+jj
            i=lab(1,ii)
            k=lab(2,ii)
            l=lab(1,jj)
            n=lab(2,jj)
            etemp=0.0d0
            do ibl=1,nb1
              etemp=etemp+prefac1(ibl)*
     &             (sm1(ibl,k,i)*(diag1(ibl,l)*csdiag1(ibl,n)
     &             +diag1(ibl,n)*csdiag1(ibl,l)
     &             +bdiag1(ibl,l)*bsdiag1(ibl,n)
     &             +bdiag1(ibl,n)*bsdiag1(ibl,l))
     &             +(diag1(ibl,i)*csdiag1(ibl,k)
     &             +diag1(ibl,k)*csdiag1(ibl,i)
     &             +bdiag1(ibl,i)*bsdiag1(ibl,k)
     &             +bdiag1(ibl,k)*bsdiag1(ibl,i))*sm1(ibl,n,l)
     &             -csdiag1(ibl,i)*(aat1(ibl,n,k)*diag1(ibl,l)
     &             +aat1(ibl,l,k)*diag1(ibl,n))
     &             -bsdiag1(ibl,i)*(tt1(ibl,n,k)*bdiag1(ibl,l)
     &             +tt1(ibl,l,k)*bdiag1(ibl,n))
     &             -csdiag1(ibl,l)*(aaa1(ibl,n,i)*diag1(ibl,k)
     &             +aaa1(ibl,n,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,l)*(tat1(ibl,n,i)*bdiag1(ibl,k)
     &             +tat1(ibl,n,k)*bdiag1(ibl,i))
     &             -csdiag1(ibl,k)*(aat1(ibl,n,i)*diag1(ibl,l)
     &             +aat1(ibl,l,i)*diag1(ibl,n))
     &             -bsdiag1(ibl,k)*(tt1(ibl,n,i)*bdiag1(ibl,l)
     &             +tt1(ibl,l,i)*bdiag1(ibl,n))
     &             -csdiag1(ibl,n)*(aaa1(ibl,l,i)*diag1(ibl,k)
     &             +aaa1(ibl,l,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,n)*(tat1(ibl,l,i)*bdiag1(ibl,k)
     &             +tat1(ibl,l,k)*bdiag1(ibl,i)))
            enddo
            etotb=etotb+g(intg)*etemp
          enddo
!$omp end parallel do
        enddo

!     Reset index into buffer
        nb1=0

        call timer_stop(35)
          
      endif

      return
      end

      
