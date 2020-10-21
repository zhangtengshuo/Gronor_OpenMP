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

!$acc data create(sm)

      call timer_start(30)
      
!$acc kernels present(aat,aaa,tt,ta,tat,sm)

      do i=1,nbas 
        do k=1,nbas
          tat(k,i)=ta(k,i)
          tt(i,k)=tat(k,i)
          aat(i,k)=aaa(k,i)
        enddo
      enddo
      
      do i=1,nbas 
        do k=1,nbas
          sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+tat(k,i)
        enddo
      enddo
!$acc end kernels

      call timer_stop(30)

      call timer_start(38)

      if(ising.le.0) then
        
        call timer_start(31)
        
        tst=ts
!$acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab)
!$acc loop reduction (+:tst)
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
!$acc end kernels
        ts=tst

        call timer_stop(31)

      else
        
        call timer_start(34)

        e2t=e2
!$acc kernels present(aat,aaa,tt,ta,tat,sm,g)
!$acc& present(lab,diag,bdiag,bsdiag,csdiag)
!$acc loop reduction (+:e2t)
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
!$acc end kernels
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

!$acc end data

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

!$acc data create(sm)

      call timer_start(30)
      
!$acc kernels present(aat,aaa,tt,ta,tat,sm)

      do i=1,nbas 
        do k=1,nbas
          tat(k,i)=ta(k,i)
          tt(i,k)=tat(k,i)
          aat(i,k)=aaa(k,i)
        enddo
      enddo
      
      do i=1,nbas 
        do k=1,nbas
          sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+tat(k,i)
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
!$acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab)
!$acc loop reduction (+:tst)
        do ii=intndx,jntndx
          do jj=ii,kl
            intg=intg+1
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
!$acc end kernels
        ts=tst

        call timer_stop(31)

      else
        
        call timer_start(34)

        e2t=e2
        kl=nbas*(nbas+1)/2
        intg=0
!$acc kernels present(aat,aaa,tt,ta,tat,sm,g)
!$acc& present(lab,diag,bdiag,bsdiag,csdiag)
!$acc loop reduction (+:e2t)
        do ii=intndx,jntndx
          do jj=ii,kl
            intg=intg+1
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
!$acc end kernels
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

!$acc end data

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

      if(ising.eq.0) then
        nb0=nb0+1
        ibl=nb0        
        prefac(ibl)=4.0d0*deta*fctr
      else
        ibl=nhc-nb1
        nb1=nb1+1
        prefac(ibl)=fctr
      endif
      factor(ibl)=fctr

      if(iamhead.eq.1) sstot=sstot+factor(ibl)*deta

      if(idbg.ge.13) write(lfndbg,600)
 600  format(' Calculating two electron matrix element',/)

      numint=mint2

      if(ising.gt.1) e1=0.0d0
      e2=0.0d0
      ts=0.0d0
      fourdet=4.0d0*deta

      call timer_start(30)

      if(ising.gt.0) then

!$acc kernels present(diag,bdiag,bsdiag,csdiag)
!$acc& present(diagl,bdiagl,bsdiagl,csdiagl) copyin(ibl)
        do i=1,nbas
          diagl(ibl,i)=diag(i)
          bdiagl(ibl,i)=bdiag(i)
          bsdiagl(ibl,i)=bsdiag(i)
          csdiagl(ibl,i)=csdiag(i)
        enddo
!$acc end kernels

      endif
!$acc kernels present(aat,aaa,tt,ta)
      do i=1,nbas 
        do k=1,nbas 
          tt(i,k)=ta(k,i)
          aat(i,k)=aaa(k,i)
        enddo
      enddo
!$acc end kernels

!     ASSEMBLE ARRAYS ON GPU
      
!$acc kernels present(aaal,aatl,ttl,tatl,tt,aat,aaa,ta) copyin(ibl)
      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo
!$acc end kernels

!     END ASSEMBLY ON GPU
      
!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
!$acc update self(ta,aaa,aat,tt)

      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo

!     END ASSEMBLY ON CPU

      call timer_stop(30)

      endif

      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if(ihc.eq.nhc) then

c!$acc update device(tatl,aaal,aatl,ttl,sml,diagl,bdiagl,bsdiagl,csdiagl)
!$acc update device(tatl,aaal,aatl,ttl,sml)

        etotb=0.0d0
        
!$acc data copyin(prefac)
        
c        if(nb0.gt.0.and.nb1.le.0) then
        if(nb0.gt.0) then
          call timer_start(31)

!$acc kernels present(lab,g,sml,aatl,aaal,ttl,tatl,prefac)
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
!$acc end kernels

          call timer_stop(31)
        endif

c        if(nb0.gt.0.and.nb1.gt.0) then
cc        if(nb0.gt.0.and.nb1.gt.0.and.1.eq.2) then
c          call timer_start(31)
c
c!$acc kernels present(lab,g,sml,aatl,aaal,ttl,tatl,prefac)
c!$acc& present(diagl,bdiagl,bsdiagl,csdiagl)
c          do intg=1,numint
c            i=lab(1,intg)
c            k=lab(2,intg)
c            l=lab(3,intg)
c            n=lab(4,intg)
c            do ibl=1,nb0
c              etotb=etotb+g(intg)*prefac(ibl)*(sml(ibl,k,i)*sml(ibl,n,l)
c     &             -aatl(ibl,n,i)*aaal(ibl,l,k)
c     &             -ttl(ibl,n,i)*tatl(ibl,l,k)
c     &             -aaal(ibl,n,i)*aatl(ibl,l,k)
c     &             -tatl(ibl,n,i)*ttl(ibl,l,k)
c     &             -aatl(ibl,l,i)*aaal(ibl,n,k)
c     &             -ttl(ibl,l,i)*tatl(ibl,n,k)
c     &             -aaal(ibl,l,i)*aatl(ibl,n,k)
c     &             -tatl(ibl,l,i)*ttl(ibl,n,k))
c            enddo
c            do ibl=nhc-nb1+1,nhc
c              etotb=etotb+g(intg)*prefac(ibl)*
c     &             (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
c     &             +diagl(ibl,n)*csdiagl(ibl,l)
c     &             +bdiagl(ibl,l)*bsdiagl(ibl,n)
c     &             +bdiagl(ibl,n)*bsdiagl(ibl,l))
c     &             +(diagl(ibl,i)*csdiagl(ibl,k)
c     &             +diagl(ibl,k)*csdiagl(ibl,i)
c     &             +bdiagl(ibl,i)*bsdiagl(ibl,k)
c     &             +bdiagl(ibl,k)*bsdiagl(ibl,i))*sml(ibl,n,l)
c     &             -csdiagl(ibl,i)*(aatl(ibl,n,k)*diagl(ibl,l)
c     &             +aatl(ibl,l,k)*diagl(ibl,n))
c     &             -bsdiagl(ibl,i)*(ttl(ibl,n,k)*bdiagl(ibl,l)
c     &             +ttl(ibl,l,k)*bdiagl(ibl,n))
c     &             -csdiagl(ibl,l)*(aaal(ibl,n,i)*diagl(ibl,k)
c     &             +aaal(ibl,n,k)*diagl(ibl,i))
c     &             -bsdiagl(ibl,l)*(tatl(ibl,n,i)*bdiagl(ibl,k)
c     &             +tatl(ibl,n,k)*bdiagl(ibl,i))
c     &             -csdiagl(ibl,k)*(aatl(ibl,n,i)*diagl(ibl,l)
c     &             +aatl(ibl,l,i)*diagl(ibl,n))
c     &             -bsdiagl(ibl,k)*(ttl(ibl,n,i)*bdiagl(ibl,l)
c     &             +ttl(ibl,l,i)*bdiagl(ibl,n))
c     &             -csdiagl(ibl,n)*(aaal(ibl,l,i)*diagl(ibl,k)
c     &             +aaal(ibl,l,k)*diagl(ibl,i))
c     &             -bsdiagl(ibl,n)*(tatl(ibl,l,i)*bdiagl(ibl,k)
c     &             +tatl(ibl,l,k)*bdiagl(ibl,i)))
c            enddo
c          enddo
c!$acc end kernels
c
c          call timer_stop(31)
c        endif
        
c        if(nb0.le.0.and.nb1.gt.0) then
        if(nb1.gt.0) then
          call timer_start(34)

!$acc kernels present(lab,g,sml,aatl,aaal,ttl,tatl,prefac)
!$acc& present(diagl,bdiagl,bsdiagl,csdiagl)
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
!$acc end kernels
    
          call timer_stop(34)
        endif
        
!$acc end data

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
      factor(ibl)=fctr

      if(iamhead.eq.1) sstot=sstot+factor(ibl)*deta

      if(idbg.ge.13) write(lfndbg,600)
 600  format(' Calculating two electron matrix element',/)

      numint=mint2

      if(ising.gt.1) e1=0.0d0
      e2=0.0d0
      ts=0.0d0
      fourdet=4.0d0*deta

      call timer_start(30)

      if(ising.gt.0) then

!$acc kernels present(diag,bdiag,bsdiag,csdiag)
!$acc& present(diagl,bdiagl,bsdiagl,csdiagl) copyin(ibl)
        do i=1,nbas
          diagl(ibl,i)=diag(i)
          bdiagl(ibl,i)=bdiag(i)
          bsdiagl(ibl,i)=bsdiag(i)
          csdiagl(ibl,i)=csdiag(i)
        enddo
!$acc end kernels

      endif
!$acc kernels present(aat,aaa,tt,ta)
      do i=1,nbas 
        do k=1,nbas 
          tt(i,k)=ta(k,i)
          aat(i,k)=aaa(k,i)
        enddo
      enddo
!$acc end kernels

!     ASSEMBLE ARRAYS ON GPU
      
!$acc kernels present(aaal,aatl,ttl,tatl,tt,aat,aaa,ta) copyin(ibl)
      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo
!$acc end kernels

!     END ASSEMBLY ON GPU
      
!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
!$acc update self(ta,aaa,aat,tt)

      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo

!     END ASSEMBLY ON CPU

      call timer_stop(30)

      endif

      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if(ihc.eq.nhc) then

c!$acc update device(tatl,aaal,aatl,ttl,sml,diagl,bdiagl,bsdiagl,csdiagl)
!$acc update device(tatl,aaal,aatl,ttl,sml)

        etotb=0.0d0
        
!$acc data copyin(prefac)
        
c        if(nb0.gt.0.and.nb1.le.0) then
        if(nb0.gt.0) then
          call timer_start(31)

          kl=nbas*(nbas+1)/2
          intg=0
!$acc kernels present(lab,g,sml,aatl,aaal,ttl,tatl,prefac)
          do ii=intndx,jntndx
            do jj=ii,kl
              intg=intg+1
              i=lab(1,ii)
              k=lab(2,ii)
              l=lab(1,jj)
              n=lab(2,jj)
              do ibl=1,nb0
                etotb=etotb+g(intg)*
     &               prefac(ibl)*(sml(ibl,k,i)*sml(ibl,n,l)
     &               -aatl(ibl,n,i)*aaal(ibl,l,k)
     &               -ttl(ibl,n,i)*tatl(ibl,l,k)
     &               -aaal(ibl,n,i)*aatl(ibl,l,k)
     &               -tatl(ibl,n,i)*ttl(ibl,l,k)
     &               -aatl(ibl,l,i)*aaal(ibl,n,k)
     &               -ttl(ibl,l,i)*tatl(ibl,n,k)
     &               -aaal(ibl,l,i)*aatl(ibl,n,k)
     &               -tatl(ibl,l,i)*ttl(ibl,n,k))
              enddo
            enddo
          enddo
!$acc end kernels

          call timer_stop(31)
        endif

c        if(nb0.gt.0.and.nb1.gt.0) then
cc        if(nb0.gt.0.and.nb1.gt.0.and.1.eq.2) then
c          call timer_start(31)
c
c       kl=nbas*(nbas+1)/2
c       intg=0
c!$acc kernels present(lab,g,sml,aatl,aaal,ttl,tatl,prefac)
c!$acc& present(diagl,bdiagl,bsdiagl,csdiagl)
c       do ii=intndx,jntndx
c         do jj=ii,kl
c           intg=intg+1
c           i=lab(1,ii)
c           k=lab(2,ii)
c           l=lab(1,jj)
c           n=lab(2,jj)
c            do ibl=1,nb0
c              etotb=etotb+g(intg)*prefac(ibl)*(sml(ibl,k,i)*sml(ibl,n,l)
c     &             -aatl(ibl,n,i)*aaal(ibl,l,k)
c     &             -ttl(ibl,n,i)*tatl(ibl,l,k)
c     &             -aaal(ibl,n,i)*aatl(ibl,l,k)
c     &             -tatl(ibl,n,i)*ttl(ibl,l,k)
c     &             -aatl(ibl,l,i)*aaal(ibl,n,k)
c     &             -ttl(ibl,l,i)*tatl(ibl,n,k)
c     &             -aaal(ibl,l,i)*aatl(ibl,n,k)
c     &             -tatl(ibl,l,i)*ttl(ibl,n,k))
c            enddo
c            do ibl=nhc-nb1+1,nhc
c              etotb=etotb+g(intg)*prefac(ibl)*
c     &             (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
c     &             +diagl(ibl,n)*csdiagl(ibl,l)
c     &             +bdiagl(ibl,l)*bsdiagl(ibl,n)
c     &             +bdiagl(ibl,n)*bsdiagl(ibl,l))
c     &             +(diagl(ibl,i)*csdiagl(ibl,k)
c     &             +diagl(ibl,k)*csdiagl(ibl,i)
c     &             +bdiagl(ibl,i)*bsdiagl(ibl,k)
c     &             +bdiagl(ibl,k)*bsdiagl(ibl,i))*sml(ibl,n,l)
c     &             -csdiagl(ibl,i)*(aatl(ibl,n,k)*diagl(ibl,l)
c     &             +aatl(ibl,l,k)*diagl(ibl,n))
c     &             -bsdiagl(ibl,i)*(ttl(ibl,n,k)*bdiagl(ibl,l)
c     &             +ttl(ibl,l,k)*bdiagl(ibl,n))
c     &             -csdiagl(ibl,l)*(aaal(ibl,n,i)*diagl(ibl,k)
c     &             +aaal(ibl,n,k)*diagl(ibl,i))
c     &             -bsdiagl(ibl,l)*(tatl(ibl,n,i)*bdiagl(ibl,k)
c     &             +tatl(ibl,n,k)*bdiagl(ibl,i))
c     &             -csdiagl(ibl,k)*(aatl(ibl,n,i)*diagl(ibl,l)
c     &             +aatl(ibl,l,i)*diagl(ibl,n))
c     &             -bsdiagl(ibl,k)*(ttl(ibl,n,i)*bdiagl(ibl,l)
c     &             +ttl(ibl,l,i)*bdiagl(ibl,n))
c     &             -csdiagl(ibl,n)*(aaal(ibl,l,i)*diagl(ibl,k)
c     &             +aaal(ibl,l,k)*diagl(ibl,i))
c     &             -bsdiagl(ibl,n)*(tatl(ibl,l,i)*bdiagl(ibl,k)
c     &             +tatl(ibl,l,k)*bdiagl(ibl,i)))
c            enddo
c          enddo
c!$acc end kernels
c
c          call timer_stop(31)
c        endif
        
c        if(nb0.le.0.and.nb1.gt.0) then
        if(nb1.gt.0) then
          call timer_start(34)

          kl=nbas*(nbas+1)/2
          intg=0
!$acc kernels present(lab,g,sml,aatl,aaal,ttl,tatl,prefac)
!$acc& present(diagl,bdiagl,bsdiagl,csdiagl)
          do ii=intndx,jntndx
            do jj=ii,kl
              intg=intg+1
              i=lab(1,ii)
              k=lab(2,ii)
              l=lab(1,jj)
              n=lab(2,jj)
              do ibl=nhc-nb1+1,nhc
                etotb=etotb+g(intg)*prefac(ibl)*
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
            enddo
          enddo
!$acc end kernels
    
          call timer_stop(34)
        endif
        
!$acc end data

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

!$OMP PARALLEL DO REDUCTION(+:ts)
!$OMP& SHARED(sm,aat,aaa,tt,ta,g,lab,numint)
!$OMP& PRIVATE(intg,i,k,l,n)
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
!$OMP END PARALLEL DO

       call timer_stop(32)

      else

       call timer_start(35)

!$OMP PARALLEL DO REDUCTION(+:e2)
!$OMP& SHARED(diag,bdiag,csdiag,bsdiag,sm,aat,tt,aaa,ta,g,lab,numint)
!$OMP& PRIVATE(intg,i,k,l,n,aai,abi,bai,bbi,aak,abk,bak,bbk)
!$OMP& PRIVATE(aaj,abj,baj,bbj,sum2,aal,abl,bal,bbl)
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
!$OMP END PARALLEL DO

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
!$OMP PARALLEL DO REDUCTION(+:ts)
!$OMP& SHARED(sm,aat,aaa,tt,ta,g,lab,numint)
!$OMP& PRIVATE(intg,i,k,l,n)
       do ii=intndx,jntndx
         do jj=ii,kl
           intg=intg+1
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
!$OMP END PARALLEL DO

       call timer_stop(32)

      else

       call timer_start(35)

       kl=nbas*(nbas+1)/2
       intg=0
!$OMP PARALLEL DO REDUCTION(+:e2)
!$OMP& SHARED(diag,bdiag,csdiag,bsdiag,sm,aat,tt,aaa,ta,g,lab,numint)
!$OMP& PRIVATE(intg,i,k,l,n,aai,abi,bai,bbi,aak,abk,bak,bbk)
!$OMP& PRIVATE(aaj,abj,baj,bbj,sum2,aal,abl,bal,bbl)
       do ii=intndx,jntndx
         do jj=ii,kl
           intg=intg+1
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
!$OMP END PARALLEL DO

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
      factor(ibl)=fctr

      if(iamhead.eq.1) sstot=sstot+factor(ibl)*deta

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

!$OMP PARALLEL DO
!$OMP& SHARED(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas)
!$OMP& PRIVATE(i,k)
      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo
!$OMP END PARALLEL DO
          
      if(ising.gt.0) then

!$OMP PARALLEL DO
!$OMP& SHARED(diagl,bdiagl,csdiagl,bsdiagl,nbas)
!$OMP& PRIVATE(i)
        do i=1,nbas
          diagl(ibl,i)=diag(i)
          bdiagl(ibl,i)=bdiag(i)
          bsdiagl(ibl,i)=bsdiag(i)
          csdiagl(ibl,i)=csdiag(i)
        enddo
!$OMP END PARALLEL DO

      endif
      
      call timer_stop(30)

      endif
      
      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if(ihc.eq.nhc) then

        call timer_start(38)
        
        if(nb0.gt.0.and.nb1.le.0) then
          call timer_start(32)

!$OMP PARALLEL DO REDUCTION(+:etotb)
!$OMP& SHARED(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,prefac)
!$OMP& PRIVATE(intg,i,k,l,n,ibl)
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
!$OMP END PARALLEL DO

          call timer_stop(32)
        endif

        if(nb0.gt.0.and.nb1.gt.0) then
          call timer_start(36)
!$OMP PARALLEL DO REDUCTION(+:etotb)
!$OMP& SHARED(diagl,bdiagl,csdiagl,bsdiagl)
!$OMP& SHARED(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,nhc,nb1,prefac)
!$OMP& PRIVATE(intg,i,k,l,n,ibl)
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
!$OMP END PARALLEL DO
          call timer_stop(36)
        endif
        
        if(nb0.le.0.and.nb1.gt.0) then
          call timer_start(35)
!$OMP PARALLEL DO REDUCTION(+:etotb)
!$OMP& SHARED(diagl,bdiagl,csdiagl,bsdiagl)
!$OMP& SHARED(sml,aatl,aaal,ttl,tatl,g,lab,numint,nb0,nhc,nb1,prefac)
!$OMP& PRIVATE(intg,i,k,l,n,ibl)
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
!$OMP END PARALLEL DO
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

      real (kind=8) :: e2n,tsn,gtemp
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
      factor(ibl)=fctr

      if(iamhead.eq.1) sstot=sstot+factor(ibl)*deta

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

      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aat(i,k)+aaa(i,k)+tt(i,k)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo
          
      if(ising.gt.0) then

        do i=1,nbas
          diagl(ibl,i)=diag(i)
          bdiagl(ibl,i)=bdiag(i)
          bsdiagl(ibl,i)=bsdiag(i)
          csdiagl(ibl,i)=csdiag(i)
        enddo

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
          intg=0
          do ii=intndx,jntndx
            do jj=ii,kl
              intg=intg+1
              i=lab(1,ii)
              k=lab(2,ii)
              l=lab(1,jj)
              n=lab(2,jj)
              do ibl=1,nb0
                etotb=etotb+g(intg)*prefac(ibl)*
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
            enddo
          enddo

          call timer_stop(32)
        endif

        if(nb0.gt.0.and.nb1.gt.0) then
          call timer_start(36)
          kl=nbas*(nbas+1)/2
          intg=0
          do ii=intndx,jntndx
            do jj=ii,kl
              intg=intg+1
              i=lab(1,ii)
              k=lab(2,ii)
              l=lab(1,jj)
              n=lab(2,jj)
              do ibl=1,nb0
                etotb=etotb+g(intg)*prefac(ibl)*
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
                etotb=etotb+g(intg)*prefac(ibl)*
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
            enddo
          enddo
          call timer_stop(36)
        endif
        
        if(nb0.le.0.and.nb1.gt.0) then
          call timer_start(35)
          kl=nbas*(nbas+1)/2
          intg=0
          do ii=intndx,jntndx
            do jj=ii,kl
              intg=intg+1
              i=lab(1,ii)
              k=lab(2,ii)
              l=lab(1,jj)
              n=lab(2,jj)
              do ibl=nhc-nb1+1,nhc
                etotb=etotb+g(intg)*prefac(ibl)*
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
            enddo
          enddo
          call timer_stop(35)
        endif
        
      endif

      return
      end

      
