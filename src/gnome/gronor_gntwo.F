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
      integer :: nmax,kbloc,noff,int

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
        do int=1,numint
          i=lab(1,int)
          k=lab(2,int)
          l=lab(3,int)
          n=lab(4,int)
          tst=tst+g(int)*(sm(k,i)*sm(n,l)
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
        if(itest.eq.12345) then
!$acc kernels present(aat,aaa,tt,ta,tat,sm,g)
!$acc& present(lab,diag,bdiag,bsdiag,csdiag)
!$acc loop reduction (+:e2t)
        do int=1,numint
          i=lab(1,int) 
          k=lab(2,int) 
          l=lab(3,int) 
          n=lab(4,int)
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
          e2t=e2t+g(int)*
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
        endif
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

#ifdef _OPENACC
      use openacc
      use cuda_functions
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: lfndbg
      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,int,ihc,nhc,ibl

      integer :: iaf0,ial0,inf0,inl0,iaf1,ial1,inf1,inl1

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet

#ifdef SINGLEP
      real (kind=4) :: e2t,e2nt,tst,tsnt,tsts
      real (kind=4) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=4) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#else
      real (kind=8) :: e2t,e2nt,tst,tsnt,tsts
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#endif

      real(kind=8), external :: timer_wall

#ifdef _OPENACC
      integer (kind=4) :: istat
      type(c_ptr) :: cpfre, cptot
#endif


!     Copying araays into batch arrays

      if(ihc.eq.1) then
        nb0=0
        nb1=0
      endif

      if(ising.eq.0) then
        nb0=nb0+1
        ibl=nb0
      else
        ibl=nhc-nb1
        nb1=nb1+1
      endif
      
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

!$acc kernels present(aat,aaa,tt,tat,sm,aatl,aaal,ttl,tatl,sml)
      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=sm(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=tat(i,k)
        enddo
      enddo
!$acc end kernels
          
      if(ising.gt.0) then
!$acc kernels present(diag,bdiag,bsdiag,csdiag)
!$acc& present(diagl,bdiagl,bsdiagl,csdiagl)
        do i=1,nbas
          diagl(ibl,i)=diag(i)
          bdiagl(ibl,i)=bdiag(i)
          bsdiagl(ibl,i)=bsdiag(i)
          csdiagl(ibl,i)=csdiag(i)
        enddo
!$acc end kernels
      endif
      
      call timer_stop(30)

      if(ihc.eq.nhc) then

        call timer_start(38)
        
        if(nb0.gt.0) then
          call timer_start(31)
          tst=ts
!$acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab)
!$acc loop reduction (+:tst)
          do int=1,numint
            i=lab(1,int)
            k=lab(2,int)
            l=lab(3,int)
            n=lab(4,int)
            do ibl=1,nb0
              tst=tst+g(int)*(sml(ibl,k,i)*sml(ibl,n,l)
     &             -aatl(ibl,n,i)*aaal(ibl,l,k)
     &             -ttl(ibl,n,i)*tatl(ibl,l,k)
     &             -aaal(ibl,n,i)*aatl(ibl,l,k)
     &             -tatl(ibl,n,i)*ttl(ibl,l,k)
     &             -aatl(ibl,l,i)*aaal(ibl,n,k)
     &             -ttl(ibl,l,i)*tatl(ibl,n,k)
     &             -aaal(ibl,l,i)*aatl(ibl,n,k)
     &             -tatl(ibl,l,i)*ttl(ibl,n,k))
            enddo
c            tst=tst+tsts*g(int)
          enddo
!$acc end kernels
          ts=tst
          call timer_stop(31)
        endif
      
        if(nb1.gt.0) then
          call timer_start(34)

          e2t=e2
        if(itest.eq.12345) then
!$acc kernels present(aat,aaa,tt,ta,tat,sm,g)
!$acc& present(lab,diag,bdiag,bsdiag,csdiag)
!$acc loop reduction (+:e2t)
          do int=1,numint
            i=lab(1,int) 
            k=lab(2,int) 
            l=lab(3,int) 
            n=lab(4,int)
c            aai=diag(i)
c            abi=bdiag(i)
c            csdiag(i)=csdiag(i)
c            bsdiag(i)=bsdiag(i)
c            diag(k)=diag(k)
c            bdiag(k)=bdiag(k)
c            csdiag(k)=csdiag(k)
c            bsdiak(k)=bsdiag(k)
c            diag(l)=diag(l)
c            bdiag(l)=bdiag(l)
c            csdiag(l)=csdiag(l)
c            bsdiag(l)=bsdiag(l)
c            diag(n)=diag(n)
c            bdiag(n)=bdiag(n)
c            csdiag(n)=csdiag(n)
c            bsdiag(n)=bsdiag(n)
            tsts=0.0d0
            do ibl=nhc-nb1+1,nhc
              sum2=diagl(ibl,i)*csdiagl(ibl,k)
     &             +diagl(ibl,k)*csdiagl(ibl,i)
     &             +bdiagl(ibl,i)*bsdiagl(ibl,k)
     &             +bdiagl(ibl,k)*bsdiagl(ibl,i)
              tsts=tsts+g(int)*
     &             (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
     &             +diagl(ibl,n)*csdiagl(ibl,l)
     &             +bdiagl(ibl,l)*bsdiagl(ibl,n)
     &             +bdiagl(ibl,n)*bsdiagl(ibl,l))+sum2*sml(ibl,n,l)
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
            e2t=e2t+tsts
          enddo
!$acc end kernels
          endif
          e2=e2t

          call timer_stop(34)
        endif
        
      endif

      
c      if(ising.le.0) then
        
c        call timer_start(31)
        
c        tst=ts
c!$acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab)
c!$acc loop reduction (+:tst)
c        do int=1,numint
c          i=lab(1,int)
c          k=lab(2,int)
c          l=lab(3,int)
c          n=lab(4,int)
c          tst=tst+g(int)*(sm(k,i)*sm(n,l)
c     &         -aat(n,i)*aaa(l,k)-tt(n,i)*tat(l,k)
c     &         -aaa(n,i)*aat(l,k)-tat(n,i)*tt(l,k)
c     &         -aat(l,i)*aaa(n,k)-tt(l,i)*tat(n,k)
c     &         -aaa(l,i)*aat(n,k)-tat(l,i)*tt(n,k))
c        enddo
c!$acc end kernels
c        ts=tst

c        call timer_stop(31)

c      else
        
c        call timer_start(34)

c        e2t=e2
c!$acc kernels present(aat,aaa,tt,ta,tat,sm,g)
c!$acc& present(lab,diag,bdiag,bsdiag,csdiag)
c!$acc loop reduction (+:e2t)
c        do int=1,numint
c          i=lab(1,int) 
c          k=lab(2,int) 
c          l=lab(3,int) 
c          n=lab(4,int)
c          aai=diag(i)
c          abi=bdiag(i)
c          bai=csdiag(i)
c          bbi=bsdiag(i)
c          aak=diag(k)
c          abk=bdiag(k)
c          bak=csdiag(k)
c          bbk=bsdiag(k)
c          sum2=aai*bak+aak*bai+abi*bbk+abk*bbi
c          aaj=diag(l)
c          abj=bdiag(l)
c          baj=csdiag(l)
c          bbj=bsdiag(l)
c          aal=diag(n)
c          abl=bdiag(n)
c          bal=csdiag(n)
c          bbl=bsdiag(n)
c          e2t=e2t+g(int)*
c     &         (sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l)
c     &         -bai*(aat(n,k)*aaj+aat(l,k)*aal)
c     &         -bbi*(tt(n,k)*abj+tt(l,k)*abl)
c     &         -baj*(aaa(n,i)*aak+aaa(n,k)*aai)
c     &         -bbj*(tat(n,i)*abk+tat(n,k)*abi)
c     &         -bak*(aat(n,i)*aaj+aat(l,i)*aal)
c     &         -bbk*(tt(n,i)*abj+tt(l,i)*abl)
c     &         -bal*(aaa(l,i)*aak+aaa(l,k)*aai)
c     &         -bbl*(tat(l,i)*abk+tat(l,k)*abi))
c        enddo
c!$acc end kernels
c        e2=e2t

c        call timer_stop(34)

c      endif

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

!      numint=ig(mod(me,mgr)+1)

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
c      write(*,'(3f12.6)') ts,e2,e1
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
      integer :: nmax,kbloc,noff,intg,ihc,nhc,ibl,int

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

c      write(*,'(a,4i5)') 'Start gntwo ',ihc,nhc,ising
!     Copying arrays into batch arrays

!     Initialize nb0 and nb1 when first matrix element in a batch

!     nb0 : number of singularity 0 elements
!     nb1 : number of singularity 1 or 2 elements
      
      if(ihc.eq.1) then
        nb0=0
        nb1=0
      endif
      ibdata(ihc,1)=-1
      ibdata(ihc,2)=-1
      
      if(ising.lt.3) then

      if(ising.eq.0) then
        nb0=nb0+1
        ibl=nb0
      else
        ibl=nhc-nb1
        nb1=nb1+1
      endif

!     Task list data arrays

!     ibdata(ihc,1) : index into task list arrays ddata and hdata
!     ibdata(ihc,2) : singularity
      
!     hdata(ihc,1)  : one-electron energy
!     hdata(ihc,2)  :
!     hdata(ihc,3)  : 

!     ddata(ibl,1)  : energy
!     ddata(ibl,2)  : deta

      
      ibdata(ihc,1)=ibl
      ibdata(ihc,2)=ising

      bdata(ibl,2)=deta
      bdata(ibl,3)=bias

      ddata(ibl,2)=deta
      
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
          tat(k,i)=ta(k,i)
          tt(i,k)=ta(k,i)
          aat(i,k)=aaa(k,i)
        enddo
      enddo

      do i=1,nbas 
       do k=1,nbas
         sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+ta(k,i)
        enddo
      enddo

      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=sm(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=tat(i,k)
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

      hdata(ihc,1)=e1
      
      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if(ihc.eq.nhc) then

        call timer_start(38)
        
        if(nb0.gt.0.and.nb1.le.0) then
          call timer_start(31)
          do i=1,nbatch
            esum(i)=0.0d0
          enddo
          do int=1,numint
            i=lab(1,int)
            k=lab(2,int)
            l=lab(3,int)
            n=lab(4,int)
            do ibl=1,nb0
              esum(ibl)=esum(ibl)+g(int)*(sml(ibl,k,i)*sml(ibl,n,l)
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
          do i=1,nb0
            ddata(i,1)=4.0d0*ddata(i,2)*esum(i)
! bdata to become obsolete
            bdata(i,1)=4.0d0*bdata(i,2)*esum(i)
          enddo
          call timer_stop(31)
        endif

        if(nb0.gt.0.and.nb1.gt.0) then
          call timer_start(31)
          do i=1,nbatch
            esum(i)=0.0d0
            fsum(i)=0.0d0
          enddo
          do int=1,numint
            i=lab(1,int)
            k=lab(2,int)
            l=lab(3,int)
            n=lab(4,int)
            do ibl=1,nb0
              esum(ibl)=esum(ibl)+g(int)*(sml(ibl,k,i)*sml(ibl,n,l)
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
              fsum(ibl)=fsum(ibl)+g(int)*
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
          do i=1,nb0
            ddata(i,1)=4.0d0*ddata(i,2)*esum(i)
! bdata to become obsolete
            bdata(i,1)=4.0d0*bdata(i,2)*esum(i)
          enddo
          do i=nhc-nb1+1,nhc
            ddata(i,1)=fsum(i)
            bdata(i,1)=fsum(i)
          enddo
          call timer_stop(31)
        endif
        
        if(nb0.le.0.and.nb1.gt.0) then
          call timer_start(31)
          do i=1,nbatch
            fsum(i)=0.0d0
          enddo
          do int=1,numint
            i=lab(1,int) 
            k=lab(2,int) 
            l=lab(3,int) 
            n=lab(4,int)
            do ibl=nhc-nb1+1,nhc
              fsum(ibl)=fsum(ibl)+g(int)*
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
          do i=nhc-nb1+1,nhc
            ddata(i,1)=fsum(i)
            bdata(i,1)=fsum(i)
          enddo
          call timer_stop(31)
        endif
        
      endif

c      print*,'Returning from gntwo'

      return
      end

      subroutine gronor_gntwo_omp_batch2(lfndbg,ihc,nhc)

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
      integer :: nmax,kbloc,noff,int,ihc,nhc,ibl

      integer :: iaf0,ial0,inf0,inl0,iaf1,ial1,inf1,inl1

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet

#ifdef SINGLEP
      real (kind=4) :: e2t,e2nt,tst,tsnt,tsts
      real (kind=4) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=4) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#else
      real (kind=8) :: e2t,e2nt,tst,tsnt,tsts
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl
#endif

      real(kind=8), external :: timer_wall

#ifdef _OPENACC
      integer (kind=4) :: istat
      type(c_ptr) :: cpfre, cptot
#endif


!     Copying araays into batch arrays

      if(ihc.eq.1) then
        nb0=0
        nb1=0
      endif

      if(ising.eq.0) then
        nb0=nb0+1
        ibl=nb0
      else
        ibl=nhc-nb1
        nb1=nb1+1
      endif
      
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

      call timer_start(30)
      
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

      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=sm(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aat(i,k)
          ttl(ibl,i,k)=tt(i,k)
          tatl(ibl,i,k)=tat(i,k)
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

      if(ihc.eq.nhc) then

c        write(*,'(a,4i5)') 'GNTWO ',nhc,nb0,nb1
        call timer_start(38)
        
        if(nb0.gt.0) then
          call timer_start(31)
          tst=ts

          do int=1,numint
            i=lab(1,int)
            k=lab(2,int)
            l=lab(3,int)
            n=lab(4,int)
            do ibl=1,nb0
              tst=tst+g(int)*(sml(ibl,k,i)*sml(ibl,n,l)
     &             -aatl(ibl,n,i)*aaal(ibl,l,k)
     &             -ttl(ibl,n,i)*tatl(ibl,l,k)
     &             -aaal(ibl,n,i)*aatl(ibl,l,k)
     &             -tatl(ibl,n,i)*ttl(ibl,l,k)
     &             -aatl(ibl,l,i)*aaal(ibl,n,k)
     &             -ttl(ibl,l,i)*tatl(ibl,n,k)
     &             -aaal(ibl,l,i)*aatl(ibl,n,k)
     &             -tatl(ibl,l,i)*ttl(ibl,n,k))
            enddo

c            tst=tst+tsts*g(int)
          enddo

          ts=tst
          call timer_stop(31)
        endif
      
        if(nb1.gt.0) then
          call timer_start(34)

          e2t=e2
        if(itest.eq.12345) then

          do int=1,numint
            i=lab(1,int) 
            k=lab(2,int) 
            l=lab(3,int) 
            n=lab(4,int)
c            aai=diag(i)
c            abi=bdiag(i)
c            csdiag(i)=csdiag(i)
c            bsdiag(i)=bsdiag(i)
c            diag(k)=diag(k)
c            bdiag(k)=bdiag(k)
c            csdiag(k)=csdiag(k)
c            bsdiak(k)=bsdiag(k)
c            diag(l)=diag(l)
c            bdiag(l)=bdiag(l)
c            csdiag(l)=csdiag(l)
c            bsdiag(l)=bsdiag(l)
c            diag(n)=diag(n)
c            bdiag(n)=bdiag(n)
c            csdiag(n)=csdiag(n)
c            bsdiag(n)=bsdiag(n)
            tsts=0.0d0
            do ibl=nhc-nb1+1,nhc
              sum2=diagl(ibl,i)*csdiagl(ibl,k)
     &             +diagl(ibl,k)*csdiagl(ibl,i)
     &             +bdiagl(ibl,i)*bsdiagl(ibl,k)
     &             +bdiagl(ibl,k)*bsdiagl(ibl,i)
              tsts=tsts+g(int)*
     &             (sml(ibl,k,i)*(diagl(ibl,l)*csdiagl(ibl,n)
     &             +diagl(ibl,n)*csdiagl(ibl,l)
     &             +bdiagl(ibl,l)*bsdiagl(ibl,n)
     &             +bdiagl(ibl,n)*bsdiagl(ibl,l))+sum2*sml(ibl,n,l)
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
            e2t=e2t+tsts
          enddo
          endif
          e2=e2t

          call timer_stop(34)
        endif
        
      endif

      
c      if(ising.le.0) then
        
c        call timer_start(31)
        
c        tst=ts
c!$acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab)
c!$acc loop reduction (+:tst)
c        do int=1,numint
c          i=lab(1,int)
c          k=lab(2,int)
c          l=lab(3,int)
c          n=lab(4,int)
c          tst=tst+g(int)*(sm(k,i)*sm(n,l)
c     &         -aat(n,i)*aaa(l,k)-tt(n,i)*tat(l,k)
c     &         -aaa(n,i)*aat(l,k)-tat(n,i)*tt(l,k)
c     &         -aat(l,i)*aaa(n,k)-tt(l,i)*tat(n,k)
c     &         -aaa(l,i)*aat(n,k)-tat(l,i)*tt(n,k))
c        enddo
c!$acc end kernels
c        ts=tst

c        call timer_stop(31)

c      else
        
c        call timer_start(34)

c        e2t=e2
c!$acc kernels present(aat,aaa,tt,ta,tat,sm,g)
c!$acc& present(lab,diag,bdiag,bsdiag,csdiag)
c!$acc loop reduction (+:e2t)
c        do int=1,numint
c          i=lab(1,int) 
c          k=lab(2,int) 
c          l=lab(3,int) 
c          n=lab(4,int)
c          aai=diag(i)
c          abi=bdiag(i)
c          bai=csdiag(i)
c          bbi=bsdiag(i)
c          aak=diag(k)
c          abk=bdiag(k)
c          bak=csdiag(k)
c          bbk=bsdiag(k)
c          sum2=aai*bak+aak*bai+abi*bbk+abk*bbi
c          aaj=diag(l)
c          abj=bdiag(l)
c          baj=csdiag(l)
c          bbj=bsdiag(l)
c          aal=diag(n)
c          abl=bdiag(n)
c          bal=csdiag(n)
c          bbl=bsdiag(n)
c          e2t=e2t+g(int)*
c     &         (sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l)
c     &         -bai*(aat(n,k)*aaj+aat(l,k)*aal)
c     &         -bbi*(tt(n,k)*abj+tt(l,k)*abl)
c     &         -baj*(aaa(n,i)*aak+aaa(n,k)*aai)
c     &         -bbj*(tat(n,i)*abk+tat(n,k)*abi)
c     &         -bak*(aat(n,i)*aaj+aat(l,i)*aal)
c     &         -bbk*(tt(n,i)*abj+tt(l,i)*abl)
c     &         -bal*(aaa(l,i)*aak+aaa(l,k)*aai)
c     &         -bbl*(tat(l,i)*abk+tat(l,k)*abi))
c        enddo
c!$acc end kernels
c        e2=e2t

c        call timer_stop(34)

c      endif

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
