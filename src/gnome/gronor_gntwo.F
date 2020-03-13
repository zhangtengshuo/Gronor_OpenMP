      subroutine gronor_gntwo(lfndbg)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENACC
      use openacc
c      use cuda_functions
#endif

      implicit none
!tps#include 'mpif.h'

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

      if(idbg.ge.13) write(lfndbg,600)
 600  format(' Calculating two electron matrix element',/)

c      numint=ig(mod(me,mgr)+1)
      numint=mint2

      if(naccel.eq.0) then
      inf0=0
      inl0=0
      inf1=0
      inl1=0
      iaf0=1
      ial0=numint
      iaf1=1
      ial1=numint
      elseif(naccel.lt.0) then
      inf0=1
      inl0=numint
      inf1=1
      inl1=numint
      iaf0=0
      ial0=0
      iaf1=0
      ial1=0
      elseif(naccel.lt.0) then
!     not yet implemented
      endif

      if(ising.gt.1) e1=0.0d0
      e2=0.0d0
      ts=0.0d0
      e2n=0.0d0
      tsn=0.0d0
      fourdet=4.0d0*deta

!$acc data create(sm)

      call timer_start(30)

      if((iaf0.gt.0).or.(iaf1.gt.0)) then
!$acc kernels present(aat,aaa,tt,ta,tat,sm)
*       do i=1,nsymf
*        do k=1,nsymf
       do i=1,nbas 
        do k=1,nbas
        tat(k,i)=ta(k,i)
        tt(i,k)=tat(k,i)
        aat(i,k)=aaa(k,i)
       enddo
      enddo
*       do i=1,nsymf
*        do k=1,nsymf
       do i=1,nbas 
        do k=1,nbas
         sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+tat(k,i)
        enddo
       enddo
!$acc end kernels
      endif

      if((inf0.gt.0).or.(inf1.gt.0)) then
!$acc update self(aaa,ta,tt,aat,tat)
*       do i=1,nsymf
*        do k=1,nsymf
       do i=1,nbas 
        do k=1,nbas
        tt(i,k)=ta(k,i)
        aat(i,k)=aaa(k,i)
       enddo
      enddo
*       do i=1,nsymf
*        do k=1,nsymf
       do i=1,nbas 
        do k=1,nbas
         sm(k,i)=aat(k,i)+aaa(k,i)+tt(k,i)+tat(k,i)
        enddo
       enddo
      endif

      call timer_stop(30)

      call timer_start(38)

      if(ising.le.0) then

      if(iaf0.gt.0) then

       call timer_start(31)

      tst=ts
!$acc kernels present(aat,aaa,tt,ta,tat,sm,g)
!$acc loop reduction (+:tst)
       do int=1,numint
        i=lab(1,int)
        k=lab(2,int)
        l=lab(3,int)
        n=lab(4,int)
        tst=tst+g(int)*(sm(k,i)*sm(n,l)
     &       -aat(n,i)*aaa(l,k)-tt(n,i)*tat(l,k)
     &       -aaa(n,i)*aat(l,k)-tat(n,i)*tt(l,k)
     &       -aat(l,i)*aaa(n,k)-tt(l,i)*tat(n,k)
     &       -aaa(l,i)*aat(n,k)-tat(l,i)*tt(n,k))
       enddo
!$acc end kernels
       ts=tst

       call timer_stop(31)
      endif
       if(inf0.gt.0) then
       call timer_start(32)

       do int=1,numint
        i=lab(1,int)
        k=lab(2,int)
        l=lab(3,int)
        n=lab(4,int)
        tsn=tsn+g(int)*(sm(k,i)*sm(n,l)
     &       -aat(n,i)*aaa(l,k)-tt(n,i)*tat(l,k)
     &       -aaa(n,i)*aat(l,k)-tat(n,i)*tt(l,k)
     &       -aat(l,i)*aaa(n,k)-tt(l,i)*tat(n,k)
     &       -aaa(l,i)*aat(n,k)-tat(l,i)*tt(n,k))
       enddo


       call timer_stop(32)
       endif
       call timer_start(33)
c!$acc wait
       call timer_stop(33)
      else
       if(iaf1.gt.0) then
       call timer_start(34)

       e2t=e2
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
     &    (sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l)
     &          -bai*(aat(n,k)*aaj+aat(l,k)*aal)
     &          -bbi*(tt(n,k)*abj+tt(l,k)*abl)
     &          -baj*(aaa(n,i)*aak+aaa(n,k)*aai)
     &          -bbj*(tat(n,i)*abk+tat(n,k)*abi)
     &          -bak*(aat(n,i)*aaj+aat(l,i)*aal)
     &          -bbk*(tt(n,i)*abj+tt(l,i)*abl)
     &          -bal*(aaa(l,i)*aak+aaa(l,k)*aai)
     &          -bbl*(tat(l,i)*abk+tat(l,k)*abi))
       enddo
!$acc end kernels
       e2=e2t

       call timer_stop(34)
       endif
       if(inf1.gt.0) then
!$acc update self(diag,bdiag,csdiag,bsdiag)
       call timer_start(35)

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
        e2n=e2n+g(int)*
     &    (sm(k,i)*(aaj*bal+aal*baj+abj*bbl+abl*bbj)+sum2*sm(n,l)
     &          -bai*(aat(n,k)*aaj+aat(l,k)*aal)
     &          -bbi*(tt(n,k)*abj+tt(l,k)*abl)
     &          -baj*(aaa(n,i)*aak+aaa(n,k)*aai)
     &          -bbj*(tat(n,i)*abk+tat(n,k)*abi)
     &          -bak*(aat(n,i)*aaj+aat(l,i)*aal)
     &          -bbk*(tt(n,i)*abj+tt(l,i)*abl)
     &          -bal*(aaa(l,i)*aak+aaa(l,k)*aai)
     &          -bbl*(tat(l,i)*abk+tat(l,k)*abi))
       enddo

       call timer_stop(35)
       endif
       call timer_start(36)
c!$acc wait
       call timer_stop(36)
      endif

#ifdef _OPENACC
      if(numdev.gt.1) then
c       cpfre=c_loc(memfre)
c       cptot=c_loc(memtot)
c       istat=cudaMemGetInfo(cpfre,cptot)
       memfre=0
       memtot=0
       memavail=min(memfre,memavail)
      endif
#endif

!$acc end data

      if(inf0.gt.0) ts=ts+tsn
      if(inf1.gt.0) e2=e2+e2n
      ts=ts*fourdet
      e2=e2+ts

      call timer_stop(38)

!     load balance between CPU and GPU

      if(naccel.gt.0.and.iload.eq.1) then

       if(ising.le.0) then
        tcalc=timer_wall(38)
        twait=timer_wall(33)
        if(twait.gt.0.01d0*tcalc) then
         nacc0=min(myints(6),nacc0+1)
        else
         nacc0=max(myints(5),nacc0-1)
        endif
       else
        tcalc=timer_wall(38)
        twait=timer_wall(36)
        if(twait.gt.0.01d0*tcalc) then
         nacc1=min(myints(6),nacc1+1)
        else
         nacc1=max(myints(5),nacc1-1)
        endif
       endif

      endif


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
!tps#include 'mpif.h'

      integer (kind=4) :: ierr
      integer :: lfndbg

      integer :: i,ii,k,l,status(MPI_STATUS_SIZE)
      integer :: im,is,ie,km,ks,ke,lm,ls,le,nm,ns,ne,n,kl,ll,nl,nn,n2
      integer :: nmax,kbloc,noff,int

      real(kind=8) :: e2n,tsn,gtemp
      real (kind=8) :: tcalc, twait

      real (kind=8) :: sum,sum1,sum2,fik,fln,fikln,gnn,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

      if(idbg.ge.13) write(lfndbg,600)
 600  format(' Calculating two electron matrix element',/)

c      numint=ig(mod(me,mgr)+1)
      numint=mint2

      if(ising.gt.1) e1=0.0d0
      e2=0.0d0
      ts=0.0d0
      fourdet=4.0d0*deta

      call timer_start(30)

*      do i=1,nsymf
*       do k=1,nsymf
      do i=1,nbas 
       do k=1,nbas 
        tt(i,k)=ta(k,i)
        aat(i,k)=aaa(k,i)
       enddo
      enddo
*      do i=1,nsymf
*       do k=1,nsymf
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
!$OMP& PRIVATE(int,i,k,l,n)
       do int=1,numint
        i=lab(1,int)
        k=lab(2,int)
        l=lab(3,int)
        n=lab(4,int)
        ts=ts+g(int)*(sm(k,i)*sm(n,l)
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
!$OMP& PRIVATE(int,i,k,l,n,aai,abi,bai,bbi,aak,abk,bak,bbk)
!$OMP& PRIVATE(aaj,abj,baj,bbj,sum2,aal,abl,bal,bbl)
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
        e2=e2+g(int)*
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
