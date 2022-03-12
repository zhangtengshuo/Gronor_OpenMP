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

      integer :: lfndbg,i,k,l,n,intg

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
      
      if(idbg.ge.50) write(lfndbg,600)
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
c        istat=cudaMemGetInfo(cpfre,cptot)
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
            tst=tst+g(intg)*(sm(i,k)*sm(l,n)
     &           -aaa(i,n)*aaa(l,k)-ta(i,n)*ta(l,k)
     &           -aat(i,n)*aat(l,k)-tt(i,n)*tt(l,k)
     &           -aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k)
     &           -aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
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
c        istat=cudaMemGetInfo(cpfre,cptot)
        memavail=min(memfre,memavail)
      endif
#endif

      ts=ts*fourdet
      e2=e2+ts
      
      call timer_stop(38)
      
      return
      end
      
      subroutine gronor_gntwo_nolab(lfndbg)

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

      integer :: lfndbg,i,ii,jj,k,l,n,kl,intg,intl,ls

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
        intg=0

!_ACCTGT_($acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab,ndx))
!_ACCTGT_($acc& copyin(intg,nbas))
!_OMPTGT_($omp taskwait)
!_OMPTGT_($omp target teams distribute parallel do reduction(+:tst))
        do k=1,nbas
          do i=1,k
            tsn=0.0d0
            do n=k,nbas
              ls=1
              if(n.eq.k) ls=i
              tsl=0.0d0
              do l=ls,n
                intg=intg+1
                tsl=tsl+g(intg)*(sm(i,k)*sm(l,n)
     &               -aaa(i,n)*aaa(l,k)-ta(i,n)*ta(l,k)
     &               -aat(i,n)*aat(l,k)-tt(i,n)*tt(l,k)
     &               -aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k)
     &               -aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
              enddo
              tsn=tsn+tsl
            enddo
            tst=tst+tsn
          enddo
        enddo
c        do ii=intndx,jntndx
c          do jj=ii,kl
c            intg=ndx(ii)+jj
c            i=lab(1,ii)
c            k=lab(2,ii)
c            l=lab(1,jj)
c            n=lab(2,jj)
c            tst=tst+g(intg)*(sm(i,k)*sm(l,n)
c     &           -aaa(i,n)*aaa(l,k)-ta(i,n)*ta(l,k)
c     &           -aat(i,n)*aat(l,k)-tt(i,n)*tt(l,k)
c     &           -aat(l,i)*aaa(n,k)-tt(l,i)*ta(n,k)
c     &           -aaa(l,i)*aat(n,k)-ta(l,i)*tt(n,k))
c          enddo
c        enddo
!_ACCTGT_($acc end kernels)
        ts=tst

        call timer_stop(31)

      else
        
        call timer_start(34)

        e2t=e2
        kl=nbas*(nbas+1)/2
        intg=0

!_ACCTGT_($acc kernels present(aat,aaa,tt,ta,tat,sm,g,lab,ndx))
!_ACCTGT_($acc& present(diag,bdiag,bsdiag,csdiag))
!_ACCTGT_($acc& copyin(intg,kl,intndx,jntndx))     
!_OMPTGT_($omp target teams distribute parallel do reduction(+:e2t))
c        do ii=intndx,jntndx
c          do jj=ii,kl
c            intg=ndx(ii)+jj
c            i=lab(1,ii)
c            k=lab(2,ii)
c            l=lab(1,jj)
c     n=lab(2,jj)
        
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
                e2l=e2l+g(intl)*
     &               (sum1*(aaj*bal+aal*baj+abj*bbl+abl*bbj)
     &               +sum2*sm(n,l)
     &               -bai*(aat(n,k)*aaj+aat(l,k)*aal)
     &               -bbi*(tt(n,k)*abj+tt(l,k)*abl)
     &               -baj*(aaa(n,i)*aak+aaa(n,k)*aai)
     &               -bbj*(tat(n,i)*abk+tat(n,k)*abi)
     &               -bak*(aat(n,i)*aaj+aat(l,i)*aal)
     &               -bbk*(tt(n,i)*abj+tt(l,i)*abl)
     &               -bal*(aaa(l,i)*aak+aaa(l,k)*aai)
     &               -bbl*(tat(l,i)*abk+tat(l,k)*abi))
              enddo
              intg=intg+n+1-ls
              e2n=e2n+e2l
            enddo
            e2t=e2t+e2n
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
c        istat=cudaMemGetInfo(cpfre,cptot)
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

      integer :: lfndbg,i,k,l,n,intg,ihc,nhc,ibl

      real (kind=8) :: ts,fourdet

      real (kind=8), external :: timer_wall

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
        
        if(idbg.ge.50) write(lfndbg,600)
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

      integer :: lfndbg,i,ii,jj,k,l,n,kl,intg,ihc,nhc,ibl

      real (kind=8) :: ts,fourdet

      real (kind=8), external :: timer_wall

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
        
        if(iamhead.eq.1) sstot=sstot+fctr*deta
        
        if(idbg.ge.50) write(lfndbg,600)
 600    format(' Calculating two electron matrix element',/)
        
        numint=mint2
        
        if(ising.gt.1) e1=0.0d0
        e2=0.0d0
        ts=0.0d0
        fourdet=4.0d0*deta
        
        if(ising.eq.0) then
          nb0=nb0+1
          ibl=nb0        
          prefac(ibl)=4.0d0*deta*fctr
        else
          ibl=nhc-nb1
          nb1=nb1+1
          prefac(ibl)=fctr
          
c!_ACCTGT_($acc kernels present(diag,bdiag,bsdiag,csdiag))
c!_ACCTGT_($acc& present(diagl,bdiagl,bsdiagl,csdiagl) copyin(ibl))
c!_OMPTGT_($omp target teams distribute parallel do)
c        do i=1,nbas
c          diagl(ibl,i)=diag(i)
c          bdiagl(ibl,i)=bdiag(i)
c          bsdiagl(ibl,i)=bsdiag(i)
c          csdiagl(ibl,i)=csdiag(i)
c        enddo
c!_ACCTGT_($acc end kernels)

        endif
        
!     ASSEMBLE ARRAYS ON GPU
      
!_ACCTGT_($acc kernels present(aaal,aatl,ttl,tatl,tt,aat,aaa,ta) copyin(ibl))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aaa(k,i)
          ttl(ibl,i,k)=ta(k,i)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
!_ACCTGT_($acc update self(ta,aaa,aat,tt))
!_OMPTGT_($omp target update from(ta,aaa,aat,tt))     
      do k=1,nbas
        do i=1,nbas
          sml(ibl,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaal(ibl,i,k)=aaa(i,k)
          aatl(ibl,i,k)=aaa(k,i)
          ttl(ibl,i,k)=ta(k,i)
          tatl(ibl,i,k)=ta(i,k)
        enddo
      enddo

!     END ASSEMBLY ON CPU

      call timer_stop(30)

      endif

      call timer_start(38)

!     After the last element of a task the list of energies is generated
      
      if(ihc.eq.nhc) then

!_ACCTGT_($acc update device(tatl,aaal,aatl,ttl,sml))
!_OMPTGT_($omp target update to(tatl,aaal,aatl,ttl,sml)) 

        etotb=0.0d0
        
!_ACCTGT_($acc data copyin(prefac))
!_OMPTGT_($omp target data map(to:prefac))
        
        if(nb0.gt.0) then
          call timer_start(31)

          kl=nbas*(nbas+1)/2
          intg=0
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
!_ACCTGT_($acc end kernels)

          call timer_stop(31)
        endif

        if(nb1.gt.0) then
          call timer_start(34)

          kl=nbas*(nbas+1)/2
          intg=0
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
!_ACCTGT_($acc end kernels)
    
          call timer_stop(34)
        endif
        
!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)
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

      integer :: lfndbg,i,k,l,n,intg,ihc,nhc,ibl

      real(kind=8) :: ts,fourdet

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
              sm0(i,k,nb0)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa0(i,k,nb0)=aaa(i,k)
              aat0(i,k,nb0)=aaa(k,i)
              ta0(i,k,nb0)=ta(i,k)
              tt0(i,k,nb0)=ta(k,i)
            enddo
          enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt)) 

c      do k=1,nbas
c        do i=1,nbas
c          sm0(i,k,nb0)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa0(i,k,nb0)=aaa(i,k)
c          aat0(i,k,nb0)=aaa(k,i)
c          tt0(i,k,nb0)=ta(k,i)
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
              sm1(i,k,nb1)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa1(i,k,nb1)=aaa(i,k)
              aat1(i,k,nb1)=aaa(k,i)
              ta1(i,k,nb1)=ta(i,k)
              tt1(i,k,nb1)=ta(k,i)
            enddo
          enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt)) 

c      do k=1,nbas
c        do i=1,nbas
c          sm1(i,k,nb1)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa1(i,k,nb1)=aaa(i,k)
c          aat1(i,k,nb1)=aaa(k,i)
c          tt1(i,k,nb1)=ta(k,i)
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
      
      if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.nbatch) then
        
        call timer_start(31)

c!_ACCTGT_($acc update device(ta0,aaa0,aat0,tt0,sm0))
c!_OMPTGT_($omp target update to(ta0,aaa0,aat0,tt0,sm0))
        
!_ACCTGT_($acc data copyin(prefac0))
!_OMPTGT_($omp target data map(to:prefac0))

!_ACCTGT_($acc kernels present(lab,g,sm0,aat0,aaa0,tt0,ta0,prefac0))
!_OMPTGT_($omp target teams distribute parallel do)
        do intg=1,numint
          i=lab(1,intg)
          k=lab(2,intg)
          l=lab(3,intg)
          n=lab(4,intg)
!_ACCTGT_($acc loop private(ibl) vector(32))
          do ibl=1,nb0
            etotb=etotb+g(intg)*prefac0(ibl)*(sm0(k,i,ibl)*sm0(n,l,ibl)
     &           -aat0(n,i,ibl)*aaa0(l,k,ibl)
     &           -tt0(n,i,ibl)*ta0(l,k,ibl)
     &           -aaa0(n,i,ibl)*aat0(l,k,ibl)
     &           -ta0(n,i,ibl)*tt0(l,k,ibl)
     &           -aat0(l,i,ibl)*aaa0(n,k,ibl)
     &           -tt0(l,i,ibl)*ta0(n,k,ibl)
     &           -aaa0(l,i,ibl)*aat0(n,k,ibl)
     &           -ta0(l,i,ibl)*tt0(n,k,ibl))
          enddo
        enddo
!_ACCTGT_($acc end kernels)
        
!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)

!     Reset index in buffer to zero
        nb0=0
        
        call timer_stop(31)
        
      endif
        
      if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.nbatch) then
        
        call timer_start(34)

c!_ACCTGT_($acc update device(ta1,aaa1,aat1,tt1,sm1))
c!_OMPTGT_($omp target update to(ta1,aaa1,aat1,tt1,sm1))
        
!_ACCTGT_($acc data copyin(prefac1))
!_OMPTGT_($omp target data map(to:prefac1))

!_ACCTGT_($acc kernels present(lab,g,sm1,aat1,aaa1,tt1,ta1,prefac1))
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
     &           (sm1(k,i,ibl)*(diag1(l,ibl)*csdiag1(n,ibl)
     &           +diag1(n,ibl)*csdiag1(l,ibl)
     &           +bdiag1(l,ibl)*bsdiag1(n,ibl)
     &           +bdiag1(n,ibl)*bsdiag1(l,ibl))
     &           +(diag1(i,ibl)*csdiag1(k,ibl)
     &           +diag1(k,ibl)*csdiag1(i,ibl)
     &           +bdiag1(i,ibl)*bsdiag1(k,ibl)
     &           +bdiag1(k,ibl)*bsdiag1(i,ibl))*sm1(n,l,ibl)
     &           -csdiag1(i,ibl)*(aat1(n,k,ibl)*diag1(l,ibl)
     &           +aat1(l,k,ibl)*diag1(n,ibl))
     &           -bsdiag1(i,ibl)*(tt1(n,k,ibl)*bdiag1(l,ibl)
     &           +tt1(l,k,ibl)*bdiag1(n,ibl))
     &           -csdiag1(l,ibl)*(aaa1(n,i,ibl)*diag1(k,ibl)
     &           +aaa1(n,k,ibl)*diag1(i,ibl))
     &           -bsdiag1(l,ibl)*(ta1(n,i,ibl)*bdiag1(k,ibl)
     &           +ta1(n,k,ibl)*bdiag1(i,ibl))
     &           -csdiag1(k,ibl)*(aat1(n,i,ibl)*diag1(l,ibl)
     &           +aat1(l,i,ibl)*diag1(n,ibl))
     &           -bsdiag1(k,ibl)*(tt1(n,i,ibl)*bdiag1(l,ibl)
     &           +tt1(l,i,ibl)*bdiag1(n,ibl))
     &           -csdiag1(n,ibl)*(aaa1(l,i,ibl)*diag1(k,ibl)
     &           +aaa1(l,k,ibl)*diag1(i,ibl))
     &           -bsdiag1(n,ibl)*(ta1(l,i,ibl)*bdiag1(k,ibl)
     &           +ta1(l,k,ibl)*bdiag1(i,ibl)))
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
      
!_ACCTGT_($acc kernels present(aaa0,aat0,tt0,ta0,tt,aat,aaa,ta,sm0))
!_ACCTGT_($acc& copyin(nb0))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
      do k=1,nbas
        do i=1,nbas
          sm0(i,k,nb0)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa0(i,k,nb0)=aaa(i,k)
          aat0(i,k,nb0)=aaa(k,i)
          ta0(i,k,nb0)=ta(i,k)
          tt0(i,k,nb0)=ta(k,i)
        enddo
      enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt))     
c      do k=1,nbas
c        do i=1,nbas
c          sm0(i,k,nb0)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa0(i,k,nb0)=aaa(i,k)
c          aat0(i,k,nb0)=aaa(k,i)
c          tt0(i,k,nb0)=ta(k,i)
c        enddo
c      enddo

c!     END ASSEMBLY ON CPU

      else
        nb1=nb1+1
        prefac1(nb1)=fctr
      
!     ASSEMBLE ARRAYS ON GPU
      
!_ACCTGT_($acc kernels present(aaa1,aat1,tt1,ta1,tt,aat,aaa,ta,sm1))
!_ACCTGT_($acc& copyin(nb1))
!_OMPTGT_($omp target teams distribute parallel do collapse(2))
      do k=1,nbas
        do i=1,nbas
          sm1(i,k,nb1)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
          aaa1(i,k,nb1)=aaa(i,k)
          aat1(i,k,nb1)=aaa(k,i)
          ta1(i,k,nb1)=ta(i,k)
          tt1(i,k,nb1)=ta(k,i)
        enddo
      enddo
!_ACCTGT_($acc end kernels)

!     END ASSEMBLY ON GPU
      
c!     ASSEMBLE ARRAYS ON CPU AFTER UPDATING ta,aaa,aat,tt
      
c!_ACCTGT_($acc update self(ta,aaa,aat,tt))
c!_OMPTGT_($omp target update from(ta,aaa,aat,tt))     
c      do k=1,nbas
c        do i=1,nbas
c          sm1(i,k,nb1)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
c          aaa1(i,k,nb1)=aaa(i,k)
c          aat1(i,k,nb1)=aaa(k,i)
c          tt1(i,k,nb1)=ta(k,i)
c        enddo
c      enddo

c!     END ASSEMBLY ON CPU

!_ACCTGT_($acc kernels present(diag,bdiag,bsdiag,csdiag))
!_ACCTGT_($acc& present(diag1,bdiag1,bsdiag1,csdiag1) copyin(nb1))
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
      
      if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.nbatch) then
        
        call timer_start(31)

c!_ACCTGT_($acc update device(ta0,aaa0,aat0,tt0,sm0))
c!_OMPTGT_($omp target update to(ta0,aaa0,aat0,tt0,sm0)) 

!_ACCTGT_($acc data copyin(prefac0))
!_OMPTGT_($omp target data map(to:prefac0))
        
        kl=nbas*(nbas+1)/2

!_ACCTGT_($acc kernels present(lab,ndx,g,sm0,aat0,aaa0,tt0,ta0,prefac0))
!_OMPTGT_($omp target teams distribute parallel do)       
        do ibl=1,nb0
          do ii=intndx,jntndx
            do jj=ii,kl
              intg=ndx(ii)+jj
              i=lab(1,ii)
              k=lab(2,ii)
              l=lab(1,jj)
              n=lab(2,jj)
              etotb=etotb+g(intg)*prefac0(ibl)*
     &             (sm0(i,k,ibl)*sm0(l,n,ibl)
     &             -aaa0(i,n,ibl)*aaa0(l,k,ibl)
     &             -ta0(i,n,ibl)*ta0(l,k,ibl)
     &             -aat0(i,n,ibl)*aat0(l,k,ibl)
     &             -tt0(i,n,ibl)*tt0(l,k,ibl)
     &             -aat0(l,i,ibl)*aaa0(n,k,ibl)
     &             -tt0(l,i,ibl)*ta0(n,k,ibl)
     &             -aaa0(l,i,ibl)*aat0(n,k,ibl)
     &             -ta0(l,i,ibl)*tt0(n,k,ibl))
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

      if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.nbatch) then
        
        call timer_start(34)
        
c!_ACCTGT_($acc update device(ta1,aaa1,aat1,tt1,sm1))
c!_OMPTGT_($omp target update to(ta1,aaa1,aat1,tt1,sm1)) 

!_ACCTGT_($acc data copyin(prefac1))
!_OMPTGT_($omp target data map(to:prefac1))

        kl=nbas*(nbas+1)/2

!_ACCTGT_($acc kernels present(lab,ndx,g,sm1,aat1,aaa1,tt1,ta1,prefac1))
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
     &             (sm1(k,i,ibl)*(diag1(l,ibl)*csdiag1(n,ibl)
     &             +diag1(n,ibl)*csdiag1(l,ibl)
     &             +bdiag1(l,ibl)*bsdiag1(n,ibl)
     &             +bdiag1(n,ibl)*bsdiag1(l,ibl))
     &             +(diag1(i,ibl)*csdiag1(k,ibl)
     &             +diag1(k,ibl)*csdiag1(i,ibl)
     &             +bdiag1(i,ibl)*bsdiag1(k,ibl)
     &             +bdiag1(k,ibl)*bsdiag1(i,ibl))*sm1(n,l,ibl)
     &             -csdiag1(i,ibl)*(aat1(n,k,ibl)*diag1(l,ibl)
     &             +aat1(l,k,ibl)*diag1(n,ibl))
     &             -bsdiag1(i,ibl)*(tt1(n,k,ibl)*bdiag1(l,ibl)
     &             +tt1(l,k,ibl)*bdiag1(n,ibl))
     &             -csdiag1(l,ibl)*(aaa1(n,i,ibl)*diag1(k,ibl)
     &             +aaa1(n,k,ibl)*diag1(i,ibl))
     &             -bsdiag1(l,ibl)*(ta1(n,i,ibl)*bdiag1(k,ibl)
     &             +ta1(n,k,ibl)*bdiag1(i,ibl))
     &             -csdiag1(k,ibl)*(aat1(n,i,ibl)*diag1(l,ibl)
     &             +aat1(l,i,ibl)*diag1(n,ibl))
     &             -bsdiag1(k,ibl)*(tt1(n,i,ibl)*bdiag1(l,ibl)
     &             +tt1(l,i,ibl)*bdiag1(n,ibl))
     &             -csdiag1(n,ibl)*(aaa1(l,i,ibl)*diag1(k,ibl)
     &             +aaa1(l,k,ibl)*diag1(i,ibl))
     &             -bsdiag1(n,ibl)*(ta1(l,i,ibl)*bdiag1(k,ibl)
     &             +ta1(l,k,ibl)*bdiag1(i,ibl)))
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

      integer :: lfndbg,i,k,l,n,intg

      real (kind=8) :: sum2,ts,fourdet
      real (kind=8) :: aai,abi,bai,bbi,aak,abk,bak,bbk
      real (kind=8) :: aaj,abj,baj,bbj,aal,abl,bal,bbl

      real(kind=8), external :: timer_wall

      if(ising.ge.3) return
      
      if(idbg.ge.50) write(lfndbg,600)
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

      integer :: lfndbg,i,k,l,n,intg,ihc,nhc,ibl

      real (kind=8) :: ts,fourdet

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

      if(idbg.ge.50) write(lfndbg,600)
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

      integer :: lfndbg,i,ii,k,l,n,kl,intg,ihc,nhc,ibl,jj

      real (kind=8) :: etemp,ts,fourdet

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

      if(idbg.ge.50) write(lfndbg,600)
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

      integer :: lfndbg

      integer :: i,k,l,n,intg,ihc,nhc,ibl

      real (kind=8) :: ts,fourdet

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

      if(idbg.ge.50) write(lfndbg,600)
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
!$omp& private(i,k) collapse(2)
          do k=1,nbas
            do i=1,nbas
              sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa0(nb0,i,k)=aaa(i,k)
c              aat0(nb0,i,k)=aaa(k,i)
              tt0(nb0,i,k)=ta(k,i)
            enddo
          enddo
!$omp end parallel do
        else
          nb1=nb1+1
          prefac1(nb1)=fctr
!$omp parallel do
!$omp& shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas)
!$omp& private(i,k) collapse(2)
          do k=1,nbas
            do i=1,nbas
              sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa1(nb1,i,k)=aaa(i,k)
c              aat1(nb1,i,k)=aaa(k,i)
              tt1(nb1,i,k)=ta(k,i)
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

      if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.nbatch) then

        call timer_start(32)

        kl=nbas*(nbas+1)/2
        do ii=intndx,jntndx
          igg=ndx(ii)
          i=lab(1,ii)
          k=lab(2,ii)
!$omp parallel
!$omp& shared(igg,i,k,ii,sm0,aaa0,tt0,g,lab,nb0,prefac)
!$omp do reduction(+:etotb)
!$omp& private(intg,jj,l,n,ibl)
          do jj=ii,kl
            intg=igg+jj
            l=lab(1,jj)
            n=lab(2,jj)
            etemp=0.0d0
            do ibl=1,nb0
              etemp=etemp+prefac0(ibl)*
     &             (sm0(ibl,k,i)*sm0(ibl,n,l)
     &             -aaa0(ibl,i,n)*aaa0(ibl,l,k)
     &             -aaa0(ibl,n,i)*aaa0(ibl,k,l)
     &             -aaa0(ibl,i,l)*aaa0(ibl,n,k)
     &             -aaa0(ibl,l,i)*aaa0(ibl,k,n)
     &             -tt0(ibl,i,n)*tt0(ibl,l,k)
     &             -tt0(ibl,n,i)*tt0(ibl,k,l)
     &             -tt0(ibl,i,l)*tt0(ibl,n,k)
     &             -tt0(ibl,l,i)*tt0(ibl,k,n))
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

      if((ihc.eq.nhc.and.nb1.gt.0).or.nb1.eq.nbatch) then
        
        call timer_start(35)
        
        kl=nbas*(nbas+1)/2
        do ii=intndx,jntndx
          igg=ndx(ii)
          i=lab(1,ii)
          k=lab(2,ii)
!$omp parallel do reduction(+:etotb)
!$omp& shared(igg,i,k,diagl,bdiagl,csdiagl,bsdiagl)
!$omp& shared(sm1,aaa1,tt1,g,lab,nb1,prefac)
!$omp& private(intg,l,n,ibl)
          do jj=ii,kl
            intg=igg+jj
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
     &             -csdiag1(ibl,i)*(aaa1(ibl,k,n)*diag1(ibl,l)
     &             +aaa1(ibl,k,l)*diag1(ibl,n))
     &             -bsdiag1(ibl,i)*(tt1(ibl,n,k)*bdiag1(ibl,l)
     &             +tt1(ibl,l,k)*bdiag1(ibl,n))
     &             -csdiag1(ibl,l)*(aaa1(ibl,n,i)*diag1(ibl,k)
     &             +aaa1(ibl,n,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,l)*(tt1(ibl,i,n)*bdiag1(ibl,k)
     &             +tt1(ibl,k,n)*bdiag1(ibl,i))
     &             -csdiag1(ibl,k)*(aaa1(ibl,i,n)*diag1(ibl,l)
     &             +aaa1(ibl,i,l)*diag1(ibl,n))
     &             -bsdiag1(ibl,k)*(tt1(ibl,n,i)*bdiag1(ibl,l)
     &             +tt1(ibl,l,i)*bdiag1(ibl,n))
     &             -csdiag1(ibl,n)*(aaa1(ibl,l,i)*diag1(ibl,k)
     &             +aaa1(ibl,l,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,n)*(tt1(ibl,i,l)*bdiag1(ibl,k)
     &             +tt1(ibl,k,l)*bdiag1(ibl,i)))
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

      
      subroutine gronor_gntwo_omp_batch2_nolab(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer :: lfndbg

      integer :: i,ii,k,l,n,kl,intg,ihc,nhc,ibl,jj,igg,ls

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
!$omp& private(i,k) collapse(2)
          do k=1,nbas
            do i=1,nbas
              sm0(nb0,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa0(nb0,i,k)=aaa(i,k)
c              aat0(nb0,i,k)=aaa(k,i)
              tt0(nb0,i,k)=ta(k,i)
            enddo
          enddo
!$omp end parallel do
        else
          nb1=nb1+1
          prefac1(nb1)=fctr
!$omp parallel do
!$omp& shared(sml,aat,aaa,tt,ta,aaal,aatl,ttl,tatl,nbas)
!$omp& private(i,k) collapse(2)
          do k=1,nbas
            do i=1,nbas
              sm1(nb1,i,k)=aaa(k,i)+aaa(i,k)+ta(k,i)+ta(i,k)
              aaa1(nb1,i,k)=aaa(i,k)
c              aat1(nb1,i,k)=aaa(k,i)
              tt1(nb1,i,k)=ta(k,i)
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

      if((ihc.eq.nhc.and.nb0.gt.0).or.nb0.eq.nbatch) then

        call timer_start(32)

        intg=0
        do k=1,nbas
          do i=1,k
            do n=k,nbas
              ls=1
              if(n.eq.k) ls=i
              do l=ls,n
                intg=intg+1
                etemp=0.0d0
                do ibl=1,nb0
                  etemp=etemp+prefac0(ibl)*
     &                 (sm0(ibl,k,i)*sm0(ibl,n,l)
     &                 -aaa0(ibl,i,n)*aaa0(ibl,l,k)
     &                 -aaa0(ibl,n,i)*aaa0(ibl,k,l)
     &                 -aaa0(ibl,i,l)*aaa0(ibl,n,k)
     &                 -aaa0(ibl,l,i)*aaa0(ibl,k,n)
     &                 -tt0(ibl,i,n)*tt0(ibl,l,k)
     &                 -tt0(ibl,n,i)*tt0(ibl,k,l)
     &                 -tt0(ibl,i,l)*tt0(ibl,n,k)
     &                 -tt0(ibl,l,i)*tt0(ibl,k,n))
                enddo
                etotb=etotb+g(intg)*etemp
              enddo
            enddo
          enddo
        enddo
     
        kl=nbas*(nbas+1)/2
c        do ii=intndx,jntndx
c          igg=ndx(ii)
c          i=lab(1,ii)
c          k=lab(2,ii)
c!$omp parallel
c!$omp& shared(igg,i,k,ii,sm0,aaa0,tt0,g,lab,nb0,prefac)
c!$omp do reduction(+:etotb)
c!$omp& private(intg,jj,l,n,ibl)
c          do jj=ii,kl
c            intg=igg+jj
c            l=lab(1,jj)
c            n=lab(2,jj)
c            etemp=0.0d0
c            do ibl=1,nb0
c              etemp=etemp+prefac0(ibl)*
c     &             (sm0(ibl,k,i)*sm0(ibl,n,l)
c     &             -aaa0(ibl,i,n)*aaa0(ibl,l,k)
c     &             -aaa0(ibl,n,i)*aaa0(ibl,k,l)
c     &             -aaa0(ibl,i,l)*aaa0(ibl,n,k)
c     &             -aaa0(ibl,l,i)*aaa0(ibl,k,n)
c     &             -tt0(ibl,i,n)*tt0(ibl,l,k)
c     &             -tt0(ibl,n,i)*tt0(ibl,k,l)
c     &             -tt0(ibl,i,l)*tt0(ibl,n,k)
c     &             -tt0(ibl,l,i)*tt0(ibl,k,n))
c            enddo
c            etotb=etotb+g(intg)*etemp
c          enddo
c!$omp end do
c!$omp end parallel
c        enddo

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
!$omp parallel do reduction(+:etotb)
!$omp& shared(igg,i,k,diagl,bdiagl,csdiagl,bsdiagl)
!$omp& shared(sm1,aaa1,tt1,g,lab,nb1,prefac)
!$omp& private(intg,l,n,ibl)
          do jj=ii,kl
            intg=igg+jj
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
     &             -csdiag1(ibl,i)*(aaa1(ibl,k,n)*diag1(ibl,l)
     &             +aaa1(ibl,k,l)*diag1(ibl,n))
     &             -bsdiag1(ibl,i)*(tt1(ibl,n,k)*bdiag1(ibl,l)
     &             +tt1(ibl,l,k)*bdiag1(ibl,n))
     &             -csdiag1(ibl,l)*(aaa1(ibl,n,i)*diag1(ibl,k)
     &             +aaa1(ibl,n,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,l)*(tt1(ibl,i,n)*bdiag1(ibl,k)
     &             +tt1(ibl,k,n)*bdiag1(ibl,i))
     &             -csdiag1(ibl,k)*(aaa1(ibl,i,n)*diag1(ibl,l)
     &             +aaa1(ibl,i,l)*diag1(ibl,n))
     &             -bsdiag1(ibl,k)*(tt1(ibl,n,i)*bdiag1(ibl,l)
     &             +tt1(ibl,l,i)*bdiag1(ibl,n))
     &             -csdiag1(ibl,n)*(aaa1(ibl,l,i)*diag1(ibl,k)
     &             +aaa1(ibl,l,k)*diag1(ibl,i))
     &             -bsdiag1(ibl,n)*(tt1(ibl,i,l)*bdiag1(ibl,k)
     &             +tt1(ibl,k,l)*bdiag1(ibl,i)))
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

      
