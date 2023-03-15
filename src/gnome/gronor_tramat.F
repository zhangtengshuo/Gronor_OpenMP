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
!! Transformed MOs matrix multiplication
!!
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!


      subroutine gronor_tramat()

!      This routine mutiplies the input-matrix (a)
!       (which is positioned in a workmatrix (aa))
!      with both transformed  m.o.'s va and vb
!      the result is positioned in matrix aa
!               a'=va(t)*a*vb
!      first an intermediate matrix is calculated
!      and put into the work-matrix

!      if ising # 0 : the array's diag and sdiag are
!      transformed also
!        diag'=va*diag    sdiag'=sdiag*vb

      use gnome_parameters
      use gnome_data
      implicit none
      integer :: i,j,k,kk,m1
      real (kind=8) :: sum

      m1=nalfa+1

#ifdef ACC
!$acc kernels present(va,vb,diag,sdiag,ta,taa,w1,w2)
#endif
#ifdef OMPTGT
!$omp target teams loop collapse(2) private(i,j)
#endif
#ifdef ACC
!$acc loop collapse(2)
#endif
      do i=1,nelecs
       do j=1,nelecs
        taa(j,i)=ta(j,i)
       enddo
      enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif

      if(ising.ne.0) then

!     transformation of diag
#ifdef OMPTGT
!$omp target teams loop private(sum)
#endif
       do j=1,nbas
        sum=0.0d0
        if(nalfa .ne. 0) then
#ifdef ACC
!$acc loop reduction(+:sum)
#endif
         do k=1,nalfa
          sum=sum+diag(k)*va(k,j)
         enddo
         if(ntcla .ne. 0) then
#ifdef ACC
!$acc loop reduction(+:sum)
#endif
          do k=1,ntcla
           kk=k+nalfa
           sum=sum+diag(kk)*va(k,j)
          enddo
         endif
        endif
        if(nalfa .ne. nveca) then
#ifdef ACC
!$acc loop reduction(+:sum)
#endif
         do k=m1,nveca
          kk=k+ntcla
          sum=sum+diag(kk)*va(k,j)
         enddo
        endif
        w1(j)=sum
      enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif
#ifdef OMPTGT
!$omp target teams loop
#endif
       do j=1,nbas
        diag(j)=w1(j)
       enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif
      endif


!     transformation of matrix a
!                first :
!     calculation of an intermediate array (w)

#ifdef OMPTGT
!$omp target teams loop private(sum)
#endif
#ifdef ACC
!$acc loop collapse(2)
#endif
      do j=1,nelecs
       do i=1,nbas
        sum=0.0d0
        if(nalfa .ne. 0) then
#ifdef ACC
!$acc loop seq reduction(+:sum)
#endif
         do k=1,nalfa
          sum=sum+va(k,i)*taa(k,j)
         enddo
         if(ntcla .ne. 0) then
#ifdef ACC
!$acc loop seq reduction(+:sum)
#endif
          do k=1,ntcla
           kk=k+nalfa
           sum=sum+va(k,i)*taa(kk,j)
          enddo
         endif
        endif
        if(nalfa .ne. nveca) then
#ifdef ACC
!$acc loop seq reduction(+:sum)
#endif
         do k=m1,nveca
          kk=k+ntcla
          sum=sum+va(k,i)*taa(kk,j)
         enddo
        endif
        w2(i,j)=sum
       enddo
      enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif

!     put w in the work-matrix (aa)
#ifdef OMPTGT
!$omp target teams loop collapse(2)
#endif
#ifdef ACC
!$acc loop collapse(2)
#endif
      do j=1,nelecs
       do i=1,nbas
        taa(i,j)=w2(i,j)
       enddo
      enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif

!     calculation of the final matrix aa
#ifdef OMPTGT
!$omp target teams loop collapse(2) private(sum)
#endif
#ifdef ACC
!$acc loop collapse(2)
#endif
      do i=1,nbas
       do j=1,nbas
        sum=0.0d0
        if(nalfa .ne. 0) then
#ifdef ACC
!$acc loop seq reduction(+:sum)
#endif
         do k=1,nalfa
          sum=sum+taa(i,k)*vb(k,j)
         enddo
         if(ntclb .ne. 0) then
#ifdef ACC
!$acc loop seq reduction(+:sum)
#endif
          do k=1,ntclb
           kk=k+nalfa
           sum=sum+taa(i,kk)*vb(k,j)
          enddo
         endif
        endif
        if(nalfa .ne. nvecb) then
#ifdef ACC
!$acc loop reduction(+:sum)
#endif
         do k=m1,nvecb
          kk=k+ntclb
          sum=sum+taa(i,kk)*vb(k,j)
         enddo
        endif
        w2(i,j)=sum
       enddo
      enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif

!     put w back in the work-matrix aa
#ifdef OMPTGT
!$omp target teams loop collapse(2)
#endif
#ifdef ACC
!$acc loop collapse(2)
#endif
      do j=1,nbas
       do i=1,nbas
        taa(i,j)=w2(i,j)
       enddo
      enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif

!        transformation of sdiag
      if(ising.ne.0) then
#ifdef OMPTGT
!$omp target teams loop private(sum)
#endif
       do j=1,nbas
        sum=0.0d0
        if(nalfa .ne. 0) then
#ifdef ACC
!$acc loop seq reduction(+:sum)
#endif
         do k=1,nalfa
          sum=sum+sdiag(k)*vb(k,j)
         enddo
         if(ntclb .ne. 0) then
#ifdef ACC
!$acc loop seq reduction(+:sum)
#endif
          do k=1,ntclb
           kk=k+nalfa
           sum=sum+sdiag(kk)*vb(k,j)
          enddo
         endif
        endif
        if(nalfa .ne. nvecb) then
#ifdef ACC
!$acc loop seq reduction(+:sum)
#endif
         do k=m1,nvecb
          kk=k+ntclb
          sum=sum+sdiag(kk)*vb(k,j)
         enddo
        endif
        w1(j)=sum
      enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif
#ifdef OMPTGT
!$omp target teams loop
#endif
       do j=1,nbas
        sdiag(j)=w1(j)
       enddo
#ifdef OMPTGT
!$omp end target teams loop
#endif
      endif

#ifdef ACC
!$acc end kernels
#endif
      return
      end

      subroutine gronor_tramat_omp()

!      This routine mutiplies the input-matrix (a)
!       (which is positioned in a workmatrix (aa))
!      with both transformed  m.o.'s va and vb
!      the result is positioned in matrix aa
!               a'=va(t)*a*vb
!      first an intermediate matrix is calculated
!      and put into the work-matrix

!      if ising # 0 : the array's diag and sdiag are
!      transformed also
!        diag'=va*diag    sdiag'=sdiag*vb

      use gnome_parameters
      use gnome_data
      implicit none
      integer :: i,j,k,kk,m1
      real (kind=8) :: sum

      m1=nalfa+1
#ifdef OMP
!$omp parallel do                                                                 &
!$omp& shared(nelecs,ta,taa) collapse(2)
#endif
      do i=1,nelecs
       do j=1,nelecs
        taa(j,i)=ta(j,i)
       enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
      if(ising.ne.0) then

!     transformation of diag

#ifdef OMP
!$omp parallel do                                                                 &
!$omp& shared(nbas,nalfa,ntcla,nveca,nelecs,diag,va,w1) private(sum)
#endif
        do j=1,nbas
          sum=0.0d0
          if(nalfa .ne. 0) then
            do k=1,nalfa
              sum=sum+diag(k)*va(k,j)
            enddo
            if(ntcla .ne. 0) then
              do k=1,ntcla
                kk=k+nalfa
                sum=sum+diag(kk)*va(k,j)
              enddo
            endif
          endif
          if(nalfa .ne. nveca) then
            do k=m1,nveca
              kk=k+ntcla
              sum=sum+diag(kk)*va(k,j)
            enddo
          endif
          w1(j)=sum
        enddo
#ifdef OMP
!$omp end parallel do
#endif
#ifdef OMP
!$omp parallel do                                                               &
!$omp& shared(diag,w1)
#endif
        do j=1,nbas
          diag(j)=w1(j)
        enddo
#ifdef OMP
!$omp end parallel do
#endif
      endif


!     transformation of matrix a
!                first :
!     calculation of an intermediate array (w)

#ifdef OMP
!$omp parallel do                                                               &
!$omp& shared(nbas,nalfa,ntcla,nveca,nelecs,diag,va,w2,taa)                     &
!$omp& private(sum,kk) collapse(2)
#endif
      do j=1,nelecs
        do i=1,nbas
          sum=0.0d0
          if(nalfa .ne. 0) then
            do k=1,nalfa
              sum=sum+va(k,i)*taa(k,j)
            enddo
            if(ntcla .ne. 0) then
              do k=1,ntcla
                kk=k+nalfa
                sum=sum+va(k,i)*taa(kk,j)
              enddo
            endif
          endif
          if(nalfa .ne. nveca) then
            do k=m1,nveca
              kk=k+ntcla
              sum=sum+va(k,i)*taa(kk,j)
            enddo
          endif
          w2(i,j)=sum
        enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
!     put w in the work-matrix (aa)

#ifdef OMP
!$omp parallel do                                                               &
!$omp& shared(nbas,nelecs,w2,taa)                                               &
!$omp& private(sum,kk) collapse(2)
#endif
      do j=1,nelecs
       do i=1,nbas
        taa(i,j)=w2(i,j)
       enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
!     calculation of the final matrix aa

#ifdef OMP
!$omp parallel do                                                               &
!$omp& shared(nbas,nelecs,nalfa,ntclb,nvecb,w2,taa,vb)                          &
!$omp& private(sum,kk) collapse(2)
#endif
      do i=1,nbas
        do j=1,nbas
          sum=0.0d0
          if(nalfa .ne. 0) then
            do k=1,nalfa
              sum=sum+taa(i,k)*vb(k,j)
            enddo
            if(ntclb .ne. 0) then
              do k=1,ntclb
                kk=k+nalfa
                sum=sum+taa(i,kk)*vb(k,j)
              enddo
            endif
          endif
          if(nalfa .ne. nvecb) then
            do k=m1,nvecb
              kk=k+ntclb
              sum=sum+taa(i,kk)*vb(k,j)
            enddo
          endif
          w2(i,j)=sum
        enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
!     put w back in the work-matrix aa

#ifdef OMP
!$omp parallel do                                                                 &
!$omp& shared(nbas,w2,taa)                                                        &
!$omp& collapse(2)
#endif
      do j=1,nbas
       do i=1,nbas
        taa(i,j)=w2(i,j)
       enddo
      enddo
#ifdef OMP
!$omp end parallel do
#endif
!        transformation of sdiag

      if(ising.ne.0) then
#ifdef OMP
!$omp parallel do                                                                &
!$omp& shared(nbas,nelecs,nalfa,ntclb,nvecb,w1,vb,sdiag)                         &
!$omp& private(sum,kk)
#endif
        do j=1,nbas
          sum=0.0d0
          if(nalfa .ne. 0) then
            do k=1,nalfa
              sum=sum+sdiag(k)*vb(k,j)
            enddo
            if(ntclb .ne. 0) then
              do k=1,ntclb
                kk=k+nalfa
                sum=sum+sdiag(kk)*vb(k,j)
              enddo
            endif
          endif
          if(nalfa .ne. nvecb) then
            do k=m1,nvecb
              kk=k+ntclb
              sum=sum+sdiag(kk)*vb(k,j)
            enddo
          endif
          w1(j)=sum
        enddo
#ifdef OMP
!$omp end parallel do
#endif
#ifdef OMP
!$omp parallel do                                                                &
!$omp& shared(nbas,sdiag,w1)
#endif
        do j=1,nbas
          sdiag(j)=w1(j)
        enddo
#ifdef OMP
!$omp end parallel do
#endif
      endif

      return
      end
