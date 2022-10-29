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
!! Processing the one-electron integrals
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
!! These routines come in two different forms
!!
!! <table>
!! <caption id="multi_row">gntwo routines</caption>
!! <tr><th> Routine name <th> OpenACC <th> OpenMP
!! <tr><td> gronor_gntwo <td> \htmlonly &#x2714;\endhtmlonly <td>
!! <tr><td> gronor_gntwo_omp <td> <td> \htmlonly &#x2714;\endhtmlonly
!! </table>
!!


      subroutine gronor_gnone(lfndbg)
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

      implicit none
      integer :: lfndbg

      integer :: j,k,ielem,jkoff,nn
      real (kind=8) :: tsum,vsum,abjk,potnuc1,tsj,vsj

      nn=0
      tsum=0.0d0
      vsum=0.0d0
      ielem=0
      jkoff=0

#ifdef ACC
!$acc kernels present(t,v,diag,bdiag,bsdiag,csdiag,ta,aaa,ndxtv)
#endif
      if(ising.eq.0) then
#ifdef ACC
!$acc loop reduction(+:tsum,vsum) private(tsj,vsj,abjk,j,k,nn)
#endif
#ifdef OMPTGT
!$omp target teams distribute parallel do simd                                  &
!$omp& private(tsj,vsj,abjk,j,k,nn)
#endif
        do j=1,nbas
         nn=ndxtv(j)
         tsj=0.0d0
         vsj=0.0d0
         do k=1,j
          abjk=(aaa(j,k)+aaa(k,j)+ta(j,k)+ta(k,j))*deta*2
!         kinetic and nuclear attraction energies
          tsj=tsj+t(nn+k)*abjk
          vsj=vsj+v(nn+k)*abjk
         enddo
         tsum=tsum+tsj
         vsum=vsum+vsj
        enddo
      else
#ifdef ACC
!$acc loop reduction(+:tsum,vsum)
#endif
#ifdef OMPTGT
!$omp target teams distribute parallel do private(tsj,vsj,abjk,j,k)
#endif
        do j=1,nbas
         nn=ndxtv(j)
         tsj=0.0d0
         vsj=0.0d0
         do k=1,j
          abjk=diag(j)*csdiag(k)+bdiag(j)*bsdiag(k)                             &
     &        +diag(k)*csdiag(j)+bdiag(k)*bsdiag(j)
!         kinetic and nuclear attraction energies
          tsj=tsj+t(nn+k)*abjk
          vsj=vsj+v(nn+k)*abjk
         enddo
         tsum=tsum+tsj
         vsum=vsum+vsj
        enddo
      endif
#ifdef ACC
!$acc end kernels
#endif
      potnuc1=potnuc*deta
      e1=tsum+vsum+potnuc1

      if(idbg.ge.12) then
       write(lfndbg,120)
       write(lfndbg,130) deta,tsum,vsum,potnuc1,e1
  110  format(1x,2i4,5(2x,e20.10))
  120  format(///,1x,'The one electron matrix elements of h are :',//)
  130  format(15x,'Total overlap              : ',f20.12,//,                    &
     &        15x,'Kinetic energy term        : ',f20.12,//,                    &
     &        15x,'Nuclear attraction term    : ',f20.12,//,                    &
     &        15x,'Nuclear repulsion term     : ',f20.12,//,                    &
     &        15x,'One-electron matrix element: ',f20.12,//)
      endif

      return
      end

      subroutine gronor_gnone_omp(lfndbg)
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

      implicit none
      integer :: lfndbg

      integer :: j,k,ielem,jkoff,nn
      real (kind=8) :: tsum,vsum,abjk,potnuc1,tsj,vsj

      nn=0
      tsum=0.0d0
      vsum=0.0d0
      ielem=0
      jkoff=0

      if(ising.eq.0) then
#ifdef OMP
!$omp parallel do simd private(tsj,vsj,abjk,j,k,nn)
#endif
        do j=1,nbas
         nn=ndxtv(j)
         tsj=0.0d0
         vsj=0.0d0
         do k=1,j
          abjk=(aaa(j,k)+aaa(k,j)+ta(j,k)+ta(k,j))*deta*2
!         kinetic and nuclear attraction energies
          tsj=tsj+t(nn+k)*abjk
          vsj=vsj+v(nn+k)*abjk
         enddo
         tsum=tsum+tsj
         vsum=vsum+vsj
        enddo
#ifdef OMP
!$omp end parallel do
#endif
      else
#ifdef OMP
!$omp parallel do simd private(tsj,vsj,abjk,j,k)
#endif
        do j=1,nbas
         nn=ndxtv(j)
         tsj=0.0d0
         vsj=0.0d0
         do k=1,j
          abjk=diag(j)*csdiag(k)+bdiag(j)*bsdiag(k)                             &
     &        +diag(k)*csdiag(j)+bdiag(k)*bsdiag(j)
!         kinetic and nuclear attraction energies
          tsj=tsj+t(nn+k)*abjk
          vsj=vsj+v(nn+k)*abjk
         enddo
         tsum=tsum+tsj
         vsum=vsum+vsj
        enddo
#ifdef OMP
!$omp end parallel do
#endif
      endif

      potnuc1=potnuc*deta
      e1=tsum+vsum+potnuc1

      if(idbg.ge.12) then
       write(lfndbg,120)
       write(lfndbg,130) deta,tsum,vsum,potnuc1,e1
  110  format(1x,2i4,5(2x,e20.10))
  120  format(///,1x,'The one electron matrix elements of h are :',//)
  130  format(15x,'Total overlap              : ',f20.12,//,                    &
     &        15x,'Kinetic energy term        : ',f20.12,//,                    &
     &        15x,'Nuclear attraction term    : ',f20.12,//,                    &
     &        15x,'Nuclear repulsion term     : ',f20.12,//,                    &
     &        15x,'One-electron matrix element: ',f20.12,//)
      endif

      return
      end
