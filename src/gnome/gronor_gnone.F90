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


subroutine gronor_gnone(lfndbg,diag,bdiag,bsdiag,csdiag,ta,aaa)
  use cidist
  use gnome_parameters
  use gnome_data
  use gnome_integrals

  implicit none
  integer :: lfndbg
  real (kind=8), intent(inout) :: diag(:),bdiag(:),bsdiag(:),csdiag(:),ta(:,:),aaa(:,:)

  integer :: j,k,ielem,jkoff,nn
  real (kind=8) :: tsum,vsum,abjk,potnuc1,tsj,vsj
  real (kind=8) :: dsum1,dsum2,dsum3
  real (kind=8) :: qsum1,qsum2,qsum3,qsum4,qsum5,qsum6
  real (kind=8) :: dsj1,dsj2,dsj3
  real (kind=8) :: qsj1,qsj2,qsj3,qsj4,qsj5,qsj6

  nn=0
  tsum=0.0d0
  vsum=0.0d0
  dsum1=0.0d0
  dsum2=0.0d0
  dsum3=0.0d0
  qsum1=0.0d0
  qsum2=0.0d0
  qsum3=0.0d0
  qsum4=0.0d0
  qsum5=0.0d0
  qsum6=0.0d0
  ielem=0
  jkoff=0

!$acc kernels present(t,v,dqm,diag,bdiag,bsdiag,csdiag,ta,aaa,ndxtv)
  if(ising.eq.0) then
!$acc loop reduction(+:tsum,vsum,dsum1,dsum2,dsum3, &
!$acc& qsum1,qsum2,qsum3,qsum4,qsum5,qsum6) private(tsj,vsj, &
!$acc& dsj1,dsj2,dsj3,qsj1,qsj2,qsj3,qsj4,qsj5,qsj6,abjk,j,k,nn)
    do j=1,nbas
      nn=ndxtv(j)
      tsj=0.0d0
      vsj=0.0d0
      dsj1=0.0d0
      dsj2=0.0d0
      dsj3=0.0d0
      qsj1=0.0d0
      qsj2=0.0d0
      qsj3=0.0d0
      qsj4=0.0d0
      qsj5=0.0d0
      qsj6=0.0d0
      do k=1,j
        abjk=(aaa(j,k)+aaa(k,j)+ta(j,k)+ta(k,j))*deta*2
        !         kinetic and nuclear attraction energies
        tsj=tsj+t(nn+k)*abjk
        vsj=vsj+v(nn+k)*abjk
        !         dipole moment     
        dsj1=dsj1-dqm(nn+k,1)*abjk
        dsj2=dsj2-dqm(nn+k,2)*abjk
        dsj3=dsj3-dqm(nn+k,3)*abjk
        !         quadrupole moment
        qsj1=qsj1-dqm(nn+k,4)*abjk
        qsj2=qsj2-dqm(nn+k,5)*abjk
        qsj3=qsj3-dqm(nn+k,6)*abjk
        qsj4=qsj4-dqm(nn+k,7)*abjk
        qsj5=qsj5-dqm(nn+k,8)*abjk
        qsj6=qsj6-dqm(nn+k,9)*abjk
      enddo
      tsum=tsum+tsj
      vsum=vsum+vsj
      dsum1=dsum1+dsj1
      dsum2=dsum2+dsj2
      dsum3=dsum3+dsj3
      qsum1=qsum1+qsj1
      qsum2=qsum2+qsj2
      qsum3=qsum3+qsj3
      qsum4=qsum4+qsj4
      qsum5=qsum5+qsj5
      qsum6=qsum6+qsj6
    enddo
  else
!$acc loop reduction(+:tsum,vsum,dsum1,dsum2,dsum3,qsum1,qsum2,qsum3,qsum4,qsum5,qsum6)
    do j=1,nbas
      nn=ndxtv(j)
      tsj=0.0d0
      vsj=0.0d0
      dsj1=0.0d0
      dsj2=0.0d0
      dsj3=0.0d0
      qsj1=0.0d0
      qsj2=0.0d0
      qsj3=0.0d0
      qsj4=0.0d0
      qsj5=0.0d0
      qsj6=0.0d0
      do k=1,j
        abjk=diag(j)*csdiag(k)+bdiag(j)*bsdiag(k)+diag(k)*csdiag(j)+bdiag(k)*bsdiag(j)
        !         kinetic and nuclear attraction energies
        tsj=tsj+t(nn+k)*abjk
        vsj=vsj+v(nn+k)*abjk
        !         dipole moment
        dsj1=dsj1-dqm(nn+k,1)*abjk
        dsj2=dsj2-dqm(nn+k,2)*abjk
        dsj3=dsj3-dqm(nn+k,3)*abjk
        !         quadrupole moment
        qsj1=qsj1-dqm(nn+k,4)*abjk
        qsj2=qsj2-dqm(nn+k,5)*abjk
        qsj3=qsj3-dqm(nn+k,6)*abjk
        qsj4=qsj4-dqm(nn+k,7)*abjk
        qsj5=qsj5-dqm(nn+k,8)*abjk
        qsj6=qsj6-dqm(nn+k,9)*abjk
      enddo
      tsum=tsum+tsj
      vsum=vsum+vsj
      dsum1=dsum1+dsj1
      dsum2=dsum2+dsj2
      dsum3=dsum3+dsj3
      qsum1=qsum1+qsj1
      qsum2=qsum2+qsj2
      qsum3=qsum3+qsj3
      qsum4=qsum4+qsj4
      qsum5=qsum5+qsj5
      qsum6=qsum6+qsj6
    enddo
  endif
!$acc end kernels
  potnuc1=potnuc*deta
  e1=tsum+vsum+potnuc1
  mpoles(1)=dsum1
  mpoles(2)=dsum2
  mpoles(3)=dsum3
  mpoles(4)=qsum1
  mpoles(5)=qsum2
  mpoles(6)=qsum3
  mpoles(7)=qsum4
  mpoles(8)=qsum5
  mpoles(9)=qsum6

  if(idbg.ge.12) then
    write(lfndbg,120)
    write(lfndbg,130) deta,tsum,vsum,potnuc1,e1,dsum1,dsum2,dsum3, &
        qsum1,qsum2,qsum3,qsum4,qsum5,qsum6
110 format(1x,2i4,5(2x,e20.10))
120 format(///,1x,'The one electron matrix elements of h are :',//)
130 format(15x,'Total overlap              : ',f20.12,//, &
        15x,'Kinetic energy term        : ',f20.12,//, &
        15x,'Nuclear attraction term    : ',f20.12,//, &
        15x,'Nuclear repulsion term     : ',f20.12,//, &
        15x,'One-electron matrix element: ',f20.12,//, &
        15x,'Dipole moment x (el)          : ',f20.12,//, &
        15x,'Dipole moment y (el)           : ',f20.12,//, &
        15x,'Dipole moment z (el)          : ',f20.12,//, &
        15x,'Quadrupole moment xx (el)      : ',f20.12,//, &
        15x,'Quadrupole moment xy (el)      : ',f20.12,//, &
        15x,'Quadrupole moment xz (el)      : ',f20.12,//, &
        15x,'Quadrupole moment yy (el)      : ',f20.12,//, &
        15x,'Quadrupole moment yz (el)      : ',f20.12,//, &
        15x,'Quadrupole moment zz (el)      : ',f20.12,//)
  endif
  return
end subroutine gronor_gnone