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
!>    Adding the correlation energy of the fragment wave functions
!>
!> @authors
!>    Kathir R.K., Remco H.A.W. Havenith (University of Groningen), Coen de Graaf (URV, Tarragona)
!>
!> @date
!>    2019
!>
!> @details
!>    After transforming the Hamiltion to an orthogonal basis, the diagonal
!>    matrix elements are shifted by the previously calculated correlation
!>    energy corrections as given by the user in the input (Ecorr). One can
!>    either directly use the shifts as given in the input or take a linear
!>    combination of shifts in line with the composition of the orthogonalized
!>    MEBFs by means of the Gallup-Norbeck weights. After applying the shift
!>    the Hamiltonian is transformed back to the original non-orthogonal basis.
!>
!> @param hbase orginal unshifted Hamiltonian
!> @param sbase overlap matrix of the non-orthogonal MEBFs
!> @param hcorr shifted Hamiltonian in non-orthogonal basis (at the end of the routine)
!> @param wgn Gallup-Norbeck weights
!> @param shft shifts applied on the diagonal of the orthogonalized Hamiltonian
!> @param slow S^-1/2 matrix


      subroutine gronor_correlation_energy(                                     &
     &     nwt,mstates,ecorr,nbase,nmol,ncombv,hbase,sbase,hcorr)
      use cidef           , only :  lfnout
      use gnome_parameters, only :  ipr
      use cidist

      implicit none


      integer                              :: i,j,info
      integer, intent(in)                  :: mstates,nbase,nmol,nwt
      integer, dimension(nbase)            :: ipiv
      integer, dimension(nmol,nbase)       :: ncombv
      real(kind=8), allocatable            :: work(:)
      real(kind=8), dimension(nbase,nbase) :: hbase,sbase,hcorr,slow,wgn
      real(kind=8), dimension(nbase)       :: shft
      real(kind=8), dimension(mstates)     :: ecorr


!     Lowdin orthogonalizing vectors
      call gronor_lowdin(lfnout,ipr,nbase,sbase,slow)

!     Lowdin orthogonalized Hamiltonian H'= (S^-1/2)' H (S^-1/2)
      hcorr=matmul( transpose(slow),matmul(hbase,slow) )

      if (ipr.ge.40) then
        write(lfnout,*)
        write(lfnout,*)
        write(lfnout,*) ' Applying a shift on the diagonal'
        write(lfnout,*)
        write(lfnout,*) ' Hamiltonian in orthogonal MEBF basis'
 671    format(6x,7(6x,i8,6x))
 672    format(i5,1x,7f20.10,/,(6x,7f20.10))
        call gronor_prtmat(lfnout,hcorr,nbase)
 677    format(/,' S^-1/2              ')
        write(lfnout,677)
        call gronor_prtmat(lfnout,slow,nbase)
      endif

      shft=0.0d0
!     Calculate shifts in diagonal,printing needs a bit of an upgrade
!     to make it more clear
      do i=1,nbase
        do j=1,nmol
          shft(i)=shft(i) + ecorr(ncombv(j,i))
        enddo
      enddo
      
      if (ipr.ge.40) then
        write(lfnout,*) ' '
        write(lfnout,*) ' Shifts applied on the diagonal of H'
        write(lfnout,*) ' MEBF    Total shift            Monomer shifts'
        write(lfnout,*) ' '
        do i=1,nbase
          write(lfnout,679) i,shft(i),                                           &
     &         (ncombv(j,i),ecorr(ncombv(j,i)),j=1,nmol)
        end do
 679    format(I4,F15.8,8x,5(I4,F15.8))
        write(lfnout,*)
      endif

      if (nwt.eq.1) then
        call gronor_gnweight(nbase,wgn,slow,sbase)
        if (ipr.ge.40) then
 678      format(/,' GN weights         ')
          write(lfnout,678)
          call gronor_prtmat(lfnout,wgn,nbase)
        endif
      else
        wgn=0.0
        do i=1,nbase
          wgn(i,i)=1.0
        end do
      endif

      do i=1,nbase
        do j=1,nbase
          hcorr(i,i)=hcorr(i,i) + (wgn(i,j)*shft(j))
        enddo
      enddo

      if (ipr.ge.40) then
        write(lfnout,*)
        write(lfnout,*) ' Shifted Hamiltonian in orthogonal MEBF basis'
        call gronor_prtmat(lfnout,hcorr,nbase)
      endif

!     Convert shifted H back to original non-orthogonal basis
      allocate( work(nbase) )
      call dgetrf(nbase,nbase,slow,nbase,ipiv,info)
      call dgetri(nbase,slow,nbase,ipiv,work,3*nbase,info)
      deallocate(work)
      hcorr=matmul(transpose(slow),matmul(hcorr,slow))

      end subroutine gronor_correlation_energy

      subroutine gronor_lowdin(lfnout,ipr,nbase,sbase,slow)

      use cidist

      implicit none

      integer, intent(in)        :: nbase,lfnout,ipr
      integer                    :: info,i,j

      real(kind=8), intent(in)   :: sbase(nbase,nbase)
      real(kind=8), intent(out)  :: slow (nbase,nbase)
      real(kind=8), allocatable  :: work(:)

      real(kind=8)               :: u(nbase,nbase)
      real(kind=8)               :: ut(nbase,nbase)
      real(kind=8)               :: ev(nbase)

!     diagonalize S matrix:
      u=sbase
      allocate(work(4*nbase))
      work=0.0
      info=0
      call dsyev('V','L',nbase,u,nbase,ev,work,4*nbase,info)
      if ( info .ne. 0 ) then
        write(lfnout,*)'Something went wrong in dsyev in gronor_lowdin '
        write(lfnout,*) 'info=',info
        call gronor_abort(330,"Error in Gronor_Lowdin dsyev")
      endif
      deallocate(work)

      if (ipr.ge.40) then
 671    format(6x,7(6x,i8,6x))
 672    format(i5,1x,7f20.10,/,(6x,7f20.10))
        write(lfnout,*) ' '
        write(lfnout,*) 'U - Eigenvectors of S: '
        call gronor_prtmat(lfnout,u,nbase)
      endif

!     S(diag)=U'SU
      ut=transpose(u)
      slow=matmul(ut,matmul(sbase,u))
      if (ipr.ge.40) then
        write(lfnout,*) ' '
        write(lfnout,*) 'U^ S U'
        call gronor_prtmat(lfnout,slow,nbase)
      endif
!     S(diag)^-1/2
      do j=1,nbase
        slow(j,j)=1.0d0/dsqrt(slow(j,j))
      enddo
!     S^-1/2 = U S(diag)^-1/2 U'
      slow=matmul(u,matmul(slow,ut))
      if(ipr.ge.40) then
        write(lfnout,*) ' '
        write(lfnout,*) 'After normalization'
        call gronor_prtmat(lfnout,slow,nbase)
      endif

      end subroutine gronor_lowdin

      subroutine gronor_gnweight(nbase,wgn,slow,sbase)
      implicit none

!     compute the Gallup-Norbeck weights
      integer, intent(in)       :: nbase
      integer                   :: i,k,info
      integer                   :: ipiv(nbase)

      real(kind=8), intent(in)  :: slow(nbase,nbase)
      real(kind=8), intent(in)  :: sbase(nbase,nbase)
      real(kind=8), intent(out) :: wgn(nbase,nbase)
      real(kind=8)              :: wsum
      real(kind=8)              :: sinv(nbase,nbase)
      real(kind=8)              :: csum(nbase)
      real(kind=8)              :: ngn(nbase)
      real(kind=8), allocatable :: work(:)

      wgn=0.0d0
!     invert S matrix
      sinv=sbase
      info=0
      call dgetrf(nbase,nbase,sinv,nbase,ipiv,info)
      allocate(work(3*nbase))
      call dgetri(nbase,sinv,nbase,ipiv,work,3*nbase,info)
      deallocate(work)

      csum=0
      wsum=0
      do k=1,nbase
        do i=1,nbase
          csum(k)=csum(k)+(slow(k,i)*slow(k,i))
        end do
        ngn(k)=1/(csum(k)/sinv(k,k))
        do i=1,nbase
          wgn(k,i)=ngn(k)*slow(k,i)*slow(k,i)/sinv(k,k)
        end do
      end do

      end subroutine gronor_gnweight

