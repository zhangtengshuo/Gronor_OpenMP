      subroutine gronor_correlation_energy(
     &     nwt,mstates,ecorr,nbase,nmol,ncombv,hbase,sbase,hcorr)
      use cidef           , only :  lfnout
      use gnome_parameters, only :  ipr

      implicit none


      integer                              :: i,j,k,l,m,n,info
      integer, intent(in)                  :: mstates,nbase,nmol,nwt
      integer, dimension(nbase)            :: ipiv
      integer, dimension(nmol,nbase)       :: ncombv
      real(kind=8), allocatable            :: work(:)
      real(kind=8), dimension(nbase,nbase) :: hbase,sbase,hcorr,slow,wgn
      real(kind=8), dimension(nbase)       :: hev
      real(kind=8), dimension(nbase)       :: shft
      real(kind=8), dimension(mstates)     :: ecorr


!     Lowdin orthogonalizing vectors
      call lowdin(lfnout,ipr,nbase,sbase,slow)

!     Lowdin orthogonalized Hamiltonian H'= (S^-1/2)' H (S^-1/2)
      hcorr = matmul( transpose(slow),matmul(hbase,slow) )

      write(lfnout,*)
      write(lfnout,*)
      write(lfnout,*) ' Applying a shift on the diagonal'
      write(lfnout,*)
      write(lfnout,*) ' Hamiltonian in orthogonal MEBF basis'
      write(lfnout,671) (i,i =1,nbase)
 671  format(6x,7(6x,i8,6x))
 672  format(i5,1x,7f20.10,/,(6x,7f20.10))
      do j = 1, nbase
       write(lfnout,672) j,(hcorr(i,j),i=1,nbase)
      enddo
      flush(lfnout)

      if ( ipr .ge. 5 ) then
 677    format(//,' S^-1/2              ',/)
        write(lfnout,677)
        write(lfnout,671) (i,i=1,nbase)
        do j = 1, nbase
          write(lfnout,672) j,(slow(i,j),i=1,nbase)
        enddo
      end if

      shft=0.0d0
!     Calculate shifts in diagonal, printing needs a bit of an upgrade
!     to make it more clear
      do i = 1, nbase
        do j = 1, nmol
          shft(i) = shft(i) + ecorr(ncombv(j,i))
        enddo
      enddo
      write(lfnout,*)
      write(lfnout,*) ' Shifts applied on the diagonal of H'
      write(lfnout,*) ' MEBF    Total shift            Monomer shifts'
      do i = 1, nbase
        write(lfnout,679) i,shft(i),
     &         (ncombv(j,i),ecorr(ncombv(j,i)),j=1,nmol)
      end do
 679  format(I4,F15.8,8x,5(I4,F15.8))
      write(lfnout,*)

      if ( nwt .eq. 1 ) then
        call gnweight(nbase,wgn,slow,sbase)
 678    format(/,' GN weights         ')
        write(lfnout,678)
        write(lfnout,671) (i,i=1,nbase)
        do j = 1, nbase
          write(lfnout,672) j,(wgn(i,j),i=1,nbase)
        enddo
      else
        wgn = 0.0
        do i = 1, nbase
          wgn(i,i) = 1.0
        end do
      end if

      do i = 1, nbase
        do j = 1, nbase
          hcorr(i,i) = hcorr(i,i) + (wgn(i,j)*shft(j))
        enddo
      enddo
      
      write(lfnout,*)
      write(lfnout,*) ' Shifted Hamiltonian in orthogonal MEBF basis'
      write(lfnout,671) (i,i =1,nbase)
      do i = 1, nbase
        write(lfnout,672)i,hcorr(i,:)
      end do      

!     Convert shifted H back to original non-orthogonal basis
      allocate( work(nbase) )
      call dgetrf(nbase,nbase,slow,nbase,ipiv,info)
      call dgetri(nbase,slow,nbase,ipiv,work,3*nbase,info)
      deallocate(work)
      hcorr=matmul(transpose(slow),matmul(hcorr,slow))

      end subroutine gronor_correlation_energy

      subroutine lowdin(lfnout,ipr,nbase,sbase,slow)
      implicit none

      integer, intent(in)        :: nbase,lfnout,ipr
      integer                    :: ipiv(nbase)
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
      work = 0.0
      info = 0
      call dsyev('V','L',nbase,u,nbase,ev,work,4*nbase,info)
      if ( info .ne. 0 ) then
        write(lfnout,*)'Something went wrong in dsyev in lowdin '
        write(lfnout,*) 'info = ',info
        stop
      end if
      deallocate(work)

      if ( ipr .ge. 5 ) then
 671    format(6x,7(6x,i8,6x))
 672    format(i5,1x,7f20.10,/,(6x,7f20.10))
        write(lfnout,*) 'U - Eigenvectors of S: '
        write(lfnout,671) (i,i=1,nbase)
        flush(lfnout)
        do j=1,nbase
          write(lfnout,672) j,(u(i,j),i=1,nbase)
        enddo
        flush(lfnout)
      end if

!     S(diag)=U'SU
      ut=transpose(u)
      slow=matmul(ut,matmul(sbase,u))
      if ( ipr .ge. 5 ) then
        write(lfnout,*) 'U^ S U'
        do j = 1, nbase
          write(lfnout,'(15F25.15)')slow(j,:)
        end do
      end if
!     S(diag)^-1/2
      do j=1,nbase
        slow(j,j)=1.0d0/dsqrt(slow(j,j))
      enddo
!     S^-1/2 = U S(diag)^-1/2 U'
      slow=matmul(u,matmul(slow,ut))
      if ( ipr .ge. 5 ) then
        write(lfnout,*) 'After normalization'
        do j = 1, nbase
          write(lfnout,'(15F25.15)')slow(j,:)
        end do
      end if

      end subroutine lowdin

      subroutine gnweight(nbase,wgn,slow,sbase)
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
      info = 0
      call dgetrf(nbase,nbase,sinv,nbase,ipiv,info)
      allocate(work(3*nbase))
      call dgetri(nbase,sinv,nbase,ipiv,work,3*nbase,info)
      deallocate(work)

      csum=0
      wsum=0
      do k = 1, nbase
        do i = 1, nbase
          csum(k) = csum(k) + (slow(k,i)*slow(k,i))
        end do
        ngn(k) = 1 / (csum(k)/sinv(k,k))
        do i = 1, nbase
          wgn(k,i) = ngn(k) * slow(k,i) * slow(k,i) / sinv(k,k)
        end do
      end do

      end subroutine gnweight
