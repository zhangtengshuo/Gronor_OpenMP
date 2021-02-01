* example for file 12 (user provided correlation energies for the fragment states)
*10
*-1.50140043 -1.53749914 -1.50507064 -1.47839925 -1.56694218 -1.50140043 -1.53749914 -1.50507064 -1.47839925 -1.56694218
*-----------------------*
* example for file 11, should become automatic (use test file?)
*6
*   1   1   2   3   4   5
*   6   7   6   8  10   9
* Hamiltonian Matrix
*
*                   1                   2                   3                   4                   5                   6
*    1      -767.4274654814        0.9608752252       -0.9608749955        2.4399522310       31.7225150798      -31.7225150527
*    2         0.9608752252     -767.2750500135       -0.0893140823       -0.8169869984      -22.4288632485        3.2446106118
*    3        -0.9608749955       -0.0893140823     -767.2750494438        0.8169870014        3.2446106370      -22.4288632348
*    4         2.4399522310       -0.8169869984        0.8169870014     -767.2410039446      -23.4137870959       23.4137871291
*    5        31.7225150798      -22.4288632485        3.2446106370      -23.4137870959     -767.2088208967        1.1691233374
*    6       -31.7225150527        3.2446106118      -22.4288632348       23.4137871291        1.1691233374     -767.2088213334
*
*
*
* Overlap Matrix
*
*                   1                   2                   3                   4                   5                   6
*    1         1.0000000000       -0.0012515804        0.0012515804       -0.0031776670       -0.0413259852        0.0413259852
*    2        -0.0012515804        1.0000000000        0.0001161449        0.0010641485        0.0292209434       -0.0042278921
*    3         0.0012515804        0.0001161449        1.0000000000       -0.0010641485       -0.0042278921        0.0292209434
*    4        -0.0031776670        0.0010641485       -0.0010641485        1.0000000000        0.0305063133       -0.0305063133
*    5        -0.0413259852        0.0292209434       -0.0042278921        0.0305063133        1.0000000000       -0.0015226440
*    6         0.0413259852       -0.0042278921        0.0292209434       -0.0305063133       -0.0015226440        1.0000000000
*-----------------------*
      program shiftH
      implicit none
      real( kind = 8 ), allocatable :: H(:,:),H_ortho(:,:),H_shift(:,:),
     &                                 H_wshift(:,:),weight(:,:)
      real( kind = 8 ), allocatable :: S(:,:),s_lowdin(:,:),ssave(:,:)
      real( kind = 8 ), allocatable :: shift(:),work(:),ecorr(:)
      real( kind = 8 )              :: wshift
      integer, allocatable          :: ncombv(:,:)
      integer                       :: nbase,i,j,k,info,nstates
      integer, allocatable          :: ipiv(:)
      character (len = 120)         :: title1,title2,title3

      open(11,file='HSdata.txt')     ! This should become a formatted version of f75
      open(12,file='shift.dat')      ! User provided shifts (correlation energies)
      read(11,*) nbase
      if ( nbase .eq. 1 ) then
        write(*,*)'nothing to do with only one MEBF'
        stop
      end if
      read(12,*) nstates
      allocate( H(nbase,nbase) )
      allocate( H_ortho(nbase,nbase) )
      allocate( H_shift(nbase,nbase) )
      allocate( H_wshift(nbase,nbase) )
      allocate( S(nbase,nbase) )
      allocate( ssave(nbase,nbase) )
      allocate( s_lowdin(nbase,nbase) )
      allocate( ncombv(2,nbase) )
      allocate( shift(nbase) )
      allocate( ecorr(nstates) )
      allocate( weight(nbase,nbase) )
      allocate( work(nbase) )
      allocate( ipiv(nbase) )


      read(11,*)(ncombv(1,j),j=1,nbase)
      read(11,*)(ncombv(2,j),j=1,nbase)
      read(11,*)
      read(11,*)
      read(11,*)
      do i = 1, nbase
        read(11,'(i5,1x,10f20.10)') k, (H(i,j),j=1,nbase)
      end do
      do i = 1, 6
        read(11,*)
      end do
      do i = 1, nbase
        read(11,'(i5,1x,10f20.10)') k, (S(i,j),j=1,nbase)
      end do
      read(12,*) (ecorr(i),i=1,nstates)


      do i = 1, nbase
        shift(i) = shift(i) + ecorr(ncombv(1,i)) + ecorr(ncombv(2,i))
      enddo
      write(*,*)
      write(*,*) ' Shifts applied on the diagonal of H'
      write(*,*) ' MEBF    Total shift            Monomer shifts'
      do i = 1, nbase
        write(*,679) i,shift(i),
     &         (ncombv(j,i),ecorr(ncombv(j,i)),j=1,2)
      end do
 679  format(I4,F15.8,8x,5(I4,F15.8))
      write(*,*)


      call lowdin(nbase,S,s_lowdin)
      H_ortho = matmul( transpose(s_lowdin) , matmul(H , s_lowdin) )
      H_shift  = H_ortho
      H_wshift = H_ortho
      do i = 1, nbase
        H_shift(i,i) = H_ortho(i,i) + shift(i)
      end do
      call gnweight(nbase,weight,s_lowdin,S)
      do i = 1, nbase
        wshift = 0.0
        do j = 1, nbase
          wshift =  wshift + weight(i,j) * shift(j)
        end do
        H_wshift(i,i) = H_ortho(i,i) + wshift
      end do


      title1 = 'Unshifted Hamiltonian'
      title2 = 'Shifted Hamiltonian'
      title3 = 'GN-weighted shifted Hamiltonian'

* PRINT the unshifted, shifted and GN-weighted shifted Hamiltonian
      write(*,*)'  *  *  Orthogonal MEBF basis  *  *'
      write(*,*)
      call print_Ham(title2,h_shift,nbase)
      call print_Ham(title3,h_wshift,nbase)

      call dgetrf(nbase,nbase,s_lowdin,nbase,ipiv,info)
      call dgetri(nbase,s_lowdin,nbase,ipiv,work,3*nbase,info)
      h_shift = matmul( transpose(s_lowdin) , matmul(h_shift,s_lowdin) )
      h_wshift = matmul( transpose(s_lowdin),matmul(h_wshift,s_lowdin) )

      write(*,*)'  *  *  Original non-orthogonal MEBF basis  *  *'
      write(*,*)
      call print_Ham(title1,H,nbase)
      call print_Ham(title2,h_shift,nbase)
      call print_Ham(title3,h_wshift,nbase)
      write(*,*)

** Calculate the couplings for the three Hamiltonians

      write(*,*) '  *  *  *  Electronic Couplings  *  *  *'
      write(*,*)
      call couplings(title1,h,s,nbase)
      call couplings(title2,h_shift,s,nbase)
      call couplings(title3,h_wshift,s,nbase)
      write(*,*)


** Print the NOCI wave function

      write(*,*) '  *  *  *  NOCI wave function  *  *  *'
      write(*,*)
      ssave = s
      call nocifunction(title1,h,s,nbase)
      s = ssave
      call nocifunction(title2,h_shift,s,nbase)
      s = ssave
      call nocifunction(title3,h_wshift,s,nbase)
      write(*,*)

      end program shiftH

      subroutine nocifunction(title,h,s,n)
      implicit none
      real ( kind = 8 ), intent(in)  :: h(n,n),s(n,n)
      real ( kind = 8 )              :: eps(n)
      real ( kind = 8 ), allocatable :: work(:)
      integer                        :: n,lwork,info,nk,ncols
      integer                        :: ii,il,i,j,k
      character ( len = 4 )          :: token
      character ( len = 120 )        :: title
      

      info = 0
      lwork = 4*n
      allocate( work(lwork) )
      call dsygv(1,'V','L',n,h,n,s,n,eps,work,lwork,info)
      deallocate(work)
 
      write(*,*) title
      ncols = 10
      nk = n / ncols
      if ( mod(n,ncols) .ne. 0 ) nk = nk + 1
      do k = 1, nk
        ii = (k-1) * ncols + 1
        il = min(n,k*ncols)

        write(*,'(a,12i20)')'    State:',(i,i=ii,il)
        write(*,678)'   Energy:    ', (eps(i),i=ii,il)
        write(*,*) ' '
 678    format((a,12f20.10),//)
        token = "MEBF"
        do j = 1, n
          write(*,679) token,j,(h(j,i),i=ii,il)
 679      format(a4,1x,i5,4x,12f20.10)
          token="    "
        enddo
      enddo
      write(*,*)
      return
      end subroutine nocifunction


      subroutine couplings(title,h,s,n)
      implicit none
      real ( kind = 8 ), intent(in)  :: h(n,n),s(n,n)
      real ( kind = 8 ), allocatable :: tc(:,:)
      integer                        :: n,i,j,k,nk,il,ii,ncols,ik
      character ( len = 120 )        :: title

      allocate( tc(n,n) )

      do i = 2, n
        do j = 1, i-1
          tc(i,j) = ( h(i,j)
     &                -0.5 * ( h(i,i) + h(j,j) ) * s(i,j) ) /
     &             ( 1.0 - s(i,j) * s(i,j) )
        end do
      end do
      write(*,*) title
      ncols = 10
      nk = n / ncols
      if ( mod(n,ncols) .ne. 0 ) nk = nk + 1
      do k = 1, nk
        ii = (k-1) * ncols + 1
        il = min(n,k*ncols)
        write(*,671) (i,i=ii,il)
        do j = 1, n
          if ( j .eq. 1 ) then
            write(*,672) 1
          else
            ik = min(j-1,il)
            write(*,672) j,(tc(j,i),i=ii,ik)
          end if
        end do
        write(*,*)
      end do
      return
 671  format(6x,7(6x,i8,6x))
 672  format(i5,1x,7f20.10,/,(6x,7f20.10))
      end subroutine couplings


      subroutine lowdin(nbase,sbase,slow)
      implicit none

      integer, intent(in)        :: nbase
!      integer                    :: ipiv(nbase)
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
        write(*,*)'Something went wrong in dsyev in lowdin '
        write(*,*) 'info = ',info
        stop
      end if
      deallocate(work)

!     S(diag)=U'SU
      ut=transpose(u)
      slow=matmul(ut,matmul(sbase,u))
!     S(diag)^-1/2
      do j=1,nbase
        slow(j,j)=1.0d0/dsqrt(slow(j,j))
      enddo
!     S^-1/2 = U S(diag)^-1/2 U'
      slow=matmul(u,matmul(slow,ut))

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
      write(*,*) 'Gallup-Norbeck weights:'
      write(*,671) (i,i =1,nbase)
      do i = 1, nbase
        write(*,672)i,wgn(i,:)
      end do
      write(*,*)
      write(*,*)
 671  format(6x,7(6x,i8,6x))
 672  format(i5,1x,7f20.10,/,(6x,7f20.10))
      return
      end subroutine gnweight


      subroutine print_Ham(title,m,n)
      implicit none
      real ( kind = 8 ), intent(in)   :: m(n,n)
      integer                         :: n,i
      character (len = 120)           :: title
      write(*,*) title
      write(*,*)
      write(*,671) (i,i =1,n)
      do i = 1, n
        write(*,672)i,m(i,:)
      end do
      write(*,*)
 671  format(6x,7(6x,i8,6x))
 672  format(i5,1x,7f20.10,/,(6x,7f20.10))
      return
      end subroutine print_Ham
