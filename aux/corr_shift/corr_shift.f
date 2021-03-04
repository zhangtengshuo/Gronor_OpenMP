      program shiftH
      use input_data2

      implicit none
      real( kind = 8 ), allocatable :: H_ortho(:,:),H_shift(:,:),
     &                                 H_wshift(:,:),weight(:,:)
      real( kind = 8 ), allocatable :: s_lowdin(:,:)
      real (kind = 8 ), allocatable :: ssave(:,:),hsave(:,:)
      real( kind = 8 ), allocatable :: shift(:),work(:)
      real( kind = 8 )              :: wshift
      integer, allocatable          :: iRoots(:,:),nRoots(:)
      integer                       :: i,j,k,info
      integer, allocatable          :: ipiv(:)
      character (len = 120)         :: title1,title2,title3
      logical                       :: dump_in_out

      call readin

      allocate( H_ortho(nbase,nbase) )
      allocate( H_shift(nbase,nbase) )
      allocate( H_wshift(nbase,nbase) )
      allocate( hsave(nbase,nbase) )
      allocate( ssave(nbase,nbase) )
      allocate( s_lowdin(nbase,nbase) )
      allocate( shift(nbase) )
      allocate( weight(nbase,nbase) )
      allocate( work(nbase) )
      allocate( ipiv(nbase) )


      title1 = 'MEBF Overlap matrix'
      call print_Ham(title1,s,nbase)

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

      dump_in_out = .true.
      write(*,*) '  *  *  *  NOCI wave function  *  *  *'
      write(*,*)
      ssave = s
      hsave = h
      call nocifunction(title1,h,s,nbase,dump_in_out)
      h = hsave
      s = ssave
      hsave = h_shift
      call nocifunction(title2,h_shift,s,nbase,dump_in_out)
      h_shift = hsave
      s = ssave
      hsave = h_wshift
      call nocifunction(title3,h_wshift,s,nbase,dump_in_out)
      h_wshift = hsave
      s = ssave
      if ( extra ) then
        write(*,*)
        write(*,*)
        write(*,*)
        write(*,*)'Extra NOCI (unshifted)'
        write(*,*)'----------------------'
        call superNOCI
        write(*,*)
        write(*,*)'Extra NOCI (shifted)'
        write(*,*)'--------------------'
        h = h_shift
        call superNOCI
        write(*,*)
        write(*,*)'Extra NOCI (GN-weight shifted)'
        write(*,*)'------------------------------'
        h = h_wshift
        call superNOCI
      end if

      end program shiftH

      subroutine nocifunction(title,h,s,n,dump_in_out)
      implicit none
      real ( kind = 8 ), intent(in)  :: h(n,n),s(n,n)
      real ( kind = 8 )              :: eps(n)
      real ( kind = 8 ), allocatable :: work(:)
      integer                        :: n,lwork,info,nk,ncols
      integer                        :: ii,il,i,j,k
      character ( len = 4 )          :: token
      character ( len = 120 )        :: title
      logical                        :: dump_in_out
      

      info = 0
      lwork = 4*n
      allocate( work(lwork) )
      call dsygv(1,'V','L',n,h,n,s,n,eps,work,lwork,info)
      deallocate(work)
 
      if ( dump_in_out ) then
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
 678      format((a,12f20.10),//)
          token = "MEBF"
          do j = 1, n
            write(*,679) token,j,(h(j,i),i=ii,il)
 679        format(a4,1x,i5,4x,12f20.10)
            token="    "
          enddo
        enddo
        write(*,*)
      end if
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
      write(*,'(a120)') title
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


      subroutine superNOCI
      use input_data2
      implicit none

      integer                      :: i,j,k,l,info,lwork,jj
      real (kind=8), allocatable   :: h_block(:,:)
      real (kind=8), allocatable   :: s_block(:,:)
      real (kind=8), allocatable   :: extraVec(:,:)
      real (kind=8), allocatable   :: extraH(:,:),extraS(:,:)
      real (kind=8), allocatable   :: eps(:)
      real (kind=8), allocatable   :: work(:)
      real (kind=8), allocatable   :: mVec(:,:),reduced_extraVec(:,:)
      real (kind=8)                :: ovlp,maxovlp
      character (len=120)          :: title

      allocate( extraVec(extraDim,nbase) )
      allocate( extraH(reduced_extraDim,reduced_extraDim) )
      allocate( extraS(reduced_extraDim,reduced_extraDim) )
      extraVec = 0.0
      extraH   = 0.0
      extraS   = 0.0

      jj = 0
      do i = 1 ,nBlocks
        lwork = 4 * dimens(i)
        allocate ( h_block(dimens(i),dimens(i)) )
        allocate ( s_block(dimens(i),dimens(i)) )
        allocate ( work(lwork) )
        allocate ( eps(dimens(i)) )
        h_block = 0.0
        s_block = 0.0
        work = 0.0
        eps = 0.0
        do j = 1, dimens(i)
          do k = 1, j
            h_block(j,k) = h(blocks(i,j),blocks(i,k))
            s_block(j,k) = s(blocks(i,j),blocks(i,k))
            if ( k .ne. j ) then
              h_block(k,j) = h_block(j,k)
              s_block(k,j) = s_block(j,k)
            end if
          end do
        end do
        info = 0
        call dsygv(1,'V','L',dimens(i),h_block,dimens(i),s_block,
     &                         dimens(i),eps,work,lwork,info)
        write(*,'(a7,i1,a)')' Block ',i,' : Energies and eigenvectors'
        write(*,'(12(2x,F14.8))')eps
        do j = 1, dimens(i)
          write(*,'(12F16.8)')h_block(j,:)
          jj = jj + 1
          do k = 1, dimens(i)
            extraVec(jj,blocks(i,k)) = h_block(k,j)
          end do
        end do
        write(*,*)
        deallocate ( work )
        deallocate ( eps )
        deallocate ( h_block )
        deallocate ( s_block )
      end do
      if (dressed_coupling ) then
        allocate( h_block(reduced_extraDim,reduced_extraDim) )
        allocate( s_block(reduced_extraDim,reduced_extraDim) )
        do i = 1, reduced_extraDim
          do j = 1, reduced_extraDim
            h_block(i,j) = h(ovlp_mebfs(i),ovlp_mebfs(j))
            s_block(i,j) = s(ovlp_mebfs(i),ovlp_mebfs(j))
          end do
        end do
        lwork = 4 * reduced_extraDim
        info = 0
        allocate( work(lwork) )
        allocate( eps(reduced_extraDim) )
        eps = 0.0
        work = 0.0
        call dsygv(1,'V','L',reduced_extraDim,h_block,reduced_extraDim,
     &              s_block,reduced_extraDim,eps,work,lwork,info)
        deallocate( work, eps )
        allocate( mVec(reduced_extraDim,nbase) )
        allocate( reduced_extraVec(reduced_extraDim,nbase) )
        mVec = 0.0
        do i = 1, reduced_extraDim
          do j = 1, nbase
            do k = 1, reduced_extraDim
              if ( j .eq. ovlp_mebfs(k) ) mVec(i,j) =
     &                                     h_block(k,i)
            end do
          end do
        end do
        do i = 1, reduced_extraDim
          maxovlp = 0.0
          do j = 1, extraDim
            ovlp = 0.0
            do k = 1, nbase
              do l = 1, nbase
                ovlp = ovlp + mVec(i,k) * extraVec(j,l) * S(k,l)
              end do
            end do
            if (abs(ovlp) .gt. maxovlp ) then
              maxovlp = abs(ovlp)
              reduced_extraVec(i,:) = extraVec(j,:)
            end if
          end do
        end do
        write(*,*) 'new MEBFs expressed in orginal MEBFs'
        do i = 1, reduced_extraDim
          write(*,'(I3,a,12F14.8)')i,':',reduced_extraVec(i,:)
        end do
        write(*,*)
        do i = 1, reduced_extraDim
          do j = 1, reduced_extraDim
            do k = 1, nbase
              do l = 1, nbase
                extraH(i,j) = extraH(i,j) + reduced_extraVec(i,k) 
     &                                * reduced_extraVec(j,l) * H(k,l)
                extraS(i,j) = extraS(i,j) + reduced_extraVec(i,k) 
     &                                * reduced_extraVec(j,l) * S(k,l)
              end do
            end do
          end do
        end do
        title = ' Hamiltonian in new MEBF basis'
        call print_Ham(title,extraH,reduced_extraDim)
        title = ' Overlap matrix new MEBFs'
        call print_Ham(title,extraS,reduced_extraDim)
        write(*,*)
        write(*,*) '  *  *  *  Electronic Couplings  *  *  *'
        title = ' '
        call couplings(title,extraH,extraS,reduced_extraDim)
        deallocate(h_block,s_block)
        deallocate(mVec,reduced_extraVec)
        deallocate(extraH,extraS)
      end if

      end subroutine superNOCI



      subroutine readin
      use input_data2
      implicit none

      integer,parameter    :: nKeys = 7
      integer              :: jj,iKey,i,j

      character(len=4), dimension(nKeys)   :: keyword
      character(len=4)                     :: key
      character(len=132)                   :: line,filename

      logical                              :: all_ok = .true.
      logical, dimension(nKeys)            :: hit = .false.


      data keyword /'MEBF','PROJ','SHIF','HAMI','OVER','BLOC','SELE'/

      extra = .false.
      dressed_coupling = .false.

      do while (all_ok)
        read(5,*,iostat=jj) line
        key = adjustl(line)
        call capitalize(key)
        do iKey = 1, nKeys
          if ( key .eq. keyword(iKey) ) hit(iKey) = .true.
        end do
        if (  jj .lt. 0 ) all_ok = .false.
      end do

      do iKey = 1, nKeys
        if ( hit(iKey) ) then
          select case(iKey)
            case(1)
              call locate(5,'MEBF')
              read(*,*) nbase
              if ( nbase .eq. 1 ) then
                write(*,*)'nothing to be done for only one MEBF'
                stop
              end if
              allocate( ncombv(2,nbase) )
              read(*,*)(ncombv(1,j),j=1,nbase)
              read(*,*)(ncombv(2,j),j=1,nbase)
              nstates = max(maxval(ncombv(1,:)),maxval(ncombv(2,:)))
            case(2)
              call locate(5,'PROJ')
              read(*,*) project
              project = trim(project)
              hit(4) = .false.                            ! deactivate reading H and S from input
              hit(5) = .false.
              write(filename,'(2a)')trim(Project),'.arx'
              write(*,'(2a)') 'H and S are read from ',filename
              open(11,file=filename,err=9000)
              allocate( h(nbase,nbase) )
              allocate( s(nbase,nbase) )
              call locate(11,'Hami')
              read(11,*) line
              do i = 1, nbase
                read(11,662)jj,(h(i,j),j=1,nbase)
              end do
 662          format(i5,1x,10e20.13)
              call locate(11,'Over')
              read(11,*) line
              do i = 1, nbase
                read(11,662)jj,(s(i,j),j=1,nbase)
              end do
            case(3)
              call locate(5,'SHIF')
              allocate( ecorr(nstates) )
              read(*,*)(ecorr(i),i=1,nstates)
            case(4)
              call locate(5,'HAMI')
              allocate( h(nbase,nbase) )
              do i = 1, nbase
                read(*,*)(h(i,j),j=1,nbase)
              end do
            case(5)
              call locate(5,'OVER')
              allocate( s(nbase,nbase) )
              do i = 1, nbase
                read(*,*)(s(i,j),j=1,nbase)
              end do
            case(6)
              call locate(5,'BLOC')
              read(*,*) nBlocks
              allocate( dimens(nBlocks) )
              allocate( blocks(nBlocks,nbase) )
              blocks = 0
              dimens = 0
              extraDim = 0
              do i = 1, nBlocks
                read(*,*) dimens(i),(blocks(i,j),j=1,dimens(i))
                extraDim = extraDim + dimens(i)
              end do
              extra = .true.
            case(7)
              call locate(5,'SELE')
              read(*,*) reduced_extraDim
              allocate( ovlp_mebfs(reduced_extraDim) )
              ovlp_mebfs = 0
              read(*,*) (ovlp_mebfs(i),i=1,reduced_extraDim)
              dressed_coupling = .true.
          end select
        end if
      end do     

      return

 9000 write(*,*)'A problem ocurred opening ',filename
      stop

      end subroutine readin
!MEBF                          :   Definition of the MEBFs (only dimers for now)
!  6                           :   number of MEBFs
!    1   1   2   3   4   5     :   monomer functions of fragmemt A
!    6   7   6   8  10   9     :   monomer functions of fragment B
!SHIFt                         :   correlation energies of the monomer functions
!  -1.50140043 -1.53749914 -1.50507064 -1.47839925 -1.56694218 -1.50140043 -1.53749914 -1.50507064 -1.47839925 -1.56694218
!HAMIltonian
!      -767.4274654814        0.9608752252       -0.9608749955        2.4399522310       31.7225150798      -31.7225150527
!         0.9608752252     -767.2750500135       -0.0893140823       -0.8169869984      -22.4288632485        3.2446106118
!        -0.9608749955       -0.0893140823     -767.2750494438        0.8169870014        3.2446106370      -22.4288632348
!         2.4399522310       -0.8169869984        0.8169870014     -767.2410039446      -23.4137870959       23.4137871291
!        31.7225150798      -22.4288632485        3.2446106370      -23.4137870959     -767.2088208967        1.1691233374
!       -31.7225150527        3.2446106118      -22.4288632348       23.4137871291        1.1691233374     -767.2088213334
!OVERlap
!         1.0000000000       -0.0012515804        0.0012515804       -0.0031776670       -0.0413259852        0.0413259852
!        -0.0012515804        1.0000000000        0.0001161449        0.0010641485        0.0292209434       -0.0042278921
!         0.0012515804        0.0001161449        1.0000000000       -0.0010641485       -0.0042278921        0.0292209434
!        -0.0031776670        0.0010641485       -0.0010641485        1.0000000000        0.0305063133       -0.0305063133
!        -0.0413259852        0.0292209434       -0.0042278921        0.0305063133        1.0000000000       -0.0015226440
!         0.0413259852       -0.0042278921        0.0292209434       -0.0305063133       -0.0015226440        1.0000000000
!BLOCk diagonalizations:  two extra groups: S0S1; S1S0; D+D-; D-D+ & T1T1; D+D-; !D-D+
! 2
! 4  2 3 5 6 
! 3  4 5 6
!SELEct  : indicate how many new MEBFs should be selected for the calculation of the coupling, selection is based on overlap with original MEBFs.
! 3  
!  2 3 4
!
! Hamiltonian and Overlap should be read from binary file in the near future


      subroutine capitalize(string)
      implicit none
      integer      :: i
      character(*) string

      do i = 1, len(string)
        if (ichar(string(i:i)).gt.96) then
          string(i:i) = char(ichar(string(i:i))-32)
        endif
      end do
      return
      end subroutine capitalize


      subroutine locate(lu,string)
      implicit none
      integer        ::  lu
      character(4)   ::  string,string2
      character(132) ::  line
      rewind(lu)
 40   read(lu,*) line
      string2=adjustl(line)
      if (lu.eq.5) call capitalize(string2)
      if (string2.ne.string) goto 40
      return
      end subroutine locate

