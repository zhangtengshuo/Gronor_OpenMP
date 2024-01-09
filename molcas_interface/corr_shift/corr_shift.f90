module corr_shift_input_data
  implicit none
  integer                       :: nbase,nstates,extraDim,nBlocks,nmol
  integer                       :: reduced_extraDim,maxLength,nUniqueExtra
  integer, allocatable          :: ncombv(:,:),dimens(:),blocks(:,:),nequiv(:)
  integer, allocatable          :: ovlp_mebfs(:),nDressed_MEBFs(:),allExtra_mebfs(:)
  real (kind=8), allocatable    :: h(:,:),s(:,:)
  real (kind=8), allocatable    :: shift(:),ecorr(:,:),ecorr2(:,:)
  character(len=80)             :: project
  character(len=6),allocatable  :: fragLabel(:,:)
  character(len=18),allocatable :: mebfLabel(:)
  logical :: extra,dressed_coupling,GN_weights,select_lowest,print_overlap,   &
             dcec,averaging
end module corr_shift_input_data


program  corr_shift_main
use corr_shift_input_data
implicit none

external  :: corr_shift_readin,corr_shift_printHam
external  :: corr_shift_lowdin,corr_shift_gnweight
external  :: dgetrf,dgetri,corr_shift_nocifunction
external  :: corr_shift_couplings,corr_shift_superNOCI

real(kind = 8), allocatable :: H_ortho(:,:),H_shift(:,:)
real(kind = 8), allocatable :: H_wshift(:,:),weight(:,:)
real(kind = 8), allocatable :: s_lowdin(:,:)
real(kind = 8), allocatable :: ssave(:,:),hsave(:,:)
real(kind = 8), allocatable :: work(:)
real(kind = 8)              :: wshift
integer                     :: i,j,info
integer, allocatable        :: ipiv(:)
character (len = 120)       :: title,title1,title2,title3
logical                     :: dump_in_out,useLabels

call corr_shift_readin

allocate(H_ortho(nbase,nbase))
allocate(H_shift(nbase,nbase))
allocate(H_wshift(nbase,nbase))
allocate(hsave(nbase,nbase))
allocate(ssave(nbase,nbase))
allocate(s_lowdin(nbase,nbase))
allocate(weight(nbase,nbase))
allocate(work(nbase))
allocate(ipiv(nbase))

write(*,*)
if (dcec) then
  write(*,*) ' Shifts applied on the diagonal of H'
  write(*,*) ' MEBF    Total shift            Monomer shifts'
  do j = 1, nbase
    write(*,679) trim(mebfLabel(j)),shift(j),(fragLabel(i,j),ecorr2(i,j),i=1,nmol)
  end do
  679  format(A4,F15.8,8x,5(A6,F15.8,3x))
  write(*,*)
endif

ssave = s
call corr_shift_lowdin(nbase,s,s_lowdin)
s = ssave
H_ortho = matmul( transpose(s_lowdin) , matmul(H , s_lowdin) )
H_shift  = H_ortho
H_wshift = H_ortho
if (dcec) then
  do i = 1, nbase
    H_shift(i,i) = H_ortho(i,i) + shift(i)
  end do
  if ( GN_weights ) then
    useLabels = .true.
    call corr_shift_gnweight(nbase,weight,s_lowdin,S,useLabels)
    do i = 1, nbase
      wshift = 0.0
      do j = 1, nbase
        wshift =  wshift + weight(i,j) * shift(j)
      end do
      H_wshift(i,i) = H_ortho(i,i) + wshift
    end do
  endif
endif

title1 = 'Unshifted Hamiltonian'
title2 = 'Shifted Hamiltonian'
title3 = 'GN-weighted shifted Hamiltonian'

!* Print the unshifted, shifted and GN-weighted shifted Hamiltonian
if (dcec) then
  write(*,*)'  *  *  Orthogonal MEBF basis  *  *'
  write(*,*)
  useLabels = .true.
  call corr_shift_printHam(title2,h_shift,nbase,useLabels)
  if (GN_weights) call corr_shift_printHam(title3,h_wshift,nbase,useLabels)
  
  call dgetrf(nbase,nbase,s_lowdin,nbase,ipiv,info)
  call dgetri(nbase,s_lowdin,nbase,ipiv,work,3*nbase,info)
  h_shift = matmul( transpose(s_lowdin) , matmul(h_shift,s_lowdin) )
  if ( GN_weights ) then
    h_wshift = matmul(transpose(s_lowdin),matmul(h_wshift,s_lowdin))
  endif
endif

write(*,*)'  *  *  Original non-orthogonal MEBF basis  *  *'
write(*,*)
useLabels = .true.
call corr_shift_printHam(title1,H,nbase,useLabels)
if (dcec) call corr_shift_printHam(title2,h_shift,nbase,useLabels)
if (GN_weights) call corr_shift_printHam(title3,h_wshift,nbase,useLabels)
write(*,*)

if (print_overlap)then
  title = 'Original MEBF overlap matrix'
  call corr_shift_printHam(title,S,nbase,useLabels)
endif

!** Calculate the couplings for the three Hamiltonians

write(*,*) '  *  *  *  Electronic Couplings  (meV) *  *  *'
write(*,*)
useLabels = .true.
call corr_shift_couplings(title1,h,s,nbase,useLabels)
if (dcec) call corr_shift_couplings(title2,h_shift,s,nbase,useLabels)
if (GN_weights) call corr_shift_couplings(title3,h_wshift,s,nbase,useLabels)
write(*,*)

!** Print the NOCI wave function

dump_in_out = .true.
write(*,*) '  *  *  *  NOCI wave function  *  *  *'
write(*,*)
ssave = s
hsave = h
call corr_shift_nocifunction(title1,h,s,nbase,dump_in_out)
h = hsave
s = ssave
if (dcec) then
  hsave = h_shift
  call corr_shift_nocifunction(title2,h_shift,s,nbase,dump_in_out)
  h_shift = hsave
  s = ssave
  if ( GN_weights ) then
    hsave = h_wshift
    call corr_shift_nocifunction(title3,h_wshift,s,nbase,dump_in_out)
    h_wshift = hsave
    s = ssave
  endif
endif
if ( extra ) then
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)'Extra NOCI (unshifted)'
  write(*,*)'----------------------'
  call corr_shift_superNOCI
  if (dcec) then
    write(*,*)
    write(*,*)'Extra NOCI (shifted)'
    write(*,*)'--------------------'
    h = h_shift
    call corr_shift_superNOCI
    if ( GN_weights ) then
      write(*,*)
      write(*,*)'Extra NOCI (GN-weight shifted)'
      write(*,*)'------------------------------'
      h = h_wshift
      call corr_shift_superNOCI
    endif
  endif
endif

end program corr_shift_main

subroutine corr_shift_nocifunction(title,h,s,n,dump_in_out)
use  corr_shift_input_data, only : mebfLabel
implicit none

external :: dsygv,corr_shift_gnweight

real ( kind = 8 ), intent(in)  :: h(n,n),s(n,n)
real ( kind = 8 )              :: eps(n)
real ( kind = 8 ), allocatable :: work(:),wGN(:,:),ssave(:,:)
integer                        :: n,lwork,info,nk,ncols
integer                        :: ii,il,i,j,k
character ( len = 4 )          :: token
character ( len = 120 )        :: title
logical                        :: dump_in_out,useLabels


allocate(wGN(n,n),ssave(n,n))
info = 0
lwork = 4*n
allocate( work(lwork) )
ssave = s
call dsygv(1,'V','L',n,h,n,ssave,n,eps,work,lwork,info)
deallocate(work)

if ( dump_in_out ) then
  write(*,*) title
  ncols = 10
  nk = n / ncols
  if ( mod(n,ncols) .ne. 0 ) nk = nk + 1
  do k = 1, nk
    ii = (k-1) * ncols + 1
    il = min(n,k*ncols)

    write(*,'(9x,a,12i20)')'    State:',(i,i=ii,il)
    write(*,678)'   Energy:    ', (eps(i),i=ii,il)
    write(*,*) ' '
    token = "MEBF"
    do j = 1, n
      write(*,679) token,adjustr(mebfLabel(j)),(h(j,i),i=ii,il)
      token="    "
    enddo
    write (*,*)
  enddo
  write(*,*)
end if
678 format((9x,a,12f20.10),//)
679 format(a4,1x,a,12f20.10)

useLabels = .true.
call corr_shift_gnweight(n,wGN,h,s,useLabels)
deallocate(wGN)

return
end subroutine corr_shift_nocifunction


subroutine corr_shift_couplings(title,h,s,n,useLabels)
use  corr_shift_input_data, only : mebfLabel,maxLength
implicit none
real ( kind = 8 ), intent(in)  :: h(n,n),s(n,n)
real ( kind = 8 ), allocatable :: tc(:,:)
integer                        :: n,i,j,k,first,last
character ( len = 120 )        :: title
logical                        :: useLabels

if (maxLength .gt. 18) useLabels = .false.

allocate( tc(n,n) )

do i = 2, n
  do j = 1, i-1
    tc(i,j) = (h(i,j) - 0.5*(h(i,i) + h(j,j))*s(i,j)) / (1.0 - s(i,j)*s(i,j))
  end do
end do
write(*,*) title
first = 1
if (n.eq.2) write(*,'(a,f20.10)') 't(1,2) = ',27211.4d0*tc(2,1)
last = min(7,n)
do while (first .lt. n-1)
  if ( useLabels ) then
    write(*,671) (adjustl(mebfLabel(j)),j=first,min(n-1,last))
    do k = first+1, n
      write(*,672) adjustr(mebfLabel(k)),(tc(k,j)*27211.4d0,j=first,min(k-1,last))
    end do
  else
    write(*,673) (j,j=first,min(n-1,last))
    do k = first+1, n
      write(*,674) k,(tc(k,j)*27211.4d0,j=first,min(k-1,last))
    end do
  end if
  write(*,*)
  first = last + 1
  last = min(last+7,n)
end do
return
671 format(22x,7(a20))
672 format(a20,1x,7f20.10)
673  format(6x,7(6x,i8,6x))
674  format(i5,1x,7f20.10)
end subroutine corr_shift_couplings


subroutine corr_shift_lowdin(nbase,sbase,slow)
implicit none

external :: dsyev

integer, intent(in)        :: nbase
integer                    :: info,j

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

! S(diag)=U'SU
ut=transpose(u)
slow=matmul(ut,matmul(sbase,u))

! S(diag)^-1/2
do j=1,nbase
  slow(j,j)=1.0d0/dsqrt(slow(j,j))
enddo

! S^-1/2 = U S(diag)^-1/2 U'
slow=matmul(u,matmul(slow,ut))


end subroutine corr_shift_lowdin




subroutine corr_shift_gnweight(nbase,wgn,slow,sbase,useLabels)
use  corr_shift_input_data, only : mebfLabel,maxLength
implicit none

external  :: dgetrf,dgetri

!  compute the Gallup-Norbeck weights
integer, intent(in)       :: nbase
integer                   :: i,k,info,first,last
integer                   :: ipiv(nbase)

!real(kind=8), intent(in)  :: slow(nbase,nbase)
real(kind=8)              :: slow(nbase,nbase)
real(kind=8), intent(in)  :: sbase(nbase,nbase)
real(kind=8), intent(out) :: wgn(nbase,nbase)
real(kind=8)              :: wsum
real(kind=8)              :: sinv(nbase,nbase)
real(kind=8)              :: csum(nbase)
real(kind=8)              :: ngn(nbase)
real(kind=8), allocatable :: work(:)

logical                   :: useLabels

wgn=0.0d0
!  invert S matrix
sinv=sbase
info = 0
call dgetrf(nbase,nbase,sinv,nbase,ipiv,info)
allocate(work(3*nbase))
call dgetri(nbase,sinv,nbase,ipiv,work,3*nbase,info)
deallocate(work)

slow = transpose(slow)

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
write(*,*)
write(*,*) 'Gallup-Norbeck weights:'
if ( maxLength .gt. 18 ) useLabels = .false.
first = 1
last = min(10,nbase)
do while (first .le. nbase)
  if (useLabels ) then
    write(*,671) (i,i = first,last)
  else
    write(*,672) (i,i = first,last)
  end if
  do i = 1, nbase
    if ( useLabels) then
      write(*,673)adjustr(mebfLabel(i)),(wgn(k,i),k=first,last)
    else
      write(*,674)i,(wgn(k,i),k=first,last)
    end if
  end do
  write(*,*)
  first = last + 1
  last = min(last + 10,nbase)
end do
write(*,*)
write(*,*)
671 format(22x,10(i8,6x))
672 format(9x,10(i8,6x))
673 format(a,1x,10f14.6)
674 format(i5,1x,10f14.6)
return
end subroutine corr_shift_gnweight


subroutine corr_shift_printHam(title,m,n,useLabels)
use  corr_shift_input_data, only : mebfLabel,maxlength
implicit none
real ( kind = 8 ), intent(in)   :: m(n,n)
integer                         :: n,i,j,first,last
character (len = 120)           :: title
logical                         :: useLabels
write(*,'(a120)') title
write(*,*)
if ( maxLength .gt. 18 ) useLabels = .false.
first = 1 ; last = min(n,7)
do while ( first .le. n )
  if (useLabels) then
    write(*,671) (adjustl(mebfLabel(i)),i=first,last)
    do i = 1, n
      write(*,672)adjustr(mebfLabel(i)),(m(i,j),j=first,last)
    end do
  else
    write(*,673) (i,i=first,last)
    do i = 1, n
      write(*,674)i,(m(i,j),j=first,last)
    end do
  end if
  first = last + 1
  last = min(last + 7,n)
  write(*,*)
end do
671 format(22x,7(a20))
672 format(a20,1x,7f20.10)
673 format(6x,7(6x,i8,6x))
674 format(i5,1x,7f20.10)
return
end subroutine corr_shift_printHam


subroutine corr_shift_superNOCI
use corr_shift_input_data
implicit none

external :: dsygv,corr_shift_printHam,corr_shift_couplings,corr_shift_gnweight

integer                      :: i,j,k,l,m,info,lwork,jj,offset,first,last
integer, allocatable         :: ncount(:)
real (kind=8), allocatable   :: h_block(:,:)
real (kind=8), allocatable   :: s_block(:,:),s_save(:,:)
real (kind=8), allocatable   :: extraVec(:,:)
real (kind=8), allocatable   :: extraH(:,:),extraS(:,:)
real (kind=8), allocatable   :: eps(:),h_unique(:),s_unique(:)
real (kind=8), allocatable   :: work(:),wGN(:,:)
real (kind=8), allocatable   :: mVec(:,:),reduced_extraVec(:,:)
real (kind=8)                :: ovlp,maxovlp
character (len=120)          :: title
logical, allocatable         :: used(:),useLabels

allocate( extraVec(extraDim,nbase) )
allocate( extraH(reduced_extraDim,reduced_extraDim) )
allocate( extraS(reduced_extraDim,reduced_extraDim) )
extraVec = 0.0
extraH   = 0.0
extraS   = 0.0


if ( averaging ) then
  write(*,*)'* *  * Warning: H and S are symmetrized upon user request * * *'
  write(*,*)
  allocate (h_unique(maxval(nequiv)))
  allocate (s_unique(maxval(nequiv))) 
  allocate (ncount(maxval(nequiv)))
  m = 0
  h_unique = 0.0d0
  s_unique = 0.0d0
  ncount = 0
  do k = 1, nUniqueExtra
    do l = 1, k
      m = m + 1
      h_unique(nequiv(m)) = h_unique(nequiv(m)) + abs(h(allExtra_mebfs(k),allExtra_mebfs(l)))
      s_unique(nequiv(m)) = s_unique(nequiv(m)) + abs(s(allExtra_mebfs(k),allExtra_mebfs(l)))
      ncount(nequiv(m)) = ncount(nequiv(m)) + 1
    end do
  end do
  do j = 1, maxval(nequiv)
    h_unique(j) = h_unique(j) / ncount(j)
    s_unique(j) = s_unique(j) / ncount(j)
  end do
  m = 0
  do k = 1, nUniqueExtra
    do l = 1, k
      m = m + 1
      h(allExtra_mebfs(k),allExtra_mebfs(l)) = h_unique(nequiv(m))*   &
               h(allExtra_mebfs(k),allExtra_mebfs(l))/abs(h(allExtra_mebfs(k),allExtra_mebfs(l)))
      h(allExtra_mebfs(l),allExtra_mebfs(k)) = h(allExtra_mebfs(k),allExtra_mebfs(l))
      s(allExtra_mebfs(k),allExtra_mebfs(l)) = s_unique(nequiv(m))*   &
               s(allExtra_mebfs(k),allExtra_mebfs(l))/abs(s(allExtra_mebfs(k),allExtra_mebfs(l)))
      s(allExtra_mebfs(l),allExtra_mebfs(k)) = s(allExtra_mebfs(k),allExtra_mebfs(l))
    end do
  end do
  title = ' (partially) symmetrized Hamiltonian'
  useLabels = .true.
  call corr_shift_printHam(title,h,nbase,useLabels)
  write(*,*)
  if (print_overlap) then
    title = ' (partially) symmetrized overlap matrix'
    call corr_shift_printHam(title,s,nbase,useLabels)
  endif
  deallocate(h_unique,s_unique,ncount)
endif  


jj = 0
do i = 1 ,nBlocks
  lwork = 4 * dimens(i)
  allocate ( h_block(dimens(i),dimens(i)) )
  allocate ( s_block(dimens(i),dimens(i)) )
  allocate ( s_save(dimens(i),dimens(i)) )
  allocate ( work(lwork) )
  allocate ( eps(dimens(i)) )
  allocate ( wGN(dimens(i),dimens(i)) )
  h_block = 0.0d0
  s_block = 0.0d0
  work = 0.0d0
  eps = 0.0d0
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
  s_save=s_block
  call dsygv(1,'V','L',dimens(i),h_block,dimens(i),s_save,dimens(i),eps,work,lwork,info)
  write(*,'(a7,i1,a)')' Block ',i,' : Energies and eigenvectors'
  write(*,'(12(2x,F14.8))')eps
  do j = 1, dimens(i)
    write(*,'(12F16.8)')h_block(j,:)
    jj = jj + 1
    do k = 1, dimens(i)
      extraVec(jj,blocks(i,k)) = h_block(k,j)
    end do
  end do
!  h_block = transpose(h_block)
  useLabels = .false.
  call corr_shift_gnweight(dimens(i),wGN,h_block,s_block,useLabels)
  write(*,*)
  deallocate ( work )
  deallocate ( eps )
  deallocate ( h_block )
  deallocate ( s_block , s_save)
  deallocate ( wGN )
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
  call dsygv(1,'V','L',reduced_extraDim,h_block,reduced_extraDim,s_block,reduced_extraDim,eps,work,lwork,info)
  deallocate( work, eps )
  allocate( mVec(reduced_extraDim,nbase) )
  allocate( reduced_extraVec(reduced_extraDim,nbase) )
  if ( select_lowest ) then
    offset = 0
    k = 0
    do i = 1, nBlocks
      do j = 1, nDressed_MEBFs(i)
        k = k + 1
        reduced_extraVec(k,:) = extraVec(k + offset,:)
      end do
      offset = offset + dimens(i) - k
    end do
  else
! keep track which vector has already been selected
    allocate( used(extraDim) )
    do j = 1, extraDim
      used(j) = .false.
    end do
    mVec = 0.0
    do i = 1, reduced_extraDim
      do j = 1, nbase
        do k = 1, reduced_extraDim
          if ( j .eq. ovlp_mebfs(k) ) mVec(i,j) = h_block(k,i)
        end do
      end do
    end do
    if ( print_overlap ) then
      write(*,*) 'Overlaps between original and dressed MEBFs'
      write(*,*) '- - - - - - - - - - - - - - - - - - - - - -'
      write(*,*) '    MEBF     dressed MEBF     overlap'
    end if
    do i = 1, reduced_extraDim
      maxovlp = 0.0
      do j = 1, extraDim
        ovlp = 0.0
        do k = 1, nbase
          do l = 1, nbase
            ovlp = ovlp + mVec(i,k) * extraVec(j,l) * S(k,l)
          end do
        end do
        if ((abs(ovlp) .gt. maxovlp) .and. (.not.used(j))) then
          maxovlp = abs(ovlp)
          jj = j
          reduced_extraVec(i,:) = extraVec(j,:)
        end if
        if ( print_overlap ) then
          write(*,'(3x,I4,7x,I4,8x,F14.8)') i,j,ovlp
        end if
      end do
      used(jj) = .true.
    end do
  end if
  write(*,*)
  write(*,*) 'new MEBFs (in rows) expressed in orginal MEBFs'
  first = 1; last = min(15,nbase)
  do while (first .le. nbase)
    if ( maxLength .le.10 ) then
      write(*,'(7x,15(A10,4x))')(mebfLabel(i),i=first,last)
    else
      write(*,'(15(10x,i4))')(i,i=first,last)
    endif
    do i = 1, reduced_extraDim
      write(*,'(I3,a,15F14.8)')i,':',(reduced_extraVec(i,j),j=first,last)
    end do
    first = last + 1
    last = min(first + 14,nbase)
    write(*,*)
  end do
  write(*,*)
  do i = 1, reduced_extraDim
    do j = 1, reduced_extraDim
      do k = 1, nbase
        do l = 1, nbase
          extraH(i,j) = extraH(i,j) + reduced_extraVec(i,k) * reduced_extraVec(j,l) * H(k,l)
          extraS(i,j) = extraS(i,j) + reduced_extraVec(i,k) * reduced_extraVec(j,l) * S(k,l)
        end do
      end do
    end do
  end do
  title = ' Hamiltonian in new MEBF basis'
  useLabels = .false.
  call corr_shift_printHam(title,extraH,reduced_extraDim,useLabels)
  title = ' Overlap matrix new MEBFs'
  call corr_shift_printHam(title,extraS,reduced_extraDim,useLabels)
  write(*,*)
  write(*,*) '  *  *  *  Electronic Couplings  (meV) *  *  *'
  title = ' '
  useLabels = .false.
  call corr_shift_couplings(title,extraH,extraS,reduced_extraDim,useLabels)
  deallocate(h_block,s_block)
  deallocate(mVec,reduced_extraVec)
  deallocate(extraH,extraS)
end if

end subroutine corr_shift_superNOCI



subroutine corr_shift_readin
use corr_shift_input_data
implicit none

external  :: corr_shift_capitalize,corr_shift_locate

integer,parameter    :: nKeys = 8
integer              :: jj,iKey,i,j,k,first, last, maxunique,iMEBF,nelem,aux
integer,allocatable  :: unique(:)

character(len=4), dimension(nKeys)   :: keyword
character(len=4)                     :: key
character(len=5)                     :: dummy
character(len=6),allocatable         :: user_label(:,:)
character(len=18),allocatable        :: blocks_label(:,:),ovlp_mebflabels(:),unique_labels(:)
character(len=132)                   :: line,filename

logical                              :: all_ok = .true.
logical                              :: new
logical, dimension(nKeys)            :: hit = .false.


data keyword /'PROJ','SHIF','BLOC','SELE','GNWE','LOWE','PROV','AVER'/

extra = .false.
dressed_coupling = .false.
GN_weights = .false.
select_lowest = .false.
print_overlap = .false.
dcec = .false.
! dcec = dynamic correlation energy correction
averaging = .false.
nUniqueExtra = 0

do while (all_ok)
  read(5,*,iostat=jj) line
  line = adjustl(line)
  key = line(1:4)
  call corr_shift_capitalize(key)
  do iKey = 1, nKeys
    if ( key .eq. keyword(iKey) ) hit(iKey) = .true.
  end do
  if (  jj .lt. 0 ) all_ok = .false.
end do

do iKey = 1, nKeys
  if ( hit(iKey) ) then
    select case(iKey)
      case(1)
        call corr_shift_locate(5,'PROJ',line)
        read(*,*) project
        project = trim(project)
        write(filename,'(2a)')trim(Project),'.arx'
        open(11,file=filename,err=9000)
        call corr_shift_locate(11,'Stat',line)
        read(line,*) dummy, nbase, nmol
        allocate(mebfLabel(nbase))
        allocate( shift(nbase) )
        maxLength = 0
        do i = 1, nbase
          read(11,'(1x,a)')mebfLabel(i)
          mebfLabel(i) = adjustl(mebfLabel(i))
          maxLength = max(maxLength,len(mebfLabel(i)))
        end do
        allocate( h(nbase,nbase) )
        allocate( s(nbase,nbase) )
        write(*,'(2a)') 'H and S are read from ',filename
        call corr_shift_locate(11,'Hami',line)
        read(11,*) line
        first = 1 
        last = min(7,nbase)
        do while ( first .lt. nbase )
          do i = 1, nbase
            read(11,662)jj,(h(i,j),j=first,last)
          end do
          read(11,*) line
          first = last + 1
          last = min(last + 7,nbase)
        end do
662     format(i5,1x,7e20.13)
        call corr_shift_locate(11,'Over',line)
        read(11,*) line
        first = 1
        last = min(7,nbase)
        do while ( first .lt. nbase )
          do i = 1, nbase
            read(11,662)jj,(s(i,j),j=first,last)
          end do
          read(11,*) line
          first = last + 1
          last = min(last + 7,nbase)
        end do
      case(2)
        dcec = .true.
        allocate( fraglabel(nmol,nbase) )
        allocate( unique(nmol) )
        call corr_shift_locate(11,'Frag',line)
        do i = 1, nmol
          read(11,*) (fragLabel(i,j),j=1,nbase)
          unique(i) = 1
          do j = 2, nbase
            new = .true.
            do k = 1, j-1
              if (trim(fraglabel(i,k)).eq.trim(fraglabel(i,j))) new = .false.
            end do
            if ( new ) unique(i) = unique(i) + 1
          end do
        end do
        maxunique = maxval(unique)
        allocate ( ecorr(nmol,maxunique) )
        allocate ( ecorr2(nmol,nbase) )
        allocate ( user_label(nmol,maxunique) )
        call corr_shift_locate(5,'SHIF',line)
        shift = 0.0d0
        do i = 1, nmol
          read(*,*) (user_label(i,k),k=1,unique(i))
          read(*,*) (ecorr(i,k),k=1,unique(i))
          jj = 1
          do j = 1, nbase
            do k = 1, unique(i)
              if (fraglabel(i,j) .eq. user_label(i,k)) then
                shift(jj) = shift(jj) + ecorr(i,k)
                ecorr2(i,j) = ecorr(i,k)
              endif
            end do
            jj = jj + 1
          end do
        end do
      case(3)
        call corr_shift_locate(5,'BLOC',line)
        read(*,*) nBlocks
        allocate( dimens(nBlocks) )
        allocate( blocks(nBlocks,nbase) )
        allocate( blocks_label(nBlocks,nbase) )
        blocks = 0
        dimens = 0
        extraDim = 0
        do i = 1, nBlocks
          read(*,*) dimens(i),(blocks_label(i,j),j=1,dimens(i))
          extraDim = extraDim + dimens(i)
          do iMEBF = 1, nbase
            do j = 1, dimens(i)
              if (trim(blocks_label(i,j)) .eq. trim(mebfLabel(iMEBF))) then
                blocks(i,j) = iMEBF
              endif
            end do
          end do
        end do
        allocate(unique_labels(sum(dimens)))
        unique_labels=''
        do i = 1, nBlocks
          do j = 1, dimens(i)
            if ( .not. any(unique_labels .eq. blocks_label(i,j)) ) then
              nUniqueExtra = nUniqueExtra + 1
              unique_labels(nUniqueExtra) = blocks_label(i,j)
            endif
          end do
        end do
        allocate(allExtra_mebfs(nUniqueExtra))
        do iMEBF = 1, nbase
          do j = 1, nUniqueExtra
            if (trim(unique_labels(j)) .eq. trim(mebfLabel(iMEBF))) allExtra_mebfs(j) = iMEBF
          end do
        end do
        do i = 1, nUniqueExtra
          do j = 1, i
            if (allExtra_mebfs(j).gt.allExtra_mebfs(i)) then
              aux = allExtra_mebfs(j)
              allExtra_mebfs(j) = allExtra_mebfs(i)
              allExtra_mebfs(i) = aux
            endif
          end do
        end do
        extra = .true.
      case(4)
        call corr_shift_locate(5,'SELE',line)
        read(*,*) reduced_extraDim
        allocate( ovlp_mebfs(reduced_extraDim) )
        allocate( ovlp_mebflabels(reduced_extraDim) )
        ovlp_mebfs = 0
        read(*,*) (ovlp_mebflabels(i),i=1,reduced_extraDim)
        do iMEBF = 1, nbase
          do j = 1, reduced_extraDim 
            if (trim(ovlp_mebflabels(j)) .eq. trim(mebfLabel(iMEBF))) ovlp_mebfs(j) = iMEBF
          end do
        end do
        dressed_coupling = .true.
      case(5)
        GN_weights = .true.
      case(6)
        call corr_shift_locate(5,'LOWE',line)
        allocate(nDressed_MEBFs(nBlocks))
        read(*,*)(nDressed_MEBFs(i),i=1,nBlocks)
        reduced_extraDim = sum(nDressed_MEBFs)
        select_lowest = .true.
        dressed_coupling = .true.
      case(7)
        print_overlap = .true.
      case(8)
        averaging = .true.
        nelem = int((nUniqueExtra*(nUniqueExtra+1))/2)
        allocate(nequiv(nelem))
        nequiv = 0
        call corr_shift_locate(5,'AVER',line)
        read(*,*)(nequiv(j),j=1,nelem)
    end select
  end if
end do     

return

9000 write(*,*)'A problem ocurred opening ',filename
stop

end subroutine corr_shift_readin


subroutine corr_shift_capitalize(string)
implicit none
integer      :: i
character(*) string

do i = 1, len(string)
  if (ichar(string(i:i)).gt.96) then
    string(i:i) = char(ichar(string(i:i))-32)
  endif
end do
return
end subroutine corr_shift_capitalize


subroutine corr_shift_locate(lu,string,line)
implicit none

external :: corr_shift_capitalize

integer        ::  lu
character(4)   ::  string
character(132) ::  line,string2
rewind(lu)
40 read(lu,'(A132)') line
string2=adjustl(line)
if (lu.eq.5) call corr_shift_capitalize(string2)
if (string2(1:4).ne.string) goto 40
return
end subroutine corr_shift_locate

!-----------------------------------------------------------!
!                                                           !
! Program for post NOCI-F analysis                          !
!    - application of correlation correction                !
!    - subblock diagonalizations                            !
!                                                           !
! Required files                                            !
!    - input file (see below)                               !
!    - arx file, generated by GronOR                        !
!                                                           !
!  written by Coen de Graaf                                 !
!             URV/ICREA, 2021-2023                          !
!                                                           !
!-----------------------------------------------------------!
!                                                           !
! Input example                                             !
!                                                           !
! Project                                                   !
!  benzene                                                  !
! Shift                                                     !
!  S0 S1 T1 D+ D-                                           !
!  0.0 -0.15 -0.11 -0.11 -0.11                              !
!  S0 S1 T1 D+ D-                                           !
!  0.0 -0.15 -0.11 -0.11 -0.11                              !
! Block                                                     !
!  2                                                        !
!  4 S0S1 S1S0 D+D- D-D+                                    !
!  3 T1T1 D+D- D-D+                                         !
! Select                                                    !
!  3                                                        !
!  S0S1 S1S0 T1T1                                           !
!                                                           !
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -!
!                                                           !
! PROJect :  Root of the arx file                           !
! SHIFt   :  Correlation correction of the fragment         !
!            states (in au)                                 !
! BLOCk   :  Number of subblock to be diagonalized,         !
!            followed by one line for each subblock.        !
!            Each line contains the dimension of the        !
!            subblock and its MEBF labels                   !
! SELEct  :  Selection of the vectors resulting from the    !
!            subblock diagonalizations to construct a       !
!            new Hamiltonian with 'dressed' MEBFs. First,   !
!            the number of vectors that will be selected,   !
!            followed by the MEBF labels that should be     !
!            dominant in the new 'dressed' MEBFs.           !
!                                                           !
!-----------------------------------------------------------!
