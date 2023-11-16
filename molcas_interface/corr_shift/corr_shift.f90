module corr_shift_input_data
  implicit none
  integer   :: nbase,nstates,extraDim,nBlocks,nmol
  integer   :: reduced_extraDim
  integer, allocatable  :: ncombv(:,:),dimens(:),blocks(:,:)
  integer, allocatable  :: ovlp_mebfs(:),nDressed_MEBFs(:)
  real (kind=8), allocatable   :: h(:,:),s(:,:)
  real (kind=8), allocatable   :: ecorr2(:,:)
  character(len=80)     :: project
  character(len=6),allocatable    :: fragLabel(:)
  character(len=18),allocatable   :: mebfLabel(:)
  logical :: extra,dressed_coupling,GN_weights,select_lowest,print_overlap
end module corr_shift_input_data


program  corr_shift_main
use corr_shift_input_data
implicit none

external  :: corr_shift_readin,corr_shift_printHam
external  :: corr_shift_lowdin,corr_shift_gnweight
external  :: dgetrf,dgetri,corr_shift_nocifunction
external  :: corr_shift_couplings,corr_shift_superNOCI

real( kind = 8 ), allocatable :: H_ortho(:,:),H_shift(:,:)
real( kind = 8 ), allocatable :: H_wshift(:,:),weight(:,:)
real( kind = 8 ), allocatable :: s_lowdin(:,:)
real (kind = 8 ), allocatable :: ssave(:,:),hsave(:,:)
real( kind = 8 ), allocatable :: shift(:),work(:)
real( kind = 8 )              :: wshift
integer                       :: i,j,info
integer, allocatable          :: ipiv(:)
character (len = 120)         :: title1,title2,title3
logical                       :: dump_in_out,useLabels

call corr_shift_readin

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


!title1 = 'MEBF Overlap matrix'
!useLabels = .true.
!call corr_shift_printHam(title1,s,nbase,useLabels)

shift = 0.0d0
do i = 1, nbase
  shift(i) = sum(ecorr2(:,i))
enddo
write(*,*)
write(*,*) ' Shifts applied on the diagonal of H'
write(*,*) ' MEBF    Total shift            Monomer shifts'
do i = 1, nbase
  write(*,679) i,shift(i),(fragLabel(i+(j-1)*nbase),ecorr2(j,i),j=1,nmol)
end do
679  format(I4,F15.8,8x,5(A6,F15.8,3x))
write(*,*)

ssave = s
call corr_shift_lowdin(nbase,s,s_lowdin)
s = ssave
H_ortho = matmul( transpose(s_lowdin) , matmul(H , s_lowdin) )
H_shift  = H_ortho
H_wshift = H_ortho
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
end if

title1 = 'Unshifted Hamiltonian'
title2 = 'Shifted Hamiltonian'
title3 = 'GN-weighted shifted Hamiltonian'

!* Print the unshifted, shifted and GN-weighted shifted Hamiltonian
write(*,*)'  *  *  Orthogonal MEBF basis  *  *'
write(*,*)
useLabels = .true.
call corr_shift_printHam(title2,h_shift,nbase,useLabels)
if ( GN_weights ) call corr_shift_printHam(title3,h_wshift,nbase,useLabels)

call dgetrf(nbase,nbase,s_lowdin,nbase,ipiv,info)
call dgetri(nbase,s_lowdin,nbase,ipiv,work,3*nbase,info)
h_shift = matmul( transpose(s_lowdin) , matmul(h_shift,s_lowdin) )
if ( GN_weights ) then
  h_wshift = matmul(transpose(s_lowdin),matmul(h_wshift,s_lowdin))
end if

write(*,*)'  *  *  Original non-orthogonal MEBF basis  *  *'
write(*,*)
call corr_shift_printHam(title1,H,nbase,useLabels)
call corr_shift_printHam(title2,h_shift,nbase,useLabels)
if ( GN_weights ) call corr_shift_printHam(title3,h_wshift,nbase,useLabels)
write(*,*)

!** Calculate the couplings for the three Hamiltonians

write(*,*) '  *  *  *  Electronic Couplings  (meV) *  *  *'
write(*,*)
useLabels = .true.
call corr_shift_couplings(title1,h,s,nbase,useLabels)
call corr_shift_couplings(title2,h_shift,s,nbase,useLabels)
if ( GN_weights ) call corr_shift_couplings(title3,h_wshift,s,nbase,useLabels)
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
hsave = h_shift
call corr_shift_nocifunction(title2,h_shift,s,nbase,dump_in_out)
h_shift = hsave
s = ssave
if ( GN_weights ) then
  hsave = h_wshift
  call corr_shift_nocifunction(title3,h_wshift,s,nbase,dump_in_out)
  h_wshift = hsave
  s = ssave
end if
if ( extra ) then
  write(*,*)
  write(*,*)
  write(*,*)
  write(*,*)'Extra NOCI (unshifted)'
  write(*,*)'----------------------'
  call corr_shift_superNOCI
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
  end if
end if

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
use  corr_shift_input_data, only : mebfLabel
implicit none
real ( kind = 8 ), intent(in)  :: h(n,n),s(n,n)
real ( kind = 8 ), allocatable :: tc(:,:)
integer                        :: n,i,j,k,first,last
character ( len = 120 )        :: title
logical                        :: useLabels

allocate( tc(n,n) )

do i = 2, n
  do j = 1, i-1
    tc(i,j) = (h(i,j) - 0.5*(h(i,i) + h(j,j))*s(i,j)) / (1.0 - s(i,j)*s(i,j))
  end do
end do
write(*,*) title
first = 1
last = min(7,n)
do while (first .lt. n-1)
  if ( useLabels ) then
    write(*,671) (adjustr(mebfLabel(j)),j=first,min(n-1,last))
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
671 format(17x,7(a20))
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
use  corr_shift_input_data, only : mebfLabel
implicit none

external  :: dgetrf,dgetri

!  compute the Gallup-Norbeck weights
integer, intent(in)       :: nbase
integer                   :: i,k,info,first,last
integer                   :: ipiv(nbase)

real(kind=8), intent(in)  :: slow(nbase,nbase)
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
      write(*,673)adjustr(mebfLabel(i)),(wgn(i,k),k=first,last)
    else
      write(*,674)i,(wgn(i,k),k=first,last)
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
use  corr_shift_input_data, only : mebfLabel
implicit none
real ( kind = 8 ), intent(in)   :: m(n,n)
integer                         :: n,i,j,first,last
character (len = 120)           :: title
logical                         :: useLabels
write(*,'(a120)') title
write(*,*)
first = 1 ; last = min(n,7)
do while ( first .le. n )
  if (useLabels) then
    write(*,671) (adjustr(mebfLabel(i)),i=first,last)
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
671 format(17x,7(a20))
672 format(a20,1x,7f20.10)
673 format(6x,7(6x,i8,6x))
674 format(i5,1x,7f20.10)
return
end subroutine corr_shift_printHam


subroutine corr_shift_superNOCI
use corr_shift_input_data
implicit none

external :: dsygv,corr_shift_printHam,corr_shift_couplings,corr_shift_gnweight

integer                      :: i,j,k,l,info,lwork,jj,offset,first,last
real (kind=8), allocatable   :: h_block(:,:)
real (kind=8), allocatable   :: s_block(:,:),s_save(:,:)
real (kind=8), allocatable   :: extraVec(:,:)
real (kind=8), allocatable   :: extraH(:,:),extraS(:,:)
real (kind=8), allocatable   :: eps(:)
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

jj = 0
do i = 1 ,nBlocks
  lwork = 4 * dimens(i)
  allocate ( h_block(dimens(i),dimens(i)) )
  allocate ( s_block(dimens(i),dimens(i)) )
  allocate ( s_save(dimens(i),dimens(i)) )
  allocate ( work(lwork) )
  allocate ( eps(dimens(i)) )
  allocate ( wGN(dimens(i),dimens(i)) )
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
  write(*,*) 'new MEBFs expressed in orginal MEBFs'
  first = 1; last = min(15,nbase)
  do while (first .le. nbase)
    write(*,'(15(8x,A6))')(mebfLabel(i),i=first,last)
    do i = 1, reduced_extraDim
      write(*,'(I3,a,15F14.8)')i,':',(reduced_extraVec(i,j),j=first,last)
    end do
    first = last + 1
    last = min(first + 15,nbase)
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

integer,parameter    :: nKeys = 9
integer              :: jj,iKey,i,j,k,first, last, maxunique
integer,allocatable  :: unique(:)

real (kind=8), allocatable   :: ecorr(:)

character(len=4), dimension(nKeys)   :: keyword
character(len=4)                     :: key
character(len=5)                     :: dummy
character(len=6),allocatable         :: user_label(:)
character(len=132)                   :: line,filename

logical                              :: all_ok = .true.
logical                              :: new
logical, dimension(nKeys)            :: hit = .false.


data keyword /'PROJ','SHIF','HAMI','OVER','BLOC','SELE','GNWE','LOWE','PROV'/

extra = .false.
dressed_coupling = .false.
GN_weights = .false.
select_lowest = .false.
print_overlap = .false.

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
!        hit(3) = .false.                            ! deactivate reading H and S from input
!        hit(4) = .false.
        write(filename,'(2a)')trim(Project),'.arx'
        open(11,file=filename,err=9000)
        call corr_shift_locate(11,'Stat',line)
        read(line,*) dummy, nbase, nmol
        allocate(mebfLabel(nbase))
        do i = 1, nbase
          read(11,'(a)')mebfLabel(i)
        end do
        allocate( h(nbase,nbase) )
        allocate( s(nbase,nbase) )
        if ( .not. hit(3) .and. .not. hit(4) ) then
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
662       format(i5,1x,7e20.13)
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
        end if
        nstates = nmol*nbase
        allocate( fraglabel(nstates) )
        allocate( unique(nmol) )
        call corr_shift_locate(11,'Frag',line)
        read(11,*)(fragLabel(i),i=1,nstates)
        do i = 1, nmol
          unique(i) = 1
          do j = 2, nbase
            new = .true.
            do k = 1, j-1
              if (fraglabel(k).eq.fraglabel(j)) new = .false.
            end do
            if ( new ) unique(i) = unique(i) + 1
          end do
        end do
        maxunique = maxval(unique)
        allocate ( ecorr(nmol*maxunique) )
        allocate ( ecorr2(nmol,nbase) )
        allocate ( user_label(nmol*maxunique) )
      case(2)
        call corr_shift_locate(5,'SHIF',line)
        j=1
        write(*,*) 'unique',unique
        do i = 1, nmol
          read(*,*) dummy
          read(*,*) (user_label(k),k=j,unique(i)+j-1)
          read(*,*) (ecorr(k),k=j,unique(i)+j-1)
          j = j + unique(i)
        end do
        do i = 1, nmol
          jj = 1
          do j = (i-1)*nbase + 1, i*nbase
            do k = 1, unique(i)
              if ( fraglabel(j).eq.user_label(k+(i-1)*unique(i)) ) ecorr2(i,jj) = ecorr(k)
            end do
            jj = jj + 1
          end do
        end do
      case(3)
        call corr_shift_locate(5,'HAMI',line)
        allocate( h(nbase,nbase) )
        do i = 1, nbase
          read(*,*)(h(i,j),j=1,nbase)
        end do
      case(4)
        call corr_shift_locate(5,'OVER',line)
        allocate( s(nbase,nbase) )
        do i = 1, nbase
          read(*,*)(s(i,j),j=1,nbase)
        end do
      case(5)
        call corr_shift_locate(5,'BLOC',line)
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
      case(6)
        call corr_shift_locate(5,'SELE',line)
        read(*,*) reduced_extraDim
        allocate( ovlp_mebfs(reduced_extraDim) )
        ovlp_mebfs = 0
        read(*,*) (ovlp_mebfs(i),i=1,reduced_extraDim)
        dressed_coupling = .true.
      case(7)
        GN_weights = .true.
      case(8)
        call corr_shift_locate(5,'LOWE',line)
        allocate(nDressed_MEBFs(nBlocks))
        read(*,*)(nDressed_MEBFs(i),i=1,nBlocks)
        reduced_extraDim = sum(nDressed_MEBFs)
        select_lowest = .true.
        dressed_coupling = .true.
      case(9)
        print_overlap = .true.
    end select
  end if
end do     

return

9000 write(*,*)'A problem ocurred opening ',filename
stop

end subroutine corr_shift_readin
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
! Hamiltonian and Overlap are normally read from the arx file
! Fragment labels need to be implemented to comply with the naming scheme in
! gcommon and GronOR


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

