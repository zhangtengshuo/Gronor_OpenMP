!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or
!     modified under
!     the Apache License version 2.0
!     (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license.
!     Unless required
!     by applicable law or agreed to in writing, software distributed
!     under the license
!     is distributed on an ‘as is’ basis, without warranties or
!     conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions
!     and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!>  @brief
!>    Construct the MEBFs from the fragment functions
!>
!>  @authors
!>    Aitor Sanchez-Mansilla, Coen de Graaf
!>    (URV, Tarragona)
!>
!>  @date
!>    November 2020, rewritten from scratch January 2023
!>
!>  @details
!>    MEBF(iMEBF) is made by multiplying the wave functions of the
!>    different fragments. First, state1 is combined with state2. 
!>    This new state is then taken as a new fragment wave function
!>    and combined with state3. This is repeated until all fragments
!>    are considered.
!>
!>  @param state1 
!>  @param state2 
!>
!>  @todo
!>    Get rid of the packing of the orbital occupations, which
!>    limits the number of active orbitals to 32
!>

subroutine gronor_make_basestate(iMEBF,first_pass)
use gnome_parameters, only : tau_CI,idbg
use gnome_data      , only : nbasis
use makebasedata
use cidef
use cidist          , only :  me,master

implicit none

integer,intent(in)              :: iMEBF
integer                         :: ndets1,ndets2,newdets
integer,allocatable             :: micro_ndets1(:),micro_ndets2(:)
integer,allocatable             :: occ_num(:)
integer                         :: spin1,spin2,target_spin
integer                         :: micro_dets1,micro_dets2
integer                         :: i,j,k,idet,jdet,iFrag,iFragWF,iAct
integer                         :: ms,ms1,ms2
integer                         :: first1,first2,last1,last2
integer                         :: jstart,kstart
integer                         :: perm_ab,isetsign
real(kind=8)                    :: prod,dnorm
real(kind=8)                    :: clebsch_gordon
real(kind=8),parameter          :: two=2.0d0
real(kind=8),allocatable        :: coef1(:),coef2(:)
real(kind=8),allocatable        :: coefmebf(:)
character(len=255)              :: dumstr
character(len=255),allocatable  :: occ1(:),occ2(:)
character(len=255),allocatable  :: occmebf(:)
logical,intent(in)              :: first_pass
logical                         :: debug


debug = .true.
!  ==== Combining the vec files of the fragments (only on the first pass)

if (first_pass) then
  if ( idbg .ge. 50 ) then
    write(lfndbg,'(a,l3,i5)') 'allocated inactb  ?',allocated(inactb),size(inactb)
    write(lfndbg,'(a,l3,i5)') 'allocated inactm  ?',allocated(inactm),size(inactm)
    write(lfndbg,'(a,l3,i5)') 'allocated nactb   ?',allocated(nactb),size(nactb)
    write(lfndbg,'(a,l3,i5)') 'allocated nactm   ?',allocated(nactm),size(nactm)
    write(lfndbg,'(a,l3,2i5)') 'allocated ncombv  ?',allocated(ncombv),size(ncombv,1),size(ncombv,2)
    write(lfndbg,'(a,20i10)') 'inactm: ',inactm
    write(lfndbg,'(a,20i10)') 'nactm : ',nactm
    flush(lfndbg)
  endif
  inactb(iMEBF)=0
  nactb(iMEBF)=0
  do i=1,nmol
    inactb(iMEBF)=inactb(iMEBF)+inactm(ncombv(i,iMEBF))
    nactb(iMEBF)=nactb(iMEBF)+nactm(ncombv(i,iMEBF))
  enddo
  if ( idbg .ge. 50 ) then
    write(lfndbg,*) 'first_pass: ',first_pass
    write(lfndbg,'(a,i10)') 'iMEBF: ',iMEBF
    write(lfndbg,'(a,20i10)') 'inactb: ',inactb
    write(lfndbg,'(a,20i10)') 'nactb : ',nactb
    write(lfndbg,'(a,20i10)') 'ncombv: ',ncombv
 
    write(lfndbg,'(a,l3,3i5)') 'allocated vecsb ?',allocated(vecsb),size(vecsb,1),size(vecsb,2),size(vecsb,3)
    write(lfndbg,'(a,l3,3i5)') 'allocated vecsm ?',allocated(vecsm),size(vecsm,1),size(vecsm,2),size(vecsm,3)
    write(lfndbg,'(a,l3,3i5)') 'allocated nbasm ?',allocated(nbasm),size(nbasm)
    write(lfndbg,'(a,20i10)') 'nbasm : ',nbasm
    flush(lfndbg)
  end if
  do i=1,nbasis
    do j=1,nbasis
      vecsb(i,j,iMEBF)=0.0d0
    enddo
  enddo
  jstart=1
  kstart=1
  do i=1,nmol
    do j=1,inactm(ncombv(i,iMEBF))
      do k=1,nbasm(ncombv(i,iMEBF))
        vecsb(kstart+k-1,jstart+j-1,iMEBF)= vecsm(k,j,ncombv(i,iMEBF))
      enddo
    enddo
    kstart=kstart+nbasm(ncombv(i,iMEBF))
    jstart=jstart+inactm(ncombv(i,iMEBF))
  enddo
  
  jstart=inactb(iMEBF)+1
  kstart=1
  
  do i=1,nmol
    do j=1,nactm(ncombv(i,iMEBF))
      do k=1,nbasm(ncombv(i,iMEBF))
        vecsb(kstart+k-1,jstart+j-1,iMEBF)=vecsm(k,j+inactm(ncombv(i,iMEBF)),ncombv(i,iMEBF))
      enddo
    enddo
    kstart=kstart+nbasm(ncombv(i,iMEBF))
    jstart=jstart+nactm(ncombv(i,iMEBF))
  enddo
endif

!  ==== Spin function ======
! First fragment, set the stage

if ( idbg .ge. 50 ) then
  write(lfndbg,'(a,l3,i5)') 'allocated idetm ?',allocated(idetm),size(idetm)
  write(lfndbg,'(a,l3,i5)') 'allocated spinm ?',allocated(spinm),size(spinm)
  write(lfndbg,'(a,20i10)') 'idetm : ',idetm
  write(lfndbg,'(a,20i10)') 'spinm : ',spinm
  write(lfndbg,'(a,l3,2i5)') 'allocated civm ?',allocated(civm),size(civm,1),size(civm,2)
  write(lfndbg,'(a,l3,2i5)') 'allocated occm_string ?',allocated(occm_string),size(occm_string,1),size(occm_string,2)
  flush(lfndbg)
endif
iFragWF = ncombv(1,iMEBF)
ndets1   = idetm(iFragWF)
spin1    = spinm(iFragWF)
spinFrag = spin1
allocate(coef(ndets1))
allocate(occ(ndets1))
allocate(micro_ndets1(spin1))
if ( idbg .ge. 50 ) then
  write(lfndbg,'(a,4i5)') 'iFragWF,ndets1,spin1,spinFrag : ',iFragWF,ndets1,spin1,spinFrag 
  write(lfndbg,'(a,l3,i5)') 'allocated coef ?',allocated(coef),size(coef)
  write(lfndbg,'(a,l3,i5)') 'allocated occ  ?',allocated(occ),size(occ)
  write(lfndbg,'(a,l3,i5)') 'allocated micro_ndets1 ?',allocated(micro_ndets1),size(micro_ndets1)
  flush(lfndbg)
endif
do idet = 1, ndets1
  coef(idet) = civm(idet,iFragWF)
  occ(idet) = occm_string(idet,iFragWF)
end do
micro_ndets1 = 0

if ( idbg .ge. 50 ) then
  write(lfndbg,'(a,i3)') 'spin for frag 1', spin1
  write(lfndbg,*) 'determinants for frag 1'
  do idet = 1, ndets1
    write(lfndbg,'(i6,f15.8,3x,a)') idet,coef(idet),trim(occ(idet))
  end do
  flush(lfndbg)
endif

if (spin1 .gt. 1) then
  call generate_microstates(ndets1,micro_dets1)
  allocate(coef1(micro_dets1))
  allocate(occ1(micro_dets1))
  do ms = 1, spin1
    micro_ndets1(ms) = micro_ndets(ms)
  end do
  do idet = 1, micro_dets1
    coef1(idet) = micro_coef(idet)
    occ1(idet)  = micro_occ(idet)
  end do
  if ( idbg .ge. 50 ) then
    write(lfndbg,*) 'micro states of frag 1'
    do idet = 1, micro_dets1
      write(lfndbg,'(I6,F15.8,3x,A)') idet,coef1(idet),trim(occ1(idet))
    end do
    flush(lfndbg)
  endif
! allocated in generate_microstates
  deallocate(micro_ndets,micro_coef,micro_occ)
else
  micro_dets1 = ndets1
  allocate(coef1(micro_dets1))
  allocate(occ1(micro_dets1))
  micro_ndets1(1) = micro_dets1
  do idet = 1, micro_dets1
    coef1(idet) = coef(idet)
    occ1(idet) = occ(idet)
  end do
endif
deallocate(coef,occ)
if (idbg .ge. 50 ) then
  write(lfndbg,*)'First fragment is initialized'
  write(lfndbg,'(a,l3,i5)') 'allocated coef1?',allocated(coef1),size(coef1)
  write(lfndbg,'(a,l3,i5)') 'allocated occ1 ?',allocated(occ1),size(occ1)
  write(lfndbg,'(a,l3,i5)') 'allocated micro_ndets1 ?',allocated(micro_ndets1),size(micro_ndets1)
  write(lfndbg,'(a,20i5)') 'micro_ndets1 : ',micro_ndets1
  write(lfndbg,'(a,l3,i5)') 'allocated coef ?',allocated(coef),size(coef)
  write(lfndbg,'(a,l3,i5)') 'allocated occ  ?',allocated(occ),size(occ)
  write(lfndbg,'(a,l3,i5)') 'allocated micro_coef ?',allocated(micro_coef),size(micro_coef)
  write(lfndbg,'(a,l3,i5)') 'allocated micro_occ  ?',allocated(micro_occ),size(micro_occ)
  write(lfndbg,'(a,l3,i5)') 'allocated micro_ndets?',allocated(micro_ndets),size(micro_ndets)
  flush(lfndbg)
endif


! When there is only one fragment, not much has to be done
if (nmol .eq. 1) then
  if (first_pass) then
    if (iMEBF .eq. 1) then
      maxcib = micro_dets1
      idetb(1) = micro_dets1
    else
      maxcib = max(maxcib,micro_dets1)
      idetb(iMEBF) = micro_dets1
    endif
  else
    if ( idbg .ge. 50 ) then
      write(lfndbg,'(a,i3)') 'civb and occupations for MEBF ',iMEBF
      flush(lfndbg)
    end if
    allocate(occ_num(nactb(iMEBF)))
    do idet = 1, micro_dets1
      dumstr = occ1(idet)
      do iAct = 1, nactb(iMEBF)
        if (dumstr(iAct:iAct) .eq. '2') occ_num(iAct) = 2
        if (dumstr(iAct:iAct) .eq. 'a') occ_num(iAct) = 1
        if (dumstr(iAct:iAct) .eq. 'b') occ_num(iAct) =-1
        if (dumstr(iAct:iAct) .eq. '0') occ_num(iAct) = 0
      end do 
      coef1(idet) = coef1(idet) * perm_ab(occ_num,nactb(iMEBF))
      coef1(idet) = coef1(idet) * isetsign(occ_num,nactb(iMEBF))
      civb(idet,iMEBF) = coef1(idet)
      call pack(ioccb(idet,iMEBF),occ_num,nactb(iMEBF))
      if ( idbg .ge. 50 ) then
        write(lfndbg,'(i6,f14.8,32i3)')idet,civb(idet,iMEBF),(occ_num(iAct),iAct=1,nactb(iMEBF))
        flush(lfndbg)
      endif
      maxcoef=max(maxcoef,abs(coef1(idet)))
    end do
  endif
endif

! Now, the other fragments
do iFrag = 2, nmol
  iFragWF = ncombv(iFrag,iMEBF)
  ndets2 = idetm(iFragWF)
  spin2 = spinm(iFragWF)
  spinFrag = spin2
  allocate(coef(ndets2))
  allocate(occ(ndets2))
  allocate(micro_ndets2(spin2))
  do idet = 1, ndets2
    coef(idet) = civm(idet,iFragWF)
    occ(idet) = occm_string(idet,iFragWF)
  end do
  micro_ndets2 = 0

  if ( idbg .ge. 50 ) then
    write(lfndbg,'(A,I3)') 'spin for frag 2', spin1
    write(lfndbg,*) 'determinants for frag 2'
    do idet = 1, ndets2
      write(lfndbg,'(I6,F15.8,3x,A)') idet,coef(idet),trim(occ(idet))
    end do
    flush(lfndbg)
  endif

  if ( spin2 .gt. 1 ) then
    call generate_microstates(ndets2,micro_dets2)
    allocate(coef2(micro_dets2))
    allocate(occ2(micro_dets2))
    do ms = 1, spin2
      micro_ndets2(ms) = micro_ndets(ms)
    end do
    do idet = 1, micro_dets2
      coef2(idet) = micro_coef(idet)
      occ2(idet)  = micro_occ(idet)
    end do
    if ( idbg .ge. 50 ) then
      write(lfndbg,*) 'micro states of frag 2'
      do idet = 1, micro_dets2
        write(lfndbg,'(I6,F15.8,3x,A)') idet,coef2(idet),trim(occ2(idet))
      end do
      flush(lfndbg)
    endif
!   allocated in generate_microstates
    deallocate(micro_ndets,micro_coef,micro_occ)
  else
    micro_dets2 = ndets2
    allocate(coef2(micro_dets2))
    allocate(occ2(micro_dets2))
    micro_ndets2(1) = micro_dets2
    do idet = 1, micro_dets2
      coef2(idet) = coef(idet)
      occ2(idet) = occ(idet)
    end do
  endif
  deallocate(coef,occ)
!
! merging occ1 and occ2 and save it into a new occ1.
! First counting the number of new determinants to correctly
! allocate the arrays
!
  if ( idbg .ge. 50 ) then
    write(lfndbg,'(A,12I5)') 'micro states 1 ',micro_ndets1
    write(lfndbg,'(A,12I5)') 'micro states 2 ',micro_ndets2
    if (nmol .gt. 2) then
      write(lfndbg,'(A,I4)') 'target spin    ',inter_couplings(iFrag-1,iMEBF)
    end if
    flush(lfndbg)
  endif
  if (nmol .gt. 2) then
    target_spin = inter_couplings(iFrag-1,iMEBF)-1
  else
    target_spin = nspin
  end if
  if ( idbg .ge. 50 ) then
    write(lfndbg,'(a,i3)') 'Target spin is ',target_spin
    flush(lfndbg)
  endif
  newdets = 0
  first1  = 1
  last1   = 0
  do i = 1, spin1
    ms1    = spin1-1 - 2*(i-1)
    last1  = last1 + micro_ndets1(i)
    first2 = 1
    last2  = 0
    do j = 1, spin2
      ms2 = spin2-1 - 2*(j-1)
      last2 = last2 + micro_ndets2(j)
      if (ms1 + ms2 .eq. target_spin) then
        do idet = first1, last1
          do jdet = first2, last2
            prod = coef1(idet) * coef2(jdet)
            if (abs(prod) .ge. tau_CI) then
              newdets = newdets + 1
            endif
          end do
        end do
      endif
      first2 = first2 + micro_ndets2(j)
    end do
    first1 = first1 + micro_ndets1(i)
  end do
  if ( idbg .ge. 50 ) then
    write(lfndbg,'(2(a,i2),a,i5)') 'number of determinants combining fragment ',    &
                         iFrag-1,' and ',iFrag,' : ',newdets
    flush(lfndbg)
  endif
! create the new intermediate occupations and coefficients
  allocate(occmebf(newdets))
  allocate(coefmebf(newdets))
  newdets = 0
  first1  = 1
  last1   = 0
  do i = 1, spin1
    ms1    = spin1-1 - 2*(i-1)
    last1  = last1 + micro_ndets1(i)
    first2 = 1
    last2  = 0
    do j = 1, spin2
      ms2   = spin2-1 - 2*(j-1)
      last2 = last2 + micro_ndets2(j)
      if (ms1 + ms2 .eq. target_spin) then
        do idet = first1, last1
          do jdet = first2, last2
            prod = coef1(idet) * coef2(jdet)
            if (abs(prod) .ge. tau_CI) then
              dnorm = clebsch_gordon((spin1-1)/two,(spin2-1)/two, &
                                      ms1/two,ms2/two,            &
                                      target_spin/two,target_spin/two)
              newdets = newdets + 1
              coefmebf(newdets) = prod * dnorm
              occmebf(newdets) = trim(occ1(idet))//trim(occ2(jdet))
            endif
          end do
        end do
      end if
      first2 = first2 + micro_ndets2(j)
    end do
    first1 = first1 + micro_ndets1(i)
  end do
  if ( idbg .ge. 50 ) then
    write(lfndbg,'(a,10i2)') 'MEBF after considering the fragments ',(i,i=1,iFrag)
    do idet = 1, newdets
      write(lfndbg,'(I6,F15.8,3x,A)') idet,coefmebf(idet),trim(occmebf(idet))
    end do
    flush(lfndbg)
  endif
  deallocate(coef1,occ1,micro_ndets1,coef2,occ2,micro_ndets2)
! Copy coefmebf into coef for further processing if more fragments are
! to be added.
! Same for occupations and spin to mimic a new "fragment 1" 
  if (iFrag .ne. nmol ) then
    spin1 = inter_couplings(iFrag-1,iMEBF)
    spinFrag = spin1
    ndets1 = newdets
    allocate(coef(ndets1))
    allocate(occ(ndets1))
    allocate(micro_ndets1(spin1))
    do idet = 1, ndets1
      coef(idet) = coefmebf(idet)
      occ(idet) = occmebf(idet)
    end do
    deallocate(coefmebf,occmebf)
    do ms = 1, spin1
      micro_ndets1(ms) = 0
    end do
    if (spin1 .gt. 1) then
      call generate_microstates(ndets1,micro_dets1)
      allocate(coef1(micro_dets1))
      allocate(occ1(micro_dets1)) 
      do ms = 1, spin1
        micro_ndets1(ms) = micro_ndets(ms)
      end do
      do idet = 1, micro_dets1
        coef1(idet) = micro_coef(idet)
        occ1(idet) = micro_occ(idet)
      end do
      if ( idbg .ge. 50 ) then
        write(lfndbg,'(a)') 'micro states for the new fragment 1'
        write(lfndbg,'(a,i3)') 'generated by merging until ',iFrag
        do idet = 1, micro_dets1
          write(lfndbg,'(I6,F15.8,3x,A)') idet,coef1(idet),trim(occ1(idet))
        end do
        flush(lfndbg)
      endif
      deallocate(micro_ndets,micro_coef,micro_occ)
!     allocated in generate_microstates
    else
      micro_dets1 = ndets1
      allocate(coef1(micro_dets1))
      allocate(occ1(micro_dets1)) 
      micro_ndets1(1) = micro_dets1
      do idet = 1, micro_dets1
        coef1(idet) = coef(idet)
        occ1(idet) = occ(idet)
      end do
    end if
    deallocate(coef,occ)
  end if    
end do
! loop over the other fragments ends here and MEBF construction
! completes
! coefficients and occupations must be saved in civb and occb
if (first_pass .and. nmol .ne. 1) then
  if (iMEBF.eq.1) then
    maxcib = newdets
    idetb(1) = micro_dets1
  else
    maxcib = max(maxcib,newdets)
    idetb(iMEBF) = micro_dets1
  endif
  if ( idbg .ge. 50 ) then
    write(lfndbg,*) 'First pass: finished with MEBF ',iMEBF
    write(lfndbg,*) 'Maximum number of determinants so far: ',maxcib
    flush(lfndbg)
  endif
endif
if ( .not. first_pass .and. nmol .ne. 1) then
  call normalize(coefmebf,newdets)
  call quicksort_number(coefmebf,occmebf,newdets)
  if (iMEBF.eq.1) then
    maxcoef = abs(coefmebf(1))
  else
    maxcoef = max(maxcoef,abs(coefmebf(1)))
  endif
  if ( idbg .ge. 50 ) then  
    write(lfndbg,*) 'Max abs coef: ',maxcoef
    write(lfndbg,'(a,i3)') 'civb and occupations for MEBF ',iMEBF
    flush(lfndbg)
  endif
  allocate(occ_num(nactb(iMEBF)))
  do idet = 1, newdets
    dumstr = occmebf(idet)
    do iAct = 1, nactb(iMEBF)
      if (dumstr(iAct:iAct) .eq. '2') occ_num(iAct) = 2
      if (dumstr(iAct:iAct) .eq. 'a') occ_num(iAct) = 1
      if (dumstr(iAct:iAct) .eq. 'b') occ_num(iAct) =-1
      if (dumstr(iAct:iAct) .eq. '0') occ_num(iAct) = 0
    end do 
    coefmebf(idet) = coefmebf(idet) * perm_ab(occ_num,nactb(iMEBF))
    coefmebf(idet) = coefmebf(idet) * isetsign(occ_num,nactb(iMEBF))
    civb(idet,iMEBF) = coefmebf(idet)
    call pack(ioccb(idet,iMEBF),occ_num,nactb(iMEBF))
    if ( idbg .ge. 50 ) then
      write(lfndbg,'(i6,f14.8,32i3)')idet,civb(idet,iMEBF),(occ_num(iAct),iAct=1,nactb(iMEBF))
      flush(lfndbg)
    endif 
  end do
  deallocate(occ_num,coefmebf,occmebf)
end if

if ( idbg .ge. 50 ) then
  write(lfndbg,*) 'Allocated in makebasestate, should all be deallocated'
  write(lfndbg,*) 'allocated micro_dets1 ?',allocated(micro_ndets1)
  write(lfndbg,*) 'allocated micro_dets2 ?',allocated(micro_ndets2)
  write(lfndbg,*) 'allocated occ_num     ?',allocated(occ_num)
  write(lfndbg,*) 'allocated coef1       ?',allocated(coef1)
  write(lfndbg,*) 'allocated coef2       ?',allocated(coef2)
  write(lfndbg,*) 'allocated occ1        ?',allocated(occ1)
  write(lfndbg,*) 'allocated occ2        ?',allocated(occ2)
  write(lfndbg,*) 'allocated coefmebf    ?',allocated(coefmebf)
  write(lfndbg,*) 'allocated occmebf     ?',allocated(occmebf)
  write(lfndbg,*) 'Allocated in gronor_main, should be allocated in the second pass'
  write(lfndbg,*) 'allocated civb        ?',allocated(civb),size(civb,1),size(civb,2)
  write(lfndbg,*) 'allocated ioccb       ?',allocated(ioccb),size(ioccb,1),size(ioccb,2)
 
  
  write(lfndbg,*) 'Allocated in generate_microstates, should be deallocated'
  write(lfndbg,*) 'Allocated micro_coef  ?',allocated(micro_coef),size(micro_coef)
  write(lfndbg,*) 'Allocated micro_occ   ?',allocated(micro_occ),size(micro_occ)
  write(lfndbg,*) 'Allocated micro_ndets ?',allocated(micro_ndets),size(micro_ndets)
  flush(lfndbg)
endif




end subroutine gronor_make_basestate


