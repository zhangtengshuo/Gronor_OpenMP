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
!! Print results to the output file
!!
!! @author  T. P. Straatsma, ORNL
!! @author  C. de Graaf, URV
!! @date    2020
!!

subroutine gronor_print_results(hb)

  use gnome_parameters, only : itest,ncols
  use cidef

  implicit none

  external :: close_tag,open_tag
  external :: writetag_matrix_real,writetag_array_real
  external :: writetag_scalar_real,writetag_scalar_integer
  external :: gronor_print_matrix,dsygv

  real(kind=8)                  :: hb(nbase,nbase)
  integer                       :: i,j,k,ii,il
  integer                       :: nk,info,lwork
  real(kind=8),allocatable      :: work(:)

  character(len=132)            :: info_cml,info_cml2,fmt_1,id
  character(len=20)             :: label
  character(len=1)              :: sep
  integer                       :: indent

  real(kind=8)                  :: debye,angstrom
  integer                       :: lfnt

  real(kind=8), allocatable :: tc(:,:),nociovlp(:,:)
  real(kind=8), allocatable :: nociwf0(:)

  allocate(tc(nbase,nbase),nociovlp(nbase,nbase))

  debye    = 2.54174644986
  angstrom = 0.529177249

  lfnt=0
  if(itest.gt.0) lfnt=lfntst

  key='Hamil'
  header='Hamiltonian Matrix'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      mebfLabels,mebfLabel,.false.,1.0d0,hb,nbase,7,6)

  key='Overl'
  header='Overlap Matrix'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      mebfLabels,mebfLabel,.false.,1.0d0,sbase,nbase,7,6)

  !      allocate( nociwf(nbase,nbase) )
  if (nbase .eq. 1) then
    hev(1)=hb(1,1)/sbase(1,1)
    nociwf(1,1) = 1.0d0
  else
    allocate(nociwf0(nbase))
    do i=1,nbase
      do j=1,nbase
        tc(i,j)=0.0d0
      enddo
    enddo
    do i=2,nbase
      do j=1,i-1
        tc(i,j)=(hb(i,j)-0.5d0*(hb(i,i)+hb(j,j))*sbase(i,j))/(1.0d0-sbase(i,j)*sbase(i,j))
        tc(j,i)=tc(i,j)
      enddo
    enddo

    key='Coupl'
    header='Electronic Couplings (meV)'
    call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
        mebfLabels,mebfLabel,.true.,27.2114d3,tc,nbase,7,3)

    !     Diagonalize the NOCI matrix and print the energies and wave functions
    !     Don't touch hb and sbase, needed again when shifting

    info=0
    lwork = 4*nbase
    do j=1,nbase
      do i=1,nbase
        nociwf(i,j)=hb(i,j)
        nociovlp(i,j)=sbase(i,j)
      enddo
    enddo
    allocate(work(lwork))
    call dsygv(1,'V','L',nbase,nociwf,nbase,nociovlp,nbase,hev,work,lwork,info)
    deallocate(work)
    if ( info .ne. 0 ) write(*,*) 'something went wrong in dyegv'

    ! Coen 2020/07/18
    ! Check if the eigenvectors come in rows or columns!
    ! a transpose(nociwf) might be needed
    ! printing takes care of it now: nociwf(j,i) instead of nociwf(i,j)

    write(lfnout,677)
677 format(//,' NOCI energies and wave functions')
    write(lfnarx,667) nbase,ncols
    write(lfnxrx,667) nbase,ncols
667 format('NOCI ',2i10)

    nk=nbase/ncols
    if(mod(nbase,ncols).ne.0) nk=nk+1
    do k=1,nk
      ii=(k-1)*ncols+1
      il=min(nbase,k*ncols)

      write(lfnout,680) (i,i=ii,il)
680   format(/,'          State:',t18,12(i15,5x))
      write(lfnarx,668) (i,i=ii,il)
      write(lfnxrx,668) (i,i=ii,il)
668   format(7(6x,i8,6x))
      write(lfnout,681) (hev(i),i=ii,il)
681   format('    Energy (Eh):',t18,12f20.10)
      write(lfnarx,669) (hev(i),i=ii,il)
      write(lfnxrx,669) (hev(i),i=ii,il)
669   format(10e20.13)
      write(lfnout,682) (27.2114d0*(hev(i)-hev(1)),i=ii,il)
682   format('  Relative (eV):',t18,12f20.10)
      write(lfnarx,669) (27.2114d0*(hev(i)-hev(1)),i=ii,il)
      write(lfnxrx,669) (27.2114d0*(hev(i)-hev(1)),i=ii,il)
      write(lfnout,683)
683   format(1x)
      do j=1,nbase
        if (.not.mebfLabels) then
          write(lfnout,684) j,(nociwf(j,i),i=ii,il)
684       format(i14,t18,12f20.10)
        else
          if(lablen.le.labmax) then
            write(lfnout,685) trim(mebfLabel(j)),(nociwf(j,i),i=ii,il)
685         format(1x,a,t26,12f20.10)
          else
            write(lfnout,686) j,(nociwf(j,i),i=ii,il)
686         format(1x,i5,2x,12f20.10)
          endif
        endif
        write(lfnarx,655) (nociwf(j,i),i=ii,il)
        write(lfnxrx,655) (nociwf(j,i),i=ii,il)
655     format(12f20.10)
        if(itest.eq.2) then
          do i=ii,il
            nociwf0(i)=nociwf(j,i)
            if(abs(nociwf0(i)).lt.1.0d-6) nociwf0(i)=0.0d0            
          enddo
          write(lfntst,687) j,(nociwf0(i),i=ii,il)
687       format(i14,t18,12f20.6)
        endif
      enddo
    enddo
    deallocate(nociwf0)
  endif

  !   Dumping the results in the cml file
  if (ncorr.eq.1) then
    if (nwt.eq.1) then
      id = 'id="GNshifted"'
    else
      id = 'id="shifted"'
    end if
  else
    id = 'id="standard"'
  end if
  ! 1. NOCI matrices
  label = 'module'
  indent = 4
  info_cml = 'dictRef="gr:matrices" id="normalized"'
  call open_tag(lfncml,label,info_cml,indent)
  label = 'propertyList'
  indent= 5
  info_cml = 'empty'
  call open_tag(lfncml,label,info_cml,indent)
  label = 'property'
  info_cml='dictRef="gr:hamiltonian" '//trim(id)
  call open_tag(lfncml,label,info_cml,6)
  info_cml='units="nonsi:hartree"'
  sep='|'
  fmt_1='f20.10'
  call writetag_matrix_real(lfncml,info_cml,7,sep,hb,size(hb,1),size(hb,2),fmt_1)
  call close_tag(lfncml,label,6)
  if (ncorr.ne.1) then
    info_cml='dictRef="gr:overlap"'
    indent = 6
    call open_tag(lfncml,label,info_cml,indent)
    info_cml='empty'
    call writetag_matrix_real(lfncml,info_cml,7,sep,sbase,size(hb,1),size(hb,2),fmt_1)
    call close_tag(lfncml,label,6)
  end if
  label = 'propertyList'
  indent= 5
  call close_tag(lfncml,label,indent)
  label='module'
  indent = 4
  call close_tag(lfncml,label,indent)
  ! 2. NOCI wf and energy (only when thre is more than one MEBF)
  if (nbase.gt.1) then
    indent = 4
    label = 'module'
    info_cml = 'dictRef="gr:nociStates"'
    call open_tag(lfncml,label,info_cml,indent)
    do i = 1, nbase
      indent= 5
      info_cml = 'dictRef="gr:nociroot"'
      call open_tag(lfncml,label,info_cml,indent)
      label = 'propertyList'
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,6)
      label = 'property'
      info_cml = 'dictRef="gr:rootNumber"'
      call open_tag(lfncml,label,info_cml,7)
      info_cml='empty'
      call writetag_scalar_integer(lfncml,info_cml,8,i)
      call close_tag(lfncml,label,7)
      info_cml = 'dictRef="gr:nociEnergy" '//trim(id)
      call open_tag(lfncml,label,info_cml,7)
      fmt_1='f16.8'
      info_cml = 'units="nonsi:hartree"'
      call writetag_scalar_real(lfncml,info_cml,8,hev(i),fmt_1)
      call close_tag(lfncml,label,7)
      info_cml = 'dictRef="gr:nociRelEnergy" '//trim(id)
      call open_tag(lfncml,label,info_cml,7)
      fmt_1='f12.4'
      info_cml2 = trim(info_cml)//' units="nonsi:electronvolt"'
      call writetag_scalar_real(lfncml,info_cml2,8,(hev(i)-hev(1))*27.2114,fmt_1)
      call close_tag(lfncml,label,7)
      info_cml = 'dictRef="gr:mebfCoeff" '//trim(id)
      call open_tag(lfncml,label,info_cml,7)
      fmt_1 = 'f20.10'
      call writetag_array_real(lfncml,info_cml,8,sep,nociwf(i,:),size(nociwf,2),fmt_1)
      call close_tag(lfncml,label,7)
      label = 'propertyList'
      call close_tag(lfncml,label,6)
      label = 'module'
      call close_tag(lfncml,label,5)
    end do
    indent = 4
    call close_tag(lfncml,label,indent)
    !3. Electronic couplings (only when there is more than one MEBF)
    do i=1,size(tc,1)
      do j=1,size(tc,2)
        tc(i,j)=tc(i,j)*27211.4d0
      end do
    end do
    indent = 4
    info_cml = 'dictRef="gr:elCoupling"'
    call open_tag(lfncml,label,info_cml,indent)
    label = 'propertyList'
    indent= 5
    info_cml = 'empty'
    call open_tag(lfncml,label,info_cml,indent)
    label = 'property'
    info_cml = 'dictRef="gr:mebfCoupling" '//trim(id)
    call open_tag(lfncml,label,info_cml,6)
    info_cml='units="nonsi2:meV"'
    fmt_1='f12.4'
    call writetag_matrix_real(lfncml,info_cml,7,sep,tc,size(tc,1),size(tc,2),fmt_1)
    label = 'property'
    call close_tag(lfncml,label,6)
    label = 'propertyList'
    indent= 5
    call close_tag(lfncml,label,indent)
    label='module'
    indent = 4
    call close_tag(lfncml,label,indent)
  end if

  flush(lfnarx)
  flush(lfnxrx)

  deallocate(tc,nociovlp)
  !      deallocate(nociwf)

  return

end subroutine gronor_print_results
