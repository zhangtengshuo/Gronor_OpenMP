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

      subroutine gronor_print_results(hbase,sbase,nbase,hev)
      
      use gnome_parameters, only : itest,ncols
      use cidef           , only : lfnout,lfntst,lfnarx,lfncml,ncorr,
     &                             nwt,mebfLabel,mebfLabels

      implicit none

      integer,intent(in)            :: nbase
      integer                       :: i,j,k,ii,ik,il
      integer                       :: nk,info,lwork
      real(kind=8),intent(in)       :: hbase(nbase,nbase)
      real(kind=8),intent(in)       :: sbase(nbase,nbase)
      real(kind=8)                  :: hbt(nbase),tc(nbase,nbase)
      real(kind=8)                  :: hsave(nbase,nbase)
      real(kind=8)                  :: ssave(nbase,nbase)
      real(kind=8)                  :: hev(nbase)
      real(kind=8),allocatable      :: work(:)
      character(len=5)              :: token

      character(len=132)            :: info_cml,info_cml2,fmt_1,id
      character(len=20)             :: label
      character(len=1)              :: sep
      integer                       :: indent
      
      write(lfnarx,660) nbase,ncols
 660  format('Hamil',2i10)
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,670)
        write(lfntst,670)
 670    format(//,' Hamiltonian Matrix',/)
        if (.not.mebfLabels) then
          write(lfnout,671) (i,i=ii,il)
        else
          write(lfnout,6711) '|',(mebfLabel(i),i=ii,il)
          write(lfnout,6712)'-------|',('--------------------',i=ii,il)
        end if
        write(lfntst,671) (i,i=ii,il)
        write(lfnarx,671) (i,i=ii,il)
 671    format(6x,7(6x,i8,6x))
 6711   format(7x,a,7(9x,a6,5x))
 6712   format(8a)
        do j=1,nbase
          if (.not.mebfLabels) then
            write(lfnout,672) j,(hbase(i,j),i=ii,il)
          else
            write(lfnout,6721) mebfLabel(j),' |',(hbase(i,j),i=ii,il)
          end if
          write(lfnarx,662) j,(hbase(i,j),i=ii,il)
 672      format(i5,1x,10f20.10)
 6721     format(a6,a2,10f20.10)
 662      format(i5,1x,10e20.13)
          do i=ii,il
            hbt(i)=hbase(i,j)
            if(dabs(hbt(i)).lt.1.0d-08) hbt(i)=0.0d0
          enddo
          write(lfntst,1672) j,(hbt(i),i=ii,il)
 1672     format(i5,1x,10f20.6)
        enddo
        write(lfnout,*) ' '
      enddo
      write(lfnarx,663) nbase,ncols
 663  format('Overl',2i10)
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,673)
        write(lfntst,673)
 673    format(//,' Overlap Matrix',/)
        if (.not.mebfLabels) then
          write(lfnout,671) (i,i=ii,il)
        else
          write(lfnout,6711) '|',(mebfLabel(i),i=ii,il)
          write(lfnout,6712)'-------|',('--------------------',i=ii,il)
        end if
        write(lfntst,671) (i,i=ii,il)
        write(lfnarx,671) (i,i=ii,il)
        do j=1,nbase
          if (.not.mebfLabels) then
            write(lfnout,672) j,(sbase(i,j),i=ii,il)
          else
            write(lfnout,6721) mebfLabel(j),' |',(sbase(i,j),i=ii,il)
          end if
          write(lfnarx,662) j,(sbase(i,j),i=ii,il)
          do i=ii,il
            hbt(i)=sbase(i,j)
            if(dabs(hbt(i)).lt.1.0d-08) hbt(i)=0.0d0
          enddo
          write(lfntst,1672) j,(hbt(i),i=ii,il)
        enddo
        write(lfnout,*) ' '
      enddo
      
      if(nbase.gt.1) then
        tc = 0.0
        do i=2,nbase
          do j=1,i-1
            tc(i,j)=(hbase(i,j)
     &           -0.5d0*(hbase(i,i)+hbase(j,j))*sbase(i,j))/
     &           (1.0d0-sbase(i,j)*sbase(i,j))
          enddo
        enddo

        write(lfnout,676)
        write(lfnarx,666) nbase,ncols
 666    format('Coupl',2i10)
 676    format(//,' Electronic Couplings (meV)',/)
        nk=nbase/ncols
        if(mod(nbase,ncols).ne.0) nk=nk+1
        do k=1,nk
          ii=(k-1)*ncols+1
          il=min(nbase,k*ncols)
          if (.not.mebfLabels) then
            write(lfnout,671) (i,i=ii,il)
          else
            write(lfnout,6711)  '|',(mebfLabel(i),i=ii,il)
            write(lfnout,6712)'-------|',
     &                     ('--------------------',i=ii,il)
          end if
          write(lfnarx,671) (i,i=ii,il)
          do j=1,nbase
            if(j.eq.1) then
              if (.not.mebfLabels) then
                write(lfnout,672) 1
              else
                write(lfnout,6721) mebfLabel(1),' |'
              end if
              write(lfnarx,662) 1
            else
              ik=min(j-1,il)
              if (.not.mebfLabels) then
                write(lfnout,672) j,(27.2114d3*tc(j,i),i=ii,ik)
              else
                write(lfnout,6721) mebfLabel(j),' |',
     &                             (27.2114d3*tc(j,i),i=ii,ik)
              end if
              write(lfnarx,662) j,(27.2114d3*tc(j,i),i=ii,ik)
            endif
          enddo
          write(lfnout,*) ' '
        enddo

!     Diagonalize the NOCI matrix and print the energies and wave functions
!     Don't touch hbase and sbase, needed again when shifting

        if(nbase.gt.1) then
          info=0
          lwork = 4 * nbase
          do j=1,nbase
            do i=1,nbase
              hsave(i,j)=hbase(i,j)
              ssave(i,j)=sbase(i,j)
            enddo
          enddo
          allocate(work(lwork))
          call dsygv(1,'V','L',nbase,hsave,nbase,
     &         ssave,nbase,hev,work,lwork,info)
          deallocate(work)
          if ( info .ne. 0 ) write(*,*) 'something went wrong in dyegv'
        else
          hev(1)=hbase(1,1)/sbase(1,1)
        endif
        
! Coen 2020/07/18
! Check if the eigenvectors come in rows or columns!
! a transpose(hsave) might be needed
! printing takes care of it now: hsave(j,i) instead of hsave(i,j)
          
        write(lfnout,677)
 677    format(//,' NOCI energies and wave functions')
        write(lfnarx,667) nbase,ncols
 667    format('NOCI ',2i10)

        nk=nbase/ncols
        if(mod(nbase,ncols).ne.0) nk=nk+1
        do k=1,nk
          ii=(k-1)*ncols+1
          il=min(nbase,k*ncols)

          write(lfnout,680) (i,i=ii,il)
 680      format(/,'          State:',t18,12(i15,5x))
          write(lfnarx,668) (i,i=ii,il)
 668      format(7(6x,i8,6x))
          write(lfnout,681) (hev(i),i=ii,il)
 681      format('    Energy (Eh):',t18,12f20.10)
          write(lfnarx,669) (hev(i),i=ii,il)
 669      format(10e20.13)
          write(lfnout,682) (27.2114d0*(hev(i)-hev(1)),i=ii,il)
 682      format('  Relative (eV):',t18,12f20.10)
          write(lfnarx,669) (27.2114d0*(hev(i)-hev(1)),i=ii,il)
          token="MEBF:"
          do j=1,nbase
             if(j.eq.1) then
                if (.not.mebfLabels) then
                   write(lfnout,683)token,j,(hsave(j,i),i=ii,il)
                else
                   write(lfnout,685)token,mebfLabel(j),
     &                  (hsave(j,i),i=ii,il)
                endif
             else
                if (.not.mebfLabels) then
                   write(lfnout,684)token,j,(hsave(j,i),i=ii,il)
                else
                   write(lfnout,686)token,mebfLabel(j),
     &                  (hsave(j,i),i=ii,il)
                endif
             endif
            write(lfnarx,655) (hsave(j,i),i=ii,il)
 655        format(12f20.10)
 683        format(/,4x,a5,i5,t18,12f20.10)
 684        format(4x,a5,i5,t18,12f20.10)
 685        format(/,4x,a5,1x,a4,t18,12f20.10)
 686        format(4x,a5,1x,a4,t18,12f20.10)
            if(itest.eq.2) then
              write(lfntst,687) token,j,(hsave(j,i),i=ii,il)
 687          format(6x,a4,i5,t18,12f20.6)
            endif
            token="    "
          enddo
        enddo
      endif

*** Dumping the results in the cml file
      if (ncorr.eq.1) then
        if (nwt.eq.1) then
          id = 'id="GNshifted"'
        else
          id = 'id="shifted"'
        end if
      else
        id = 'id="standard"'
      end if
* 1. NOCI matrices
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
      call writetag_matrix_real(lfncml,info_cml,7,sep,hbase,
     &   size(hbase,1),size(hbase,2),fmt_1)
      call close_tag(lfncml,label,6)
      if (ncorr.ne.1) then
        info_cml='dictRef="gr:overlap"'
        indent = 6
        call open_tag(lfncml,label,info_cml,indent) 
        info_cml='empty'
        call writetag_matrix_real(lfncml,info_cml,7,sep,sbase,
     &     size(hbase,1),size(hbase,2),fmt_1)
        call close_tag(lfncml,label,6)
      end if
      label = 'propertyList'
      indent= 5
      call close_tag(lfncml,label,indent)
      label='module'
      indent = 4
      call close_tag(lfncml,label,indent)
* 2. NOCI wf and energy (only when thre is more than one MEBF)
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
          call writetag_scalar_real(lfncml,info_cml2,8,(hev(i)-hev(1))
     &                                                   *27.2114,fmt_1)
          call close_tag(lfncml,label,7)
          info_cml = 'dictRef="gr:mebfCoeff" '//trim(id)
          call open_tag(lfncml,label,info_cml,7)
          fmt_1 = 'f20.10'
          call writetag_array_real(lfncml,info_cml,8,sep,hsave(i,:),
     &              size(hsave,2),fmt_1)
          call close_tag(lfncml,label,7)
          label = 'propertyList'
          call close_tag(lfncml,label,6)
          label = 'module'
          call close_tag(lfncml,label,5)
        end do
        indent = 4
        call close_tag(lfncml,label,indent)
*3. Electronic couplings (only when there is more than one MEBF)
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
        call writetag_matrix_real(lfncml,info_cml,7,sep,tc,
     &              size(tc,1),size(tc,2),fmt_1)
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
      
      return
      
      end subroutine gronor_print_results
