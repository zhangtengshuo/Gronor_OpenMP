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

      subroutine gronor_print_results(hbase,sbase,nbase)
      
      use gnome_parameters, only : itest,ncols
      use cidef           , only : lfnout,lfntst,lfnarx,lfncml,ncorr,
     &                             mebfLabel,mebfLabels

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
*
      character(len=132)            :: info_cml
      character(len=20)             :: label
      character(len=1)              :: sep
      integer                       :: indent
      logical                       :: scientific_format

      write(lfnarx,660) nbase,ncols
 660  format('Hamil',2i10)
      write(*,670)
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,670)
        write(lfntst,670)
 670    format(//,' Hamiltonian Matrix',/)
        write(*,671) (i,i=ii,il)
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
          write(*,672) j,(hbase(i,j),i=ii,il)
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
 1672     format(i5,1x,10f20.8)
        enddo
        write(lfnout,*) ' '
        write(*,*) ' '
      enddo
      write(lfnarx,663) nbase,ncols
 663  format('Overl',2i10)
      write(*,673)
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,673)
        write(lfntst,673)
 673    format(//,' Overlap Matrix',/)
        write(*,671) (i,i=ii,il)
        if (.not.mebfLabels) then
          write(lfnout,671) (i,i=ii,il)
        else
          write(lfnout,6711) '|',(mebfLabel(i),i=ii,il)
          write(lfnout,6712)'-------|',('--------------------',i=ii,il)
        end if
        write(lfntst,671) (i,i=ii,il)
        write(lfnarx,671) (i,i=ii,il)
        do j=1,nbase
          write(*,672) j,(sbase(i,j),i=ii,il)
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
        write(*,*) ' '
      enddo
      
      if(nbase.gt.1) then
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
        write(*,676)
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
          write(*,671) (i,i=ii,il)
          do j=1,nbase
            if(j.eq.1) then
              if (.not.mebfLabels) then
                write(lfnout,672) 1
              else
                write(lfnout,6721) mebfLabel(1),' |'
              end if
              write(lfnarx,662) 1
              write(*,672) 1
            else
              ik=min(j-1,il)
              if (.not.mebfLabels) then
                write(lfnout,672) j,(27.2114d3*tc(j,i),i=ii,ik)
              else
                write(lfnout,6721) mebfLabel(j),' |',
     &                             (27.2114d3*tc(j,i),i=ii,ik)
              end if
              write(lfnarx,662) j,(27.2114d3*tc(j,i),i=ii,ik)
              write(*,672) j,(27.2114d3*tc(j,i),i=ii,ik)
            endif
          enddo
          write(lfnout,*) ' '
          write(*,*) ' '
        enddo

!     Diagonalize the NOCI matrix and print the energies and wave functions
!     Don't touch hbase and sbase, needed again when shifting
          
        info=0
        lwork = 4 * nbase
        hsave = hbase
        ssave = sbase
        allocate(work(lwork))
        call dsygv(1,'V','L',nbase,hsave,nbase,
     &                   ssave,nbase,hev,work,lwork,info)
        deallocate(work)
        if ( info .ne. 0 ) write(*,*) 'something went wrong in dyegv'
        
! Coen 2020/07/18
! Check if the eigenvectors come in rows or columns!
! a transpose(hsave) might be needed
! printing takes care of it now: hsave(j,i) instead of hsave(i,j)
          
        write(*,677)
        write(lfnout,677)
 677    format(//,' NOCI energies and wave functions (columns)',/)

        nk=nbase/ncols
        if(mod(nbase,ncols).ne.0) nk=nk+1
        do k=1,nk
          ii=(k-1)*ncols+1
          il=min(nbase,k*ncols)

          write(lfnout,680) (i,i=ii,il)
          write(*,680) (i,i=ii,il)
 680      format('          State:',t18,12(i15,5x))
          write(lfnout,681) (hev(i),i=ii,il)
          write(*,681) (hev(i),i=ii,il)
 681      format('    Energy (Eh):',t18,12f20.10)
          write(lfnout,682) (27.2114d0*(hev(i)-hev(1)),i=ii,il)
          write(*,682) (27.2114d0*(hev(i)-hev(1)),i=ii,il)
 682      format('  Relative (eV):',t18,12f20.10)
          token="MEBF:"
          do j=1,nbase
            if(j.eq.1) then
              write(*,683)token,j,(hsave(j,i),i=ii,il)
              if (.not.mebfLabels) then
                write(lfnout,683)token,j,(hsave(j,i),i=ii,il)
              else
                write(lfnout,685)token,mebfLabel(j),(hsave(j,i),i=ii,il)
              end if
            else
              write(*,684) token,j,(hsave(j,i),i=ii,il)
              if (.not.mebfLabels) then
                write(lfnout,684)token,j,(hsave(j,i),i=ii,il)
              else
                write(lfnout,686)token,mebfLabel(j),(hsave(j,i),i=ii,il)
              end if
            endif
 683        format(/,4x,a5,i5,t18,12f20.10)
 684        format(4x,a5,i5,t18,12f20.10)
 685        format(/,4x,a5,1x,a4,t18,12f20.10)
 686        format(4x,a5,1x,a4,t18,12f20.10)
            if(itest.eq.2) then
              write(lfntst,687) token,j,(hsave(j,i),i=ii,il)
 687          format(6x,a4,i5,t18,12f20.8)
            endif
            token="    "
          enddo
        enddo
      endif

*** Dumping the results in the cml file
      label = 'propertyList'
      indent= 3
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,indent)
      label = 'property'
      info_cml='dictRef="gr:Hamiltonian"'
      call open_tag(lfncml,label,info_cml,4)
      info_cml='units="nonsi:hartree"'
      sep='|'
      scientific_format=.true.
      call writetag_matrix_real(lfncml,info_cml,5,sep,hbase,
     &   size(hbase,1),size(hbase,2),scientific_format)
      call close_tag(lfncml,label,4)
 
      return
      
      end subroutine gronor_print_results
