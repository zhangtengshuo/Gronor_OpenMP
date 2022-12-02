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

      subroutine gronor_print_results(hbase,sbase,mbase,mnuc,nbase,hev)

      use gnome_parameters, only : itest,ncols
      use gnome_data      , only : CoM
      use cidef           , only : lfnout,lfntst,lfnarx,lfncml,ncorr,           &
     &                             nwt,mebfLabel,mebfLabels

      implicit none

      integer,intent(in)            :: nbase
      integer                       :: i,j,k,ii,ik,il,l,m
      integer                       :: nk,info,lwork
      real(kind=8),intent(in)       :: hbase(nbase,nbase)
      real(kind=8),intent(in)       :: sbase(nbase,nbase)
      real(kind=8),intent(in)       :: mnuc(9)
      real(kind=8),intent(inout)    :: mbase(nbase,nbase,9)
      real(kind=8),allocatable      :: noci(:,:)
      real(kind=8),allocatable      :: mnoci(:,:,:)
      real(kind=8),allocatable      :: obase(:,:)
      real(kind=8),allocatable      :: onoci(:,:)
      real(kind=8),allocatable      :: qmtraceless(:,:,:)
      real(kind=8)                  :: Debye,Angstrom,trace      
      
      real(kind=8)                  :: hbt(nbase),tc(nbase,nbase)
      real(kind=8)                  :: hsave(nbase,nbase)
      real(kind=8)                  :: ssave(nbase,nbase)
      real(kind=8)                  :: hev(nbase)
      real(kind=8),allocatable      :: work(:)
      real(kind=8)                  :: Two3rds
      character(len=5)              :: token

      character(len=132)            :: info_cml,info_cml2,fmt_1,id
      character(len=20)             :: label
      character(len=1)              :: sep
      integer                       :: indent

!     ------- Constants -------
      Debye    = 2.54174644986
      Angstrom = 0.529177249
      Two3rds  = 2.0d0/3.0d0
!     ------------------------
      
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
            tc(i,j)=(hbase(i,j)                                                 &
     &           -0.5d0*(hbase(i,i)+hbase(j,j))*sbase(i,j))/                    &
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
            write(lfnout,6712)'-------|',                                       &
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
                write(lfnout,6721) mebfLabel(j),' |',                           &
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
          call dsygv(1,'V','L',nbase,hsave,nbase,                               &
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
                   write(lfnout,685)token,mebfLabel(j),                         &
     &                  (hsave(j,i),i=ii,il)
                endif
             else
                if (.not.mebfLabels) then
                   write(lfnout,684)token,j,(hsave(j,i),i=ii,il)
                else
                   write(lfnout,686)token,mebfLabel(j),                         &
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
      
      allocate(noci(nbase,nbase))
      allocate(mnoci(nbase,nbase,9))
      allocate(obase(nbase,nbase))
      allocate(onoci(nbase,nbase))
      allocate(qmtraceless(nbase,nbase,6))
      
      do i=1,nbase
        do j=1,nbase
          noci(i,j)=0.0d0
          obase(i,j)=0.0d0
          onoci(i,j)=0.0d0
          do k=1,9
            mnoci(i,j,k)=0.0d0
          enddo
          do k=1,6
            qmtraceless(i,j,k)=0.0d0
          enddo
        enddo
      enddo
      
      if(nbase.eq.1)then
        noci(1,1)=1.0d0
      else
        do i=1,nbase
          do j=1,nbase
            noci(j,i)=hsave(j,i)            
          enddo
        enddo
      endif
      
      do i=1,9
        do j=1,nbase
          do k=1,nbase            
            do l=1,nbase
              do m=1,nbase
                mnoci(k,j,i)=mnoci(k,j,i)+                              &
     &               noci(l,j)*noci(m,k)*mbase(m,l,i)                
              enddo              
            enddo
          enddo
        enddo
      enddo

      do j=1,nbase
        do k=1,j          
          obase(k,j)=obase(k,j)+Two3rds*(hev(j)-hev(k))*                &
     &         (mbase(k,j,1)**2+mbase(k,j,2)**2+                        &
     &         mbase(k,j,3)**2)
        enddo        
      enddo
      
      do j=1,nbase
        do k=1,j          
          onoci(k,j)=onoci(k,j)+Two3rds*(hev(j)-hev(k))*                &
     &         (mnoci(k,j,1)**2+mnoci(k,j,2)**2+                        &
     &         mnoci(k,j,3)**2)
        enddo        
      enddo
      
      do j=1,nbase
        do k=1,j
          obase(j,k)=-obase(k,j)
        enddo
      enddo

      do j=1,nbase
        do k=1,j
          onoci(j,k)=-onoci(k,j)
        enddo
      enddo

      do i=1,9
        do j=1,nbase
          do k=1,nbase
            mbase(k,j,i)=mbase(k,j,i)+mnuc(i)
            mnoci(k,j,i)=mnoci(k,j,i)+mnuc(i)
          enddo
        enddo        
      enddo
 688  format(//,'      ------- Multipoles -------',//)
 689  format(//,' Origin of the dipole operator (x,y,z)',//,7x,3f10.4,//        &
     &     ,' Origin of the quadrupole operator (x,y,z)',//,7x,3f10.4)          
 690  format(//,' Dipole moment of MEBFs (Debye)',/)
 691  format(7x,i5,2x,' X:',2x,f10.4,4x,' Y:',2x,f10.4,4x,                      &
     &     ' Z:',2x,f10.4,4x,' Total:',2x,f10.4,/)
 692  format(7x,a5,2x,' X:',2x,f10.4,4x,' Y:',2x,f10.4,4x,                      &
     &     ' Z:',2x,f10.4,4x,' Total:',2x,f10.4,/)
      write(lfnout,688)
      write(lfnout,689) 0.0d0,0.0d0,0.0d0,(CoM(j),j=1,3)
      write(lfnout,690)
      do j=1,nbase
        if(.not.mebfLabels)then
          write(lfnout,691) j,(Debye*mbase(j,j,k),k=1,3),Debye*         &
     &         sqrt(mbase(j,j,1)**2+mbase(j,j,2)**2+mbase(j,j,3)**2)
        else
          write(lfnout,692) mebfLabel(j),                               &
     &         (Debye*mbase(j,j,k),k=1,3),Debye*                        &
     &         sqrt(mbase(j,j,1)**2+mbase(j,j,2)**2+mbase(j,j,3)**2)
        endif        
      enddo
 694  format(//,' Dipole moment of NOCI states (Debye)',/)
      write(lfnout,694)
      do j=1,nbase
        write(lfnout,691) j,(Debye*mnoci(j,j,k),k=1,3),                 &
     &       Debye*sqrt(mnoci(j,j,1)**2+mnoci(j,j,2)**2+mnoci(j,j,3)**2)        
      enddo
      
 695  format(//,' Quadrupole moment of MEBFs (Debye*Angstrom)',/)
 696  format(//,' In traceless form',/)
      write(lfnout,695)
 697  format(7x,i5,2x,                                                  &
     &     ' XX:',2x,f10.4,4x,' XY:',2x,f10.4,4x,' XZ:',2x,f10.4)      
 698  format(7x,a5,2x,                                                  &
     &     ' XX:',2x,f10.4,4x,' XY:',2x,f10.4,4x,' XZ:',2x,f10.4)       
 699  format(14x,                                                       &
     &     ' YY:',2x,f10.4,4x,' YZ:',2x,f10.4,4x,' ZZ:',2x,f10.4,/)       
      do j=1,nbase        
        if(.not.mebfLabels)then
          write(lfnout,697) j,(Debye*Angstrom*mbase(j,j,k),k=4,6)
        else
          write(lfnout,698) mebfLabel(j),(Debye*Angstrom*mbase(j,j,k),  &
     &         k=4,6)
        endif        
        write(lfnout,699) (Debye * Angstrom * mbase(j,j,k),k=7,9)
      enddo      
      write(lfnout,696)
      do j=1,nbase        
        trace=(mbase(j,j,4)+mbase(j,j,7)+mbase(j,j,9))/3
        do k=1,6
          if(k.eq.1.or.k.eq.4.or.k.eq.6)then
            qmtraceless(j,j,k)=mbase(j,j,k+3)-trace
          else
            qmtraceless(j,j,k)=mbase(j,j,k)
          endif
        enddo
        
        if(.not.mebfLabels)then
          write(lfnout,697) j,(Debye*Angstrom*(qmtraceless(j,j,k)),     &
     &         k=1,3)
        else
          write(lfnout,698) mebfLabel(j),(Debye*Angstrom*               &
     &         (qmtraceless(j,j,k)),k=1,3)
        endif
        write(lfnout,699) (Debye*Angstrom*(qmtraceless(j,j,k)),         &
     &       k=4,6)
      enddo
     
 700  format(//,' Quadrupole moment of NOCI states (Debye*Angstrom)',/)
      write(lfnout,700)
      do j=1,nbase
        write(lfnout,697) j,(Debye * Angstrom * mnoci(j,j,k),k=4,6)
        write(lfnout,699) (Debye * Angstrom * mnoci(j,j,k),k=7,9)
      enddo
      write(lfnout,696)      
      do j=1,nbase
        trace=(mbase(j,j,4)+mbase(j,j,7)+mbase(j,j,9))/3        
        do k=1,6
          if(k.eq.1.or.k.eq.4.or.k.eq.6)then
            qmtraceless(j,j,k)=mnoci(j,j,k+3)-trace
          else
            qmtraceless(j,j,k)=mnoci(j,j,k)
          endif
        enddo
        write(lfnout,697) j,(Debye * Angstrom * (qmtraceless(j,j,k)),    &
     &       k=1,3)
        write(lfnout,699) (Debye * Angstrom * (qmtraceless(j,j,k)),      &
     &       k=4,6)
      enddo
 704  format(//,' Transition dipole moments of MEBFs (Debye)',/)
 705  format(//,' Transition dipole moments of NOCI states (Debye)',/)     
 703  format(//,' Component',a2,//)
      
      write(lfnout,704)
      write(lfnout,703)'X'
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        if (.not.mebfLabels) then
          write(lfnout,671) (i,i=ii,il)
        else
          write(lfnout,6711)  '|',(mebfLabel(i),i=ii,il)
          write(lfnout,6712)'-------|',                                 &
     &         ('--------------------',i=ii,il)
        end if        
        do j=1,nbase
          if(j.eq.1) then
            if (.not.mebfLabels) then
              write(lfnout,672) 1,(mbase(j,i,1),i=1,nbase)
            else
              write(lfnout,6721) mebfLabel(1),' |',                     &
     &             (mbase(j,i,1),i=1,nbase)
            end if            
          else
            ik=min(j-1,il)
            if (.not.mebfLabels) then
              write(lfnout,672) j,(mbase(j,i,1),i=1,nbase)
            else
              write(lfnout,6721) mebfLabel(j),' |',                     &
     &             (mbase(j,i,1),i=1,nbase)
            end if
          endif
        enddo
        write(lfnout,*) ' '
      enddo
      write(lfnout,703)'Y'
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        if (.not.mebfLabels) then
          write(lfnout,671) (i,i=ii,il)
        else
          write(lfnout,6711)  '|',(mebfLabel(i),i=ii,il)
          write(lfnout,6712)'-------|',                                 &
     &         ('--------------------',i=ii,il)
        end if        
        do j=1,nbase
          if(j.eq.1) then
            if (.not.mebfLabels) then
              write(lfnout,672) 1,(mbase(j,i,2),i=1,nbase)
            else
              write(lfnout,6721) mebfLabel(1),' |',                     &
     &             (mbase(j,i,2),i=1,nbase)
            end if            
          else
            ik=min(j-1,il)
            if (.not.mebfLabels) then
              write(lfnout,672) j,(mbase(j,i,2),i=1,nbase)
            else
              write(lfnout,6721) mebfLabel(j),' |',                     &
     &             (mbase(j,i,2),i=1,nbase)
            end if
          endif
        enddo
        write(lfnout,*) ' '
      enddo            
      write(lfnout,703)'Z'
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        if (.not.mebfLabels) then
          write(lfnout,671) (i,i=ii,il)
        else
          write(lfnout,6711)  '|',(mebfLabel(i),i=ii,il)
          write(lfnout,6712)'-------|',                                 &
     &         ('--------------------',i=ii,il)
        end if        
        do j=1,nbase
          if(j.eq.1) then
            if (.not.mebfLabels) then
              write(lfnout,672) 1,(mbase(j,i,3),i=1,nbase)
            else
              write(lfnout,6721) mebfLabel(1),' |',                     &
     &             (mbase(j,i,3),i=1,nbase)
            end if            
          else
            ik=min(j-1,il)
            if (.not.mebfLabels) then
              write(lfnout,672) j,(mbase(j,i,3),i=1,nbase)
            else
              write(lfnout,6721) mebfLabel(j),' |',                     &
     &             (mbase(j,i,3),i=1,nbase)
            end if
          endif
        enddo
        write(lfnout,*) ' '
      enddo      
      write(lfnout,705)
      write(lfnout,703)'X'
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,671) (i,i=ii,il)
        do j=1,nbase
          if(j.eq.1) then            
            write(lfnout,672) 1,(mnoci(j,i,1),i=1,nbase)
          else
            ik=min(j-1,il)
            write(lfnout,672) j,(mnoci(j,i,1),i=1,nbase)
          endif
        enddo
        write(lfnout,*) ' '
      enddo
      write(lfnout,703)'Y'
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,671) (i,i=ii,il)
        do j=1,nbase
          if(j.eq.1) then            
            write(lfnout,672) 1,(mnoci(j,i,2),i=1,nbase)
          else
            ik=min(j-1,il)
            write(lfnout,672) j,(mnoci(j,i,2),i=1,nbase)
          endif
        enddo
        write(lfnout,*) ' '
      enddo
      write(lfnout,703)'Z'
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,671) (i,i=ii,il)
        do j=1,nbase
          if(j.eq.1) then            
            write(lfnout,672) 1,(mnoci(j,i,3),i=1,nbase)
          else
            ik=min(j-1,il)
            write(lfnout,672) j,(mnoci(j,i,3),i=1,nbase)
          endif
        enddo
        write(lfnout,*) ' '
      enddo


 701  format(//,' Oscillator strength of MEBFs',/)
      write(lfnout,701)
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        if (.not.mebfLabels) then
          write(lfnout,671) (i,i=ii,il)
        else
          write(lfnout,6711)  '|',(mebfLabel(i),i=ii,il)
          write(lfnout,6712)'-------|',                                 &
     &         ('--------------------',i=ii,il)
        end if        
        do j=1,nbase
          if(j.eq.1) then
            if (.not.mebfLabels) then
              write(lfnout,672) 1,(obase(j,i),i=1,nbase)
            else
              write(lfnout,6721) mebfLabel(1),' |',                     &
     &             (obase(j,i),i=1,nbase)
            end if            
          else
            ik=min(j-1,il)
            if (.not.mebfLabels) then
              write(lfnout,672) j,(obase(j,i),i=1,nbase)
            else
              write(lfnout,6721) mebfLabel(j),' |',                     &
     &             (obase(j,i),i=1,nbase)
            end if
          endif
        enddo
        write(lfnout,*) ' '
      enddo      
      
 702  format(//,' Oscillator strength of NOCI states',/)
      write(lfnout,702)
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,671) (i,i=ii,il)
        do j=1,nbase
          if(j.eq.1) then            
            write(lfnout,672) 1,(onoci(j,i),i=1,nbase)
          else
            ik=min(j-1,il)
            write(lfnout,672) j,(onoci(j,i),i=1,nbase)
          endif
        enddo
        write(lfnout,*) ' '
      enddo

!** Dumping the results in the cml file
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
      call writetag_matrix_real(lfncml,info_cml,7,sep,hbase,                    &
     &   size(hbase,1),size(hbase,2),fmt_1)
      call close_tag(lfncml,label,6)
      if (ncorr.ne.1) then
        info_cml='dictRef="gr:overlap"'
        indent = 6
        call open_tag(lfncml,label,info_cml,indent)
        info_cml='empty'
        call writetag_matrix_real(lfncml,info_cml,7,sep,sbase,                  &
     &     size(hbase,1),size(hbase,2),fmt_1)
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
          call writetag_scalar_real(lfncml,info_cml2,8,(hev(i)-hev(1))          &
     &                                                   *27.2114,fmt_1)
          call close_tag(lfncml,label,7)
          info_cml = 'dictRef="gr:mebfCoeff" '//trim(id)
          call open_tag(lfncml,label,info_cml,7)
          fmt_1 = 'f20.10'
          call writetag_array_real(lfncml,info_cml,8,sep,hsave(i,:),            &
     &              size(hsave,2),fmt_1)
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
        call writetag_matrix_real(lfncml,info_cml,7,sep,tc,                     &
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

      deallocate(noci)
      deallocate(mnoci)
      deallocate(obase)
      deallocate(onoci)
      deallocate(qmtraceless)

      return

      end subroutine gronor_print_results
