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

      subroutine gronor_print_dipole_moments(dqbase,mnuc,nbase,hev)

      use gnome_parameters, only : itest,ncols
      use gnome_data      , only : com
      use cidef           , only : lfnout,lfntst,lfnarx,lfncml,ncorr,           &
     &                             nwt,mebfLabel,mebfLabels,header,key,         & 
     &                             nociwf

      implicit none

      integer,intent(in)            :: nbase
      integer                       :: i,j,k,ii,ik,il,l,m
      integer                       :: nk,info,lwork
      real(kind=8)                  :: hbt(nbase),tc(nbase,nbase)
!      real(kind=8)                  :: hsave(nbase,nbase)
!      real(kind=8)                  :: ssave(nbase,nbase)
      real(kind=8)                  :: hev(nbase)
      real(kind=8),allocatable      :: work(:)
      character(len=5)              :: token

      character(len=132)            :: info_cml,info_cml2,fmt_1,id
      character(len=20)             :: label
      character(len=1)              :: sep
      integer                       :: indent
      
      real(kind=8),intent(in)       :: mnuc(9)
      real(kind=8),intent(inout)    :: dqbase(nbase,nbase,9)
!      real(kind=8),allocatable      :: noci(:,:)
      real(kind=8),allocatable      :: dqnoci(:,:,:)
      real(kind=8),allocatable      :: obase(:,:)
      real(kind=8),allocatable      :: onoci(:,:)
      real(kind=8),allocatable      :: qtraceless(:,:,:)
      real(kind=8)                  :: debye,angstrom,trace
      real(kind=8)                  :: two3rds

      debye    = 2.54174644986
      angstrom = 0.529177249
      two3rds  = 2.0d0/3.0d0

!      allocate(noci(nbase,nbase))
      allocate(dqnoci(nbase,nbase,9))
      allocate(obase(nbase,nbase))
      allocate(onoci(nbase,nbase))
      allocate(qtraceless(nbase,nbase,6))
      
      do i=1,nbase
        do j=1,nbase
!          noci(i,j)=0.0d0
          obase(i,j)=0.0d0
          onoci(i,j)=0.0d0
          do k=1,9
            dqnoci(i,j,k)=0.0d0
          enddo
          do k=1,6
            qtraceless(i,j,k)=0.0d0
          enddo
        enddo
      enddo

!      if(nbase.eq.1) then
!        noci(1,1)=1.0d0
!      else
!        do i=1,nbase
!          do j=1,nbase
!            noci(j,i)=hsave(j,i)            
!          enddo
!        enddo
!      endif
      
      do i=1,9
        do j=1,nbase
          do k=1,nbase            
            do l=1,nbase
              do m=1,nbase
                dqnoci(k,j,i)=dqnoci(k,j,i)+                              &
     &               nociwf(l,j)*nociwf(m,k)*dqbase(m,l,i)                
              enddo              
            enddo
          enddo
        enddo
      enddo

      do j=1,nbase
        do k=1,j          
          obase(k,j)=obase(k,j)+two3rds*(hev(j)-hev(k))*                  &
     &         (dqbase(k,j,1)**2+dqbase(k,j,2)**2+                        &
     &         dqbase(k,j,3)**2)
        enddo        
      enddo
      
      do j=1,nbase
        do k=1,j          
          onoci(k,j)=onoci(k,j)+two3rds*(hev(j)-hev(k))*                  &
     &         (dqnoci(k,j,1)**2+dqnoci(k,j,2)**2+                        &
     &         dqnoci(k,j,3)**2)
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
            dqbase(k,j,i)=dqbase(k,j,i)+mnuc(i)
            dqnoci(k,j,i)=dqnoci(k,j,i)+mnuc(i)
          enddo
        enddo        
      enddo
 688  format(//,' Multipole moments')
 689  format(/,' Origin of the dipole operator (x,y,z)',t50,3f10.4,/            &
     &     ,' Origin of the quadrupole operator (x,y,z)',t50,3f10.4)          
 690  format(/,' Dipole moment of MEBFs (Debye)',/)
 691  format(7x,i5,2x,' X:',2x,f10.4,4x,' Y:',2x,f10.4,4x,                      &
     &     ' Z:',2x,f10.4,4x,' Total:',2x,f10.4)
 692  format(7x,a5,2x,' X:',2x,f10.4,4x,' Y:',2x,f10.4,4x,                      &
     &     ' Z:',2x,f10.4,4x,' Total:',2x,f10.4)
      write(lfnout,688)
      write(lfnout,689) 0.0d0,0.0d0,0.0d0,(com(j),j=1,3)
      write(lfnout,690)
      do j=1,nbase
        if(.not.mebfLabels) then
          write(lfnout,691) j,(debye*dqbase(j,j,k),k=1,3),debye*                &
     &         sqrt(dqbase(j,j,1)**2+dqbase(j,j,2)**2+dqbase(j,j,3)**2)
        else
          write(lfnout,692) mebfLabel(j),                                       &
     &         (debye*dqbase(j,j,k),k=1,3),debye*                               &
     &         sqrt(dqbase(j,j,1)**2+dqbase(j,j,2)**2+dqbase(j,j,3)**2)
        endif        
      enddo
 694  format(/,' Dipole moment of NOCI states (Debye)',/)
      write(lfnout,694)
      do j=1,nbase
        write(lfnout,691) j,(debye*dqnoci(j,j,k),k=1,3),                 &
     &       debye*sqrt(dqnoci(j,j,1)**2+dqnoci(j,j,2)**2+               &
     &       dqnoci(j,j,3)**2)        
      enddo
      
 695  format(/,' Quadrupole moment of MEBFs (Debye*Angstrom)',/)
 696  format(/,' In traceless form',/)
      write(lfnout,695)
 697  format(7x,i5,2x,                                                  &
     &     ' XX:',2x,f10.4,4x,' XY:',2x,f10.4,4x,' XZ:',2x,f10.4)      
 698  format(7x,a5,2x,                                                  &
     &     ' XX:',2x,f10.4,4x,' XY:',2x,f10.4,4x,' XZ:',2x,f10.4)       
 699  format(14x,                                                       &
     &     ' YY:',2x,f10.4,4x,' YZ:',2x,f10.4,4x,' ZZ:',2x,f10.4)       
      do j=1,nbase        
        if(.not.mebfLabels) then
          write(lfnout,697) j,(debye*angstrom*dqbase(j,j,k),k=4,6)
        else
          write(lfnout,698) mebfLabel(j),(debye*angstrom*dqbase(j,j,k),  &
     &         k=4,6)
        endif        
        write(lfnout,699) (debye*angstrom*dqbase(j,j,k),k=7,9)
      enddo      
      write(lfnout,696)
      do j=1,nbase        
        trace=(dqbase(j,j,4)+dqbase(j,j,7)+dqbase(j,j,9))/3
        do k=1,6
          if(k.eq.1.or.k.eq.4.or.k.eq.6) then
            qtraceless(j,j,k)=dqbase(j,j,k+3)-trace
          else
            qtraceless(j,j,k)=dqbase(j,j,k)
          endif
        enddo
        
        if(.not.mebfLabels) then
          write(lfnout,697) j,(debye*angstrom*(qtraceless(j,j,k)),      &
     &         k=1,3)
        else
          write(lfnout,698) mebfLabel(j),(debye*angstrom*               &
     &         (qtraceless(j,j,k)),k=1,3)
        endif
        write(lfnout,699) (debye*angstrom*(qtraceless(j,j,k)),          &
     &       k=4,6)
      enddo
     
 700  format(/,' Quadrupole moment of NOCI states (Debye*Angstrom)',/)
      write(lfnout,700)
      do j=1,nbase
        write(lfnout,697) j,(debye*angstrom*dqnoci(j,j,k),k=4,6)
        write(lfnout,699) (debye*angstrom*dqnoci(j,j,k),k=7,9)
      enddo
      write(lfnout,696)      
      do j=1,nbase
        trace=(dqbase(j,j,4)+dqbase(j,j,7)+dqbase(j,j,9))/3        
        do k=1,6
          if(k.eq.1.or.k.eq.4.or.k.eq.6) then
            qtraceless(j,j,k)=dqnoci(j,j,k+3)-trace
          else
            qtraceless(j,j,k)=dqnoci(j,j,k)
          endif
        enddo
        write(lfnout,697) j,(debye*angstrom*(qtraceless(j,j,k)),        &
     &       k=1,3)
        write(lfnout,699) (debye*angstrom*(qtraceless(j,j,k)),          &
     &       k=4,6)
      enddo     
 703  format(' Component',a2,/)

      write(lfnout,704)
 704  format(/,' Transition dipole moments of MEBFs (Debye)')
      key='     '
      header='Component X'
      call gronor_print_matrix(lfnout,0,0,header,key,                   &
     &     mebfLabels,mebfLabel,.false.,1.0d0,dqbase(1,1,1),nbase,7)
      header='Component Y'
      call gronor_print_matrix(lfnout,0,0,header,key,                   &
     &     mebfLabels,mebfLabel,.false.,1.0d0,dqbase(1,1,2),nbase,7)
      header='Component Z'
      call gronor_print_matrix(lfnout,0,0,header,key,                   &
     &     mebfLabels,mebfLabel,.false.,1.0d0,dqbase(1,1,3),nbase,7)
      header='Oscillator strength'
      call gronor_print_matrix(lfnout,0,0,header,key,                   &
     &     mebfLabels,mebfLabel,.false.,1.0d0,obase,nbase,7)
      
      
      write(lfnout,705)
 705  format(/,' Transition dipole moments of NOCI states (Debye)',/)
      key='     '
      header='Component X'
      call gronor_print_matrix(lfnout,0,0,header,key,                   &
     &     .false.,mebfLabel,.false.,1.0d0,dqnoci(1,1,1),nbase,7)
      header='Component Y'
      call gronor_print_matrix(lfnout,0,0,header,key,                   &
     &     .false.,mebfLabel,.false.,1.0d0,dqnoci(1,1,2),nbase,7)
      header='Component Z'
      call gronor_print_matrix(lfnout,0,0,header,key,                   &
     &     .false.,mebfLabel,.false.,1.0d0,dqnoci(1,1,3),nbase,7)
      header='Oscillator strength'
      call gronor_print_matrix(lfnout,0,0,header,key,                   &
     &     .false.,mebfLabel,.false.,1.0d0,onoci,nbase,7)

!      deallocate(nociwf)
      deallocate(dqnoci)
      deallocate(obase)
      deallocate(onoci)
      deallocate(qtraceless)

      return

      end subroutine gronor_print_dipole_moments
