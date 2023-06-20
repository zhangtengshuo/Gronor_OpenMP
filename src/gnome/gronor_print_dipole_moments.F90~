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

subroutine gronor_print_dipole_moments(dqbase,mnuc,nbase,hev,hbase)

  use gnome_parameters, only : itest,idbg
  use gnome_data, only : com
  use cidef, only : lfnout,lfntst,lfnarx,mebfLabel,mebfLabels,header,key,nociwf,lfnxrx,lfndbg

  implicit none

  external :: gronor_print_matrix

  integer,intent(in)            :: nbase
  integer                       :: i,j,k,l,m
  real(kind=8)                  :: hev(nbase)

  integer                       :: lfnt

  real(kind=8),intent(in)       :: mnuc(9),hbase(nbase,nbase)
  real(kind=8),intent(inout)    :: dqbase(nbase,nbase,9)
  real(kind=8),allocatable      :: dqnoci(:,:,:)
  real(kind=8),allocatable      :: obase(:,:)
  real(kind=8),allocatable      :: onoci(:,:)
  real(kind=8),allocatable      :: qtraceless(:,:,:)
  real(kind=8)                  :: debye,angstrom
  real(kind=8)                  :: two3rds

  debye    = 2.54174644986
  angstrom = 0.529177249
  two3rds  = 2.0d0/3.0d0

  lfnt=0
  if(itest.gt.0) lfnt=lfntst

  allocate(dqnoci(nbase,nbase,9))
  allocate(obase(nbase,nbase))
  allocate(onoci(nbase,nbase))
  allocate(qtraceless(nbase,nbase,6))

  do i=1,nbase
    do j=1,nbase
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


  do i=1,9
    do j=1,nbase
      do k=1,nbase            
        do l=1,nbase
          do m=1,nbase
            dqnoci(k,j,i)=dqnoci(k,j,i)+nociwf(l,j)*nociwf(m,k)*dqbase(m,l,i)                
          enddo
        enddo
      enddo
    enddo
  enddo

  do j=1,nbase
    do k=1,j          
      obase(k,j)=two3rds*(hbase(j,j)-hbase(k,k))*(dqbase(k,j,1)**2+ &
          dqbase(k,j,2)**2+dqbase(k,j,3)**2)    
    enddo
  enddo

  do j=1,nbase
    do k=1,j          
      onoci(k,j)=two3rds*(hev(j)-hev(k))*(dqnoci(k,j,1)**2+dqnoci(k,j,2)**2+dqnoci(k,j,3)**2) 
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

  if ( idbg .ge. 5 ) then
    write(lfndbg,*)'Nuclear contribution'
    write(lfndbg,'(9f18.10)')(mnuc(i),i=1,9)
  endif
  do i=1,9
    do j=1,nbase
      dqbase(j,j,i)=dqbase(j,j,i)+mnuc(i)
      dqnoci(j,j,i)=dqnoci(j,j,i)+mnuc(i)
    enddo
  enddo

688 format(//,' Multipole moments')
689 format(/,' Origin of the dipole operator (x,y,z)',t50,3f10.4,/, &
        ' Origin of the quadrupole operator (x,y,z)',t50,3f10.4)          
690 format(/,' Dipole moment of MEBFs (Debye)',/)
693 format(t27,'X',t41,'Y',t55,'Z',t67,'Total',/)
699 format(t27,'XX',t41,'XY',t55,'XZ',t69,'YY',t83,'YZ',t97,'ZZ',/)
691 format(1x,i15,t18,4f14.4)
692 format(1x,a,t18,4f14.4)
  write(lfnout,688)
  write(lfnout,689) 0.0d0,0.0d0,0.0d0,(com(j),j=1,3)
  write(lfnout,690)
  write(lfnout,693)
  do j=1,nbase
    if(.not.mebfLabels) then
      write(lfnout,691) j,(debye*dqbase(j,j,k),k=1,3),debye* &
          sqrt(dqbase(j,j,1)**2+dqbase(j,j,2)**2+dqbase(j,j,3)**2)
    else
      write(lfnout,692) trim(mebfLabel(j)),(debye*dqbase(j,j,k),k=1,3),debye* &
          sqrt(dqbase(j,j,1)**2+dqbase(j,j,2)**2+dqbase(j,j,3)**2)
    endif
  enddo
694 format(/,' Dipole moment of NOCI states (Debye)',/)
  write(lfnout,694)
  write(lfnout,693)
  do j=1,nbase
    write(lfnout,691) j,(debye*dqnoci(j,j,k),k=1,3), &
        debye*sqrt(dqnoci(j,j,1)**2+dqnoci(j,j,2)**2+dqnoci(j,j,3)**2)        
  enddo

695 format(/,' Quadrupole moment of MEBFs (Debye*Angstrom)',/)
696 format(/,' In traceless form',/)
  write(lfnout,695)
  write(lfnout,699)
697 format(1x,i15,t18,6f14.4)
698 format(1x,a,t18,6f14.4)
  do j=1,nbase        
    if(.not.mebfLabels) then
      write(lfnout,697) j,(debye*angstrom*dqbase(j,j,k),k=4,9)
    else
      write(lfnout,698) trim(mebfLabel(j)),(debye*angstrom*dqbase(j,j,k),k=4,9)
    endif
  enddo
  write(lfnout,696)
  write(lfnout,699)

  !     Quadrupole moment traceless transformation
  !     
  !     1    0    0    -0.5  0   -0.5
  !     0    1.5  0     0    0     0
  !     0    0    1.5   0    0     0
  !    -0.5  0    0     1    0   -0.5
  !     0    0    0     0    1.5   0
  !    -0.5  0    0    -0.5  0     1

  do j=1,nbase        
    qtraceless(j,j,1) = dqbase(j,j,4)-0.5d0*dqbase(j,j,7)-0.5d0*dqbase(j,j,9)
    qtraceless(j,j,2) = 1.5d0*dqbase(j,j,5)
    qtraceless(j,j,3) = 1.5d0*dqbase(j,j,6)
    qtraceless(j,j,4) = -0.5d0*dqbase(j,j,4)+dqbase(j,j,7)-0.5d0*dqbase(j,j,9)
    qtraceless(j,j,5) = 1.5d0*dqbase(j,j,8)
    qtraceless(j,j,6) = -0.5d0*dqbase(j,j,4)-0.5d0*dqbase(j,j,7)+dqbase(j,j,9)
    if(.not.mebfLabels) then
      write(lfnout,697) j,(debye*angstrom*(qtraceless(j,j,k)),k=1,6)
    else
      write(lfnout,698) trim(mebfLabel(j)),(debye*angstrom*(qtraceless(j,j,k)),k=1,6)
    endif
  enddo

700 format(/,' Quadrupole moment of NOCI states (Debye*Angstrom)',/)
  write(lfnout,700)
  write(lfnout,699)
  do j=1,nbase
    write(lfnout,697) j,(debye*angstrom*dqnoci(j,j,k),k=4,9)
  enddo
  write(lfnout,696)
  write(lfnout,699)
  do j=1,nbase
    qtraceless(j,j,1) = dqnoci(j,j,4)-0.5d0*dqnoci(j,j,7)-0.5d0*dqnoci(j,j,9)
    qtraceless(j,j,2) = 1.5d0*dqnoci(j,j,5)
    qtraceless(j,j,3) = 1.5d0*dqnoci(j,j,6)
    qtraceless(j,j,4) = -0.5d0*dqnoci(j,j,4)+dqnoci(j,j,7)-0.5d0*dqnoci(j,j,9)
    qtraceless(j,j,5) = 1.5d0*dqnoci(j,j,8)
    qtraceless(j,j,6) = -0.5d0*dqnoci(j,j,4)-0.5d0*dqnoci(j,j,7)+dqnoci(j,j,9)        
    write(lfnout,697) j,(debye*angstrom*(qtraceless(j,j,k)),k=1,6)
  enddo
703 format(' Component',a2,/)


  write(lfnout,704)
704 format(/,' Transition dipole moments of MEBFs (a.u.)')
  key='dipMx'
  header='Component X'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      mebfLabels,mebfLabel,.false.,1.0d0,dqbase(1,1,1),nbase,7,6)
  key='dipMy'
  header='Component Y'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      mebfLabels,mebfLabel,.false.,1.0d0,dqbase(1,1,2),nbase,7,6)
  key='dipMz'
  header='Component Z'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      mebfLabels,mebfLabel,.false.,1.0d0,dqbase(1,1,3),nbase,7,6)
  key='oscM '
  header='Oscillator strength'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      mebfLabels,mebfLabel,.false.,1.0d0,obase,nbase,7,6)

  write(lfnout,705)
705 format(/,' Transition dipole moments of NOCI states (a.u.)',/)
  key='dipNx'
  header='Component X'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      .false.,mebfLabel,.false.,1.0d0,dqnoci(1,1,1),nbase,7,6)
  key='dipNy'
  header='Component Y'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      .false.,mebfLabel,.false.,1.0d0,dqnoci(1,1,2),nbase,7,6)
  key='dipNz'
  header='Component Z'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      .false.,mebfLabel,.false.,1.0d0,dqnoci(1,1,3),nbase,7,6)
  key='oscN '
  header='Oscillator strength'
  call gronor_print_matrix(lfnout,lfnarx,lfnxrx,lfnt,header,key, &
      .false.,mebfLabel,.false.,1.0d0,onoci,nbase,7,6)

  deallocate(dqnoci)
  deallocate(obase)
  deallocate(onoci)
  deallocate(qtraceless)

  return

end subroutine gronor_print_dipole_moments
