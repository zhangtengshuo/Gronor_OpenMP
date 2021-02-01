
!> @brief
!! GNOME input of sym file
!!     
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  

      subroutine gronor_gnome_molcas_input
      use gnome_data
      use cidef
      implicit none
*      integer               :: lfnsym, lfnout, lfnarx, ipr
      integer               :: i,l,iAtom
      integer, allocatable  :: nCntr(:,:)
      character (len = 19)  :: lAO

      read(lfnsym,*) text
!     first line of $Project.sym is a dummy
      read(lfnsym,'(a72)') name(1)
      read(lfnsym,'(a72)') name(2)
      write(lfnout,*) ' '
      write(lfnout,*) ' '
      write(lfnout,*) ' OpenMolcas header from Seward (integral module)'
      write(lfnout,*) name(1)
      write(lfnout,*) name(2)
      write(lfnout,*) ' '
      flush(lfnout) 
      read(lfnsym,'(a17,i4)') text,nNucl

      allocate( centn(nNucl)   )
      allocate( xCord(nNucl)   )
      allocate( yCord(nNucl)   )
      allocate( zCord(nNucl)   )
      allocate(  zNuc(nNucl)   )
      allocate( nCntr(nNucl,7) )

      read(lfnsym,*) text
      lAO = 's  p  d  f  g  h  i'
      write(lfnout,625)
     & 'Atom','x             y             z',
     &     '          Znuc',lAO
      write(lfnarx,624) nNucl
 624  format('Atoms',i10)
      write(lfnarx,625)
     & 'Atom','x             y             z',
     &     '          Znuc',lAO
 625  format(6x,A,10x,2A,5x,A)
      write(lfnout,*) ' '
      do iAtom = 1, nNucl
        read(lfnsym,*) i,centn(iAtom),xCord(iAtom),yCord(iAtom),
     &                   zCord(iAtom),zNuc(iAtom),(nCntr(iAtom,l),
     &                   l=1,7)
        write(lfnout,626)i,centn(iAtom),xCord(iAtom),yCord(iAtom),
     &                   zCord(iAtom),zNuc(iAtom),(nCntr(iAtom,l),
     &                   l=1,7)
        write(lfnarx,626)i,centn(iAtom),xCord(iAtom),yCord(iAtom),
     &                   zCord(iAtom),zNuc(iAtom),(nCntr(iAtom,l),
     &                   l=1,7)
      end do
      flush(lfnout)
 626  format (I4,2x,A4,2x,3F14.8,2x,F7.2,3x,7I3)
      deallocate( nCntr )
      return
      end
