

      subroutine gronor_gnome_molcas_input(lfnsym,lfnout,ipr)
      use gnome_data
      implicit none
      integer               :: lfnsym, lfnout, ipr
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
      write(lfnout,'(6x,A,10x,2A,5x,A)')
     & 'Atom','x             y             z',
     &         '          Znuc',lAO
      do iAtom = 1, nNucl
        read(lfnsym,*) i,centn(iAtom),xCord(iAtom),yCord(iAtom),
     &                   zCord(iAtom),zNuc(iAtom),(nCntr(iAtom,l),
     &                   l=1,7)
        write(lfnout,626)i,centn(iAtom),xCord(iAtom),yCord(iAtom),
     &                   zCord(iAtom),zNuc(iAtom),(nCntr(iAtom,l),
     &                   l=1,7)
      end do
      flush(lfnout)
 626  format (I4,2x,A4,2x,3F14.8,2x,F7.2,3x,7I3)
      deallocate( nCntr )
      return
      end
