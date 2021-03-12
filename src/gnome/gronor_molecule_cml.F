      subroutine gronor_molecule_cml(indent)
      use cidef
      use gnome_data

      implicit none
      character(len=2), external :: element
      real (kind=8), external :: atom_radius
      integer :: indent,i,j,l,numi,numj,nDiffEl
      integer,allocatable :: stoichio(:)
      character(len=132)  :: info_cml,concise,info_cml2,info_cml3
      character(len=80)   :: label,fmt_1,fmt_2
      character(len=4)    :: id,id1,id2,dummy
      character(len=2)    :: element_name
      character(len=2),allocatable  :: molFormula(:)
      character(len=1)    :: bond_order
      real(kind=8)        :: distance
      real(kind=8),parameter  :: c=0.529177249d0
      logical :: new
      character (len=1), dimension(7)  :: lMoment
      character (len=4), dimension(7)  :: dumstr
      character (len=28)               :: basis_str

      data lMoment /'s','p','d','f','g','h','i'/


      label = 'molecule'
      info_cml = 'id="molgeom"'
      call open_tag(lfncml,label,info_cml,indent)
      indent = indent + 1
      label = 'atomArray'
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent = indent + 1
      do i = 1, nnucl
        id = trim(centn(i))
        element_name = element(int(znuc(i)))
        call writetag_atom(lfncml,indent,id,element_name,c*xcord(i),
     &                        c*ycord(i),c*zcord(i))
      end do 
      indent = indent - 1
      label = 'atomArray'
      call close_tag(lfncml,label,indent)
      allocate(molFormula(nnucl))
      allocate(stoichio(nnucl))
      molFormula = ''
      stoichio = 0
      nDiffEl = 1
      molFormula(nDiffEl) = element(int(znuc(1)))
      stoichio(nDiffEl) = 1
      do i = 2, nnucl
        new = .true.
        j = 1
        do while ((new).and.(j.le.nDiffEl))
          if (element(int(znuc(i))).eq.molFormula(j)) then
            stoichio(j) = stoichio(j) + 1
            new = .false.
          end if
          j = j + 1 
        end do
        if (new) then
          nDiffEl = nDiffEl + 1
          molFormula(nDiffEl) = element(int(znuc(i)))
          stoichio(nDiffEl) = 1
        end if
      end do
      write(dummy,'(i4)')nDiffEl
      dummy=adjustl(dummy)
      write(fmt_1,'(2a)')trim(dummy),'(a,i3)'
      write(fmt_2,'(3a)')'(',trim(adjustl(fmt_1)),')'
      write(concise,fmt=fmt_2)(molFormula(i),stoichio(i),i=1,nDiffEl)
      label = 'formula'
      write(info_cml,'(3a)')'concise="',trim(concise),'"'
      call open_tag(lfncml,label,info_cml,indent) 
      write(fmt_1,'(6a)')'(a,',trim(dummy),'(a2,1x)',
     &    ',a,',trim(dummy),'(f5.1),a)'
      write(info_cml,fmt=fmt_1)'elementType="',
     &  (molFormula(i),i=1,nDiffEl),
     &  '" count="',(float(stoichio(i)),i=1,nDiffEl),'"'
      label = 'atomArray'
      call write_simpletag(lfncml,label,info_cml,indent+1)
      label = 'formula'
      call close_tag(lfncml,label,indent)
      label = 'bondArray'
      info_cml='empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent = indent + 1
      do i = 1, nnucl-1
        numi = int(znuc(i))
        do j = i+1, nnucl
          numj = int(znuc(j))
          distance = c * sqrt((xcord(i)-xcord(j))**2 +
     &          (ycord(i)-ycord(j))**2 + (zcord(i)-zcord(j))**2) * 0.1
          if (distance .lt. atom_radius(numi)+atom_radius(numj)) then
            bond_order='S'
            id1 = trim(centn(i))
            id2 = trim(centn(j))
            call writetag_bond(lfncml,indent,id1,id2,bond_order)
          endif
        enddo
      enddo
      indent = indent - 1
      call close_tag(lfncml,label,indent)
      indent = indent - 1
      label = 'molecule'
      call close_tag(lfncml,label,indent)
      label = 'molecule'
      info_cml = 'id="basis_set"'
      call open_tag(lfncml,label,info_cml,indent)
      indent= indent + 1
      label = 'atomArray'
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent = indent + 1
      label = 'atom'
      do i = 1, nnucl
        info_cml = 'id="'//trim(centn(i))//'"'
        info_cml2 = info_cml//' array dataType="xsd:integer"'
        info_cml3 = info_cml2//' dictRef="gr:contracted"'
        call write_simpletag(lfncml,label,info_cml3,indent)
      end do
      indent= indent - 1
      label = 'atomArray'
      info_cml = 'empty'
      call close_tag(lfncml,label,indent)
      indent = indent - 1
      label = 'molecule'
      call close_tag(lfncml,label,indent)

      deallocate(molFormula)
      deallocate(stoichio)
      return
      end subroutine gronor_molecule_cml

