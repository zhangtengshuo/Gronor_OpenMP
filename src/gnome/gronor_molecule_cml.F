      subroutine gronor_molecule_cml()
      use cidef
      use gnome_data


      implicit none
      character(len=2), external :: element
      real (kind=8), external :: atom_radius
      integer :: indent,i,j,numi,numj
      character(len=132)  :: info_cml
      character(len=20)   :: label
      character(len=4)    :: id,id1,id2
      character(len=2)    :: element_name
      character(len=1)    :: bond_order
      real(kind=8)        :: distance
      real(kind=8),parameter  :: c=0.529177249d0

      indent = 3
      label = 'molecule'
      info_cml = 'id="molgeom"'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 4
      label = 'atomArray'
      info_cml = 'empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 5
      do i = 1, nnucl
        id = trim(centn(i))
        element_name = element(int(znuc(i)))
        call writetag_atom(lfncml,indent,id,element_name,c*xcord(i),
     &                        c*ycord(i),c*zcord(i))
      end do 
      indent = 4
      label = 'atomArray'
      call close_tag(lfncml,label,indent)
      label = 'bondArray'
      info_cml='empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent = 5
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
      indent=4
      call close_tag(lfncml,label,indent)
      indent = 3
      label = 'molecule'
      call close_tag(lfncml,label,indent)

      return
      end subroutine gronor_molecule_cml

