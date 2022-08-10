      subroutine gronor_molecule_cml(indent)
      use cidef
      use gnome_data

      implicit none
      character(len=2), external :: element
      real (kind=8), external :: atom_radius
      integer :: indent,i,j,l,numi,numj,nDiffEl,lmax
      integer,allocatable :: stoichio(:)
      character(len=132)  :: info_cml,concise,info_cml2,info_cml3
      character(len=80)   :: label,fmt_1,fmt_2
      character(len=4)    :: id,id1,id2,dummy
      character(len=2)    :: element_name
      character(len=2),allocatable  :: molFormula(:)
      character(len=1)    :: bond_order,sep
      real(kind=8)        :: distance
      real(kind=8),parameter  :: c=0.529177249d0
      logical :: new
      character (len=1), dimension(7)  :: lMoment

      data lMoment /'s','p','d','f','g','h','i'/

      label='molecule'
      info_cml='id="molgeom"'
      call open_tag(lfncml,label,info_cml,indent)
      indent=indent+1
      label='atomArray'
      info_cml='empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent=indent+1
      do i=1, nnucl
        id=trim(centn(i))
        element_name=element(int(znuc(i)))
        call writetag_atom(lfncml,indent,id,element_name,c*xcord(i),
     &                        c*ycord(i),c*zcord(i))
      enddo 
      indent=indent-1
      label='atomArray'
      call close_tag(lfncml,label,indent)
      allocate(molFormula(nnucl))
      allocate(stoichio(nnucl))
      molFormula=''
      stoichio=0
      nDiffEl=1
      molFormula(nDiffEl)=element(int(znuc(1)))
      stoichio(nDiffEl)=1
      do i=2, nnucl
        new=.true.
        j=1
        do while ((new).and.(j.le.nDiffEl))
          if (element(int(znuc(i))).eq.molFormula(j)) then
            stoichio(j)=stoichio(j)+1
            new=.false.
          endif
          j=j+1 
        enddo
        if (new) then
          nDiffEl=nDiffEl+1
          molFormula(nDiffEl)=element(int(znuc(i)))
          stoichio(nDiffEl)=1
        endif
      enddo
      write(dummy,'(i4)')nDiffEl
      dummy=adjustl(dummy)
      write(fmt_1,'(2a)')trim(dummy),'(a,i3)'
      write(fmt_2,'(3a)')'(',trim(adjustl(fmt_1)),')'
      write(concise,fmt=fmt_2)(molFormula(i),stoichio(i),i=1,nDiffEl)
      label='formula'
      write(info_cml,'(3a)')'concise="',trim(concise),'"'
      call open_tag(lfncml,label,info_cml,indent) 
      write(fmt_1,'(6a)')'(a,',trim(dummy),'(a2,1x)',
     &    ',a,',trim(dummy),'(f5.1),a)'
      write(info_cml,fmt=fmt_1)'elementType="',
     &  (molFormula(i),i=1,nDiffEl),
     &  '" count="',(float(stoichio(i)),i=1,nDiffEl),'"'
      label='atomArray'
      call write_simpletag(lfncml,label,info_cml,indent+1)
      label='formula'
      call close_tag(lfncml,label,indent)
      label='bondArray'
      info_cml='empty'
      call open_tag(lfncml,label,info_cml,indent)
      indent=indent+1
      do i=1, nnucl-1
        numi=int(znuc(i))
        do j=i+1, nnucl
          numj=int(znuc(j))
          distance=c*sqrt((xcord(i)-xcord(j))**2 +
     &          (ycord(i)-ycord(j))**2+(zcord(i)-zcord(j))**2)*0.1d0
          if (distance .lt. atom_radius(numi)+atom_radius(numj)) then
            bond_order='S'
            id1=trim(centn(i))
            id2=trim(centn(j))
            call writetag_bond(lfncml,indent,id1,id2,bond_order)
          endif
        enddo
      enddo
      indent=indent-1
      call close_tag(lfncml,label,indent)
      indent=indent-1
      label='molecule'
      call close_tag(lfncml,label,indent)
      label='module'
      info_cml='dictRef="gr:atomicBasisSet"'
      call open_tag(lfncml,label,info_cml,indent)
      indent=indent+1
      do i=1, nnucl
        label='list'
        info_cml='dictRef="gr:contracted"'
        call open_tag(lfncml,label,info_cml,indent)
        indent=indent+1
        id=trim(centn(i))
        info_cml='dictRef="gr:atomLabel"'
        call writetag_scalar_string(lfncml,info_cml,indent,id)
        do l=1, 7
          if (nCntr(i,l).ne.0) lmax=l
        enddo
        sep='|'
        info_cml='dictRef="gr:shells"'
        call writetag_array_string(lfncml,info_cml,indent,sep,
     &                  lMoment,lmax)
        sep=' '
        info_cml='dictRef="gr:contractedFunctions"'
        call writetag_array_integer(lfncml,info_cml,indent,sep,
     &                  nCntr(i,:),lmax)
        indent=indent-1
        call close_tag(lfncml,label,indent)
      enddo
      indent=indent -1
      label='module'
      call close_tag(lfncml,label,indent)
      call flush(lfncml)

      deallocate(molFormula)
      deallocate(stoichio)
      return
      end subroutine gronor_molecule_cml

