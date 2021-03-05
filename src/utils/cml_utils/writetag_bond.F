      subroutine writetag_bond(lfn,indent,id1,id2,bond_order)
      implicit none
      integer  :: lfn,indent
      character(len=132) fmt_1
      character(len=*) :: id1,id2
      character(len=4)  :: ts
      character(len=1)  :: bond_order
      character(len=4), external :: t_spaces

      id1 = adjustl(id1)
      id2 = adjustl(id2)
      ts = t_spaces(indent)
      if (indent.eq.0) then
        fmt_1 = '(7a)'
      else
        write(fmt_1,'(3a)')'(',trim(ts),',7a)'
      end if
      write(lfn,fmt=fmt_1)'<bond atomRefs2="',trim(id1),' ',trim(id2),
     &  '" order="',bond_order,'"/>'
      
      end subroutine writetag_bond
