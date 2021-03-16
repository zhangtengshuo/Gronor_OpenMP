!>  @brief
!>    Writes a tag to specify an atom
!>    <atom id=" " elementType=" " x3=" " y3=" " z3=" "/>
!>
!>  @authors
!>    Coen de Graaf (URV, Tarragona)
!>    
!>  @date
!>    March 2021
!>
      subroutine writetag_atom(lfn,indent,id,element,x,y,z)
      implicit none
      integer  :: lfn,indent
      real(kind=8) :: x,y,z
      character(len=132) fmt_1
      character(len=80) :: x3,y3,z3
      character(len=*) :: id
      character(len=4)  :: ts
      character(len=2)  :: element
      character(len=4), external :: t_spaces

      write(x3,'(f15.8)') x
      write(y3,'(f15.8)') y
      write(z3,'(f15.8)') z
      x3 = adjustl(x3)
      y3 = adjustl(y3)
      z3 = adjustl(z3)
      id = adjustl(id)
      element = adjustl(element)
      ts = t_spaces(indent)
      if (indent.eq.0) then
        fmt_1 = '(11a)'
      else
        write(fmt_1,'(3a)')'(',trim(ts),',11a)'
      end if
      write(lfn,fmt=fmt_1)'<atom id="',trim(id),'" elementType="',
     &  trim(element),'" x3="',trim(x3),'" y3="',trim(y3),'" z3="',
     &  trim(z3),'"/>'
      
      end subroutine writetag_atom
