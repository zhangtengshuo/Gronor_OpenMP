      subroutine writetag_scalar_integer(lfn,info,indent,num)
      implicit none
      integer  :: lfn,num,indent
      character(len=4)  :: ts
      character(len=80) :: num2
      character(len=4), external :: t_spaces
      character(len=132)   :: info,fmt_1,fmt_2

      ts = t_spaces(indent)
      write(num2,'(i10)') num
      num2=adjustl(num2)
      if (indent.eq.0) then
        fmt_1 = "(4a)"
        fmt_2 = "(5a)"
      else
        write(fmt_1,'(3a)')'(',trim(ts),',4a)'
        write(fmt_2,'(3a)')'(',trim(ts),',5a)'
      end if
      if (info(1:5).eq.'empty') then
        write(lfn,fmt=fmt_1)'<scalar dataType="xsd:integer"','>',
     &       trim(num2),'</scalar>'
      else
        write(lfn,fmt=fmt_2)'<scalar dataType="xsd:string" ',
     &       trim(info),'>',trim(num2),'</scalar>'
      end if

      end subroutine writetag_scalar_integer
