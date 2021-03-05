      subroutine writetag_scalar_real(lfn,info,indent,num,
     &                      numfmt)
      implicit none
      integer  :: lfn,indent
      real(kind=8) :: num
      character(len=80)  :: num2
      character(len=4)  :: ts
      character(len=4), external :: t_spaces
      character(len=132)   :: info,fmt_1,fmt_2,fmt_3,numfmt
      logical  :: scientific_format

      write(fmt_3,'(3a)')'(',trim(adjustl(numfmt)),')'
      write(num2,fmt=fmt_3)num
      num2=adjustl(num2)
      ts = t_spaces(indent)
      if (indent.eq.0) then
        fmt_1 = "(4a)"
        fmt_2 = "(5a)"
      else
        write(fmt_1,'(3a)')'(',trim(ts),',4a)'
        write(fmt_2,'(3a)')'(',trim(ts),',5a)'
      end if
      if (info(1:5).eq.'empty') then
        write(lfn,fmt=fmt_1)'<scalar dataType="xsd:double" ',
     &   '>',trim(num2),'</scalar>'
      else
        write(lfn,fmt=fmt_2)'<scalar dataType="xsd:double" ',trim(info),
     &   '>',trim(num2),'</scalar>'
      end if

      end subroutine writetag_scalar_real
