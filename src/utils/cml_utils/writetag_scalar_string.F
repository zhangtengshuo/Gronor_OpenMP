!>  @brief
!>    Writes an opening and closing tag for a string
!>    <scalar dataType="xsd:string" info > string </scalar>
!>
!>  @authors
!>    Coen de Graaf (URV, Tarragona)
!>    
!>  @date
!>    March 2021
!> 
      subroutine writetag_scalar_string(lfn,info,indent,charstring)
      implicit none
      integer  :: lfn,indent
      character(len=4)  :: ts
      character(len=4), external :: t_spaces
      character(len=132)   :: info,fmt_1,fmt_2
      character(len=*)     :: charstring

      ts = t_spaces(indent)
      charstring=adjustl(charstring)
      if (indent.eq.0) then
        fmt_1 = '(4a)'
        fmt_2 = '(5a)'
      else
        write(fmt_1,'(3a)') '(',trim(ts),',4a)'
        write(fmt_2,'(3a)') '(',trim(ts),',5a)'
      end if
      if (info(1:5).eq.'empty') then
        write(lfn,fmt=fmt_1)'<scalar dataType="xsd:string"',
     &  '>',trim(charstring),'</scalar>'
      else
        write(lfn,fmt=fmt_2)'<scalar dataType="xsd:string" ',trim(info),
     &  '>',trim(charstring),'</scalar>'
      end if

      end subroutine writetag_scalar_string
