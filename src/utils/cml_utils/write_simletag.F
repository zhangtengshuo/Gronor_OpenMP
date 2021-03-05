!>  @brief
!>    Writes a 'standalone' tag <label info />
!>
!>  @authors
!>    Coen de Graaf (URV, Tarragona)
!>
!>  @date
!>    March 2021
!>
      subroutine write_simpletag(lfn,label,info,indent)
      implicit none
      integer  :: lfn,indent
      character(len=4)     :: ts
      character(len=4), external :: t_spaces
      character(len=20)    :: label
      character(len=132)   :: info,fmt_1,fmt_2

      ts = t_spaces(indent)
      if (indent.eq.0) then
        fmt_1 = "(3a)"
        fmt_2 = "(2a,1x,2a)"
      else
        write(fmt_1,'(3a)')'(',trim(ts),',3a)'
        write(fmt_2,'(3a)')'(',trim(ts),',2a,1x,2a)'
      end if
      if (info(1:5).eq.'empty') then
        write(lfn,fmt=fmt_1) '<',trim(label),'/>'
      else
        write(lfn,fmt=fmt_2) '<',trim(label),trim(info),'/>'
      end if

      end subroutine write_simpletag
