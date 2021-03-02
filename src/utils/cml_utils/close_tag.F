      subroutine close_tag(lfn,label,indent)
      implicit none
      integer  :: lfn,indent
      character(len=4)     :: ts
      character(len=4), external :: t_spaces
      character(len=20)    :: label
      character(len=132)   :: fmt_1

      ts = t_spaces(indent)
      if (indent.eq.0) then
        fmt_1 = "(3a)"
      else
        write(fmt_1,'(3a)')'(',trim(ts),',3a)'
      end if
      write(lfn,fmt=fmt_1) '</',trim(label),'>'

      end subroutine close_tag
