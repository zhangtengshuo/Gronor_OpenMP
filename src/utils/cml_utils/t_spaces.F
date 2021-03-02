      character(len=4) function t_spaces(indent)
      implicit none
      integer :: indent,spaces,i
      character(len=4) s

      spaces = 0
      do i = 1, indent
        spaces = spaces + 4
      end do
      write(s,'(i4)') spaces
      t_spaces = 't'//trim(adjustl(s))
      return
      end 
      
