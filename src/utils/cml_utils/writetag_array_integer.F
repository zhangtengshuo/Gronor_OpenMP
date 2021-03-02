      subroutine writetag_array_integer(lfn,info,indent,sep,num,length)
      implicit none
      integer  :: lfn,indent,length,j
      integer  :: num(length)
      character(len=80) :: num2(length)
      character(len=1)     :: sep
      character(len=4)     :: ts,l,l1
      character(len=4), external :: t_spaces
      character(len=132)   :: info,fmt_1

      do j = 1, length
        write(num2(j),'(i10)')num(j)
        num2(j) = adjustl(num2(j))
      end do
      write(l,'(i4)') length
      write(l1,'(i4)') 2*length
      ts = t_spaces(indent)
      if (info(1:5).eq.'empty') info=''
      if (sep.ne." ") then
        info=trim(info)//' delimiter="'//sep//'" '
      end if
      if (indent.eq.0) then
        write(fmt_1,'(3a)')'(5a,',trim(adjustl(l1)),'a)'
      else
        write(fmt_1,'(5a)')'(',trim(ts),',5a,',trim(adjustl(l1)),'a)'
      end if
      write(lfn,fmt=fmt_1)'<array dataType="xsd:integer"',trim(info),
     &  ' size="',trim(adjustl(l)),'">',
     &  (trim(num2(j)),sep,j=1,length-1),trim(num2(length)),'</array>'
      
      end subroutine writetag_array_integer
