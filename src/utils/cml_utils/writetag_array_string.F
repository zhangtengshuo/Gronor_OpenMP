!>  @brief
!>    Writes an opening and closing tag for an array of strings
!>    of length 'length', values are separated by 'sep'
!>    <array dataType="xsd:string" info> array of strings </array>
!>
!>  @authors
!>    Coen de Graaf (URV, Tarragona)
!>    
!>  @date
!>    March 2021
!>
      subroutine writetag_array_string(lfn,info,indent,sep,num,length)
      implicit none
      integer  :: lfn,indent,length,j
      character(len=*) :: num(length)
      character(len=1)     :: sep
      character(len=4)     :: ts,l,l1
      character(len=4), external :: t_spaces
      character(len=132)   :: info,fmt_1

      do j = 1, length
        num(j) = adjustl(num(j))
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
      write(lfn,fmt=fmt_1)'<array dataType="xsd:string" ',trim(info),
     &     ' size="',trim(adjustl(l)),'">',
     &     (trim(num(j)),sep,j=1,length-1),trim(num(length)),'</array>'
      
      end subroutine writetag_array_string
