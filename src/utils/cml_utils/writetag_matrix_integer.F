!>  @brief
!>    Writes an opening and closing tag for a matrix of integers
!>    of dimension 'row' x 'col', values are separated by 'sep'
!>    <matrix dataType="xsd:integer" info columns=" " rows=" "> matrix of integers </matrix>
!>
!>  @authors
!>    Coen de Graaf (URV, Tarragona)
!>
!>  @date
!>    March 2021
!>
      subroutine writetag_matrix_integer(lfn,info,indent,sep,                   &
     &                                               num,row,col)
      implicit none
      integer  :: lfn,indent,i,j,k,col,row
      integer  :: num(row,col)
      character(len=80)    ::num2(row*col)
      character(len=1)     :: sep
      character(len=4)     :: rc,ts,r,c
      character(len=4), external :: t_spaces
      character(len=132)   :: info,fmt_1

      num2 = ''
      k = 0
      do i = 1, row
        do j = 1, col
          k = k + 1
          write(num2(k),'(i10)') num(i,j)
          num2(k) = adjustl(num2(k))
        end do
      end do
      write(r,'(i4)') row
      write(c,'(i4)') col
      write(rc,'(i4)') 5 + 2*row*col + 2
      if (info(1:5).eq.'empty') info =''
      if (sep.ne.' ') then
        info=trim(info)//' delimiter="'//sep//'" '
      endif
      ts = t_spaces(indent)
      if (indent.eq.0) then
        write(fmt_1,'(3a)') '(',trim(adjustl(rc)),'a)'
      else
        write(fmt_1,'(5a)') '(',trim(ts),',',trim(adjustl(rc)),'a)'
      end if
      write(lfn,fmt=fmt_1) '<matrix dataType="xsd:integer" ',trim(info),        &
     & ' columns="',trim(adjustl(c)),'" rows="',trim(adjustl(r)),'">',          &
     &   (trim(num2(k)),sep,k=1,row*col-1),trim(num2(row*col)),                 &
     &                     '</matrix>'

      end subroutine writetag_matrix_integer
