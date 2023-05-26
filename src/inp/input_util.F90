      module inp

      implicit none

!     Maximum input width of an input line
      integer (kind=8), parameter:: max_width=1024
!     Maximum no. of fields in an input line
      integer (kind=8), parameter :: max_field = max_width/2 + 1
      character (len=1024) :: ja, ia ! Input buffers ... MUST match max_width
      character (len=1024) :: tmp    ! Same size work space
      character (len=80) :: errmsg   ! Error message
      character (len=1) :: xcomm     ! Comment character
      character (len=1) :: xsplit    ! Character to split physical input lines
      character (len=1) :: xback     ! Backslash for concatenation and quoting
      character (len=1) :: xquote    ! Quotation marks for strings
      character (len=1) :: xblnk     ! Space
      character (len=1) :: xtab      ! Tab

      integer (kind=8) :: iread, iwrite
      integer (kind=8) :: jrec       ! No. of current field
      integer (kind=8) :: jump       ! No. of fields in current line
      integer (kind=8) :: istrt(max_field) ! Start of fields
      integer (kind=8) :: inumb(max_field) ! Length of fields
      integer (kind=8) :: nend(max_field)  ! End of fields
      integer (kind=8) :: iwidth     ! Length of current logical input line
      integer (kind=8) :: nline      ! Current logical line inside physical line
      integer (kind=8) :: noline     ! No. of logical lines inside physical line
      integer (kind=8) :: input_line ! No. of current physical input line
      integer (kind=8) :: nerr       ! ?
      logical :: oswit                  ! True if EOF has beeen detected
      integer (kind=8) :: ierrpos    ! Input char position where error detected
      integer (kind=8) :: nstart(max_field)

      integer (kind=8), parameter :: max_include_level=3
      integer (kind=8) :: include_level
      character (len=255) :: include_file_name(max_include_level)

      data iread   /5/
      data iwrite  /6/
      data jrec    /-1/
      data jump    /0/
      data oswit   /.false. /
      data nerr    /999/
      data nline   /0/
      data noline  /0/
      data ierrpos /-1/
      data errmsg  /' '/
      data input_line /0/
      data xblnk /' '/
#if defined(DECOSF) || defined(SGITFP) || defined(SGI_N32) || defined(CRAY) || defined(HPUX) || defined(WIN32) || defined(PSCALE) ||defined(GCC4)
      data xtab  /'	'/      ! Tab ... no backslash necessary
#elif (defined(LINUX) || defined(MACX)) && !defined(PGLINUX) && !defined(XLFLINUX) &&!defined(GCC4)
      data xtab  /'9'/            ! Tab ... g77 has trouble with escape sequence
#else
      data xtab  /'\	'/      ! Tab ... note backslash for cpp
#endif
      data xsplit/';'/
      data xcomm /'#'/
      data xback /'\\'/         ! Backslash ... note backslash for cpp
      data xquote/'"'/

      data include_level /0/

      contains

      subroutine inp_init(ir,iw)
      implicit none

      integer (kind=8) :: ir
      integer (kind=8) :: iw

      iread = ir
      iwrite = iw
      jrec = -1
      jump = 0
      oswit = .false.
      nerr = 999
      nline = 0
      noline = 0
      ierrpos = -1
      errmsg = ' '
      input_line = 0

      end
      integer (kind=8) function inp_n_field()
      implicit none

!     return no. of fields in the input line ... 0 = EOF

      inp_n_field = jump

      end

      integer (kind=8) function inp_cur_field()
      implicit none

!     return no. of fields processed so far (0,...,inp_n_field())

      inp_cur_field = jrec

      end

      subroutine inp_set_field(ivalue)
      implicit none

      integer (kind=8) :: ivalue

!     set field to be read next (ivalue=0,...,inp_n_field())

      if (ivalue.lt.0 .or. ivalue.gt.inp_n_field()) &
       call inpquit('inp_set_field: stupid field value',ivalue)
      jrec = ivalue

      ierrpos = -1
      errmsg = ' '

      end

      logical function inp_line(z)
      implicit none

      character (len=*) :: z

!     set the variable z to be as much of the current input line
!     that it can hold

      if (jump .gt. 0) then
         z = ia
         inp_line = .true.
         ierrpos = -1
         errmsg = ' '
      else
         errmsg = 'no input line available'
         ierrpos = -1
         inp_line = .false.
      endif

      end
      logical function inp_read_physical_line(buf)
      implicit none
!#include "inpP.fh"
      character (len=*) :: buf

!     Read a physical line of input into buf() WITHOUT any
!     tokenizing, scanning, etc.

!     This routine is intended for programms embedded into NWChem
!     that need to read from the NWChem input file bypassing any
!     processing that NWChem does on input data, but while still
!     maintaining a correct count of the input line number for
!     error reporting.

!     First void out any info about the current input line

      ierrpos = -1
      jump = 0
      jrec = 0

      buf = ' '
      read(5,'(a)',end=10,err=10) buf
      input_line = input_line + 1

      inp_read_physical_line = .true.
      return

 10   inp_read_physical_line = .false.
      errmsg  = 'unexpected end of data file'

      end

      logical function inp_read()
      implicit none

!     this routine reads a data card and scans it for non - space fields
!     the number of fields is stored in jump, the starting point of a
!     field in istrt(i) and the number of characters in that field
!     in inumb(i).

      character (len=1) :: xprev
      integer (kind=8) :: lenja, i, k, mark, j, nbegin, nfini, iopt
      integer (kind=8) :: jwidth
      integer (kind=8) :: ncol(max_field)
      logical :: ois_ws            ! Inline funtion
      character (len=1) :: xtest
      integer (kind=8) :: ivalue

      ois_ws(xtest) = (iachar(xtest).eq.32 .or. iachar(xtest).eq.9)

      inp_read = .false.  ! First assume things will not work

 1    nline=nline+1
      if(nline.le.noline)go to 150

      if (oswit) then
         ierrpos = -1
         errmsg = 'unexpected end of data file'
         inp_read = .false.
         jump=0
         jrec=0
         return
      else
         ierrpos = -1
         errmsg = ' '
      endif

!     read next physical input line

 101  lenja = 0
 100  read(iread,'(a)',end=300)ja(lenja+1:max_width)
      input_line = input_line + 1
      lenja=len(trim(ja))

!     Check for . * eof at beginning of line to indicate EOF

      if (lenja.eq.1 .and. (ja(1:1).eq.'.' .or. ja(1:1).eq.'*')) goto 300
      if (lenja.eq.3 .and. inp_compare(.false., 'eof', ja(1:3))) goto 300

!     Handle include statement

      if (inp_compare(.false.,ja(1:7),'include')) then
         if (include_level .eq. max_include_level) & 
             call inpquit('inp_read: include nested too deep ', include_level)
         include_level = include_level + 1
         include_file_name(include_level) = ja(9:)
         write(6,*) ' include: start of ',trim(include_file_name(include_level))
         open(80+include_level,file=trim(include_file_name(include_level)), &
             form='formatted', status='old', err=105)
!     call inp_save_state()
         iopt=6
         call inp_init(80+include_level,iopt)
         goto 1
      endif

!     handle blank lines and concatenation using backslash

      if (lenja.eq.0) then
         goto 100
      else
         if (ja(lenja:lenja) .eq. xback) then
            ja(lenja:lenja) = xblnk
            goto 100
         endif
      endif
      jwidth=len(trim(ja))

!     handle comments from # to eol ... allow for backslash quoting

      xprev = xblnk
      do i=1, jwidth
 91      if (ja(i:i) .eq. xcomm .and. xprev.ne.xback) then
            lenja=len(trim(ja))
            ja(i:max_width) = xblnk
            goto 80
         else if (ja(i:i) .eq. xcomm .and. xprev.eq.xback) then
!     Shuffle string down to overwrite quoting backslash
            tmp = ja
            ja(i-1:jwidth) = tmp(i:jwidth)
            xprev = xblnk
            goto 91
         endif
         xprev = ja(i:i)
      enddo

 80   if (len(trim(ja(1:i))) .eq. 0) goto 101 ! All line comments

!     figure out where ; splits physical line into multiple logical lines
!     again handling quoted backslash

      k=jwidth
      mark=0
      xprev = xblnk
      do i=1,jwidth
 81      if(ja(i:i).eq.xsplit .and. xprev.ne.xback) then
            mark=mark+1
            ncol(mark)=i
         else if(ja(i:i).eq.xsplit .and. xprev.eq.xback) then
!     Shuffle string down to overwrite quoting backslash
            tmp = ja
            ja(i-1:jwidth) = tmp(i:jwidth)
            xprev = xblnk
            goto 81
         endif
         xprev = ja(i:i)
      enddo
      noline=1
      if(mark.eq.0) then
         nstart(noline)=1
         nend(noline)=jwidth
      else
         i=ncol(mark)+1
         if(i.le.jwidth) then
            do j=i,jwidth
               if(.not. ois_ws(ja(j:j))) go to 170
            enddo
         endif
         k=ncol(mark)-1
         mark=mark-1

 170     noline=mark+1
         nstart(1)=1
         do i=1,mark
            j=ncol(i)
            nend(i)=j-1
            nstart(i+1)=j+1
         enddo
         nend(noline)=k
      endif
      nline=1

!     Start processing next logical input line (put into ia(1:iwidth))

 150  jump=0
      jrec=0
      nbegin = nstart(nline)
      nfini  = nend(nline)
      iwidth = nfini-nbegin+1
      ia = xblnk
      ia(1:iwidth)=ja(nbegin:nfini)

!     partition input line into strings inside double quotes or
!     white space separated fields

      i = 1
 151  continue
      do j = i, iwidth
         if (.not. ois_ws(ia(j:j))) goto 152
       enddo
       ! Done
      goto 155
 152  i = j

      jump = jump + 1
      istrt(jump) = i
      if (ia(i:i) .eq. xquote) then

!     Quoted string ... look for closing quote

         do j = i+1, iwidth
 154        if (ia(j:j).eq.xquote .and. ia(j-1:j-1).ne.xback) then
               goto 153
            else if (ia(j:j).eq.xquote .and. ia(j-1:j-1).eq.xback) then
               tmp = ia
               ia(j-1:max_width) = tmp(j:max_width)
               goto 154
            endif
         enddo
         ierrpos = j
         errmsg = 'no terminating quote for string'
         inp_read = .false.
         oswit = .false.
         return
 153     continue
      else

!     Simple field ... look for next ws

         do j = i+1, iwidth
            if (ois_ws(ia(j:j))) goto 157
         enddo
 157     j = j - 1
      endif

      inumb(jump) = j - istrt(jump) + 1
      i = j + 1
      goto 151

 155  continue                  ! Finished tokenizing
      if (jump .gt. 0) iwidth = istrt(jump)+inumb(jump)-1
      inp_read = .true.
      return

 300  if (include_level .gt. 0) then ! End of file detected
         close(80+include_level)
         write(6,*) ' include: end of ',trim(include_file_name(include_level))
         include_level = include_level - 1
!         call inp_restore_state()
         goto 1
      else
         oswit = .true.
         ierrpos = -1
         errmsg = 'unexpected end of data file'
         inp_read = .false.
         jump=0
         jrec=0
      endif
      return

 105  continue

      ivalue=0
      call inpquit('inp_read: failed to open include file',ivalue)

      end
      logical function inp_eof()
      implicit none

      inp_eof = oswit

      end
      subroutine inp_clear_err()
      implicit none

!     Clear error conditions and messages

      ierrpos = -1
      errmsg = ' '

      end
      subroutine inp_errout
      implicit none

      integer (kind=8) :: length, i
      character (len=1) :: xpt, xstp
      data xpt,xstp/'*', '.'/

!     If an error has occured print out the error message
!     and the position in the current input line

      if (include_level .gt. 0) then
         write(6,*) ' Include file stack '
         do i = 1, include_level
            write(6,321) i, trim(include_file_name(i))
 321        format(1x,i5,2x,a)
         enddo
      endif
      if (errmsg .ne. ' ') then
         length=len(trim(errmsg))
         write(iwrite, 40) input_line, errmsg(1:length)
 40      format(' input error at line', i5,': ', a)

         write(iwrite,50)ia(1:iwidth)
 50      format(1x,a)
         if (ierrpos .gt. 0) then
            do 60 i=1,iwidth
               tmp(i:i)=xstp
 60         continue
            tmp(ierrpos:ierrpos)=xpt
            write(iwrite,50)tmp(1:iwidth)
         endif
      endif

      end
      subroutine inp_outrec
      implicit none

!     Write out the current input line

      write(iwrite,50) input_line, ia(1:iwidth)
 50   format(1x,i5,': ',a)

      end
      logical function inp_a(a)
      implicit none

      integer (kind=8) :: i1, i2, length
      character (len=*) :: a

!     Return field as character string, minus any enclosing quotes
!     with an error if it does not fit

      ierrpos = -1
      errmsg = ' '
      if(jrec .ge. jump) then
         a = xblnk
         inp_a = .false.
         ierrpos = 0
         errmsg = 'at end of line looking for character string'
         return
      endif
      i1 = istrt(jrec+1)
      i2 = istrt(jrec+1)+inumb(jrec+1)-1
      if (ia(i1:i1).eq.xquote .and. ia(i2:i2).eq.xquote) then
         i1 = i1+1
         length = inumb(jrec+1)-2
      else
         length = inumb(jrec+1)
      endif
      if (len(a) .lt. length) then
         a = xblnk
         inp_a = .false.
         ierrpos = 0
         errmsg = 'inp_a: string is too large for argument'
         return
      else
         jrec = jrec + 1
         a = ia(i1:i1+length-1)
         inp_a = .true.
         return
      endif

      end
      logical function inp_a_trunc(a)
      implicit none

      integer (kind=8) :: i1, i2, length
      character (len=*) :: a

!     Return field as character string, minus any enclosing quotes
!     quietly truncating if it does not fit

      ierrpos = -1
      errmsg = ' '
      if(jrec .ge. jump) then
         a = xblnk
         inp_a_trunc = .false.
         ierrpos = 0
         errmsg = 'at end of line looking for character string'
         return
      endif
      i1 = istrt(jrec+1)
      i2 = istrt(jrec+1)+inumb(jrec+1)-1
      if (ia(i1:i1).eq.xquote .and. ia(i2:i2).eq.xquote) then
         i1 = i1+1
         length = inumb(jrec+1)-2
      else
         length = inumb(jrec+1)
      endif
      jrec = jrec + 1
      a = ia(i1:i1+length-1)
      inp_a_trunc = .true.
      return

      end
      logical function inp_f (buf)
      implicit none

      double precision ten, buf, dtmp
      integer (kind=8) :: i1, i2, ie2, isign, ie, iexp, ie1, itmp, i, j
      logical :: orep
      character (len=1) :: xchar(17)
      data xchar /'0','1','2','3','4','5','6','7','8','9','+','-','.','e','d','E','D'/
      data ten/10.0d0/

      ierrpos = -1
      errmsg = ' '
      dtmp=0.0d0
      if (jrec.ge.jump) then
         inp_f = .false.
         errmsg = 'at end of line looking for floating point number'
         ierrpos=-1
         return
      endif
      jrec=jrec+1
      i1=istrt(jrec)
      i2=i1+inumb(jrec)-1
      ie2=i2
!...  sign
      isign=1
      if (ia(i1:i1).eq.xchar(12))isign=-1
      if (ia(i1:i1).eq.xchar(12).or.ia(i1:i1).eq.xchar(11)) i1=i1+1
!...  exponent
      do ie=i1+1,i2
         if (ia(ie:ie).eq.xchar(14) .or. ia(ie:ie).eq.xchar(15) .OR. &
             ia(ie:ie).eq.xchar(16) .or. ia(ie:ie).eq.xchar(17)) goto 20
      enddo
      iexp=0
      go to 50
 20   i2=ie-1
      iexp=1
      ie1=ie+1
      if (ia(ie1:ie1).eq.xchar(12))iexp=-1
      if (ia(ie1:ie1).eq.xchar(12).or.ia(ie1:ie1).eq.xchar(11)) ie1=ie1+1
      itmp=0
      do i=ie1,ie2
         do j=1,10
            if (ia(i:i).eq.xchar(j)) go to 41
         enddo
         goto 100
 41      itmp=itmp*10+j-1
      enddo
      iexp=iexp*itmp
!.... the number itself
 50   orep=.false.
      do i=i1,i2
         if(ia(i:i).ne.xchar(13)) then
            do j=1,10
               if (ia(i:i).eq.xchar(j)) go to 70
            enddo
            goto 100
 70         dtmp=dtmp*ten+ dble(j-1)
         else
            if(orep)go to 100
            iexp=iexp+i-i2
            orep=.true.
         endif
      enddo
      dtmp = dtmp * dble(isign) * ten**iexp
      inp_f = .true.
      buf = dtmp
      return

 100  inp_f = .false.
      jrec = jrec-1             ! Position to re-read the field
      ierrpos = i
      errmsg = 'illegal character reading floating point number'

      end
      subroutine inp_mark_err(message)
      implicit none

      character (len=*) :: message

!     Mark an input error at the beginning of the current input field

      ierrpos = istrt(min(max_field,max(1,jrec)))
      errmsg  = message

      end
      logical function inp_i(jbuf)
      implicit none

      character (len=1) :: xchar(12)
      integer (kind=8) :: n, ifact, ist, nstrt, i, j
      character (len=1) :: xtemp
      integer :: jbuf, jtmp
      data xchar /'0','1','2','3','4','5','6','7','8','9','+','-'/

!     subroutine for reading integers from the array ia,
!     starting at ia(istrt(jrec)) and going on for inumb(jrec))
!     elements. plus signs are ignored, the answer is accumulated
!     in jtmp

      ierrpos = -1
      errmsg = ' '
      jtmp = 0
      if(jrec.ge.jump) then
         inp_i = .false.
         ierrpos = -1
         errmsg = 'at end of line looking for integer'
         return
      endif
      jrec = jrec + 1
      n = inumb(jrec)
      ifact = 1
      ist=istrt(jrec)
      nstrt = ist + n - 1
      do i = 1,n
         xtemp = ia(nstrt:nstrt)
         do j=1,12
            if(xchar(j).eq.xtemp)go to 130
         enddo
         goto 120

 130     if(j.ge.11) then
            if(nstrt.ne.ist)go to 120
            if(j.ge.12)jtmp=-jtmp
            go to 160
         endif
         jtmp=jtmp+(j-1)*ifact
         ifact = ifact * 10
         nstrt=nstrt-1
      enddo
 160  continue
      inp_i = .true.
      jbuf = jtmp
      return

 120  ierrpos = nstrt
      errmsg  = 'illegal character when reading integer'
      inp_i = .false.
      jrec = jrec-1
      return

      end
      logical function inp_compare(ocase, a, b)
      implicit none
      logical :: ocase
      character (len=*) :: a, b
      integer (kind=8) :: la, lb, i
      character (len=1) :: atest, btest

      inp_compare = .false.
      la=len(trim(a))
      lb=len(trim(b))
      if (la .ne. lb) then      ! use .gt. for short match
         return
      else if (ocase) then
         inp_compare = a(1:la) .eq. b(1:lb)
         return
      else
         do i = 1, la
            atest = a(i:i)
            btest = b(i:i)
            call inp_lcase(atest)
            call inp_lcase(btest)
            if (atest.ne.btest) return
         enddo
         inp_compare = .true.
         return
      endif

      end
      logical function inp_match(nrec, ocase, test, array, ind)
      implicit none
      integer (kind=8) :: nrec, ind
      logical :: ocase
      character (len=*) :: test, array(*)
      integer (kind=8) :: i, j, l

      l=len(trim(test))
      inp_match = .false.
      ind = -1

      do i=1,nrec
         if (inp_compare(ocase, test(1:l), array(i))) then
            if (inp_match) then
               inp_match = .false. ! Ambiguity
               ind = 0
               write(6,1) test(1:l), (array(j),j=1,nrec)
 1             format('inp: ambiguous match for ', a,', in:',/,100(1x,a/))
               return
            else
               inp_match = .true. ! First match
               ind = i
            endif
         endif
      enddo

      end
      logical function inp_contains(ocase, a, b,ipos)

!     check if string a is contained in b. return starting
!     location of string a in b.

      implicit none
      logical :: ocase
      character (len=*) :: a, b
      integer (kind=8) :: la, lb, i, j, ipos
      character (len=1) :: atest, btest

      ipos = -1
      inp_contains = .false.
      la=len(trim(a))
      lb=len(trim(b))
      if (la .gt. lb) then
         return
      else if (ocase) then
         do i = 1, lb - la + 1
            inp_contains = a(1:la) .eq. b(i:i+la)
            if (inp_contains) then
               ipos = i
               return
            endif
         enddo
         return
      else
         do j = 0, lb - la
            do i = 1, la
               atest = a(i:i)
               btest = b(j+i:j+i)
               call inp_lcase(atest)
               call inp_lcase(btest)
               if (atest.ne.btest) goto 11
            enddo
            inp_contains = .true.
            return
   11    continue
         enddo
         return
      endif
!
      end
      subroutine inp_prev_field()
      implicit none
      integer (kind=8) :: ivalue

      ivalue=max(0,inp_cur_field()-1)
!      call inp_set_field(max(0,inp_cur_field()-1))
      call inp_set_field(ivalue)

      end

      integer (kind=8) function inp_strlen(a)
      implicit none

      character (len=*) :: a
      integer (kind=8) :: i,j
      logical :: ois_ws
      intrinsic len
      character (len=1) :: xtest
      ois_ws(xtest) = (iachar(xtest).eq.32 .or. iachar(xtest).eq.9)

      do i = len(a),1,-1
         j=i
         if (.not. ois_ws(a(i:i))) goto 10
      enddo

      i=j
 10   inp_strlen = i

      end

      subroutine inp_lcase(string)
      implicit none

      character (len=*) :: string
      intrinsic ichar, len
      integer (kind=8) :: i, length, uca, ucz, lca, shift, test
      integer (kind=8) :: ivalue

      uca = ichar('A')          ! MUST be uppercase A
      ucz = ichar('Z')          ! MUST be uppercase Z
      lca = ichar('a')          ! MUST be lowercase a
      shift = lca - uca

      ivalue=0
      if (shift .eq. 0) &
          call inpquit('inp_lcase: check case of program source',ivalue)

      length = len(string)
      do i = 1, length
         test = ichar(string(i:i))
         if (test.ge.uca .and. test.le.ucz) string(i:i) = char(test+shift)
      enddo

      end
      subroutine inp_adjustl(a)
      implicit none

      character (len=*) :: a
      intrinsic  len
      logical :: ois_ws
      character (len=1) :: xtest
      integer (kind=8) :: i, length,s

      ois_ws(xtest) = (iachar(xtest).eq.32 .or. iachar(xtest).eq.9)

      length = len(a)
      do i = 1, length
        if (.not. ois_ws(a(i:i))) goto 10
      end do
 10    continue
      s=i-1
      do i=s+1,length
       a(i-s:i-s) = a(i:i)
       a(i:i)=""
      enddo

      end
      subroutine inp_ucase(string)
      implicit none

      character (len=*) :: string
      intrinsic ichar, len
      integer (kind=8) :: i, length, lca, lcz, uca, shift, test
      integer (kind=8) :: ivalue

      lca = ichar('a')          ! MUST be lowercase A
      lcz = ichar('z')          ! MUST be lowercase Z
      uca = ichar('A')          ! MUST be uppercase a
      shift = uca - lca

      ivalue=0
      if (shift .eq. 0) call inpquit('inp_ucase: check case of program source',ivalue)

      length = len(string)
      do i = 1, length
         test = ichar(string(i:i))
         if (test.ge.lca .and. test.le.lcz) string(i:i) = char(test+shift)
      enddo

      end
      logical function inp_search(ocase, z, nz)
      implicit none

      integer (kind=8) :: nz
      character (len=*) :: z(nz)
      logical :: ocase
      character (len=1024) :: tmp

      integer (kind=8) :: i
      integer (kind=8) :: maxz
      parameter (maxz = 100)
      integer (kind=8) :: length(maxz)

      if (maxz .lt. nz) call inpquit('inp_search: hard dim fail',nz)
      do i = 1, nz
         length(i)=len(trim(z(i)))
      enddo

 10   if (inp_read()) then
         if (inp_a(tmp)) then
            do i = 1, nz
               if (inp_compare(ocase, z(i)(1:length(i)), tmp)) then
                  call inp_prev_field()
                  inp_search = .true.
                  return
               endif
            enddo
         endif
         goto 10
      endif

      inp_search = .false.

      end
      logical function inp_search_fast(z)
      implicit none

      character (len=*) :: z
      integer (kind=8) :: length
      integer (kind=8) :: ivalue

!     Quicker search that

!     1) matches case
!     2) assumes the token being searched for is at the beginning of the line
!     3) ignores continutation lines, comments, quotes etc.
!     4) Still attempts to track line numbers EOF.

!     Only called from inside the basis set input routine?

      inp_search_fast = .false.
      length = max(3,len(trim(z))) ! 3 for EOF/eof detection

 10   read(iread,'(a)',end=300) ja(1:length)
      input_line = input_line + 1
      if (ja(1:3).eq.'EOF' .or. ja(1:3).eq.'eof') goto 300
      if (z(1:length) .eq. ja(1:length)) then
         backspace(iread)       ! Re-read line with full input routine
         input_line = input_line - 1
         ivalue=0
         if (.not. inp_read()) call inpquit('inp_search_fast: inp?',ivalue)
         inp_search_fast = .true.
         return
      endif
      goto 10

 300  oswit = .true.            ! EOF code copied from inp_read
      ierrpos = -1
      errmsg = 'unexpected end of data file'
      jump=0
      jrec=0

      end

      logical function inp_strtok(z, sep, istart, iend)
      implicit none
      character (len=*) :: z           ! [input] string to parse
      character (len=*) :: sep         ! [input] token separators
      integer (kind=8) :: istart, iend      ! [output] start/end of next token

!     Returns the number of the start and end character of the
!     next token in the character string.  Tokens are separated
!     by one of the characters in sep.  Note that all characters
!     in sep are used including any trailing blanks.

!     Before the first call initialize istart to zero, and leave
!     istart and iend UNCHANGED for subsequent calls.
!     Repeated calls return the next token and true, or false if
!     there are no more tokens.  The separators may be changed
!     between calls.

!     No internal state is maintained (which is istart and iend
!     must not be modified between calls) so multiple strings
!     may be parsed simultaneously.

!     E.g., to split list = 'robert:rick:jeff' into tokens separated
!     by ':'. You execute

!     istart = 0
!  10 if (inp_strtok(list, ':', istart, iend)) then
!     write(6,*) list(istart:iend)
!     goto 10
!     endif

      integer (kind=8) :: i, k, length, nsep

      if (istart .eq. 0) then
         istart = 1
      else
         istart = iend + 1
      endif

!     Scan start forward to next non-separator

      length = len(z)
      nsep   = len(sep)

      do i = istart, length
         do k = 1, nsep
            if (z(i:i) .eq. sep(k:k)) goto 10
         enddo
         goto 20
 10      continue
      enddo
      inp_strtok = .false.      ! No more tokens
      return

 20   istart = i                ! Beginning of next token

!     Scan end forward to one-before next separator

      do i = istart+1, length
         do k = 1, nsep
            if (z(i:i) .eq. sep(k:k)) goto 30
         enddo
      enddo
 30   iend = i - 1

      inp_strtok = .true.

      return
      end function inp_strtok

      logical Function Inp_IList( MaxList, List, N)
      Implicit NONE
      integer (kind=8) :: MaxList           ! [in]    Size of List
      integer (kind=8) :: List( MaxList )   ! [inout] Contents of list
      integer (kind=8) :: N                 ! [out]   Number of elements in list

!     DESCRIPTION
!     Reads the line for a list of integers and puts
!     the results in List.  Ranges of integers may be compactly using
!     F90-style triplet notation. The number of elements set from the
!     input is returned in N.

!     It USED to read the entire remainder of the input line. BUT now
!     reads until it finds something that is not an integer range.

!     IF (the remainder of the input line was read as a valid list)
!     .  inp_list returns .true. with N as the number of elements read
!     ELSE
!     .  inp_list returns .false. with N as the number of elements
!     .  and also sets the INP internal error message appropriately
!     .  which can be cleared with inp_clear_err

!     N may be returned in either case as zero.

!     If N > MaxList, it indicates that there is too much data on
!     the line to fit in List. NOTE that inp_ilist now returns
!     true in this scenario so the value of N must be checked
!     for this condition.  A recovery action would be to reposition
!     the input to the beginning of the list, allocate a sufficiently
!     large array and then call inp_ilist again.  The input position may
!     be recorded with inp_cur_field() and reset with inp_set_field().
!     The INP error message is also set in this circumstance.

!     LOCAL VARIABLES
      integer (kind=8) :: First, Last, Stride, I

      N = 0

!     Read as many fields as possible as integer ranges

 10   If ( Inp_IRange(First, Last, Stride) ) then

!     Expand the triplet.  If there are too many elements for the List,
!     we count them but do not overwrite the array bounds.

         Do I = First, Last, Stride
            N = N + 1
            If ( N .le. MaxList) List(N) = I
            If ( N .eq. MaxList+1) call inp_mark_err('insufficient space for integer list')
         EndDo
         goto 10
      EndIf

!     OK if we are we are the end of the input line

      inp_ilist = (inp_cur_field() .eq. inp_n_field())

      End

      logical function inp_ilist_size(n)
      implicit none
      integer (kind=8) :: n                 ! [out]   number of elements in list

!     description
!     reads the line of a list of integers and
!     and returns total number of elements
!     designed to be use a helper routine for inp_ilist
!     on exit from this routine the input is repositioned back
!     to the original value

!     local variables
      integer (kind=8) :: first, last, stride, i
      integer (kind=8) :: nc

      nc = inp_cur_field()
      n = 0

!     read as many fields as possible as integer ranges

 10   if ( Inp_IRange(first, last, stride) ) then

!     expand the triplet.  if there are too many elements for the list,
!     we count them but do not overwrite the array bounds.

         do i = first, last, stride
            n = n + 1
         enddo
         goto 10
      endif

      inp_ilist_size = (inp_cur_field() .eq. inp_n_field())
      call inp_set_field(nc)
      end

      logical function inp_irange(JFirst, JLast, JStride)

      implicit none

      integer (kind=8) :: n, ifact, ist, nstrt, i, j
      character (len=1) :: xtemp
      integer (kind=8) :: jFirst, JLast, JStride, jtmp(3), term
      logical :: Expect_Sep, Expect_Dig
      character (len=1) :: xchar(13)
      data xchar /'0','1','2','3','4','5','6','7','8','9','+','-',':'/

!     subroutine for reading an integer range specification (in Fortran90-
!     style triplet notation, stride optional) from the array ia,
!     starting at ia(istrt(jrec)) and going on for inumb(jrec))
!     elements. plus signs are ignored, the answer is accumulated
!     in jtmp.

!     F90 triple notation:  <first>[:<last>[:<stride>]]

!     Note that a simple integer "<first>" will be accepted and interpreted
!     as "<first>:<first>:<stride>"

      ierrpos = -1
      errmsg = ' '

      term = 1
      jtmp(1) = 0
      jtmp(2) = 0
      jtmp(3) = 0

      if(jrec.ge.jump) then
         inp_irange = .false.
         ierrpos = -1
         errmsg = 'at end of line looking for integer range'
         return
      endif
      jrec = jrec + 1
      n = inumb(jrec)
      ifact = 1
      ist=istrt(jrec)
      nstrt = ist + n - 1
      Expect_Sep = .FALSE.
      Expect_Dig = .TRUE.   ! Should be a digit before anything else appears
      do i = 1,n
         xtemp = ia(nstrt:nstrt)
         do j=1,13
            if(xchar(j).eq.xtemp)go to 130
         enddo
         goto 120

 130     Continue
         If ( J .eq. 13) then                    ! ":" separating terms
            If ( Expect_Dig ) Goto 120
            Term = Term + 1
            IFact = 1
            Expect_Dig = .TRUE.   ! Must have a digit next (or EOR)
            Expect_Sep = .FALSE.

         ElseIf(j.eq.11 .OR. J .eq. 12) then     ! +/- signs
            If ( Expect_Sep .OR. Expect_Dig) Goto 120
            if(j.eq.12)jtmp(term) =-jtmp(term)
            Expect_sep = .TRUE. ! Need separator (or EOR) next

         Else                                    ! Add new digit
            If ( Expect_Sep) Goto 120
            jtmp(term) =jtmp(term) +(j-1)*ifact
            ifact = ifact * 10
            Expect_Dig = .FALSE.
         EndIf
         nstrt=nstrt-1
      enddo

!     Finished the loop over the character in this field.  See if
!     we are missing anything (like ending with a : at the front).

      If ( Expect_Dig ) Goto 120

!     All done.  Now we must sort out how much of the triplet we
!     were given.  Remember, the terms are back to front.

 160  continue
      inp_irange = .true.

      If (term .eq. 1) then       ! Just an integer
         JFirst = JTmp(1)
         JLast  = JTmp(1)
         JStride = 1
      ElseIf (Term .eq. 2) then   ! Start and end, stride defaults to 1
         If ( JTmp(2) .gt. JTmp(1) ) Goto 125 ! Stride 1 will not work here
         JFirst = JTmp(2)
         JLast  = JTmp(1)
         JStride = 1
      ElseIf (Term .eq. 3) then   ! Full triplet
         If ( JTmp(1) .eq. 0) Goto 126 ! Stride 0 not allowed
         If ( JTmp(2) .ne. JTmp(3) ) then
            If ( (JTmp(2) - JTmp(3)) * JTmp(1) .lt. 0) Goto 127
         EndIf
         JFirst = JTmp(3)
         JLast  = JTmp(2)
         JStride = JTmp(1)
      EndIf

      return

 120  ierrpos = nstrt
      errmsg  = 'illegal character when reading integer range'
      inp_irange = .false.
      jrec = jrec-1
      return

 125  ierrpos = ist
      errmsg = 'invalid range spec. -- first > last'
      Inp_IRange = .FALSE.
      JRec = JRec - 1
      Return

 126  ierrpos = ist
      errmsg = 'invalid range spec. -- stride can not be 0'
      Inp_IRange = .FALSE.
      JRec = JRec - 1
      Return

 127  ierrpos = ist
      errmsg = 'invalid range spec. -- stride incompatable with range'
      Inp_IRange = .FALSE.
      JRec = JRec - 1
      Return

      end

      subroutine inpquit(string, icode)
      implicit none
      character (len=*), optional :: string
      integer (kind=8), optional :: icode

!     error termination

!     string = error message printed to stdout
!     icode  = informative value printed to stdout

!     try to minimize the amount of output if all nodes
!     detect an error

      write(*,999) icode,string
 999  format("Error code ",i2,": ",a)
      stop
      return
      end

      end module inp
