      subroutine errquit(irank,icode,string)
      implicit none
      character*(*) :: string
      integer :: icode
      integer (kind=4) :: irank
!
!     error termination
!
!     string = error message printed to stdout
!     icode  = informative value printed to stdout
!
      write(*,'(i6,a,i6,a,a)') irank," error ",icode,": ",string
!
      stop
      return
      end
