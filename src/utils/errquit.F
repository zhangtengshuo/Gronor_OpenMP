      subroutine errquit(irank,icode,string)
      implicit none
      character*(*) :: string
      integer :: icode
      integer (kind=4) :: irank
c     
c     error termination 
c     
c     string = error message printed to stdout
c     icode  = informative value printed to stdout
c
      write(*,'(i6,a,i6,a,a)') irank," error ",icode,": ",string
c
      stop
      return
      end
