      subroutine errquit(string, icode, errcode)
      implicit none
      character*(*) string
      integer icode
      integer errcode
c     
c     error termination 
c     
c     string = error message printed to stdout
c     icode  = informative value printed to stdout
c
c     try to minimize the amount of output if all nodes 
c     detect an error
c
      stop
      return
      end
