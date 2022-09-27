      subroutine gronor_abort(icode,string)

      use mpi
      use cidef
      use cidist
      
      implicit none

      character*(*) :: string
      integer :: icode
      integer (kind=4) :: ierror,ierr
      integer (kind=4) :: irank
      character (len=255) :: filabt
!     
!     error termination 
!     
!     string = error message printed to filabt
!     icode  = informative value printed to filabt
!
      write(filabt,'(a,i5.5,a)') "GronOR_",me,".abt"
      open(unit=lfnabt,file=trim(filabt),form="formatted")

      select case(icode)

      case(0:99)
         write(lfnabt,'(a,i6,a,a)')
     &        "Aborting with code ",icode,": ",string

      case(101)
         write(lfnabt,'(a)')
     &        "Input: Groups : Number of groups not specified"
         
      case(102)
         write(lfnabt,'(a)')
     &        "Input: Size : Number of ranks per groups not specified"
         
      case(103)
         write(lfnabt,'(a,a)')
     &        "Input: Size : Number of ranks per groups not supported ",
     &        "by the total number of available ranks"
         
      case(104)
         write(lfnabt,'(a)')
     &        "Input: Fragments : Number of fragments not specified"
         
      case(105)
         write(lfnabt,'(a)')
     &        "Input: MEBFs : Root name of MEBFs not specified"
         
      case(106)
         write(lfnabt,'(a)')
     &        "Input: MEBFs : Number of MEBF spin states not specified"
         
      case(107)
         write(lfnabt,'(a)')
     &        "Input: MEBFs : Line MEBF spin states could not be read"
         
      case(108)
         write(lfnabt,'(a)')
     &        "Input: MEBFs : MEBF spin state could not be read"
         
      case(110)
         write(lfnabt,'(a)')
     &        "Input: Spin : Total spin not specified"
         
      case(111)
         write(lfnabt,'(a)')
     &        "Input: Couplings : line with couplings could not be read"
         
      case(112)
         write(lfnabt,'(a)')
     &        "Input: Threshold : CI threshold not specified"




         
      case DEFAULT
         write(lfnabt,'(a,i6,a,a)')
     &        "Execution aborted with error code ",icode,": ",string
         
      end select
      
      close(unit=lfnabt)
!
      ierror=0
      call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
!
      return
      end
