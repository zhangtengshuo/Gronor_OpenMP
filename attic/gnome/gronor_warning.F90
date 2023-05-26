subroutine gronor_warning(icode,string)

  use mpi
  use cidef
  use cidist

  implicit none

  character*(*) :: string
  integer :: icode
  character (len=255) :: filwrn
  !
  !     error termination
  !
  !     string = error message printed to filwrn
  !     icode  = informative value printed to filwrn
  !
  write(filwrn,'(a,i5.5,a)') "GronOR_",me,".wrn"
  open(unit=lfnwrn,file=trim(filwrn),form="formatted",position="append")

  select case(icode)

  case(0:99)
    write(lfnwrn,'(a)') "Input: Groups : Number of groups not specified"

  case(102)
    write(lfnwrn,'(a,a)')  "Input: Size : Number of ranks per groups not supported ", &
        "by the total number of available ranks"

  case(104)
    write(lfnwrn,'(a)') "Input: MEBFs : Root name of MEBFs not specified"

  case(106)
    write(lfnwrn,'(a)') "Input: MEBFs : Number of MEBF spin states not specified"

  case(107)
    write(lfnwrn,'(a)') "Input: MEBFs : Line with MEBF spin states could not be read"

  case(108)
    write(lfnwrn,'(a)') "Input: MEBFs : MEBF spin state could not be read"

  case(110)
    write(lfnwrn,'(a)') "Input: Spin : Total spin not specified"

  case(111)
    write(lfnwrn,'(a,i4,a,a)') "Execution aborted with error code ",icode,": ",string

  end select

  close(unit=lfnwrn)

  return
end subroutine gronor_warning
