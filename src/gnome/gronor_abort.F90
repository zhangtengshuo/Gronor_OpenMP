      subroutine gronor_abort(icode,string)

      use mpi
      use cidef
      use cidist

      implicit none

      character*(*) :: string
      integer :: icode
      integer (kind=4) :: ierror,ierr
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
         write(lfnabt,'(a,i6,a,a)') "Aborting with code ",icode,": ",string

      case(101)
         write(lfnabt,'(a)') "Input: Groups : Number of groups not specified"

      case(102)
         write(lfnabt,'(a)') "Input: Size : Number of ranks per groups not specified"

      case(103)
         write(lfnabt,'(a,a)') "Input: Size : Number of ranks per groups not supported ", &
             "by the total number of available ranks"

      case(104)
         write(lfnabt,'(a)') "Input: Fragments : Number of fragments not specified"

      case(105)
         write(lfnabt,'(a)') "Input: MEBFs : Root name of MEBFs not specified"

      case(106)
         write(lfnabt,'(a)') "Input: MEBFs : Number of MEBF spin states not specified"

      case(107)
         write(lfnabt,'(a)') "Input: MEBFs : Line MEBF spin states could not be read"

      case(108)
         write(lfnabt,'(a)') "Input: MEBFs : MEBF spin state could not be read"

      case(109)
         write(lfnabt,'(a)') "Input: MEBFs : MEBF 1st dim fragname/state too small"

      case(110)
         write(lfnabt,'(a)') "Input: Spin : Total spin not specified"

      case(111)
         write(lfnabt,'(a)') "Input: Couplings : line with couplings could not be read"

      case(112)
         write(lfnabt,'(a)') "Input: Threshold : CI threshold not specified"

      case(113)
         write(lfnabt,'(a)') "Input: Couplings : no couplings specified"

      case(114)
         write(lfnabt,'(a)') "Input: Expert : expert option not specified"

      case(115)
         write(lfnabt,'(a)') "Input: Abort : abort option not specified"

      case(116)
         write(lfnabt,'(a)') "Input: Thresh_SIN : singularoty threshold not specified"

      case(117)
         write(lfnabt,'(a)') "Input: Task : number of tasks not specified"

      case(118)
         write(lfnabt,'(a)') "Input: Load : load balance option not specified"

      case(119)
         write(lfnabt,'(a)') "Input: Batch : Size of batch not specified"

      case(120)
         write(lfnabt,'(a)') "Input: Print : print option not properly specified"

      case(121)
        write(lfnabt,'(a)') "Input: Label : label length maximum not properly specified"
        
      case(122)
        write(lfnabt,'(a)') "Input: Managers : number of managers not properly specified"
        
      case(123)
        write(lfnabt,'(a)') "Input: Solver : solver option not specified"

      case(130)
         write(lfnabt,'(a)') "Input: MPIbuffer : MPI buffer size not specified"

      case(140)
         write(lfnabt,'(a)') "Input: Ecorr : Correlation energies not specified"

      case(150)
         write(lfnabt,'(a,a,a)') "Input: Input file ",trim(string)," could not be opened"

      case(200)
         write(lfnabt,'(a)') "gronor_main: Groups size incompatible with available ranks"

      case(201)
         write(lfnabt,'(a)') "gronor_main: number of groups is zero"

      case(202)
         write(lfnabt,'(a)') "gronor_main: number of tasks is zero"

      case(203)
         write(lfnabt,'(a)') "gronor_main: number of rasks exceeds available"

      case(204)
         write(lfnabt,'(a,a)') "gronor_main: error opening timing file ",trim(string)

      case(205)
         write(lfnabt,'(a,a)') "gronor_main: error opening log file ",trim(string)

      case(206)
         write(lfnabt,'(a,a)') "gronor_main: error opening dbg file ",trim(string)

      case(207)
         write(lfnabt,'(a,a)') "gronor_main: error opening sys file ",trim(string)

      case(208)
         write(lfnabt,'(a,a)') "gronor_main: error opening output file ",trim(string)

      case(209)
         write(lfnabt,'(a,a)') "gronor_main: error opening vector file ",trim(string)

      case(210)
         write(lfnabt,'(a,a)') "gronor_main: error opening civ file ",trim(string)

      case(220)
         write(lfnabt,'(a)')  "gronor_makebasestate: too many active orbitals"

      case(221)
         write(lfnabt,'(a)')  "gronor_makebasestate: incompatible spins"

      case(230)
         write(lfnabt,'(a,a)') "gronor_number_integrals: error reading ",trim(string)

      case(231)
         write(lfnabt,'(a,a)') "gronor_number_integrals: error opening ",trim(string)

      case(240)
         write(lfnabt,'(a,a)') "gronor_parallel_integral_input: error reading", &
             " one-electron integrals from ",trim(string)

      case(241)
         write(lfnabt,'(a,a)') "gronor_parallel_integral_input: error reading", &
             " one-electron integrals from ",trim(string)

      case(242)
         write(lfnabt,'(a,a)') "gronor_parallel_integral_input: error reading", &
             " two-electron integrals from ",trim(string)

      case(243)
         write(lfnabt,'(a,a)') "gronor_parallel_integral_input: error reading", &
             " two-electron integrals from ",trim(string)

      case(250)
         write(lfnabt,'(a,a)') "gronor_read_integrals: error reading two-electron", &
             " integrals from ",trim(string)

      case(251)
         write(lfnabt,'(a,a)') "gronor_read_integrals: error reading two-electron", &
             " integrals from ",trim(string)
         
      case(252)
         write(lfnabt,'(a,a)') "gronor_read_integrals: integrals not in canonical order on", &
             trim(string)

      case(260)
         write(lfnabt,'(a,a)') "gronor_read_vectors_and_determinants: ", &
             "inconsistent occupation"

      case(261)
         write(lfnabt,'(a,a,a)') "gronor_read_vectors_and_determinants: ", &
             "error reading system file ",trim(string)

      case(262)
         write(lfnabt,'(a,a,a)') "gronor_read_vectors_and_determinants: ", &
             "error reading vector file ",trim(string)

      case(263)
         write(lfnabt,'(a,a,a)') "gronor_read_vectors_and_determinants: ", &
             "error reading determinant file ",trim(string)

      case(270)
         write(lfnabt,'(a)')  "gronor_detemine_maxcib: incompatible spin state"
      case(271)
         write(lfnabt,'(a)')  "gronor_detemine_maxcib: incompatible spin state"
      case(272)
         write(lfnabt,'(a)')  "gronor_detemine_maxcib: incompatible spin state"

      case(290)
         write(lfnabt,'(a)')  "gronor_leftovers: problem in seticlosed"

      case(300)
         write(lfnabt,'(a)')  "gronor_calculate: inconsistent number of electrons"
      case(305)
         write(lfnabt,'(a)')  "gronor_gnome: inconsistent number of electrons"
      case(306)
         write(lfnabt,'(a)')  "gronor_gnome: incompatible nveca"
      case(307)
         write(lfnabt,'(a)')  "gronor_gnome: incompatibler nvecb"
      case(310)
         write(lfnabt,'(a)')  "gronor_cofac1: no overlap between MOs"
      case(320)
         write(lfnabt,'(a)')  "gronor_moover: inconsistent number of electron spins"

      case(330)
         write(lfnabt,'(a)')  "gronor_lowdin: error in dsyevs"

      case(400)
         write(lfnabt,'(a,a)') 'Incompatible spin coupling: ',string

      case DEFAULT
         write(lfnabt,'(a,i6,a,a)') "Execution aborted with error code ",icode,": ",string

      end select

      close(unit=lfnabt)
!
      ierror=0
      call MPI_Abort(MPI_COMM_WORLD,ierror,ierr)
!
      return
      end
