!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ bases, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!> @brief
!! Setup the parallel computing environment
!!
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!

subroutine gronor_environment()

  use mpi
  use inp
  use cidist

#ifdef _OPENACC
  use openacc
#ifdef CUDA
  use cuda_functions
#endif
  !      use cuda_rt
#endif
#ifdef _OPENMP
  use omp_lib
#endif

  implicit none

  external :: hostnm
!  external :: MPI_AllReduce

  integer :: i,j,k,node
  integer (kind=4) :: length, ierr, ncount
  !      integer (kind=4) :: istat

  integer :: getcpucount
  external :: getcpucount

  character (len=MPI_MAX_PROCESSOR_NAME) :: nodename
  character (len=20) :: string
  character (len=40) :: numeric
  character (len=128) :: value

  integer :: lenv,statv

  logical ohost

#ifdef CUDA
  type(c_ptr) :: cpfre, cptot
#endif

  string="HOME"
  call get_environment_variable(string,value,lenv,statv)

  !     me          : MPI global rank i
  !     np          : number of MPI ranks
  !     mstr        : id of master MPI rank
  !     numdev      : number of accelerator devices on current node
  !     num_threads : number of OMP threads defined through OMP_NUM_THREAD

  me=0
  np=1

  call mpi_init(ierr)
  !     call mpi_init_thread(MPI_THREAD_SINGLE,iout,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,me,ierr)
  call mpi_comm_size(MPI_COMM_WORLD,np,ierr)

  !     master process is last in the list to enable more effective
  !     allocation of the worker processors, starting at rank 0

  mstr=np-1

  role=idle

  if(me.eq.mstr) role=master

  if(me.eq.mstr) num_threads=1

  allocate(map1(np,7),map2(np,9))

  ! Retrieve the hostname from which to extract the specific computer used
  ! This allows to set some variables to optimal values
  ! Currently the code sets these for:
  !     Juwels Booster (JFZ)
  !     Piz Daint (CSCS)
  !     Snellius (SARA)
  !     Crusher (OLCF)
  !     Frontier (OLCF)
  !     Summit (OLCF)
  
  !     Using function "hostnm()" in stead of
  !     call MPI_get_processor_name(nodename, length, ierr)
  
#ifdef IBM
  call hostnm_(nodename)
#else
  call hostnm(nodename)
#endif
  length=len(trim(nodename))

  machine='            '

  ohost=.false.

#ifdef FRONTIER
  machine='Frontier    '
  ohost=.true.
#endif
#ifdef DISCOVERY
  machine='Discovery   '
  ohost=.true.
#endif
#ifdef LEONARDO
  machine='Leonardo    '
  ohost=.true.
#endif
#ifdef CRUSHER
  machine='Crusher     '
  ohost=.true.
#endif
#ifdef SUMMIT
  machine='Summit      '
  ohost=.true.
#endif
#ifdef LUMI
  machine='Lumi        '
  ohost=.true.
#endif
#ifdef JUWELS
  machine='Juwels      '
  ohost=.true.
#endif
#ifdef PIZDAINT
  machine='Piz Daint   '
  ohost=.true.
#endif
#ifdef LISA
  machine='Lisa        '
  ohost=.true.
#endif
#ifdef SNELLIUS
  machine='Snellius    '
  ohost=.true.
#endif
  
  if(.not.ohost) then
    
    !     Trim the JFZ Juwels Booster node name
    if(nodename(1:3).eq.'jwb') then
      nodename(1:4)=nodename(4:7)
      length=4
      machine='Juwels      '
    endif
    
    !     Trim the CSCS Piz Daint node name
    if(nodename(1:3).eq.'nid') then
      nodename(1:5)=nodename(4:8)
      length=5
      machine='Piz Daint   '
    endif
    
    !     Trim the SurfSara Lisa node name
    if(index(nodename,'.lisa').gt.0) then
      length=index(nodename,'.lisa')-1
      nodename(length+1:length+1)=' '
      machine='Lisa        '
    endif
    
    !     Trim the Sara Snellius node name
    if(nodename(1:3).eq.'gcn') then
      nodename(1:2)=nodename(4:5)
      length=2
      machine='Snellius    '
    endif
    
    !     Trim the Crusher node name
    if(nodename(1:7).eq.'crusher') then
      nodename(1:3)=nodename(8:10)
      length=3
      machine='Crusher     '
    endif
    
    !     Trim the Frontier node name
    if(nodename(1:8).eq.'frontier') then
      nodename(1:5)=nodename(9:13)
      length=5
      machine='Frontier    '
    endif
    
    !     Trim the OLCF Summit node name
    if(ichar(nodename(1:1)).ge.97.and. &
        ichar(nodename(1:1)).le.122.and. &
        ichar(nodename(4:4)).ge.97.and. &
        ichar(nodename(4:4)).le.122.and. &
        ichar(nodename(2:2)).ge.48.and. &
        ichar(nodename(2:2)).le.57.and. &
        ichar(nodename(3:3)).ge.48.and. &
        ichar(nodename(3:3)).le.57.and. &
        ichar(nodename(5:5)).ge.48.and. &
        ichar(nodename(5:5)).le.57.and. &
        ichar(nodename(6:6)).ge.48.and. &
        ichar(nodename(6:6)).le.57) then
      length=6
      machine='Summit      '
    endif
    
  endif

  !     Set the device number

  numdev=-1
  num_threads=1

  memfre=0
  memtot=0

#ifdef _OPENACC
#ifdef GPUAMD
  numdev=acc_get_num_devices(ACC_DEVICE_AMD)
  if(numdev.gt.1) then
    if(machine.ne.'Juwels      ') then
      mydev=mod(me,numdev)
      call acc_set_device_num(mydev,ACC_DEVICE_AMD)
    endif
#else
    numdev=acc_get_num_devices(ACC_DEVICE_NVIDIA)
    if(numdev.gt.1) then
      if(machine.ne.'Juwels      ') then
        mydev=mod(me,numdev)
        call acc_set_device_num(mydev,ACC_DEVICE_NVIDIA)
      endif
#endif
    endif
#endif

#ifdef CUDA
    if(numdev.gt.0) then
      call gronor_update_device_info()
    endif
#endif
    
#ifdef OMPTGT
    numdev=omp_get_num_devices()
#endif

    ! Convert the nodename into a nodenumber stored in variable node

    j=1
    i=1

    do while(i.le.length.and.nodename(i:i).ne.' ')
      if(nodename(i:i).eq.'0'.or.nodename(i:i).eq.'1'.or. &
          nodename(i:i).eq.'2'.or.nodename(i:i).eq.'3'.or. &
          nodename(i:i).eq.'4'.or.nodename(i:i).eq.'5'.or. &
          nodename(i:i).eq.'6'.or.nodename(i:i).eq.'7'.or. &
          nodename(i:i).eq.'8'.or.nodename(i:i).eq.'9') then
        numeric(j:j)=nodename(i:i)
        j=j+1
      elseif(nodename(i:i).eq.'a') then
        numeric(j:j+1)='01'
        j=j+2
      elseif(nodename(i:i).eq.'b') then
        numeric(j:j+1)='02'
        j=j+2
      elseif(nodename(i:i).eq.'c') then
        numeric(j:j+1)='03'
        j=j+2
      elseif(nodename(i:i).eq.'d') then
        numeric(j:j+1)='04'
        j=j+2
      elseif(nodename(i:i).eq.'e') then
        numeric(j:j+1)='05'
        j=j+2
      elseif(nodename(i:i).eq.'f') then
        numeric(j:j+1)='06'
        j=j+2
      elseif(nodename(i:i).eq.'g') then
        numeric(j:j+1)='07'
        j=j+2
      elseif(nodename(i:i).eq.'h') then
        numeric(j:j+1)='08'
        j=j+2
      elseif(nodename(i:i).eq.'i') then
        numeric(j:j+1)='09'
        j=j+2
      elseif(nodename(i:i).eq.'j') then
        numeric(j:j+1)='10'
        j=j+2
      elseif(nodename(i:i).eq.'k') then
        numeric(j:j+1)='11'
        j=j+2
      elseif(nodename(i:i).eq.'l') then
        numeric(j:j+1)='12'
        j=j+2
      elseif(nodename(i:i).eq.'m') then
        numeric(j:j+1)='13'
        j=j+2
      elseif(nodename(i:i).eq.'n') then
        numeric(j:j+1)='14'
        j=j+2
      elseif(nodename(i:i).eq.'o') then
        numeric(j:j+1)='15'
        j=j+2
      elseif(nodename(i:i).eq.'p') then
        numeric(j:j+1)='16'
        j=j+2
      elseif(nodename(i:i).eq.'q') then
        numeric(j:j+1)='17'
        j=j+2
      elseif(nodename(i:i).eq.'r') then
        numeric(j:j+1)='18'
        j=j+2
      elseif(nodename(i:i).eq.'s') then
        numeric(j:j+1)='19'
        j=j+2
      elseif(nodename(i:i).eq.'t') then
        numeric(j:j+1)='20'
        j=j+2
      elseif(nodename(i:i).eq.'u') then
        numeric(j:j+1)='21'
        j=j+2
      elseif(nodename(i:i).eq.'v') then
        numeric(j:j+1)='22'
        j=j+2
      elseif(nodename(i:i).eq.'w') then
        numeric(j:j+1)='23'
        j=j+2
      elseif(nodename(i:i).eq.'x') then
        numeric(j:j+1)='24'
        j=j+2
      elseif(nodename(i:i).eq.'y') then
        numeric(j:j+1)='25'
        j=j+2
      elseif(nodename(i:i).eq.'z') then
        numeric(j:j+1)='26'
        j=j+2
      elseif(nodename(i:i).eq.'A') then
        numeric(j:j+1)='01'
        j=j+2
      elseif(nodename(i:i).eq.'B') then
        numeric(j:j+1)='02'
        j=j+2
      elseif(nodename(i:i).eq.'C') then
        numeric(j:j+1)='03'
        j=j+2
      elseif(nodename(i:i).eq.'D') then
        numeric(j:j+1)='04'
        j=j+2
      elseif(nodename(i:i).eq.'E') then
        numeric(j:j+1)='05'
        j=j+2
      elseif(nodename(i:i).eq.'F') then
        numeric(j:j+1)='06'
        j=j+2
      elseif(nodename(i:i).eq.'G') then
        numeric(j:j+1)='07'
        j=j+2
      elseif(nodename(i:i).eq.'H') then
        numeric(j:j+1)='08'
        j=j+2
      elseif(nodename(i:i).eq.'I') then
        numeric(j:j+1)='09'
        j=j+2
      elseif(nodename(i:i).eq.'J') then
        numeric(j:j+1)='10'
        j=j+2
      elseif(nodename(i:i).eq.'K') then
        numeric(j:j+1)='11'
        j=j+2
      elseif(nodename(i:i).eq.'L') then
        numeric(j:j+1)='12'
        j=j+2
      elseif(nodename(i:i).eq.'M') then
        numeric(j:j+1)='13'
        j=j+2
      elseif(nodename(i:i).eq.'N') then
        numeric(j:j+1)='14'
        j=j+2
      elseif(nodename(i:i).eq.'O') then
        numeric(j:j+1)='15'
        j=j+2
      elseif(nodename(i:i).eq.'P') then
        numeric(j:j+1)='16'
        j=j+2
      elseif(nodename(i:i).eq.'Q') then
        numeric(j:j+1)='17'
        j=j+2
      elseif(nodename(i:i).eq.'R') then
        numeric(j:j+1)='18'
        j=j+2
      elseif(nodename(i:i).eq.'S') then
        numeric(j:j+1)='19'
        j=j+2
      elseif(nodename(i:i).eq.'T') then
        numeric(j:j+1)='20'
        j=j+2
      elseif(nodename(i:i).eq.'U') then
        numeric(j:j+1)='21'
        j=j+2
      elseif(nodename(i:i).eq.'V') then
        numeric(j:j+1)='22'
        j=j+2
      elseif(nodename(i:i).eq.'W') then
        numeric(j:j+1)='23'
        j=j+2
      elseif(nodename(i:i).eq.'X') then
        numeric(j:j+1)='24'
        j=j+2
      elseif(nodename(i:i).eq.'Y') then
        numeric(j:j+1)='25'
        j=j+2
      elseif(nodename(i:i).eq.'Z') then
        numeric(j:j+1)='26'
        j=j+2
      endif
      i=i+1
    enddo

    j=j-1
    if(j.eq.0) then
      node=-1
    else
      string='                    '
      if(j.le.20) then
        string(21-j:20)=numeric(1:j)
      else
        string(1:20)=numeric(1:20)
      endif
      if(j.eq.0) string='                   0'
      read(string,'(i20)') node
    endif

    ! Setup a mapping of the ranks
    
    do j=1,5
      do i=1,np
        map1(i,j)=0
      enddo
    enddo

    !     map1(rank,1) : number of accelerator devices on the node
    !     map1(rank,2) : number of OpenMP threads
    !     map1(rank,3) : group id
    !     map1(rank,4) : node id
    !     map1(rank,5) :

    map1(me+1,1)=numdev
    map1(me+1,2)=num_threads
    map1(me+1,3)=0
    map1(me+1,4)=node
    map1(me+1,5)=0
    map1(me+1,6)=0
    ncount=4*np
    call MPI_AllReduce(map1,map2,ncount,MPI_INTEGER4,MPI_SUM,                 &
        & MPI_COMM_WORLD,ierr)

    deallocate(map1)


    !     node counts the number ranks consecutive node-ids
    !     usually this is the number of nodes in a system
    !     on systems (Summit) with round-robin assignment to multiple
    !     resource sets per node this is the total number of resource sets
    !     nranks counts the number or ranks per resource set

    node=1
    map2(1,5)=1
    nranks=1
    j=1
    do i=2,np
      if(map2(i,4).ne.map2(i-1,4)) node=node+1
      map2(i,5)=node
      if(map2(i,4).eq.map2(i-1,4)) then
        j=j+1
        nranks=max(nranks,j)
      else
        j=1
      endif
    enddo

    !     nnodes counts the number of nodes

    nnodes=1
    map2(1,6)=1
    do i=2,np
      k=0
      do j=1,i-1
        if(map2(i,4).eq.map2(j,4)) then
          k=k+1
          map2(i,6)=map2(j,6)
        endif
      enddo
      if(k.eq.0) then
        nnodes=nnodes+1
        map2(i,6)=nnodes
      endif
    enddo

    !     ncycls is the number of times blocks of consecutive ranks
    !     are assigned to the same node

    !     nrsets is the number of ranks in each block of consecutively
    !     assigned ranks, which on Summit would be the number of resource sets

    !     For example, for 48 ranks on three nodes with 4 resource sets per node,
    !     assigned in round-robin fashion:

    !     Node:         0              1              2
    !     RSet:    0  1  2  3     0  1  2  3     0  1  2  3

    !     Ranks:   0  1  2  3     4  5  6  7     8  9 10 11
    !     Ranks:  12 13 14 15    16 17 18 19    20 21 22 23
    !     Ranks:  24 25 26 27    28 29 30 31    32 33 34 35
    !     Ranks:  36 37 38 39    40 41 42 43    44 45 46 47

    !     np is the total number of ranks: 48
    !     node counts the number of blocks with consecutive ranks: 12
    !     nnodes counts the number nodes with assigned ranks: 3
    !     ncycls is the number of 'Ranks' rows: 12 / 3 = 4
    !     nrsets is the total number of resource sets: 48 / 4 = 12
    !     nranks is the number of ranks in concecutive block: 4
    !     nranks is reset to the total number of ranks per node: 4 * 4 = 16
    !     nrnsets is the number of resource sets per node: 12 / 3 = 4

    !     For example, for 12 ranks on a two node system without resource sets:

    !     Node:           0                      1

    !     Ranks:  0  1  2  3  4  5       6  7  8  9 10 11

    !     np is the total number of ranks: 12
    !     node counts the number of blocks with consecutive ranks: 2
    !     nnodes counts the number nodes with assigned ranks: 2
    !     ncycls is the number of 'Ranks' rows: 2 / 2 = 1
    !     nrsets is the total number of resource sets: 12 / 1 = 12
    !     nranks is the number of ranks in concecutive block: 6
    !     nranks is reset to the total number of ranks per node: 1 * 6 = 6
    !     nrnsets is the number of resource sets per node: 12 / 2 = 6

    !     For example, for 8 ranks on single node system such as a workstation or laptop:

    !     Node:           0

    !     Ranks:  0  1  2  3  4  5  6  7

    !     np is the total number of ranks: 8
    !     node counts the number of blocks with consecutive ranks: 1
    !     nnodes counts the number nodes with assigned ranks: 1
    !     ncycls is the number of 'Ranks' rows: 1 / 1 = 1
    !     nrsets is the total number of resource sets: 8 / 1 = 8
    !     nranks is the number of ranks in concecutive block: 8
    !     nranks is reset to the total number of ranks per node: 1 * 8 = 8
    !     nrnsets is the number of resource sets per node: 8 / 1 = 8

    ncycls=node/nnodes
    nrsets=np/ncycls
    nranks=ncycls*nranks
    nrnsets=nrsets/nnodes

    do i=1,ncycls
      do j=1,nrsets
        map2((i-1)*nrsets+j,3)=i
      enddo
    enddo

    do i=1,np
      map2(i,4)=map2(i,5)
    enddo

    if(ncycls.eq.1) then
      node=map2(1,4)
      j=map2(1,1)
      do i=1,np
        if(map2(i,4).eq.node) then
          if(j.gt.0) then
            map2(i,5)=map2(i,1)
            j=j-1
          else
            map2(i,5)=-map2(i,2)
          endif
        else
          node=map2(i,4)
          j=map2(i,1)
          if(j.gt.0) then
            map2(i,5)=map2(i,1)
            j=j-1
          else
            map2(i,5)=-map2(i,2)
          endif
        endif
      enddo
    else
      do i=1,np
        if(map2(i,1).gt.0.and.map2(i,3).le.map2(i,1))then
          map2(i,5)=map2(i,1)
        else
          map2(i,5)=-map2(i,2)
        endif
      enddo
    endif

    
#ifdef _OPENACC
    if(numdev.ge.1) then
#ifdef GPUAMD
      mydev=acc_get_device_num(ACC_DEVICE_AMD)
#else
      mydev=acc_get_device_num(ACC_DEVICE_NVIDIA)
#endif
      map2(me+1,7)=mydev
    endif
#endif

#ifdef OMPTGT
    if(numdev.ge.1) then
#ifdef IBM
      mydev=omp_get_default_device()
#else
#ifdef GPUAMD
      mydev=omp_get_default_device()
#else
      mydev=omp_get_device_num()
#endif
#endif
      map2(me+1,7)=mydev
    endif
#endif

    do i=1,np-1
      map2(i,8)=worker
    enddo
    map2(np,8)=master
    
    ! Upon exit the array map2 contains the following rank specific values
    
    !      map2(me+1,1) : number of acc devices
    !      map2(me+1,2) : number of omp threads
    !      map2(me+1,3) : group id
    !      map2(me+1,4) : resource set
    !      map2(me+1,5) : acc device or -(omp thread)
    !      map2(me+1,6) : node
    !      map2(me+1,7) : acc device id
    !      map2(me+1,8) : role
    !      map2(me+1,9) :
    
    return
end subroutine gronor_environment
