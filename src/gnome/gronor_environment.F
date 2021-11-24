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
      use cuda_functions
c      use cuda_rt
#endif
#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer :: i,j,k,node
      integer (kind=4) :: len, ierr, istat, ncount

      integer :: getcpucount
      external :: getcpucount

      character (len=MPI_MAX_PROCESSOR_NAME) :: nodename
      character (len=20) :: string
      character (len=40) :: numeric
      character (len=128) :: value

      integer (kind=8) :: lenv,statv

#ifdef CUDA
      type(c_ptr) :: cpfre, cptot
#endif

      string="HOME"
      call get_environment_variable(string,value,lenv,statv)

!     me          : MPI global rank id
!     np          : number of MPI ranks
!     master      : id of master MPI rank
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

      master=np-1
      numdev=0
      num_threads=1

      memfre=0
      memtot=0
      
#ifdef _OPENACC
      numdev=acc_get_num_devices(ACC_DEVICE_NVIDIA)
      if(numdev.gt.1) then
c     mydev=mod(me,numdev)
c     call acc_set_device_num(mydev,ACC_DEVICE_NVIDIA)
#ifdef CUDA
        cpfre=c_loc(memfre)
        cptot=c_loc(memtot)
c        istat=cudaMemGetInfo(cpfre,cptot)
#endif
        memavail=memfre
      endif
#endif
#ifdef _OL_OMP_
      numdev=omp_get_num_devices()
#endif
      
      num_threads=1

      if(me.eq.master) num_threads=0

      allocate(map1(np,7),map2(np,7))

!     call MPI_get_processor_name(nodename, len, ierr)
      call hostnm(nodename)
      len=inp_strlen(nodename)

      machine='      '
      
!     Trim the JFZ Juwels Booster node name
      if(nodename(1:3).eq.'jwb') then
        nodename(1:4)=nodename(4:7)
        len=4
        machine='Juwels'
      endif

!     Trim the CSCS Piz Daint node name
      if(nodename(1:3).eq.'nid') then
        nodename(1:5)=nodename(4:8)
        len=5
        machine='Daint '
      endif
      
      j=1
      i=1

      do while(i.le.len.and.nodename(i:i).ne.' ')
        if(nodename(i:i).eq.'0'.or.nodename(i:i).eq.'1'.or.
     &       nodename(i:i).eq.'2'.or.nodename(i:i).eq.'3'.or.
     &       nodename(i:i).eq.'4'.or.nodename(i:i).eq.'5'.or.
     &       nodename(i:i).eq.'6'.or.nodename(i:i).eq.'7'.or.
     &       nodename(i:i).eq.'8'.or.nodename(i:i).eq.'9') then
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
        read(string,'(i20)') node
      endif
      
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
      call MPI_AllReduce(map1,map2,ncount,MPI_INTEGER8,MPI_SUM,
     & MPI_COMM_WORLD,ierr)

      deallocate(map1)

!      map2(me+1,1) : number of acc devices
!      map2(me+1,2) : number of omp threads
!      map2(me+1,3) : group
!      map2(me+1,4) : resource set
!      map2(me+1,5) : 0
!      map2(me+1,6) : node
!      map2(me+1,7) : acc device id
      
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
        mydev=acc_get_device_num(ACC_DEVICE_NVIDIA)
        map2(me+1,7)=mydev
      endif
#endif
      
      return
      end
