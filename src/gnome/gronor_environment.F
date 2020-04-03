      subroutine gronor_environment()
 
      use mpi
      use cidist

#ifdef _OPENACC
      use openacc
      use cuda_functions
#endif
#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer :: i,j,node
      integer (kind=4) :: len, ierr
      integer (kind=4) :: istat

      integer :: getcpucount
      external :: getcpucount

      character (len=MPI_MAX_PROCESSOR_NAME) :: nodename
      character (len=20) :: string
      integer::  cudaGetDeviceCount,jstat,iout

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
c      if(numdev.gt.1) then
c       mydev=mod(me,numdev)
c       call acc_set_device_num(mydev,ACC_DEVICE_NVIDIA)
c       cpfre=c_loc(memfre)
c       cptot=c_loc(memtot)
c       istat=cudaMemGetInfo(cpfre,cptot)
c       memavail=memfre
c      endif
#endif

      num_threads=1

      if(me.eq.master) num_threads=0

      allocate(map1(np,5),map2(np,5))

      call MPI_get_processor_name(nodename, len, ierr)

      j=len
      i=1
      do while(i.le.len.and.j.gt.0)
       if(nodename(i:i).ne.'0'.and.nodename(i:i).ne.'1'.and.
     &    nodename(i:i).ne.'2'.and.nodename(i:i).ne.'3'.and.
     &    nodename(i:i).ne.'4'.and.nodename(i:i).ne.'5'.and.
     &    nodename(i:i).ne.'6'.and.nodename(i:i).ne.'7'.and.
     &    nodename(i:i).ne.'8'.and.nodename(i:i).ne.'9') then
         j=j-1
         nodename(i:j)=nodename(i+1:j+1)
       else
        i=i+1
       endif
      enddo
      i=i-1
      if(i.eq.0) then
       node=-1
      else
       string='                    '
       if(i.le.20) then
        string(21-i:20)=nodename(1:i)
       else
        string(1:20)=nodename(1:20)
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

      call MPI_AllReduce(map1,map2,5*np,MPI_INTEGER8,MPI_SUM,
     & MPI_COMM_WORLD,ierr)

      deallocate(map1)

      node=1
      map2(1,5)=1
      do i=2,np
       if(map2(i,4).ne.map2(i-1,4)) node=node+1
       map2(i,5)=node
      enddo

      do i=1,np
       map2(i,4)=map2(i,5)
      enddo

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

      return
      end
