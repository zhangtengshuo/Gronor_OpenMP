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
!! Read integral file
!!     
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  
      subroutine gronor_read_integrals()

      use mpi
      use inp
      use cidef
      use cidist
      use gnome_parameters
      use gnome_data
      use gnome_integrals

#ifdef _OPENMP
      use omp_lib
#endif

      implicit none

      integer (kind=4) :: ierr
      integer :: i,ig0,j,kl,k,n,nn,igr,ngi,mdum,ielem,jkoff
      integer, allocatable :: idum(:)
      real (kind=8) :: b(510)
      integer :: status(MPI_STATUS_SIZE)

      integer :: kint, kintbatch, mint,nint,jint,lint
      integer :: nacc
      integer :: ilast,ii,jj
      integer(kind=2)::ibuf(1020)

      allocate(ig(mgr))
      allocate(ntarget(mgr))

      call daytag(10)
      
      call MPI_Comm_Group(MPI_COMM_WORLD,group_world,ierr)

      call daytag(11)
      
      allocate(new_comm(mgr),new_comm1(mgr),new_comm2(mgr,nnodes))

      allocate(ranks_list(numgrp))

!     Create communicators for each first, second, etc rank from each group
      
      do i=1,mgr
       do j=1,numgrp
         ranks_list(j)=allgroups(j,i+1)
       enddo
       call MPI_Group_Incl(group_world,int(numgrp,kind=4),
     & ranks_list,new_group,ierr)
       call MPI_Comm_Create(MPI_COMM_WORLD,new_group,
     &                      new_comm(i),ierr)
      enddo

      call daytag(12)
!     Create communicators for each first, second, etc rank from each first group per node

      lcomm1=.false.
      do i=1,mgr
        n=0
        do j=1,numgrp
          if(map2(allgroups(j,i+1)+1,6).eq.n+1) then
            n=n+1
            ranks_list(n)=allgroups(j,i+1)
            if(me.eq.allgroups(j,i+1)) lcomm1=.true.
          endif
        enddo
        call MPI_Group_Incl(group_world,int(n,kind=4),
     &       ranks_list,new_group,ierr)
        call MPI_Comm_Create(MPI_COMM_WORLD,new_group,
     &       new_comm1(i),ierr)
      enddo

      call daytag(13)
!     Create communicators for each first, second, etc rank from groups on the same node
      
      do i=1,mgr
        do n=1,nnodes
          k=0
          do j=1,numgrp
            if(map2(allgroups(j,i+1)+1,6).eq.n) then
              k=k+1
              ranks_list(k)=allgroups(j,i+1)
            endif
          enddo
          call MPI_Group_Incl(group_world,int(k,kind=4),
     &         ranks_list,new_group,ierr)
          call MPI_Comm_Create(MPI_COMM_WORLD,new_group,
     &         new_comm2(i,n),ierr)
        enddo
      enddo
      deallocate(ranks_list)

      call daytag(14)
      allocate(ranks_list(mgr))
      allocate(grp_comm(numgrp))

      do i=1,numgrp
       do j=1,mgr
        ranks_list(j)=allgroups(i,j+1)
       enddo

       call MPI_Group_Incl(group_world,int(mgr,kind=4),
     &  ranks_list,grp_group,ierr)
       call MPI_Comm_Create(MPI_COMM_WORLD,grp_group,
     &                      grp_comm(i),ierr)
      enddo

      deallocate(ranks_list)
      call daytag(15)

      call gronor_parallel_integral_input()

      call daytag(16)
      return
 993  write(lfnout,983) filint(1:inp_strlen(filint))
      call errquit(me,5,"Error in reading integrals")
 994  write(lfnout,984) filint(1:inp_strlen(filint))
      call errquit(me,5,"Error in reading integrals")
 983  format('Error reading integral file ',a)
 984  format('Unable to open integral file ',a)
      end
