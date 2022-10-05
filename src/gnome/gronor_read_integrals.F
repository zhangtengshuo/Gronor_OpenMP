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
      integer :: ilast,ii,jj,myigr,lc,nc,mc,numc,nnc
      integer(kind=2)::ibuf(1020)

      allocate(ig(mgr))
      allocate(ntarget(mgr))

      mynode=map2(me+1,6)
      
      if(me.eq.master) then
        if(ipr.ge.20) then
          numc=16
          lc=numgrp/numc
          mc=mod(numgrp,numc)
          if(mc.gt.0) lc=lc+1
          nc=numc
          write(lfnrnk,110)
 110      format(//,' Allgroups G=Group A=Accel R=Rank')
          do i=1,lc
            if(mc.gt.0.and.i.eq.lc) nc=mc
            write(lfnrnk,111) ((i-1)*numc+j,j=1,nc)
 111        format(/,'       G:',20i8)
            write(lfnrnk,112) (allgroups((i-1)*numc+j,1),j=1,nc)
 112        format('       A:',20i8)
            do k=2,mgr+1
              write(lfnrnk,113) (allgroups((i-1)*numc+j,k),j=1,nc)
 113          format('       R:',20i8)
            enddo
          enddo
        endif
        if(ipr.ge.40) then
          lc=np/3
          mc=mod(np,3)
          if(mc.gt.0) lc=lc+1
          nc=3
          write(lfnrnk,114)
 114      format(//,' Rank map2 ND=NumDev DI=DevId Nt=NumThr',
     &         ' Ac=Accel',//,
     &         '    Rank  ND NT   Group    RSet Ac    Node DI  ',
     &         '    Rank  ND NT   Group    RSet Ac    Node DI  ',
     &         '    Rank  ND NT   Group    RSet Ac    Node DI',/)
          do i=1,lc
            if(mc.gt.0.and.i.eq.lc) nc=mc
            nnc=nc
            if(i+(nc-1)*lc.gt.np) nnc=nnc-1
            write(lfnrnk,115) (i+(j-1)*lc-1,
     &           (map2(i+(j-1)*lc,k),k=1,7),j=1,nnc)
 115        format(3(i8,':',2i3,2i8,i3,i8,i3,2x))
          enddo
        endif
      endif

      myigr=0
      do i=1,numgrp
        do j=1,mgr
          if(allgroups(i,j+1).eq.me) myigr=j
        enddo
      enddo
      
      call MPI_Comm_Group(MPI_COMM_WORLD,group_world,ierr)

      allocate(new_comm(mgr),new_comm1(mgr),new_comm2(mgr,nnodes))

      allocate(ranks_list(numgrp))

!     Create communicators for each first, second, etc rank from each group

      if(idist.eq.0) then
        do i=1,mgr
          do j=1,numgrp
            ranks_list(j)=allgroups(j,i+1)
          enddo
          call MPI_Group_Incl(group_world,int(numgrp,kind=4),
     &         ranks_list,new_group,ierr)
          call MPI_Comm_Create(MPI_COMM_WORLD,new_group,
     &         new_comm(i),ierr)
        enddo
      endif

!     Create communicators for each first, second, etc rank from each first group per node

      if(idist.ne.0) then
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
     &         ranks_list,new_group,ierr)
          call MPI_Comm_Create(MPI_COMM_WORLD,new_group,
     &         new_comm1(i),ierr)
        enddo
      endif

!     Create communicators for each first, second, etc rank from groups on the same node
      
      if(idist.ne.0) then

        if(idist.eq.1) then
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
     &             ranks_list,new_group,ierr)
              call MPI_Comm_Create(MPI_COMM_WORLD,new_group,
     &             new_comm2(i,n),ierr)
            enddo
          enddo
          
        else
          numgrn=0
          if(mygroup.gt.0) then
            do i=1,numgrp
              if(map2(allgroups(i,myigr+1)+1,6).eq.mynode)
     &             numgrn=numgrn+1
            enddo
          endif
          allocate(igrn(numgrn))
          numgrn=0
          do i=1,numgrp
            if(map2(allgroups(i,myigr+1)+1,6).eq.mynode) then
              numgrn=numgrn+1
              igrn(numgrn)=allgroups(i,myigr+1)
            endif
          enddo
        endif

      endif
      deallocate(ranks_list)

      call gronor_parallel_integral_input()

      return
 993  write(lfnout,983) trim(filint)
      call gronor_abort(250,trim(filint))
 994  write(lfnout,984) trim(filint)
      call gronor_abort(251,trim(filint))
 983  format('Error reading integral file ',a)
 984  format('Unable to open integral file ',a)
      end
