!     This file is part of the GronOR software

!     GronOR is free software, and can be used, re-distributed and/or modified under
!     the Apache License version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
!     Any use of the software has to be in compliance with this license. Unless required
!     by applicable law or agreed to in writing, software distributed under the license
!     is distributed on an ‘as is’ basis, without warranties or conditions of any kind,
!     either express or implied.
!     See the license for the specific language governing permissions and limitations
!     under the license.

!     GronOR is copyright of the University of Groningen

!> @brief
!! Manager rank assignments
!! @author  Tjerk P. Straatsma, ORNL
!! @author  Coen de Graaf, URV
!! @date    2016
!>

subroutine gronor_assign_managers()

  use cidist

  integer :: i,j,m,n
  
  do i=1,np
    map2(i,8)=0
    map2(i,9)=-1
  enddo

  ! The last rank is already assigned to be the overall master rank
  map2(mstr+1,8)=master

  ! If the number of nodes is one, this is most likely a workstation
  if(nnodes.eq.1) then
    map2(mstr,8)=manager
    if(managers.gt.1) then
      m=mod(np-1,managers)
      n=(np-1-m)/managers
      do i=1,managers-1
        map2(i*n,8)=manager
      enddo
    endif
  endif

  ! If the number of nodes is larger than one, assign one manager per node
  if(nnodes.gt.1) then
  endif
  

  ! Set the worker ranks and assign the manager rank
  m=0
  do i=1,np
    if(m.eq.0.and.map2(i,8).ne.manager) m=i
    if(m.ne.0.and.map2(i,8).eq.manager) then
      do j=m,i-1
        map2(j,8)=worker
        map2(j,9)=i-1
      enddo
      m=0
    endif
  enddo

  ! Set all idle ranks

  do i=1,np
    if(map2(i-1,8).eq.0) map2(i-1,8)=idle
  enddo


  ! Set role of current rank

  role=map2(me+1,8)
  
  
  if(me.eq.mstr) then
    write(*,'(i4,a,i10)') me," Number of nodes is ",nnodes
    write(*,'(i4,a,i10)') me," Number of managers is ",managers
  endif

  write(*,'(16i8)') me,role,(map2(me+1,j),j=1,9)

  
end subroutine gronor_assign_managers
