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
!! Read records from integral file
!!     
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  
      subroutine gronor_read0(lfn,n)

      implicit none
      integer :: lfn,n
      integer :: i
      real (kind=8) :: b(510)

      do i=1,n,510
       read(lfn) b
      enddo

      return
      end

      subroutine gronor_read1(lfn,a,n)

      implicit none
      integer :: lfn,n
      real (kind=8) :: a(n)
      integer :: i,j,mxj
      real (kind=8) :: b(510)

      do i=1,n,510
       read(lfn) b
       mxj=min0(510,(n-i+1))
       do j=1,mxj
        a(i+j-1)=b(j)
       enddo
      enddo

      return
      end
