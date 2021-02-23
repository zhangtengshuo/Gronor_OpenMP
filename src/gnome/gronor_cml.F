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
!! GNOME input of sym file
!!     
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  

      subroutine gronor_cml()
      use gnome_data
      use cidef
      use inp
      use cidist

      implicit none
      character(len=2), external :: element
      real (kind=8), external :: atom_radius

      integer :: i,j,numi,numj,icnt(105)
      real (kind=8) :: c,d

      c=0.529177249d0
      
      do i=1,105
        icnt(i)=0
      enddo
      
      open(unit=lfncml,file=filcml(1:inp_strlen(filcml)),
     &     form='formatted',status='unknown',err=996)

      write(lfncml,2401) trim(root)
 2401 format('<molecule title="',a,'">')
      write(lfncml,2402)
 2402 format("   <atomarray>")
      do i=1,nnucl
        numi=int(znuc(i))
        icnt(numi)=icnt(numi)+1
        write(lfncml,2403) trim(centn(i)),
     &       trim(adjustl(element(numi))),
     &       c*xcord(i),c*ycord(i),c*zcord(i)
 2403   format('      <atom id="',a,'" elementType="',a,
     &       '" x3="',f14.8,'" y3="',f14.8,'" z3="',f14.8,'"/>')
      enddo
      write(lfncml,2404)
 2404 format("   </atomarray>")
      write(lfncml,2405)
 2405 format("   <bondarray>")
      do i=1,nnucl-1
        numi=int(znuc(i))
        do j=i+1,nnucl
          numj=int(znuc(i))
          d=c*sqrt((xcord(i)-xcord(j))**2+(ycord(i)-ycord(j))**2+
     &         (zcord(i)-zcord(j))**2)
          if(1.0d-1*d.lt.atom_radius(numi)+atom_radius(numj)) then
            write(lfncml,2406) trim(centn(i)),trim(centn(j))
 2406       format('      <bond atomrefs2="',a,1x,a'" order="S"/>')
          endif
        enddo
      enddo
      write(lfncml,2407)
 2407 format("   </bondarray>")

      write(lfncml,2499)
 2499 format("</molecule>")
      
      close(unit=lfncml,status='keep')
      
      return
      
 996  write(lfnout,986) filcml(1:inp_strlen(filcml))
 986  format('Unable to open file ',a)
      
      call errquit(me,0,"Unable to open cml file")
      
      end subroutine gronor_cml
