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
!! GNOME input of sym file in Gamess format
!!     
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  
      subroutine gronor_gnome_symmetry_input_gamess
     &     (lfnsym,lfnout,ipr)
      
      use cidist
      use gnome_data
      
      implicit none

      integer :: lfnsym, lfnout, ipr

      integer :: i,ii,j
! rkk basis
      integer :: nsh(750)
      integer,parameter :: maxorb=4096,maxat=750,mxshel=2048
      
      integer :: nat,nshell,i3,nsf,npf,nspf,ndf,nff,ngf
      integer :: ks(maxat),iat,ns1,ns2,i1,i2,
     &     ntype,intyp(maxorb),ityp,ish,jshell(6),
     &     kstart(mxshel),kng(mxshel),kmax(mxshel),kmin(mxshel)
      
      integer :: nfac
      
      character*8 ylabel(6),ytemp(mxshel)
      integer :: nynum(5)
      
      data ntype/6/
      data ylabel/'s','p','d','f','g','sp'/
      data jshell/ 1 , 3 , 6 ,10 ,15,  4 /
! ---

      write(lfnout,601)
 601  format(/,' Molecular coordinates and symmetry')

      read(lfnsym,100,end=999) text
      read(lfnsym,100,end=999) name(1)
      read(lfnsym,100,end=999) name(2)
 100  format(a80)

      allocate(centn(nnucl))
      allocate(xcord(nnucl),ycord(nnucl),zcord(nnucl),znuc(nnucl))

      write(lfnout,602)

 602  format(/,' Atom      x-coordinate   y-coordinate',
     & '   z-coordinate   znuc   basis',/)
! rkk basis
      do i=1,nnucl
       read(lfnsym,120,end=999) centn(i),xcord(i),ycord(i),zcord(i),
     &  znuc(i),nsh(i)
 120   format(2x,a4,3f20.14,f8.4,1i4)
      enddo

        nat=nnucl
        read(lfnsym,'(1i4)') nshell
        read(lfnsym,'(18i4)') (ks(i), i=1,nat)
        read(lfnsym,'(18i4)') (kstart(i), i=1,nshell)
        read(lfnsym,'(18i4)') (kng(i), i=1,nshell)
        read(lfnsym,'(18i4)') (kmax(i), i=1,nshell)
        read(lfnsym,'(18i4)') (kmin(i), i=1,nshell)

        do i = 1 , nshell
         ii = kmax(i) - kmin(i) + 1
         do j = 1 , ntype
          if (ii.eq.jshell(j)) then
           intyp(i) = j
          end if
         end do
        end do
        do iat=1,nat
        ns2=nsh(iat)
        if(ns2.eq.0)go to 1180
        i=(ns2-1)/3+1
        i3=3*i
        do 7401 j=1,i3
 7401   ytemp(j)=' '
        i=0
        ns1 = ks(iat)
        ns2 = ns1+ns2-1
        do 7402 ish=ns1,ns2
        i=i+1
 7402   ytemp(i)=ylabel(intyp(ish))
        nsf=0
        npf=0
        nspf=0
        ndf=0
        nff=0
        ngf=0

        nynum=0
        do ish=ns1,ns2
         if((intyp(ish).eq.1).or.(intyp(ish).eq.6)) nynum(1)=nynum(1)+1
         if((intyp(ish).eq.2).or.(intyp(ish).eq.6)) nynum(2)=nynum(2)+1
         if(intyp(ish).eq.3)  nynum(3)=nynum(3)+1 
         if(intyp(ish).eq.4)  nynum(4)=nynum(4)+1 
         if(intyp(ish).eq.5)  nynum(5)=nynum(5)+1 
        end do

        
        write(lfnout,603,advance='no') centn(iat),xcord(iat),ycord(iat),
     &  zcord(iat), znuc(iat)
  603   format(2x,a4,3f15.6,3x,f6.2,2x)
        
        do i=1,5
         if(nynum(i).gt.0) then 
          write(lfnout,604,advance='no') nynum(i),ylabel(i)
         end if
        end do
  604   format(5(i3,a1,x))
        write(lfnout,'(a)') ' '

 1180   continue

        enddo

        if(ipr.ge.6) then
        write(lfnout,6050)

        i1 = 1
        do iat = 1 , nat
         ns2 = nsh(iat)
         if (ns2.ne.0) then
          write (lfnout,6040) centn(iat)
          ns1 = ks(iat)
          ns2 = ns1 + ns2 - 1
          do ish = ns1 , ns2
           write (lfnout,6020)
           i2 = i1 + kng(ish) - 1
           ityp = intyp(ish)
           nfac=1
           i1 = i1 + nfac*kng(ish)
          end do
         end if
        end do
        write (lfnout,6060)
 6010   format (15x,i3,3x,a4,3x,i3,1x,2f15.6)
 6020   format (/)
 6030   format (15x,i3,3x,a4,3x,i3,1x,3f15.6)
 6040   format (/1x,a8)
 6050   format (/1x,'atom',8x,'shell',3x,'type',2x,'prim',7x,
     &        'exponents',7x,'contraction coefficients'/,87('*'))
 6060   format (/,87('*'))
        end if

      read(lfnsym,100,end=999) text

      if(me.eq.master.and.ipr.gt.2) write(lfnout,611)
 611  format(/,' Reading symmetry for INS2')

      return

 999  write(lfnout,9999)
 9999 format('Error reading symmetry file')
      call errquit(me,9,"Error reading symmetry file")
      end

