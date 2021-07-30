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
!! Input routine for the GronOR application
!! @author  Tjerk P. Straatsma, ORNL
!! @author  Coen de Graaf, URV
!! @date    2020

      subroutine gronor_input()

      use inp
      use cidef
      use cidist
      use gnome_parameters
      use gnome_integrals
      
      implicit none

      integer :: i,j,im,ib
      character (len=255) :: item
      integer, allocatable :: itemp(:,:)
      
      open(unit=lfninp,file=filinp,form='formatted',status='old',
     &     err=999)
      
 1    continue

      if(inp_read()) then

 2      continue

        if(.not.inp_a(item)) goto 1

 3      continue

        if(inp_compare(.false.,'Groups',item)) then
          if(.not.inp_i(npg)) call errquit(me,1,"Input error")
          goto 2
        endif

        if(inp_compare(.false.,'Size',item)) then
          if(.not.inp_i(mgr)) call errquit(me,1,"Input error")
          npg=(np-1)/mgr
          if(np.lt.mgr+1) then
            write(*,'(a,i5,a,i5,a)') 'Error: Size specified (',mgr,
     &           ') cannot be supported by available ranks (',np,')'
            call errquit(me,1,"Input error")
          endif
          goto 2
        endif
        
        if(inp_compare(.false.,'Fragments',item)) then
          if(.not.inp_i(nmol)) call errquit(me,1,"Input error")
          goto 2
        endif
        
        if(inp_compare(.false.,'MEBFs',item)) then
          if(nmol.gt.0) then
            if(.not.inp_i(nbase)) call errquit(me,1,"Input error")
            if(nmol.eq.0) call errquit(me,1,
     &           "Fragments should precede MEBFs in input")
            allocate(ncombv(nmol,nbase))
            do i=1,nmol
              if(.not.inp_read()) call errquit(me,1,"Input error")
              do j=1,nbase
                if(.not.inp_i(ncombv(i,j)))
     &               call errquit(me,1,"Input error")
                if(ncombv(i,j).eq.0) ncombv(i,j)=9
                mstates=max(mstates,ncombv(i,j))
              enddo
            enddo
            goto 2
          else
            if(.not.inp_i(nbase)) call errquit(me,1,"Input error")
            allocate(itemp(25,nbase))
            if(.not.inp_read()) call errquit(me,1,"Input error")
 4          continue
            nmol=nmol+1
            do j=1,nbase
              if(.not.inp_i(itemp(nmol,j))) then
                if(j.eq.1) then
                  nmol=nmol-1
                  allocate(ncombv(nmol,nbase))
                  do im=1,nmol
                    do ib=1,nbase
                      ncombv(im,ib)=itemp(im,ib)
                      if(ncombv(im,ib).eq.0) ncombv(im,ib)=9
                      mstates=max(mstates,ncombv(im,ib))
                    enddo
                  enddo
                  deallocate(itemp)
                  goto 2
                endif
                call errquit(me,1,"Input error")
              endif
            enddo
            if(inp_read()) goto 4
            allocate(ncombv(nmol,nbase))
            do im=1,nmol
              do ib=1,nbase
                ncombv(im,ib)=itemp(im,ib)
                if(ncombv(im,ib).eq.0) ncombv(im,ib)=9
                mstates=max(mstates,ncombv(im,ib))
              enddo
            enddo
            deallocate(itemp)
            goto 2            
          endif
        endif
        
        if(inp_compare(.false.,'States',item)) then
          if(.not.inp_i(mstates)) mstates=1
          goto 2
        endif
        
        if(inp_compare(.false.,'Spin',item)) then
          if(.not.inp_i(nspin)) call errquit(me,1,"Input error")
          nspin = nspin - 1
          goto 2
        endif

        if(inp_compare(.false.,'Threshold',item)) then
          if(.not.inp_f(tau_CI)) call errquit(me,1,"Input error")
          goto 2
        endif

        if(inp_compare(.false.,'Labels',item)) then
          if(.not.inp_i(labels)) call errquit(me,1,"Input error")
          goto 2
        endif
        
        if(inp_compare(.false.,'Expert',item)) then
          if(.not.inp_i(ixpert)) ixpert=1
          goto 2
        endif
        
        if(inp_compare(.false.,'Thresh_SIN',item)) then
          if(.not.inp_f(tau_SIN)) tau_SIN=1.0d-12
          goto 2
        endif
        
        if(inp_compare(.false.,'Task',item)) then
          if(.not.inp_i(ntask)) ntask=1
          goto 2
        endif
        
        if(inp_compare(.false.,'Batch',item)) then
          if(.not.inp_i(nbatch)) call errquit(me,1,"Input error Batch")
          nbatcha=nbatch
          if(.not.inp_i(nbatch)) then
          endif
          goto 2
        endif
        
        if(inp_compare(.false.,'Print',item)) then
          ipr=20
          if(inp_i(ipr)) goto 2
          if(.not.inp_a(item)) goto 1
          if(inp_compare(.false.,'Low',item)) then
            ipr=10
            goto 2
          elseif(inp_compare(.false.,'Medium',item)) then
            ipr=20
            goto 2
          elseif(inp_compare(.false.,'High',item)) then
            ipr=30
            goto 2
          elseif(inp_compare(.false.,'Debug',item)) then
            ipr=40
            goto 2
          elseif(inp_compare(.false.,'None',item)) then
            ipr=0
            goto 2
          else
            call errquit(me,1,"Input error in Print")
          endif
          goto 2
        endif
        
        if(inp_compare(.false.,'Timings',item)) then
          if(.not.inp_i(itim)) itim=1
          goto 2
        endif
        
        if(inp_compare(.false.,'Fault',item)) then
          if(.not.inp_i(ifault)) ifault=1
          goto 2
        endif
        
        if(inp_compare(.false.,'Solver',item).or.
     &       inp_compare(.false.,'Solvers',item)) then
          if(.not.inp_i(iaslvr)) then
            iaslvr=-1
            jaslvr=-1
            inslvr=-1
            jnslvr=-1
          elseif(.not.inp_i(jaslvr)) then
            jaslvr=-1
            inslvr=-1
            jnslvr=-1
          elseif(.not.inp_i(inslvr)) then
            inslvr=-1
            jnslvr=-1
          elseif(.not.inp_i(jnslvr)) then
            jnslvr=-1
          endif
          goto 2
        endif
        
        if(inp_compare(.false.,'Tolerance',item)) then
          if(.not.inp_f(tolsvj)) then
            tolsvj=1.0d-07
            tolevj=1.0d-07
          elseif(.not.inp_f(tolevj)) then
            tolevj=1.0d-07
          endif
          goto 2
        endif
        
        if(inp_compare(.false.,'Sweep',item)) then
          if(.not.inp_i(iswsvj)) then
            iswsvj=15
            iswevj=15
          elseif(.not.inp_i(iswevj)) then
            iswevj=15
          endif
          goto 2
        endif
        
        if(inp_compare(.false.,'Debug',item)) then
          if(.not.inp_i(idbg)) idbg=1
          goto 2
        endif
        
        if(inp_compare(.false.,'Dayfile',item)) then
          if(.not.inp_i(iday)) iday=10
          goto 2
        endif
        
        if(inp_compare(.false.,'Progress',item)) then
          if(.not.inp_i(ipro)) ipro=1
          goto 2
        endif
        
        if(inp_compare(.false.,'Test',item)) then
          if(.not.inp_i(itest)) itest=0
          goto 2
        endif
        
        if(inp_compare(.false.,'Distribution',item)) then
          if(.not.inp_i(idist)) idist=0
          goto 2
        endif
        
        if(inp_compare(.false.,'Development',item)) then
          if(.not.inp_i(idevel)) idevel=0
          goto 2
        endif
        
        if(inp_compare(.false.,'Maxcib',item)) then
          if(.not.inp_i(inpcib)) inpcib=0
          goto 2
        endif
        
        if(inp_compare(.false.,'Columns',item)) then
          if(.not.inp_i(ncols)) ncols=7
          goto 2
        endif
        
        if(inp_compare(.false.,'MPIbuffer',item)) then
          if(.not.inp_i(mpibuf)) call errquit(me,1,"Input error")
          goto 2
        endif
        
        if(inp_compare(.false.,'Accelerate',item)) then
          if(.not.inp_i(naccel)) naccel=0
          goto 2
        endif
        
        if(inp_compare(.false.,'MPS',item)) then
          if(.not.inp_i(nummps)) nummps=1
          goto 2
        endif
        
        if(inp_compare(.false.,'OMP',item)) then
          if(.not.inp_i(numomp)) numomp=1
          goto 2
        endif
        
        if(inp_compare(.false.,'cpr',item)) then
          lcpr=.true.
          goto 2
        endif

        if (inp_compare(.false.,'Ecorr',item)) then
          ncorr = 1
          allocate(ecorr(mstates))
          do i=1,mstates
            if (.not.inp_f(ecorr(i))) ecorr(i) = 0.0d0
          end do
          goto 2
        end if
        
*     Use Gallup-Norbeck weights? Yes, nwt=1; No, nwt .ne. 1
        if (inp_compare(.false.,'GNWeights',item)) then
          if(.not.inp_i(nwt)) nwt = 0
          goto 2
        endif
        
        goto 1
      endif

      close(unit=lfninp)
      
      return
      
 999  write(lfnout,989) filinp(1:inp_strlen(filinp))
 989  format('Unable to open input file ',a)
      
      call errquit(me,1,"Input error")
      
      end subroutine gronor_input
