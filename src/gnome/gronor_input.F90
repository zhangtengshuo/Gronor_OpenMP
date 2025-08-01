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
  use gnome_solvers

  implicit none

  external :: gronor_abort

  integer :: i,j,m
  character (len=255) :: item,item2,fil
  real (kind=8) value
  logical :: last,exist,lcouple

  lcouple=.false.

  open(unit=lfninp,file=filinp,form='formatted',status='old',err=999)

1 continue

  if(inp_read()) then

2   continue

    if(.not.inp_a(item)) goto 1

3   continue

    if(inp_compare(.false.,'Groups',item)) then
      if(.not.inp_i(npg)) call gronor_abort(101,"Groups")
      goto 2
    endif

    if(inp_compare(.false.,'Size',item)) then
      if(.not.inp_i(mgr)) call gronor_abort(102,"Size")
      npg=(np-1)/mgr
      if(np.lt.mgr+1) then
        write(*,'(a,i5,a,i5,a)') 'Error: Size specified (',mgr, &
            ') cannot be supported by available ranks (',np,')'
        call gronor_abort(103,"Input error")
      endif
      goto 2
    endif

    if(inp_compare(.false.,'Fragments',item)) then
      if(.not.inp_i(nmol)) call gronor_abort(104,"Input error")
      goto 2
    endif

    if(inp_compare(.false.,'Label',item)) then
      if(.not.inp_i(labmax)) call gronor_abort(121,"Input error")
      goto 2
    endif

    if(inp_compare(.false.,'MEBFs',item)) then
      if(.not.inp_a(mebfroot)) call gronor_abort(105,"Input error")
      if(.not.inp_i(nbase)) call gronor_abort(106,"Input error")
      if(.not.inp_a(combas)) combas=""
      allocate(fragname(36))
      allocate(fragstate(36,nbase))
      nmol=0
      do while(.true.)
        if(.not.inp_read()) call gronor_abort(107,"Input error")
        last=.not.inp_a(item)
        if(last.or.len(trim(item)).gt.3) then
          allocate(ncombv(nmol,nbase))
          mstates=nmol*nbase
          allocate(fragfile(mstates))
          allocate(vecfile(mstates))
          allocate(detfile(mstates))
          if(len(trim(combas)).eq.0) then
            do i=1,nmol
              write(combas,'(a,a)') trim(combas),trim(fragname(i))
            enddo
            write(fil,'(a,a,a)') trim(mebfroot),trim(combas),".sys"
            inquire(file=trim(fil),exist=EXIST)
            if(.not.exist) combas=""
          endif
          do i=1,nmol
            do j=1,nbase
              ncombv(i,j)=(i-1)*nbase+j
              write(fragfile((i-1)*nbase+j),'(a,a,a,a)') &
                  trim(mebfroot),trim(fragname(i)),"_",trim(fragstate(i,j))
              if(len(trim(combas)).gt.0) then
                write(vecfile((i-1)*nbase+j),'(a,a,a,a,a,a)') &
                    trim(mebfroot),trim(combas),trim(fragname(i)),"_",trim(fragstate(i,j)),'.vec'
                write(detfile((i-1)*nbase+j),'(a,a,a,a,a,a)') &
                    trim(mebfroot),trim(combas),trim(fragname(i)),"_",trim(fragstate(i,j)),'.det'
              else
                vecfile((i-1)*nbase+j)=fragfile((i-1)*nbase+j)//'.vec'
                detfile((i-1)*nbase+j)=fragfile((i-1)*nbase+j)//'.det'
              endif
              write(fil,'(a,a,a)') trim(fragfile((i-1)*nbase+j)),".vec"
              inquire(file=trim(fil),exist=EXIST)
              if(exist) vecfile((i-1)*nbase+j)=fil
              write(fil,'(a,a,a)') trim(fragfile((i-1)*nbase+j)),".det"
              inquire(file=trim(fil),exist=EXIST)
              if(exist) detfile((i-1)*nbase+j)=fil
            enddo
          enddo
          if(last) goto 1
          goto 3
        endif
        nmol=nmol+1
        if(nmol.gt.36) call gronor_abort(109,"Dimension error")
        fragname(nmol)=trim(item)
        do i=1,nbase
          if(.not.inp_a(item)) call gronor_abort(108,"Input error")
          fragstate(nmol,i)=trim(item)
        enddo
      enddo
    endif

    if(inp_compare(.false.,'States',item)) then
      if(.not.inp_i(mstates)) mstates=1
      goto 2
    endif

    if(inp_compare(.false.,'Spin',item)) then
      if(.not.inp_i(nspin)) call gronor_abort(110,"Input error")
      nspin = nspin - 1
      goto 2
    endif

    if(inp_compare(.false.,'Managers',item)) then
      if(.not.inp_i(managers)) managers=1
      goto 2
    endif

    if(inp_compare(.false.,'Couplings',item)) then
      lcouple=.true.
      allocate(inter_couplings(nmol-1,nbase))
      do i=1,nmol-2
        if(.not.inp_read()) call gronor_abort(111,"Input error")
        do j=1,nbase
          if(.not.inp_i(inter_couplings(i,j))) call gronor_abort(111,"Input error")
        enddo
        do j=1,nbase
          inter_couplings(i,j)=inter_couplings(i,j)-1
        enddo
      enddo
      do j=1,nbase
        inter_couplings(nmol-1,j)=nspin
      enddo
      goto 2
    endif

    if(inp_compare(.false.,'Threshold',item)) then
      if(.not.inp_f(tau_CI)) call gronor_abort(112,"Input error")
      if(.not.inp_f(tau_CI_off)) tau_CI_off=tau_CI
      goto 2
    endif

    if(inp_compare(.false.,'Expert',item)) then
      if(.not.inp_i(ixpert)) ixpert=1
      goto 2
    endif

    if(inp_compare(.false.,'Abort',item)) then
      if(.not.inp_i(nabort)) call gronor_abort(115,"Input error")
      goto 2
    endif

    if(inp_compare(.false.,'Thresh_SIN',item)) then
      if(.not.inp_f(tau_SIN)) call gronor_abort(116,"Input error")
      goto 2
    endif

    if(inp_compare(.false.,'Task',item)) then
      ntaska=1
      if(.not.inp_i(ntaska)) call gronor_abort(117,"Input error")
      if(.not.inp_i(ntask)) ntask=ntaska
      goto 2
    endif

    if(inp_compare(.false.,'Load',item)) then
      loada=2
      if(.not.inp_i(loada)) call gronor_abort(118,"Input error")
      if(.not.inp_i(load)) load=loada
      goto 2
    endif

    if(inp_compare(.false.,'Batch',item)) then
      if(.not.inp_i(nbatch)) call gronor_abort(119,"Input error Batch")
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
        call gronor_abort(120,"Input error in Print")
      endif
      goto 2
    endif

    if(inp_compare(.false.,'Timings',item)) then
      if(.not.inp_i(itim)) itim=12
      goto 2
    endif

    if(inp_compare(.false.,'Fault',item)) then
      if(.not.inp_i(ifault)) ifault=1
      goto 2
    endif

    if(inp_compare(.false.,'Solver',item).or.inp_compare(.false.,'Solvers',item)) then
4     continue
        if(.not.inp_a(item)) goto 1
        5 continue
      if(inp_compare(.false.,'SVD',item)) then
        if(.not.inp_a(item))  call gronor_abort(123,"Input error Solver")
        if(inp_compare(.false.,item,'EISPACK')) then
          iaslvr=SOLVER_EISPACK
        elseif(inp_compare(.false.,item,'MKL')) then
          iaslvr=SOLVER_MKL
        elseif(inp_compare(.false.,item,'MKLD')) then
          iaslvr=SOLVER_MKLD
        elseif(inp_compare(.false.,item,'MKLJ')) then
          iaslvr=SOLVER_MKLJ
        else
          call gronor_abort(123,"Input error Solver")
        endif
        if(.not.inp_a(item)) then
          inslvr=iaslvr
          goto 1
        endif
        if(inp_compare(.false.,item,'EISPACK')) then
          inslvr=SOLVER_EISPACK
        elseif(inp_compare(.false.,item,'MKL')) then
          inslvr=SOLVER_MKL
        elseif(inp_compare(.false.,item,'MKLD')) then
          inslvr=SOLVER_MKLD
        elseif(inp_compare(.false.,item,'MKLJ')) then
          inslvr=SOLVER_MKLJ
        else
          inslvr=iaslvr
          goto 5
        endif
        goto 4
      elseif(inp_compare(.false.,'EV',item).or.inp_compare(.false.,'EVD',item) &
        .or.inp_compare(.false.,'SYEV',item).or.inp_compare(.false.,'SYEVD',item)) then
        if(.not.inp_a(item))  call gronor_abort(123,"Input error Solver")
        if(inp_compare(.false.,item,'EISPACK')) then
          jaslvr=SOLVER_EISPACK
        elseif(inp_compare(.false.,item,'MKL')) then
          jaslvr=SOLVER_MKL
        elseif(inp_compare(.false.,item,'MKLD')) then
          jaslvr=SOLVER_MKLD
        elseif(inp_compare(.false.,item,'MKLJ')) then
          jaslvr=SOLVER_MKLJ
        else
          call gronor_abort(123,"Input error Solver")
        endif
        if(.not.inp_a(item)) then
          jnslvr=jaslvr
          goto 1
        endif
        if(inp_compare(.false.,item,'EISPACK')) then
          jnslvr=SOLVER_EISPACK
        elseif(inp_compare(.false.,item,'MKL')) then
          jnslvr=SOLVER_MKL
        elseif(inp_compare(.false.,item,'MKLD')) then
          jnslvr=SOLVER_MKLD
        elseif(inp_compare(.false.,item,'MKLJ')) then
          jnslvr=SOLVER_MKLJ
        else
          jnslvr=jaslvr
          goto 5
        endif
        goto 4
      else
        call gronor_abort(123,"Input error Solver")
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

    if(inp_compare(.false.,'Resilient',item)) then
      if(.not.inp_i(ires)) ires=1
      goto 2
    endif

    if(inp_compare(.false.,'Interrupt',item)) then
      if(.not.inp_i(iint)) iint=1
      goto 2
    endif

    if(inp_compare(.false.,'Temp',item)) then
      if(.not.inp_i(itmp)) itmp=1
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
      if(.not.inp_i(mpibuf)) call gronor_abort(130,"Input error")
      goto 2
    endif

    if(inp_compare(.false.,'Accelerate',item)) then
      if(.not.inp_i(naccel)) naccel=0
      goto 2
    endif

    if(inp_compare(.false.,'Idle',item)) then
      if(.not.inp_i(nidle)) nidle=0
      goto 2
    endif

    if(inp_compare(.false.,'NumGPUs',item)) then
      if(.not.inp_i(ngpus)) ngpus=1
      goto 2
    endif

    if(inp_compare(.false.,'MPService',item)) then
      if(.not.inp_i(nummps)) nummps=1
      goto 2
    endif

    if(inp_compare(.false.,'Threads',item)) then
      if(.not.inp_i(num_threads)) num_threads=1
      goto 2
    endif

    if(inp_compare(.false.,'Checkpoint',item)) then
      lcpr=.true.
      goto 2
    endif

    if (inp_compare(.false.,'Ecorrelation',item)) then
      ncorr = 1
      allocate(ecorr(mstates))
      do i=1,mstates
        ecorr(i)=0.0d0
      enddo
      !          do i=1,mstates
      !            if (.not.inp_f(ecorr(i))) ecorr(i) = 0.0d0
      !          end do
      do while(.true.)
        if(.not.inp_read()) call gronor_abort(140,"Input error")
        if(.not.inp_a(item)) goto 1
        if(len(trim(item)).gt.1) goto 3
        if(.not.inp_a(item2)) goto 1
        if(.not.inp_f(value)) goto 1
        m=0
        do i=1,nmol
          do j=1,nbase
            m=m+1
            if(fragname(i).eq.trim(item).and.fragstate(i,j).eq.trim(item2)) then
              ecorr(m)=value
              write(*,'(2i4,1x,a,1x,a,1x,f20.6)') i,j,trim(item),trim(item2),value
            endif
          enddo
        enddo
      enddo
      goto 2
    endif

    !     Use Gallup-Norbeck weights? Yes, nwt=1; No, nwt .ne. 1
    if (inp_compare(.false.,'GNWeights',item)) then
      if(.not.inp_i(nwt)) nwt = 0
      goto 2
    endif

    goto 1
  endif

  close(unit=lfninp)

  deallocate(fragstate,fragname)

  if(nmol.ge.3.and..not.lcouple) call gronor_abort(113,"Input error")

  return

999 write(lfnout,989) trim(filinp)
989 format('Unable to open input file ',a)

  call gronor_abort(150,trim(filinp))

end subroutine gronor_input
