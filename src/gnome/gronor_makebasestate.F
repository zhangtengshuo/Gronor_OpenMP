!>  @brief
!>    Construct the MEBFs from the fragment functions
!>
!>  @authors
!>    Aitor Sanchez-Mansilla, Coen de Graaf 
!>    (URV, Tarragona)
!>
!>  @date
!>    November 2020
!>
!>  @details
!>    MEBF(ibase) is made by multplying state1 and state2.
!>    First the vectors of the two state are combined, then
!>    three scenarios are foreseen: \n 
!>    (i) there is only one fragment \n 
!>    (ii) state1 and state2 are spin singlets \n 
!>    (iii) all other cases \n
!>    In the first case the fragment wave function is copied directly on the MEBF.
!>    In the second case, the coefficients of state1 and state2 are multiplied and
!>    copied to the MEBF if the product is above a certain threshold.  Case (iii) 
!>    starts with the generation of the M_S components of state 1 and 2 by sminop.
!>    By looping of the ms values of state 1 and state2, determinants of the fragments
!>    are multiplied when ms1+ms2=total spin. To ensure correct spin coupling the
!>    the product is multiplied by the corresponding Clebsch-Gordon coeefficient
!>    calculated by the function cg. \n 
!>    All products are properly multplied by (-1)^p, where p is the number of permutations
!>    to set the alpha-beta-alpha-beta-... order (function perm_ab) and to move all the
!>    doubly occupied orbitals to the left side of the determinant.
!>
!>  @param state1 many-electron state of fragment 1
!>  @param state2 many-electron state of fragment 2
!>
!>  @todo
!>    Get rid of the packing of the orbital occupations, which
!>    limits the number of active orbitals to 32
!>
!>  @todo
!>    Design a strategy for more than two fragments. There are two possibilities.\n 
!>    i) Fully general coupling, only defining the total spin and combine the fragments
!>    in all possbile ways with the total spin being the only restriction.\n 
!>    ii) Applying an intermediate coupling scheme. Fragment 1 and 2 are coupled to an
!>    intermediate spin, and coupling with fragment 3 is then performed to obtain the
!>    final total spin.\n
!>    The first is very costly (excessively?), while the second can possibly simulate
!>    the first by using several MEBFs with different coupling schemes. Taking three
!>    doublets as example: In MEBF 1 we use D1 x D2=singlet as intermediate coupling,
!>    followed by coupling with D3 to the final doublet. Adding a second MEBF with D1 x
!>    D2=triplet, followed by coupling with D3 to doublet, we mimic the fully general
!>    case of option (i).
!>
!>

      subroutine gronor_make_basestate(ibase)
      use cidef
      use cidist
      use gnome_data
      use gnome_parameters
      implicit none
      

      integer              :: ibase,state1,state2,icoupling
      integer              :: i,j,k,kk,m,iact,jact,l
      integer              :: jstart,kstart
      integer              :: ndets
      integer              :: maxcib_prev
      integer              :: target_spin
      integer              :: ms1,ms2,det,det1,det2     
      integer              :: inter_nact(2),inter_ndet(2)
      integer              :: spincomb(2)
      integer, external    :: perm_ab,isetsign4
      integer, parameter   :: maxact=32 ! because of the packing, maxact is limited to 32 (I guess)
      integer, allocatable :: spin(:)
      integer, allocatable :: ndet_ms(:),ndet1_ms(:),ndet2_ms(:)
      integer, allocatable :: occ1_ms(:,:,:),occ2_ms(:,:,:)
      integer, allocatable :: iocc_tmp(:,:)
      integer, allocatable :: inter_occ(:,:,:)
      integer, allocatable :: ioccb_copy(:,:)
      logical              :: all_singlet

      real (kind=8)              :: prod,dnorm
      real (kind=8), external    :: clebsch_gordon
      real (kind=8), allocatable :: ci_ms(:,:)
      real (kind=8), allocatable :: ci1_ms(:,:)
      real (kind=8), allocatable :: ci2_ms(:,:)
      real (kind=8), allocatable :: cicoef_tmp(:)
      real (kind=8), allocatable :: inter_civ(:,:)
      real (kind=8), allocatable :: civb_copy(:,:)

      if(idbg.ge.1.and.me.eq.master) idbg=100

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8),
     &       ' makebasestate for nmol ',nmol
        flush(lfndbg)
      endif
      
      inactb(ibase)=0
      nactb(ibase)=0
      do i=1,nmol
        inactb(ibase)=inactb(ibase)+inactm(ncombv(i,ibase))
        nactb(ibase)=nactb(ibase)+nactm(ncombv(i,ibase))
      enddo

      if ( nactb(ibase) .gt. maxact) then
        call gronor_abort(220,"Too many active orbitals")
!     write(lfnout,*)' * * *  Error  * * *'
!     write(lfnout,*)'Max. number of active orbitals: ',maxact
!     write(lfnout,*)'Actual active orbitals        : ',nactb(ibase)
!     stop 'Too many active orbitals'
      endif

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i4,a,i6)') date(1:8),time(1:8),
     &       ' number of inactive orbitals for state',ibase,
     &       ' is ',inactb(ibase)
        write(lfndbg,'(a,1x,a,1x,a,i4,a,i6)') date(1:8),time(1:8),
     &       ' number of active orbitals for state  ',ibase,
     &       ' is ',nactb(ibase)
        flush(lfndbg)
      endif
      
!     Combining the vectors of the fragments

      do i=1,nbasis
        do j=1,nbasis
          vecsb(i,j,ibase)=0.0d0
        enddo
      enddo
      jstart=1
      kstart=1
      do i=1,nmol
        do j=1,inactm(ncombv(i,ibase))
          do k=1,nbasm(ncombv(i,ibase))
            vecsb(kstart+k-1,jstart+j-1,ibase)=
     &           vecsm(k,j,ncombv(i,ibase))
          enddo
        enddo
        kstart=kstart+nbasm(ncombv(i,ibase))
        jstart=jstart+inactm(ncombv(i,ibase))
      enddo

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' vecsb step 1 completed'
        flush(lfndbg)
      endif
      
      jstart=inactb(ibase)+1
      kstart=1

      do i=1,nmol
        do j=1,nactm(ncombv(i,ibase))
          do k=1,nbasm(ncombv(i,ibase))
            vecsb(kstart+k-1,jstart+j-1,ibase)=vecsm(k,j+
     &           inactm(ncombv(i,ibase)),ncombv(i,ibase))
          enddo
        enddo
        kstart=kstart+nbasm(ncombv(i,ibase))
        jstart=jstart+nactm(ncombv(i,ibase))
      enddo

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' vecsb step 2 completed'
        flush(lfndbg)
      endif

!     Generation of the spin functions

      allocate ( spin(nmol) )
      all_singlet=.true.
      do i=1,nmol
        spin(i)=spinm(ncombv(i,ibase))
        if ( spin(i) .ne. 0 ) all_singlet=.false.
      enddo

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,l1)') date(1:8),time(1:8),
     &       ' spin completed, all_singlet is ',all_singlet
        flush(lfndbg)
      endif

      do i=1,maxcib
        civb(i,ibase)=0.0d0
        ioccb(i,ibase)=0
      enddo

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' civb/ioccb initialized'
        flush(lfndbg)
      endif

!     Just one molecule
      if ( nmol .eq. 1 ) then
        do i=1,idetm(ibase)
          do jact=1,nactm(ibase)
            iocch(jact)=ioccm(jact,i,ibase)
          enddo
          civb(i,ibase)=civm(i,ibase) * perm_ab(iocch,nactm(ibase))
          civb(i,ibase)=civb(i,ibase) * isetsign4(iocch,nactm(ibase))
          call pack(ioccb(i,ibase),iocch,nactm(ibase))
        enddo
        ndets=idetm(ibase)

        if(idbg.ge.50) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8),
     &         ' vecsb civb/ioccb completed for nmol ',nmol
          flush(lfndbg)
        endif
        goto 1581
      endif 

      state1=ncombv(1,ibase)
      state2=ncombv(2,ibase)

*     Two singlet fragments
      if ( all_singlet .and. nmol.eq.2 ) then
        k=0
        do i=1,idetm(state1)
          do j=1,idetm(state2)
            prod=civm(i,state1) * civm(j,state2)
            if ( abs(prod) .gt. tau_CI ) then
              k=k+1
              civb(k,ibase)=prod
              kk=0
              do iact=1,nactm(state1)
                kk=kk+1
                iocch(kk)=ioccm(iact,i,state1)
              enddo
              do jact=1,nactm(state2)
                kk=kk+1
                iocch(kk)=ioccm(jact,j,state2)
              enddo
              call pack(ioccb(k,ibase),iocch,nactb(ibase))
              civb(k,ibase)=civb(k,ibase) * 
     &             perm_ab(iocch,nactb(ibase))
              civb(k,ibase)=civb(k,ibase) * 
     &             isetsign4(iocch,nactb(ibase))
            endif
          enddo
        enddo
        ndets=k
        if(idbg.ge.50) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8),
     &         ' vecsb civb/ioccb completed for nmol ',nmol
          flush(lfndbg)
        endif
        goto 1581
      endif
      
      if (nmol.gt.2) then
        if(idbg.ge.50) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8),
     &         ' intercouplings for nmol ',nmol
          flush(lfndbg)
        endif
!     write(lfnout,*)'INTER COUPLINGS'
!     do i=1,nbase
!     write(lfnout,'(200I2)')(inter_couplings(j,i),j=1,nmol-1)
!     enddo
!     flush(lfnout)
      endif
      
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' allocations spin(1)+1 ',spin(1)+1
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' allocations spin(2)+1 ',spin(2)+1
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' allocations ncombv(1) ',ncombv(1,ibase)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' allocations ncombv(2) ',ncombv(2,ibase)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' allocations nci       ',nci
        flush(lfndbg)
      endif
      
*     All other cases
      
      spincomb(1)=spin(1)
      spincomb(2)=spin(2)

      allocate(ndet1_ms(spin(1)+1))
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' ndet1_ms allocated'
        flush(lfndbg)
      endif
      allocate(ci1_ms(nci,spin(1)+1))
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' ci1_ms allocated'
        flush(lfndbg)
      endif
      allocate(occ1_ms(nactm(ncombv(1,ibase)),nci,spin(1)+1))
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' occ1_ms allocated'
        flush(lfndbg)
      endif      
      allocate(ndet2_ms(spin(2)+1))
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' ndet2_ms allocated'
        flush(lfndbg)
      endif      
      allocate(ci2_ms(nci,spin(2)+1))
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' ci2_ms allocated'
        flush(lfndbg)
      endif      
      allocate(occ2_ms(nactm(ncombv(2,ibase)),nci,spin(2)+1))
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' occ2_ms allocated'
        flush(lfndbg)
      endif

      do i=1,spin(1)+1
        ndet1_ms(i)=0
      enddo
      do i=1,spin(1)+1
        do j=1,nci
          ci1_ms(j,i)=0.0d0
        enddo
      enddo
      do i=1,spin(1)+1
        do j=1,nci
          do k=1,nactm(ncombv(1,ibase))
            occ1_ms(k,j,i)=0
          enddo
        enddo
      enddo
      do i=1,spin(2)+1
        ndet2_ms(i)=0
      enddo
      do i=1,spin(2)+1
        do j=1,nci
          ci2_ms(j,i)=0.0d0
        enddo
      enddo
      do i=1,spin(2)+1
        do j=1,nci
          do k=1,nactm(ncombv(2,ibase))
            occ2_ms(k,j,i)=0
          enddo
        enddo
      enddo
      
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' allocations maxcib    ',maxcib
        write(lfndbg,'(a,1x,a,1x,a,2i10)') date(1:8),time(1:8),
     &       ' allocations nactb(ibase)   ',nactb(ibase),maxact
        flush(lfndbg)
      endif
      
*     Allocate all intermediate arrays for multi-fragment mebfs      

      allocate(inter_civ(maxcib,2))
      allocate(inter_occ(nactb(ibase),maxcib,2))
      
      do i=1,maxcib
        inter_civ(i,1)=0.0d0
        inter_civ(i,2)=0.0d0
      enddo
      do i=1,maxcib
        do j=1,nactb(ibase)
          inter_occ(j,i,1)=0
          inter_occ(j,i,2)=0
        enddo
      enddo
      
      inter_nact(1)=nactm(ncombv(1,ibase))
      inter_nact(2)=nactm(ncombv(2,ibase))
      inter_ndet(1)=idetm(ncombv(1,ibase))
      inter_ndet(2)=idetm(ncombv(2,ibase))
      
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' allocations completed'
        flush(lfndbg)
      endif
      
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,6i10)') date(1:8),time(1:8),
     &       ' civ/occ ',maxcib,nactb(ibase),
     &       ncombv(1,ibase),ncombv(2,ibase),mstates
        do i=1,mstates
          write(lfndbg,'(a,1x,a,1x,a,6i10)') date(1:8),time(1:8),
     &         ' idetm ',i,idetm(i)
        enddo
        flush(lfndbg)
        do i=1,10
          write(lfndbg,'(a,i5,2e20.10)') 'CIVM',i,
     &         civm(i,1),civm(i,2)
          flush(lfndbg)
        enddo
      endif
      
      do i=1,inter_ndet(1)
        inter_civ(i,1)=civm(i,ncombv(1,ibase))
        do j=1,inter_nact(1)
          inter_occ(j,i,1)=ioccm(j,i,ncombv(1,ibase))
        enddo
      enddo
      do i=1,inter_ndet(2)
        inter_civ(i,2)=civm(i,ncombv(2,ibase))
        do j=1,inter_nact(2)
          inter_occ(j,i,2)=ioccm(j,i,ncombv(2,ibase))
        enddo
      enddo
      
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' inter_civ, inter_occ initialized'
        flush(lfndbg)
      endif
      
!     Make the compound product of all couplings

      do icoupling=1,nmol-1
        if(idbg.ge.100) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &         ' coupling for ',icoupling
          flush(lfndbg)
        endif
        allocate(iocc_tmp(inter_nact(1),inter_ndet(1)))
        allocate(cicoef_tmp(inter_ndet(1)))        
        iocc_tmp=0
        cicoef_tmp=0
        do i=1,inter_ndet(1)
          cicoef_tmp(i)=inter_civ(i,1)
          do j=1,inter_nact(1)
            iocc_tmp(j,i)=inter_occ(j,i,1)
          enddo
        enddo
        call sminop(cicoef_tmp,ci1_ms,iocc_tmp,occ1_ms,inter_ndet(1),
     &       ndet1_ms,spincomb(1),inter_nact(1),nci)        

        deallocate(cicoef_tmp)
        deallocate(iocc_tmp)
        
        if(idbg.ge.50) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a,10i10)') date(1:8),time(1:8),
     &         ' Frag 1 ndet1_ms ',ndet1_ms
          do i=1,spincomb(1)+1
            write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &           ' ms state ',i
            do j=1,ndet1_ms(i)
              write(lfndbg,'(2i5,f12.8,2x,100i2)') j,i,ci1_ms(j,i),
     &             (occ1_ms(k,j,i),k=1,inter_nact(1))
            enddo
          enddo
          flush(lfndbg)
        endif
        
!     M_S of fragment 2
        
        allocate(iocc_tmp(inter_nact(2),inter_ndet(2)))
        allocate(cicoef_tmp(inter_ndet(2)))
        iocc_tmp=0
        cicoef_tmp=0
        do i=1,inter_ndet(2)
          cicoef_tmp(i)=inter_civ(i,2)
          do j=1,inter_nact(2)
            iocc_tmp(j,i)=inter_occ(j,i,2)
          enddo
        enddo

        call sminop(cicoef_tmp,ci2_ms,iocc_tmp,occ2_ms,inter_ndet(2),
     &       ndet2_ms,spincomb(2),inter_nact(2),nci)
        deallocate(cicoef_tmp)
        deallocate(iocc_tmp)
        
        if(idbg.ge.50) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a,10i10)') date(1:8),time(1:8),
     &         ' Frag 2 ndet2_ms ',ndet2_ms
          do i=1,spincomb(2)+1
            write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &           ' ms state ',i
            do j=1,ndet2_ms(i)
              write(lfndbg,'(2i5,f12.8,2x,100i2)') j,i,ci2_ms(j,i),
     &             (occ2_ms(k,j,i),k=1,inter_nact(2))
            enddo
          enddo
          flush(lfndbg)
        endif
        
!     Product

        if(nmol.gt.2)then
          target_spin=inter_couplings(icoupling,ibase) ! Remove nspin in the future
        else
          target_spin=nspin
        endif
        
!     Determine new maxcib for icoupling > 1 and reallocate the product arrays
        
        maxcib_prev=maxcib
        if (icoupling.gt.1)then
          do i=1,spincomb(1)+1
            ms1=spincomb(1)-2*(i-1)
            do j=1,spincomb(2)+1
              ms2=spincomb(2)-2*(j-1)
              if ( ms1+ms2 .eq. target_spin ) then
                do det1=1,ndet1_ms(i)
                  do det2=1,ndet2_ms(j)                  
                    prod=ci1_ms(det1,i)*ci2_ms(det2,j)
                    if (abs(prod) .gt. tau_CI)k=k+1
                  enddo
                enddo
              endif
            enddo
          enddo

          if(k .eq. 0) then
            if(idbg.ge.50) then
              call swatch(date,time)
              write(lfndbg,'(a,1x,a,1x,a,3i4)') date(1:8),time(1:8),
     &             ' Incompatible spins ',
     &             spincomb(1)+1,spincomb(2)+1,target_spin+1
              flush(lfndbg)
            endif
            call gronor_abort(221,"Incompatible spins")
          endif
          
          if(idbg.ge.50) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &           ' Old maxcib ',maxcib
            flush(lfndbg)
          endif

          maxcib=max(k,maxcib)

          if(idbg.ge.50) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &           ' New maxcib ',maxcib
            flush(lfndbg)
          endif

!     Reallocate intermediate civ and occ if couplings loop is not ended
!     otherwise reallocate civb and ioccb in order to store the mebfs

          deallocate(inter_civ,inter_occ)
          allocate(inter_civ(maxcib,2),inter_occ(nactb(ibase),maxcib,2))
          
          do i=1,maxcib
            inter_civ(i,1)=0.0d0
            inter_civ(i,2)=0.0d0
          enddo
          do i=1,maxcib
            do j=1,nactb(ibase)
              inter_occ(j,i,1)=0
              inter_occ(j,i,2)=0
            enddo
          enddo

          if (icoupling.eq.(nmol-1)) then

            allocate(civb_copy(maxcib,nbase))
            allocate(ioccb_copy(maxcib,nbase))
            do i=1,nbase
              do j=1,idetb(i)
                civb_copy(j,i)=civb(j,i)
                ioccb_copy(j,i)=ioccb(j,i)
              enddo
            enddo
            deallocate(civb)
            deallocate(ioccb)
            allocate(civb(maxcib,nbase))
            allocate(ioccb(maxcib,nbase))            
            civb=0.0d0
            ioccb=0
            do i=1,nbase
              do j=1,idetb(i)
                civb(j,i)=civb_copy(j,i)
                ioccb(j,i)=ioccb_copy(j,i)
              enddo
            enddo
            deallocate(civb_copy,ioccb_copy)
          endif
          
          if(idbg.ge.50) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &           ' Reallocated civb and ioccb'
            flush(lfndbg)
          endif

        endif
        
        if(idbg.ge.50) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &         ' Product state result'
          flush(lfndbg)
        endif

        k=0
        do i=1,spincomb(1)+1
          ms1=spincomb(1)-2*(i-1)
          do j=1,spincomb(2)+1
            ms2=spincomb(2)-2*(j-1)
            if ( ms1+ms2 .eq. target_spin ) then
              do det1=1,ndet1_ms(i)
                do det2=1,ndet2_ms(j)
                  dnorm=clebsch_gordon(spincomb(1)/2.d0,
     &                 spincomb(2)/2.d0,ms1/2.d0,ms2/2.d0,
     &                 target_spin/2.d0,target_spin/2.d0)
                  prod=ci1_ms(det1,i)*ci2_ms(det2,j)                                                      
                  if (abs(prod) .gt. tau_CI) then
                    iocch=0
                    k=k+1                    
                    inter_civ(k,1)=prod * dnorm                    
                    kk=0
                    do iact=1,inter_nact(1)
                      kk=kk+1
                      iocch(kk)=occ1_ms(iact,det1,i)
                    enddo
                    do jact=1,inter_nact(2)
                      kk=kk+1
                      iocch(kk)=occ2_ms(jact,det2,j)
                    enddo
                    if (icoupling.eq.(nmol-1)) then ! if loop is ended, save in civb and ioccb 

                      call pack(ioccb(k,ibase),
     &                     iocch,nactb(ibase))
                      civb(k,ibase)=inter_civ(k,1) * 
     &                     perm_ab(iocch,nactb(ibase))
                      civb(k,ibase)=civb(k,ibase) * 
     &                     isetsign4(iocch,nactb(ibase))
                      
                      if(idbg.ge.50) then
                        write(lfndbg,'(2i5,f12.8,2x,100i2)')
     &                       k,ibase,civb(k,ibase),
     &                       (iocch(m),m=1,nactb(ibase))
                        flush(lfndbg)
                      endif
                    else
!     if loop is not ended, prepare for the next iteration                      
                      do iact=1,kk
                        inter_occ(iact,k,1)=iocch(iact)
                      enddo                      
                    endif
                  endif
                enddo
              enddo
            endif
          enddo
        enddo
        
!     After the product the result is in inter_civ(:,1) or in civb(:,ibase) if the loop is ended
!     Prepare the next iteration if the loop is not ended
        if (icoupling.lt.(nmol-1)) then
          
          if(idbg.ge.50) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &           ' Next coupling ',icoupling
            flush(lfndbg)
          endif
          
          inter_nact(1)=inter_nact(1)+inter_nact(2)
          inter_nact(2)=nactm(ncombv(icoupling+2,ibase))
          inter_ndet(1)=k
          inter_ndet(2)=idetm(ncombv(icoupling+2,ibase))          
          
          do i=1,inter_ndet(2)
            inter_civ(i,2)=civm(i,ncombv(icoupling+2,ibase))
            do j=1,inter_nact(2)
              inter_occ(j,i,2)=ioccm(j,i,ncombv(icoupling+2,ibase))
            enddo
          enddo

          if(idbg.ge.50) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &           ' Recalculating nci'
            flush(lfndbg)
          endif
          
          spincomb(1)=inter_couplings(icoupling,ibase)
          spincomb(2)=spin(icoupling+2)
          call gronor_determine_nci
     &         (2,inter_nact,
     &         inter_ndet,inter_occ,nactb(ibase),maxcib_prev)
          if(idbg.ge.50) then
            call swatch(date,time)
            write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &           ' New nci ',nci
            write(lfndbg,'(2i5,f12.8,2x,100i2)') i,1,inter_civ(i,1),
     &           (inter_occ(j,i,1),j=1,inter_nact(1))
            write(lfndbg,'(2i5,f12.8,2x,100i2)') i,2,inter_civ(i,2),
     &           (inter_occ(j,i,2),j=1,inter_nact(2))
            flush(lfndbg)
          endif                   

          deallocate(ndet1_ms,occ1_ms,ci1_ms)            
          deallocate(ndet2_ms,occ2_ms,ci2_ms)
          allocate(ndet1_ms(spincomb(1)+1),ndet2_ms(spincomb(2)+1))          
          allocate(ci1_ms(nci,spincomb(1)+1))
          allocate(ci2_ms(nci,spincomb(2)+1))
          allocate(occ1_ms(inter_nact(1),nci,spincomb(1)+1))
          allocate(occ2_ms(inter_nact(2),nci,spincomb(2)+1))
          
          do i=1,spincomb(1)+1
            ndet1_ms(i)=0
          enddo
          do i=1,spincomb(1)+1
            do j=1,nci
              ci1_ms(j,i)=0.0d0
            enddo
          enddo
          do i=1,spincomb(1)+1
            do j=1,nci
              do k=1,inter_nact(1)
                occ1_ms(k,j,i)=0
              enddo
            enddo
          enddo
          do i=1,spincomb(2)+1
            ndet2_ms(i)=0
          enddo
          do i=1,spincomb(2)+1
            do j=1,nci
              ci2_ms(j,i)=0.0d0
            enddo
          enddo
          do i=1,spincomb(2)+1
            do j=1,nci
              do k=1,inter_nact(2)
                occ2_ms(k,j,i)=0
              enddo
            enddo
          enddo
          
        endif
        
        if(idbg.ge.50) then
          call swatch(date,time)
          write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &         ' couplings loop completed'
          flush(lfndbg)
        endif
      enddo                     ! end couplings loop

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' deallocating arrays'
        flush(lfndbg)
      endif
      
      deallocate(ndet1_ms,occ1_ms,ci1_ms)
      deallocate(ndet2_ms,occ2_ms,ci2_ms)
      if(allocated(inter_civ)) deallocate(inter_civ)
      if(allocated(inter_occ)) deallocate(inter_occ)

      ndets=k
      
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' deallocated arrays, ndets is ',ndets
        flush(lfndbg)
      endif
      
 1581 continue
      
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i4,a,i10,a)') date(1:8),time(1:8),
     &       ' completed for ',nmol,' with ',ndets,' determinants'
        flush(lfndbg)
      endif
      
      dnorm=0.0d0
      idetb(ibase)=ndets
      do i=1,idetb(ibase)
        dnorm=dnorm+civb(i,ibase)**2
      enddo
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' dnorm completed'
        flush(lfndbg)
      endif
      dnorm=1/dsqrt(dnorm)
      do i=1,idetb(ibase)
        civb(i,ibase)=civb(i,ibase) * dnorm
      enddo
      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' civb completed'
        flush(lfndbg)
      endif
      deallocate(spin)

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i10)') date(1:8),time(1:8),
     &       ' make_basestate completed ',nactb(ibase)
        flush(lfndbg)
      endif
      
      return
      end 
