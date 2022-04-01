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
!>    are multiplied when ms1 + ms2 = total spin. To ensure correct spin coupling the
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
!>    doublets as example: In MEBF 1 we use D1 x D2 = singlet as intermediate coupling,
!>    followed by coupling with D3 to the final doublet. Adding a second MEBF with D1 x
!>    D2 = triplet, followed by coupling with D3 to doublet, we mimic the fully general
!>    case of option (i).
!>
!>

      subroutine gronor_make_basestate(ibase)
      use cidef
      use gnome_data
      use gnome_parameters
      implicit none
      

      integer              :: ibase,state1,state2,icoupling
      integer              :: i,j,k,kk,iact,jact,l
      integer              :: jstart,kstart
      integer              :: ndets
      integer              :: maxcib_prev
      integer              :: target_spin
      integer              :: ms1,ms2,det,det1,det2
      integer, external    :: perm_ab,isetsign4
      integer, parameter   :: maxact = 32 ! because of the packing, maxact is limited to 32 (I guess)
      integer, allocatable :: spin(:),spincomb(:)
      integer, allocatable :: ndet_ms(:),ndet1_ms(:),ndet2_ms(:)
      integer, allocatable :: occ1_ms(:,:,:),occ2_ms(:,:,:)
      integer, allocatable :: iocc_tmp(:,:)
      integer, allocatable :: inter_nact(:),inter_ndet(:)
      integer, allocatable :: inter_occ(:,:,:)
      integer, allocatable :: ioccb_copy(:,:)
      logical              :: all_singlet

      real (kind = 8)              :: prod,dnorm
      real (kind = 8), external    :: clebsch_gordon
      real (kind = 8), allocatable :: ci_ms(:,:)
      real (kind = 8), allocatable :: ci1_ms(:,:)
      real (kind = 8), allocatable :: ci2_ms(:,:)
      real (kind = 8), allocatable :: cicoef_tmp(:)
      real (kind = 8), allocatable :: inter_civ(:,:)
      real (kind = 8), allocatable :: civb_copy(:,:)


      write(lfnout,*)'start gronor_makebasestate'
      flush(lfnout)
      inactb(ibase) = 0
      nactb(ibase)  = 0
      do i = 1, nmol
        inactb(ibase) = inactb(ibase) + inactm(ncombv(i,ibase))
        nactb(ibase)  = nactb(ibase) + nactm(ncombv(i,ibase))
      end do
      write(lfnout,*)'NACTB makebasestate',ibase,nactb(ibase)
      write(lfnout,*)'combination',ncombv(:,ibase)
      if ( nactb(ibase) .gt. maxact) then
        write(lfnout,*)' * * *  Error  * * *'
        write(lfnout,*)'Max. number of active orbitals: ',maxact
        write(lfnout,*)'Actual active orbitals        : ',nactb(ibase)
        stop 'Too many active orbitals'
      endif

*     Combining the vectors of the fragments

      do i = 1, nbasis
        do j = 1, nbasis
          vecsb(i,j,ibase) = 0.0d0
        enddo
      enddo
      jstart = 1
      kstart = 1
      do i = 1, nmol
        do j = 1, inactm(ncombv(i,ibase))
          do k = 1, nbasm(ncombv(i,ibase))
            vecsb(kstart+k-1,jstart+j-1,ibase) =
     &           vecsm(k,j,ncombv(i,ibase))
          end do
        end do
        kstart = kstart +  nbasm(ncombv(i,ibase))
        jstart = jstart + inactm(ncombv(i,ibase))
      end do

      jstart = inactb(ibase) + 1
      kstart = 1

      do i = 1, nmol
        do j = 1, nactm(ncombv(i,ibase))
          do k = 1, nbasm(ncombv(i,ibase))
            vecsb(kstart + k-1, jstart + j-1,ibase) = vecsm(k, j +
     &           inactm(ncombv(i,ibase)),ncombv(i,ibase))
          end do
        end do
        kstart = kstart + nbasm(ncombv(i,ibase))
        jstart = jstart + nactm(ncombv(i,ibase))
      end do

*     Generation of the spin functions

      allocate ( spin(nmol) )
      all_singlet = .true.
      do i = 1, nmol
        spin(i) = spinm(ncombv(i,ibase))
        if ( spin(i) .ne. 0 ) all_singlet = .false.
      end do

      do i = 1, maxcib
        civb(i,ibase) = 0.0d0
        ioccb(i,ibase)= 0
      enddo

*     Just one molecule
      if ( nmol .eq. 1 ) then
        do i = 1, idetm(ibase)
          do jact = 1, nactm(ibase)
            iocch(jact) = ioccm(jact,i,ibase)
          end do
          civb(i,ibase) = civm(i,ibase) * perm_ab(iocch,nactm(ibase))
          civb(i,ibase) = civb(i,ibase) * isetsign4(iocch,nactm(ibase))
          call pack(ioccb(i,ibase),iocch,nactm(ibase))
        end do
        ndets = idetm(ibase)
        goto 1581
      end if 

      state1 = ncombv(1,ibase)
      state2 = ncombv(2,ibase)

*     Two singlet fragments
      if ( all_singlet .and. nmol.eq.2 ) then
        k = 0
        do i = 1, idetm(state1)
          do j = 1, idetm(state2)
            prod = civm(i,state1) * civm(j,state2)
            if ( abs(prod) .gt. tau_CI ) then
              k = k + 1
              civb(k,ibase) = prod
              kk = 0
              do iact = 1, nactm(state1)
                kk = kk + 1
                iocch(kk) = ioccm(iact,i,state1)
              end do
              do jact = 1, nactm(state2)
                kk = kk + 1
                iocch(kk) = ioccm(jact,j,state2)
              end do
              call pack(ioccb(k,ibase),iocch,nactb(ibase))
              civb(k,ibase) = civb(k,ibase) * 
     &             perm_ab(iocch,nactb(ibase))
              civb(k,ibase) = civb(k,ibase) * 
     &             isetsign4(iocch,nactb(ibase))
            end if
          end do
        end do
        ndets = k
        goto 1581
      end if
      
      if (nmol.gt.2) then
        write(lfnout,*)'INTER COUPLINGS'
        do i=1,nbase
          write(lfnout,*)(inter_couplings(j,i),j=1,nmol-1)
        end do
      end if
*     All other cases
      allocate(spincomb(2))
      
      spincomb(1)=spin(1)
      spincomb(2)=spin(2)
      write(lfnout,*)'SPINCOMB'
      write(lfnout,*)spincomb
      allocate(ndet1_ms(spin(1)+1))
      allocate(ci1_ms(nci,spin(1)+1))
      allocate(occ1_ms(nactm(ncombv(1,ibase)),nci,spin(1)+1))      
      allocate(ndet2_ms(spin(2)+1))
      allocate(ci2_ms(nci,spin(2)+1))
      allocate(occ2_ms(nactm(ncombv(2,ibase)),nci,spin(2)+1))
      ndet1_ms = 0
      ci1_ms   = 0.0d0
      occ1_ms  = 0
      ndet2_ms = 0
      ci2_ms   = 0.0d0
      occ2_ms  = 0
      
*     Allocate all intermediate arrays for multi-fragment mebfs      
      allocate(inter_nact(2),inter_ndet(2))
      allocate(inter_civ(maxcib,2))
      allocate(inter_occ(maxact,maxcib,2))
      inter_nact = 0
      inter_ndet = 0
      inter_civ = 0.0d0
      inter_occ = 0
      inter_nact(1)=nactm(ncombv(1,ibase))
      inter_nact(2)=nactm(ncombv(2,ibase))
      inter_ndet(1)=idetm(ncombv(1,ibase))
      inter_ndet(2)=idetm(ncombv(2,ibase))
      
      do i=1,maxcib
        inter_civ(i,1)=civm(i,ncombv(1,ibase))
        inter_civ(i,2)=civm(i,ncombv(2,ibase))
        do j=1,maxnact
          inter_occ(j,i,1)=ioccm(j,i,ncombv(1,ibase))
          inter_occ(j,i,2)=ioccm(j,i,ncombv(2,ibase))
        enddo
      enddo
*     Make the compound product of all couplings

      do icoupling=1,nmol-1
*       M_S of fragment 1
        write(lfnout,*)'+++++++++++++++++++++++++++++'
        write(lfnout,*)'ICOUPLING',icoupling
        write(lfnout,*)'+++++++++++++++++++++++++++++'
        write(lfnout,*)
        write(lfnout,*)'nci: ',nci
        write(lfnout,*)'maxcib: ',maxcib
        write(lfnout,*)'spin combination',spincomb
        allocate(iocc_tmp(inter_nact(1),inter_ndet(1)))
        allocate(cicoef_tmp(inter_ndet(1)))        
        iocc_tmp = 0
        cicoef_tmp = 0
        do i = 1, inter_ndet(1)
          cicoef_tmp(i) = inter_civ(i,1)
          do j = 1, inter_nact(1)
            iocc_tmp(j,i) = inter_occ(j,i,1)
          end do
        end do
        call sminop2(cicoef_tmp,ci1_ms,iocc_tmp,occ1_ms,inter_ndet(1),
     &       ndet1_ms,spincomb(1),inter_nact(1),nci)        

        deallocate(cicoef_tmp)
        deallocate(iocc_tmp)        
        write(lfnout,*)'FRAGMENT 1'
        write(lfnout,*)'NDETMS'
        write(lfnout,*)ndet1_ms
        do i=1,spincomb(1)+1
          write(lfnout,*)'MS STATE',i
          do j=1,ndet1_ms(i)
            write(lfnout,*)ci1_ms(j,i),
     &           (occ1_ms(k,j,i),k=1,inter_nact(1))
          end do
          flush(lfnout)
        end do
        
*       M_S of fragment 2
        
        allocate(iocc_tmp(inter_nact(2),inter_ndet(2)))
        allocate(cicoef_tmp(inter_ndet(2)))
        iocc_tmp = 0
        cicoef_tmp = 0
        do i = 1, inter_ndet(2)
          cicoef_tmp(i) = inter_civ(i,2)
          do j = 1, inter_nact(2)
            iocc_tmp(j,i) = inter_occ(j,i,2)
          end do
        end do
        call sminop2(cicoef_tmp,ci2_ms,iocc_tmp,occ2_ms,inter_ndet(2),
     &       ndet2_ms,spincomb(2),inter_nact(2),nci)
        deallocate(cicoef_tmp)
        deallocate(iocc_tmp)
        
        
        write(lfnout,*)'FRAGMENT 2'
        write(lfnout,*)'NDETMS'
        write(lfnout,*)ndet2_ms
        do i=1,spincomb(2)+1
          write(lfnout,*)'MS STATE',i
          do j=1,ndet2_ms(i)
            write(lfnout,*)ci2_ms(j,i),
     &           (occ2_ms(k,j,i),k=1,inter_nact(2))
          end do
          flush(lfnout)
        end do
*       Product
        if(nmol.gt.2)then
          target_spin=inter_couplings(icoupling,ibase) ! Remove nspin in the future
        else
          target_spin=nspin
        end if
*       Determine new maxcib for icoupling > 1 and reallocate the product arrays
        maxcib_prev=maxcib
        if (icoupling.gt.1)then
          do i = 1, spincomb(1) + 1
            ms1  = spincomb(1) - 2*(i-1)
            do j = 1, spincomb(2) + 1
              ms2 = spincomb(2) - 2*(j-1)
              if ( ms1 + ms2 .eq. target_spin ) then
                do det1 = 1, ndet1_ms(i)
                  do det2 = 1, ndet2_ms(j)                  
                    prod = ci1_ms(det1,i)*ci2_ms(det2,j)
                    if (abs(prod) .gt. tau_CI)k=k+1
                  end do
                end do
              end if
            end do
          end do
          if(k .eq. 0) then
            write(lfnout,1574) '   spin(1) = ',spincomb(1)+1
            write(lfnout,1574) '   spin(2) = ',spincomb(2)+1
            write(lfnout,1574) 'total spin = ',target_spin+1
            write(lfnout,*) 'not possible'
            flush(lfnout)
            stop
 1574       format(a,i4)
          end if
          write(lfnout,*)'old maxcib',maxcib
          maxcib=max(k,maxcib)
          write(lfnout,*)'new maxcib',maxcib
*         Reallocate intermediate civ and occ if couplings loop is not ended
*         otherwise reallocate civb and ioccb in order to store the mebfs          
          deallocate(inter_civ,inter_occ)
          allocate(inter_civ(maxcib,2),inter_occ(maxact,maxcib,2))          
          inter_civ=0.0d0
          inter_occ=0
          if (icoupling.eq.(nmol-1)) then
            write(lfnout,*)'Last coupling, re-allocate maxcib'
            allocate(civb_copy(maxcib,nbase))
            allocate(ioccb_copy(maxcib,nbase))
            do i=1,nbase
              do j=1,idetb(i)
                civb_copy(j,i)=civb(j,i)
                ioccb_copy(j,i)=ioccb(j,i)
              end do
            end do
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
              end do
            end do
            deallocate(civb_copy,ioccb_copy)
          end if

        end if
        write(lfnout,*)'Do the product'
        k = 0
        do i = 1, spincomb(1) + 1
          ms1  = spincomb(1) - 2*(i-1)
          do j = 1, spincomb(2) + 1
            ms2 = spincomb(2) - 2*(j-1)
            if ( ms1 + ms2 .eq. target_spin ) then
              do det1 = 1, ndet1_ms(i)
                do det2 = 1, ndet2_ms(j)
                  dnorm = clebsch_gordon(spincomb(1)/2.d0,
     &                 spincomb(2)/2.d0,ms1/2.d0,ms2/2.d0,
     &                 target_spin/2.d0,target_spin/2.d0)
                  prod = ci1_ms(det1,i)*ci2_ms(det2,j)                                                      
                  if (abs(prod) .gt. tau_CI) then
                    iocch = 0
                    k = k + 1                    
                    inter_civ(k,1) = prod * dnorm                    
                    kk = 0
                    do iact = 1, inter_nact(1)
                      kk = kk + 1
                      iocch(kk) = occ1_ms(iact,det1,i)
                    end do
                    do jact = 1, inter_nact(2)
                      kk = kk + 1
                      iocch(kk) = occ2_ms(jact,det2,j)
                    end do
                    if (icoupling.eq.(nmol-1)) then ! if loop is ended, save in civb and ioccb 
                      call pack(ioccb(k,ibase),iocch,nactb(ibase))
                      civb(k,ibase) = inter_civ(k,1) * 
     &                     perm_ab(iocch,nactb(ibase))
                      civb(k,ibase) = civb(k,ibase) * 
     &                     isetsign4(iocch,nactb(ibase))
                    else ! if loop is not ended, prepare for the next iteration                      
                      do iact=1,kk
                        inter_occ(iact,k,1)=iocch(iact)
                      end do                      
                    end if
                  end if
                end do
              end do
            end if
          end do
        end do
        
*       After the product the result is in inter_civ(:,1) or in civb(:,ibase) if the loop is ended
*       Prepare the next iteration if the loop is not ended
        if (icoupling.lt.(nmol-1)) then
          write(lfnout,*)'not last coupling',icoupling
          
          inter_nact(1)=inter_nact(1)+inter_nact(2)
          inter_nact(2)=nactm(ncombv(icoupling+2,ibase))
          inter_ndet(1)=k
          inter_ndet(2)=idetm(ncombv(icoupling+2,ibase))          
          
          do i=1,inter_ndet(2)
            inter_civ(i,2) = civm(i,ncombv(icoupling+2,ibase))
            do j=1,inter_nact(2)
              inter_occ(j,i,2) = ioccm(j,i,ncombv(icoupling+2,ibase))
            enddo
          enddo

          write(lfnout,*)'re-calculate nci'
          spincomb(1)=inter_couplings(icoupling,ibase)
          spincomb(2)=spin(icoupling+2)
          call gronor_determine_nci
     &         (2,inter_nact,
     &         inter_ndet,inter_occ,maxact,maxcib_prev)
          write(lfnout,*)'new nci:',nci
          write(lfnout,*)'NEW inter 1'
          do i=1,inter_ndet(1)
            write(lfnout,*)inter_civ(i,1),
     &           (inter_occ(j,i,1),j=1,inter_nact(1))
          enddo
          write(lfnout,*)'NEW inter 2'
          do i=1,inter_ndet(2)
            write(lfnout,*)inter_civ(i,2),
     &           (inter_occ(j,i,2),j=1,inter_nact(2))
          enddo
          deallocate(ndet1_ms,occ1_ms,ci1_ms)            
          deallocate(ndet2_ms,occ2_ms,ci2_ms)
          allocate(ndet1_ms(spincomb(1)+1),ndet2_ms(spincomb(2)+1))          
          allocate(ci1_ms(nci,spincomb(1)+1),ci2_ms(nci,spincomb(2)+1))
          allocate(occ1_ms(inter_nact(1),nci,spincomb(1)+1))
          allocate(occ2_ms(inter_nact(2),nci,spincomb(2)+1))
          ndet1_ms = 0
          ndet2_ms = 0
          ci1_ms = 0.0d0
          ci2_ms = 0.0d0
          occ1_ms = 0
          occ2_ms = 0
        end if
      end do ! end couplings loop
*     Deallocate all
      deallocate(ndet1_ms,occ1_ms,ci1_ms)            
      deallocate(ndet2_ms,occ2_ms,ci2_ms)
      deallocate(spincomb)
      deallocate(inter_nact,inter_ndet)
      deallocate(inter_civ,inter_occ)
      ndets = k
 1581 dnorm = 0.0d0
      idetb(ibase) = ndets
      do i = 1, idetb(ibase)
        dnorm = dnorm + civb(i,ibase)**2
      end do
      dnorm = 1/dsqrt(dnorm)
      do i = 1, idetb(ibase)
        civb(i,ibase) = civb(i,ibase) * dnorm
      end do
      deallocate(spin)
      write(lfnout,*)'nactb ibase final makebasestate',nactb(ibase)
      return
      end 
