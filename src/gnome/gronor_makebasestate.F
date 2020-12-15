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

      subroutine gronor_make_basestate(ibase,maxcib0)
      use cidef
      use gnome_data
      use gnome_parameters
      implicit none
      

      integer              :: ibase,state1,state2
      integer              :: i,j,k,l,kk,iact,jact,icoupling
      integer              :: jstart,kstart      
      integer              :: ndets,ms_states,maxcib0
      integer              :: ms1,ms2,det,det1,det2
      integer              :: inter_nactb,inter_maxnact
      integer, external    :: perm_ab,isetsign4
      integer, parameter   :: maxact = 32 ! because of the packing, maxact is limited to 32 (I guess)
      integer              :: ioccb_copy(maxcib,nbase)      
      integer, allocatable :: spin(:)
      integer, allocatable :: ndet_ms(:),ndet1_ms(:),ndet2_ms(:)
      integer, allocatable :: occ_ms(:,:,:),occ1_ms(:,:,:),
     &     occ2_ms(:,:,:)
      integer, allocatable :: spin_comb(:),inter_idetm(:),
     &     inter_ioccm(:,:,:),inter_ioccb(:),inter_nactm(:)

      logical              :: all_singlet

      real (kind = 8)              :: prod,dnorm
      real (kind = 8), external    :: clebsch_gordon
      real (kind=8)                :: civb_copy(maxcib,nbase)
      real (kind = 8), allocatable :: ci_ms(:,:)
      real (kind = 8), allocatable :: ci1_ms(:,:)
      real (kind = 8), allocatable :: ci2_ms(:,:)
      real (kind = 8), allocatable :: inter_civm(:,:)
      real (kind = 8), allocatable :: inter_civb(:)      


      inactb(ibase) = 0
      nactb(ibase)  = 0
      do i = 1, nmol
         inactb(ibase) = inactb(ibase) + inactm(ncombv(i,ibase))
         nactb(ibase)  = nactb(ibase) + nactm(ncombv(i,ibase))
      end do

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
     &              vecsm(k,j,ncombv(i,ibase))
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
     &              inactm(ncombv(i,ibase)),ncombv(i,ibase))
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
      write(lfnout,*)'spin'
      write(lfnout,*)spin
      flush(lfnout)
*     Just one molecule
      if ( nmol .eq. 1 ) then
         do i = 1, idetm(ibase)
            do jact = 1, nactm(ibase)
               iocch(jact) = ioccm(jact,i,ibase)
            end do
            civb(i,ibase) = civm(i,ibase)*perm_ab(iocch,nactm(ibase))
            civb(i,ibase) = civb(i,ibase)*isetsign4(iocch,nactm(ibase))
            call pack(ioccb(i,ibase),iocch,nactm(ibase))
         end do
         ndets = idetm(ibase)
         goto 1581
      end if 

      state1 = ncombv(1,ibase)
      state2 = ncombv(2,ibase)

*     Two fragments, both singlet functions
      if ( all_singlet ) then
         k = 0
         do i = 1, idetm(state1)
            do j = 1, idetm(state2)
               prod = civm(i,state1) * civm(j,state2)
               if ( abs(prod) .gt. thresh_CI ) then
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
     &                 perm_ab(iocch,nactb(ibase))
                  civb(k,ibase) = civb(k,ibase) * 
     &                 isetsign4(iocch,nactb(ibase))
               end if
            end do
         end do
         ndets = k
         goto 1581
      end if

*     nmol=2
      if(nmol.eq.2)then
         ms_states = max(spin(1)+1 , spin(2)+1)
         allocate(ndet_ms(ms_states))
         allocate(ci_ms(nci,ms_states))
         allocate(occ_ms(maxnact,nci,ms_states))
         ndet_ms = 0
         ci_ms   = 0.0
         occ_ms  = 0
         allocate(ndet1_ms(ms_states))
         allocate(ci1_ms(nci,ms_states))
         allocate(occ1_ms(maxnact,nci,ms_states))
         allocate(ndet2_ms(ms_states))
         allocate(ci2_ms(nci,ms_states))
         allocate(occ2_ms(maxnact,nci,ms_states))
         ndet1_ms = 0
         ci1_ms   = 0.0
         occ1_ms  = 0
         ndet2_ms = 0
         ci2_ms   = 0.0
         occ2_ms  = 0
*     M_S of fragment 1
         call sminop2(civm(1,ncombv(1,ibase)),ci_ms,
     &        ioccm(1,1,ncombv(1,ibase)),
     &        occ_ms(1:nactm(ncombv(1,ibase)),:,:),
     &        idetm(ncombv(1:,ibase)),ndet_ms,
     &        spin(1),nactm(ncombv(1,ibase)),nci,lfnout)
         do j = 1, spin(1) + 1
            ndet1_ms(j) = ndet_ms(j)
            do det = 1, ndet_ms(j)
               ci1_ms(det,j) = ci_ms(det,j)
               occ1_ms(:,det,j) = occ_ms(:,det,j)
            end do
         end do
*     M_S of fragment 2
         ndet_ms = 0
         ci_ms   = 0  
         occ_ms  = 0
         call sminop2(civm(1,ncombv(2,ibase)),ci_ms,
     &        ioccm(1,1,ncombv(2,ibase)),
     &        occ_ms(1:nactm(ncombv(2,ibase)),:,:),
     &        idetm(ncombv(2:,ibase)),ndet_ms,
     &        spin(2),nactm(ncombv(2,ibase)),nci,lfnout)
         do j = 1, spin(2) + 1
            ndet2_ms(j) = ndet_ms(j)
            do det = 1, ndet_ms(j)
               ci2_ms(det,j) = ci_ms(det,j)
               occ2_ms(:,det,j) = occ_ms(:,det,j)
            end do
         end do
         deallocate(ndet_ms,ci_ms,occ_ms)         
         k = 0
         do i = 1, spin(1) + 1
            ms1  = spin(1) - 2*(i-1)
            do j = 1, spin(2) + 1
               ms2 = spin(2) - 2*(j-1)
               if ( ms1 + ms2 .eq. nspin ) then
                  do det1 = 1, ndet1_ms(i)
                     do det2 = 1, ndet2_ms(j)
                        dnorm=clebsch_gordon(spin(1)/2.d0,spin(2)/2.d0,
     &                       ms1/2.d0,ms2/2.d0,
     &                       nspin/2.d0,nspin/2.d0)
                        prod = ci1_ms(det1,i)*ci2_ms(det2,j)
                        if (abs(prod) .gt. thresh_CI) then
                           iocch = 0
                           k = k + 1
                           civb(k,ibase) = prod * dnorm
                           kk = 0
                           do iact = 1, nactm(state1)
                              kk = kk + 1
                              iocch(kk) = occ1_ms(iact,det1,i)
                           end do
                           do jact = 1, nactm(state2)
                              kk = kk + 1
                              iocch(kk) = occ2_ms(jact,det2,j)
                           end do
                           call pack(ioccb(k,ibase),iocch,nactb(ibase))
                           civb(k,ibase) = civb(k,ibase) * 
     &                          perm_ab(iocch,nactb(ibase))
                           civb(k,ibase) = civb(k,ibase) * 
     &                          isetsign4(iocch,nactb(ibase))
                        end if
                     end do
                  end do
               end if
            end do
         end do
         deallocate(ndet1_ms,occ1_ms)            
         deallocate(ndet2_ms,occ2_ms)            
         deallocate(ci1_ms)
         deallocate(ci2_ms)
         ndets = k
 1581    dnorm = 0.0
         idetb(ibase) = ndets
         do i = 1, idetb(ibase)
            dnorm = dnorm + civb(i,ibase)**2
         end do
         dnorm = 1/dsqrt(dnorm)
         do i = 1, idetb(ibase)
            civb(i,ibase) = civb(i,ibase) * dnorm
         end do         
      endif
      
*     3 or more molecules
      if(nmol.ge.3)then
         allocate(spin_comb(2))
         spin_comb(1)=spin(1)
         spin_comb(2)=spin(2)
         inter_maxnact=maxnact
         write(lfnout,*)'INTER MAXNACT',inter_maxnact
         
         write(lfnout,*)'initial nci',nci
         write(lfnout,*)'maxci',maxci
         write(lfnout,*)'maxcib',maxcib
         ms_states = max(spin_comb(1)+1,spin_comb(2)+1)
         allocate(ndet_ms(ms_states))
         allocate(ci_ms(nci,ms_states))
         allocate(occ_ms(inter_maxnact,nci,ms_states))
         ndet_ms = 0
         ci_ms   = 0.0
         occ_ms  = 0
         allocate(ndet1_ms(ms_states))
         allocate(ci1_ms(nci,ms_states))
         allocate(occ1_ms(inter_maxnact,nci,ms_states))
         allocate(ndet2_ms(ms_states))
         allocate(ci2_ms(nci,ms_states))
         allocate(occ2_ms(inter_maxnact,nci,ms_states))
         ndet1_ms = 0
         ci1_ms   = 0.0
         occ1_ms  = 0
         ndet2_ms = 0
         ci2_ms   = 0.0
         occ2_ms  = 0
         allocate(inter_civm(maxcib,2))         
         allocate(inter_ioccm(inter_maxnact,maxcib,2))
         allocate(inter_civb(maxcib))         
         allocate(inter_ioccb(maxcib))
         allocate(inter_nactm(2))
         allocate(inter_idetm(2))
         
         inter_civm=0.0
         inter_nactm=0         
         do i=1,maxci
            inter_civm(i,1)=civm(i,ncombv(1,ibase))
            inter_civm(i,2)=civm(i,ncombv(2,ibase))
            do j=1,inter_maxnact
               inter_ioccm(j,1,1)=ioccm(j,i,ncombv(1,ibase))
               inter_ioccm(j,1,2)=ioccm(j,i,ncombv(2,ibase))
            enddo
         enddo
         inter_nactm(1)=nactm(ncombv(1,ibase))
         inter_nactm(2)=nactm(ncombv(2,ibase))
         inter_idetm(1)=idetm(ncombv(1,ibase))
         inter_idetm(2)=idetm(ncombv(2,ibase))
         
         inter_nactb=inter_nactm(1)+inter_nactm(2)
         
         do icoupling=1,nmol-1
            write(lfnout,*)
            write(lfnout,*)
            write(lfnout,*)
            write(lfnout,*)'-------COUPLING-----',icoupling
            write(lfnout,*)'spin comb'
            write(lfnout,*)spin_comb
            write(lfnout,*)'N ACTIVE ORBITALS 1 AND 2'
            write(lfnout,*)inter_nactm(1),inter_nactm(2)
            flush(lfnout)
*     M_S of fragment 1            
            call sminop2(inter_civm(1,1),ci_ms,
     &           inter_ioccm(1,1,1),occ_ms(1:inter_nactm(1),:,:),
     &           inter_idetm(1),ndet_ms,
     &           spin_comb(1),inter_nactm(1),nci,lfnout)
            write(lfnout,*)'***ndet1_ms***'
            write(lfnout,*)ndet_ms
            write(lfnout,*)'***occ1_ms***'
            do i=1,spin_comb(1)+1
               write(lfnout,*)'state',i
               do j=1,ndet_ms(i)           
                  write(lfnout,*)(occ_ms(k,j,i),k=1,inter_nactm(1))
               enddo
            enddo
            do j = 1, spin_comb(1) + 1
               ndet1_ms(j) = ndet_ms(j)
               do det = 1, ndet1_ms(j)
                  ci1_ms(det,j) = ci_ms(det,j)
                  do k=1,inter_nactm(1)
                     occ1_ms(k,det,j) = occ_ms(k,det,j)
                  enddo
               end do
            end do            
            write(lfnout,*)'---CONFIRM SMINOP 1'
            do j=1,spin_comb(1)+1
               write(lfnout,*)
               write(lfnout,*)'MS STATE CONFIRM',j
               write(lfnout,*)
               do k=1,ndet1_ms(j)
                  write(lfnout,*)ci1_ms(k,j),
     &                 (occ1_ms(l,k,j),l=1,inter_nactm(1))
               enddo
            enddo
*     M_S of fragment 2
            ndet_ms = 0
            ci_ms   = 0
            occ_ms  = 0            

            call sminop2(inter_civm(1,2),ci_ms,
     &           inter_ioccm(1,1,2),occ_ms(1:inter_nactm(2),:,:),
     &           inter_idetm(2),ndet_ms,
     &           spin_comb(2),inter_nactm(2),nci,lfnout)
            write(lfnout,*)'***ndet2_ms***'
            write(lfnout,*)ndet_ms            
            write(lfnout,*)'***occ2_ms***'
            do i=1,spin_comb(2)+1
               write(lfnout,*)'state',i
               do j=1,ndet_ms(i)      
                  write(lfnout,*)(occ_ms(k,j,i),k=1,inter_nactm(2))
               enddo
            enddo            
            do j = 1, spin_comb(2) + 1
               ndet2_ms(j) = ndet_ms(j)
               do det = 1, ndet2_ms(j)
                  ci2_ms(det,j) = ci_ms(det,j)
                  do k=1,inter_nactm(2)
                     occ2_ms(k,det,j) = occ_ms(k,det,j)
                  enddo
               end do                              
            end do
            write(lfnout,*)'---CONFIRM SMINOP 2'
            do j=1,spin_comb(1)+1
               write(lfnout,*)
               write(lfnout,*)'MS STATE CONFIRM',j
               write(lfnout,*)
               do k=1,ndet2_ms(j)
                  write(lfnout,*)ci2_ms(k,j),
     &                 (occ2_ms(l,k,j),l=1,inter_nactm(2))
               enddo
            enddo
                                  
            k=0

            !recalculate maxcib for icoupling>1
            if (icoupling>1)then
               do i=1, spin_comb(1)+1
                  ms1  = spin_comb(1) - 2*(i-1)               
                  do j = 1, spin_comb(2) + 1
                     ms2 = spin_comb(2) - 2*(j-1)                  
                     if (ms1+ms2.eq.
     &                    inter_couplings(icoupling,ibase)) then
                        do det1 = 1, ndet1_ms(i)
                           do det2 = 1, ndet2_ms(j)
                              prod = ci1_ms(det1,i)*ci2_ms(det2,j)
                              if (abs(prod) .gt. thresh_CI) k=k+1
                           enddo
                        enddo
                     endif
                  enddo
               enddo            
               if(k .eq. 0) then
                  write(lfnout,1574) '   spin(1) = ',spin(1)+1
                  write(lfnout,1574) '   spin(2) = ',spin(2)+1
                  write(lfnout,1574) 'total spin = ',
     &                 inter_couplings(icoupling,ibase)+1
                  write(lfnout,*) 'not possible'
                  flush(lfnout)
                  stop
 1574             format(a,i4)
               endif               
               maxcib=max(k,maxcib)
               write(lfnout,*)'New maxcib',maxcib
               deallocate(inter_civb)
               deallocate(inter_ioccb)               
               allocate(inter_civb(maxcib))
               allocate(inter_ioccb(maxcib))
            endif            
            
            k = 0
            
            do i = 1, spin_comb(1) + 1
               ms1  = spin_comb(1) - 2*(i-1)               
               do j = 1, spin_comb(2) + 1
                  ms2 = spin_comb(2) - 2*(j-1)                  
                  if (ms1+ms2.eq.inter_couplings(icoupling,ibase)) then
                     do det1 = 1, ndet1_ms(i)
                        do det2 = 1, ndet2_ms(j)
                           dnorm = clebsch_gordon
     &                          (spin_comb(1)/2.d0,spin_comb(2)/2.d0,
     &                          ms1/2.d0,ms2/2.d0,
     &                          inter_couplings(icoupling,ibase)/2.d0,
     &                          inter_couplings(icoupling,ibase)/2.d0)
                           prod = ci1_ms(det1,i)*ci2_ms(det2,j)
                           if (abs(prod) .gt. thresh_CI) then
                              iocch = 0
                              k = k + 1
                              inter_civb(k) = prod * dnorm
                              kk = 0
                              write(lfnout,*)'act 1'
                              write(lfnout,*)occ1_ms(:,det1,i)
                              write(lfnout,*)'act 2'
                              write(lfnout,*)occ2_ms(:,det2,j)                              
                              do iact = 1, inter_nactm(1)
                                 kk = kk + 1
                                 iocch(kk) = occ1_ms(iact,det1,i)
                              end do
                              
                              do jact = 1, inter_nactm(2)
                                 kk = kk + 1
                                 iocch(kk) = occ2_ms(jact,det2,j)
                              end do
                              
                              call pack(inter_ioccb(k),iocch,
     &                             inter_nactb)
                              if(icoupling.eq.(nmol-1))then
                                 write(lfnout,*)'Permute final MEBF',k
                                 inter_civb(k) = inter_civb(k) * 
     &                                perm_ab(iocch,inter_nactb)
                                 inter_civb(k) = inter_civb(k) * 
     &                                isetsign4(iocch,inter_nactb)
                                 flush(lfnout)
                              endif
                           end if
                        end do
                     end do
                  end if
               end do               
            end do
                      
            write(lfnout,*)'deallocate and set the next iteration'
            flush(lfnout)
            deallocate(ndet_ms)
            deallocate(ci_ms)
            deallocate(occ_ms)
            deallocate(ndet1_ms)
            deallocate(occ1_ms)            
            deallocate(ndet2_ms)
            deallocate(occ2_ms)            
            deallocate(ci1_ms)
            deallocate(ci2_ms)
            write(lfnout,*)'deallocated'
            flush(lfnout)
            ndets = k
 1582       dnorm = 0.0
            idetb(ibase) = ndets
            do i = 1, idetb(ibase)
               dnorm = dnorm + inter_civb(i)**2
            end do
            dnorm = 1/dsqrt(dnorm)
            do i = 1, idetb(ibase)
               inter_civb(i) = inter_civb(i) * dnorm
            end do                        
            if(icoupling.lt.(nmol-1))then
*     update values for the next coupling
               write(lfnout,*)'update values for the next coupling'
               flush(lfnout)
               spin_comb(1)=inter_couplings(icoupling,ibase)
               spin_comb(2)=spin(icoupling+2)
               ms_states = max(spin_comb(1)+1 , spin_comb(2)+1) 
               write(lfnout,*)'new ms_states'
               write(lfnout,*)ms_states               
               inter_civm=0.0
               inter_ioccm=0
               
               do i=1,ndets
                  inter_civm(i,1)=inter_civb(i)
               enddo
               do i=1,maxci
                  inter_civm(i,2)=civm(i,ncombv(icoupling+2,ibase))
               enddo
               write(lfnout,*)'INTER NACTB'               
               inter_idetm(1)=ndets
               inter_idetm(2)=idetm(ncombv(icoupling+2,ibase))
               inter_nactm(1)=inter_nactb
               inter_maxnact=inter_nactb
               

               deallocate(inter_ioccm)               
               allocate(inter_ioccm(inter_maxnact,maxcib,2))
               inter_ioccm=0
               do i=1,ndets
                  call unpack(inter_ioccb(i),inter_ioccm(:,i,1),
     &                 inter_nactb)                  
               enddo
               call gronor_determine_nci
     &              (2,inter_nactm,inter_idetm,inter_ioccm)
               write(lfnout,*)'NEW NCI',nci 
               
               do i=1, ndets                 
                  write(lfnout,*)inter_civb(i),
     &                 (inter_ioccm(j,i,1),j=1,inter_nactb)
               enddo

               flush(lfnout)
               inter_nactm(2)=nactm(ncombv(icoupling+2,ibase))
               inter_nactb=inter_nactm(1)+inter_nactm(2)
               
               do i=1,maxci
                  do j=1,maxnact
                     inter_ioccm(j,i,2)=
     &                    ioccm(j,i,ncombv(icoupling+2,ibase))
                  enddo
               enddo                
              
               write(lfnout,*)'Reallocate arrays'
               flush(lfnout)               
               allocate(ndet_ms(ms_states))
               allocate(ci_ms(nci,ms_states))
               allocate(occ_ms(inter_maxnact,nci,ms_states))
               allocate(ndet1_ms(ms_states))
               allocate(ci1_ms(nci,ms_states))
               allocate(occ1_ms(inter_maxnact,nci,ms_states))
               allocate(ndet2_ms(ms_states))
               allocate(ci2_ms(nci,ms_states))
               allocate(occ2_ms(maxnact,nci,ms_states))               

               
               ndet_ms = 0
               ci_ms   = 0.0
               occ_ms  = 0
               ndet1_ms = 0
               ci1_ms   = 0.0
               occ1_ms  = 0
               ndet2_ms = 0
               ci2_ms   = 0.0
               occ2_ms  = 0
            else               
               if(maxcib>maxcib0)then ! we need to re-allocate civb and ioccb
!    without deleting the previous MEBFs
                  do i=1,ibase-1                     
                     do j=1,idetb(ibase-1)
                        civb_copy(j,i)=civb(j,i)
                        ioccb_copy(j,i)=ioccb(j,i)
                     enddo
                  enddo
                  deallocate(civb,ioccb)
                  allocate(civb(maxcib,nbase))
                  allocate(ioccb(maxcib,nbase))
                  civb=0.0
                  ioccb=0
                  do i=1,ibase-1                     
                     do j=1,idetb(ibase-1)
                        civb(j,i)=civb_copy(j,i)
                        ioccb(j,i)=ioccb_copy(j,i)
                     enddo
                  enddo
               endif
               do i=1,maxcib
                  civb(i,ibase)=inter_civb(i)
                  ioccb(i,ibase)=inter_ioccb(i)
               enddo
            endif                        
         enddo        
         deallocate(spin_comb)
         deallocate(inter_civm)
         deallocate(inter_ioccm)
         deallocate(inter_idetm)
         deallocate(inter_nactm)
         deallocate(inter_civb)
         deallocate(inter_ioccb)
      endif
      deallocate(spin)
      return
      end 

      subroutine gronor_makebasestate_remove                                     
     &     (lfnout,civm,vecsm,ioccm,nbasm,nactm,                            
     &     inactm,idetm,maxci,maxvec,maxnact,mstates,                     
     &     civb,vecsc,ioccb,inactb,nactb,idetb,                           
     &     maxcib,nbasis,mactb,                        
     &     nmol,ncomb,nspin,thresh_CI)

      implicit none
      integer :: maxact
      parameter (maxact=31)
      integer :: nspin,nmol,mstates,maxci,maxvec,maxnact
      integer :: inactb,nactb,idetb,maxcib,nbasis,mactb
      integer :: ioccm(maxnact,maxci,mstates)
      integer :: nbasm(mstates),nactm(mstates),inactm(mstates)
      integer :: idetm(mstates)
      integer :: ioccb(maxcib),ncomb(nmol)
      integer :: nmstates,nci,mhighs,jjstate,is,iprod,ivecs
      integer :: nspintot,ioff,idetoutoutd,llspin,kocc,lfnout,ioff2
      integer :: addspin,spin_sum,ns,opp_jjstate
      integer, allocatable :: intocc(:,:,:),intoccout(:,:,:),           
     &     jspin(:),idetout(:),kspin(:),jcomb(:),kcomb(:),
     &     nactms(:),jjspin(:),spin_counter(:),idetoutout(:),
     &     multiplicity(:) 
      integer :: i,j,k,l,m,jj,kk,ll,kstart,jstart
      real (kind=8) :: dnorm,thresh_CI
      real (kind=8) :: civb(maxcib),vecsc(nbasis,nbasis)
      real (kind=8),allocatable :: ciint(:,:),ciout(:,:)
      real (kind=8) :: civm(maxci,mstates),vecsm(maxvec,maxvec,mstates)
      real (kind=8), external :: ddot,clebsch_gordon
      logical :: oalls,onew,ocomparocc

      inactb=0
      nactb=0
      idetb=1
      do i=1,nmol
         inactb=inactb+inactm(ncomb(i))
         nactb=nactb+nactm(ncomb(i))
         idetb=idetb*idetm(ncomb(i))
      enddo

      if(nactb.gt.maxact) stop 'too many active orbitals > maxact'
      if(idetb.gt.maxcib) stop 'too many determinants > maxcib'

      do i=1,nbasis
         do j=1,nbasis
            vecsc(i,j)=0.0d0
         enddo
      enddo
      jstart=1
      kstart=1
      do i=1,nmol
         do j=1,inactm(ncomb(i))
            do k=1,nbasm(ncomb(i))
               vecsc(kstart+k-1,jstart+j-1)=vecsm(k,j,ncomb(i))
            enddo
         enddo
         kstart=kstart+nbasm(ncomb(i))
         jstart=jstart+inactm(ncomb(i))
      enddo
      
      jstart=inactb+1
      kstart=1
      
      do i=1,nmol
         do j=1,nactm(ncomb(i))
            do k=1,nbasm(ncomb(i))
               vecsc(kstart+k-1,jstart+j-1)=vecsm(k,j+inactm(ncomb(i)),
     &              ncomb(i))
            enddo
         enddo
         kstart=kstart+nbasm(ncomb(i))
         jstart=jstart+nactm(ncomb(i))
      enddo

      allocate(jspin(nmol))
      
      do i=1,nmol
        jspin(i)=0
        do j=1,nactm(ncomb(i))
          if(abs((ioccm(j,1,ncomb(i)))).eq.1)
     &         jspin(i)=jspin(i)+ioccm(j,1,ncomb(i))
        enddo
      enddo

      oalls=.true.
      
      do i=2,nmol
         if(jspin(i).ne.jspin(1))oalls=.false.
      enddo
      
      if(jspin(1).ne.0)oalls=.false.
      
!     Singlet-singlet couple
      
      if(oalls.or.nmol.eq.1) then
        call productf(civm,ioccm,nactm,idetm,nmol,ncomb,civb,ioccb,
     &       maxci,mstates,maxnact,mactb,maxcib,idetb,thresh_CI,lfnout)
      else
        nmstates=0
        nci=0
        mhighs=0
        spin_sum=0
        do i=1,nmol
          nmstates=nmstates+jspin(i)+1
          nci=max(nci,idetm(ncomb(i)))
          mhighs=max(mhighs,jspin(i))
          spin_sum=spin_sum+jspin(i)+1
        enddo
        
!     Construction of ms states

        nci=(spin_sum+2)*nci
!     safe also for > 2 fragments 

!     nci=(mhighs+2)*nci
!     nci=int(sqrt(maxcib)/(2*mhighs))

        allocate(ciint(nci,nmstates),intocc(maxnact,nci,nmstates),
     &       idetout(nmstates),kspin(nmstates),nactms(nmstates))
        jjstate=0
        do i=1,nmol
          jjstate=jjstate+1
          do jj=1,idetm(ncomb(i))
            ciint(jj,jjstate)=civm(jj,ncomb(i)) ! first ms state from molcas -> civec
            do kk=1,nactm(ncomb(i))
              intocc(kk,jj,jjstate)=ioccm(kk,jj,ncomb(i))
            enddo
          enddo
          
          idetout(jjstate)=idetm(ncomb(i))
          kspin(jjstate)=jspin(i)
          nactms(jjstate)=nactm(ncomb(i))

          if(jspin(i).eq.0) goto 10
          
          if (mod(jspin(i)+1,2).eq.0)  then
!     two cases, odd or even multiplicity
            ns=((jspin(i)+1)/2)-1
!     multiplicity/2 - 1
          else
            ns=jspin(i)/2
!     (multiplicity-1)/2
          endif
          
          is=jspin(i)           ! ms spin of the current state

          do j=1,ns
!     do is=jspin(i)-2,-jspin(i)+2,-2 ! build all intermediate ms states with sminop
            jjstate=jjstate+1
            is=is-2
               call sminop(ciint(1,jjstate-1),intocc(1,1,jjstate-1),
     &              ciint(1,jjstate),intocc(1,1,jjstate),
     &              idetout(jjstate-1),nactm(ncomb(i)),
     &              idetout(jjstate),maxnact,lfnout)
               kspin(jjstate)=is ! ms spin of jjstate
               nactms(jjstate)=nactm(ncomb(i))
               if (idetout(jjstate).gt.nci)  stop 'error: nci too small'
               dnorm=ddot(idetout(jjstate),ciint(1,jjstate),1,
     &              ciint(1,jjstate),1)
               dnorm=1.0d0/dsqrt(dnorm)
               call dscal(idetout(jjstate),dnorm,ciint(1,jjstate),1)
            enddo
            
            if(mod(jspin(i)+1,2).eq.0)then
               opp_jjstate=jjstate
               ns=(jspin(i)+1)/2
            else
               opp_jjstate=jjstate-1
               ns=jspin(i)/2
            endif

            do j=1,ns
               jjstate=jjstate+1
               do jj=1,idetout(opp_jjstate)
                  ciint(jj,jjstate)=ciint(jj,opp_jjstate)
                  do kk=1,nactm(ncomb(i))
                     if(abs(intocc(kk,jj,opp_jjstate)).eq.1) then
                        intocc(kk,jj,jjstate)=-intocc(kk,jj,opp_jjstate)
                     else
                        intocc(kk,jj,jjstate)=intocc(kk,jj,opp_jjstate)
                     endif
                  enddo
               enddo
               idetout(jjstate)=idetout(opp_jjstate)
               kspin(jjstate)=-kspin(opp_jjstate)
               nactms(jjstate)=nactm(ncomb(i))
               opp_jjstate=opp_jjstate-1
            enddo

 10         continue
         enddo
         
!     Product of the states

         iprod=1
         allocate(multiplicity(nmol))
         do i=1,nmol
            iprod=iprod*(jspin(i)+1)
            multiplicity(i)=jspin(i)+1
         enddo

         allocate(jcomb(nmol),kcomb(nmol),jjspin(nmol),
     &        spin_counter(nmol))
         if(nci**2.gt.maxcib) then
           write(*,*)' Warning maxcib too small',maxcib,nci**2
           stop
         endif
         allocate(idetoutout(iprod))
         ivecs=0
         ioff2=1

         do i=1,iprod
            call calciloop(i,jcomb,nmol,multiplicity)
            nspintot=0
            ioff=0
            idetoutoutd=1
            do j=1,nmol
               nspintot=nspintot+kspin(ioff+jcomb(j))
               spin_counter(j)=nspintot
               jjspin(j)=kspin(ioff+jcomb(j))
               kcomb(j)=ioff+jcomb(j)
               ioff=ioff+multiplicity(j)
               idetoutoutd=idetoutoutd*idetout(kcomb(j))
            enddo
            if(nspintot.eq.nspin) then
               ivecs=ivecs+1
               if(nmol.le.2)then ! TWO MOLECULES CASE
                  call productf(ciint,intocc,nactms,idetout,nmol,kcomb,
     &                 civb(ioff2),ioccb(ioff2),nci,
     &                 nmstates,maxnact,mactb,maxcib-ioff2+1,
     &                 idetoutoutd,thresh_CI,lfnout)   
                  idetoutout(ivecs)=idetoutoutd
                  dnorm=clebsch_gordon(jspin(1)/2.0d0,jspin(2)/2.0d0,
     &                 jjspin(1)/2.0d0,jjspin(2)/2.0d0,
     &                 nspin/2.0d0,nspin/2.0d0)
                  call dscal(idetoutout(ivecs),dnorm,civb(ioff2),1)
                else
                  
!     Three or more molecules
                  write(*,*)'work in progress'
                  
               endif
               ioff2=ioff2+idetoutoutd
            endif
         enddo

         idetb=0
         llspin=ivecs-1
         ioff2=1
         do i=1,ivecs
            ioff2=ioff2+idetoutout(i)
         enddo
         idetb=ioff2-1
         deallocate(ciint)
         deallocate(nactms)
         deallocate(intocc)
         deallocate(idetout)
         deallocate(kspin)
         deallocate(jjspin)
         deallocate(spin_counter)
         deallocate(jcomb)
         deallocate(kcomb)
         deallocate(idetoutout)
         deallocate(multiplicity)
      endif
      deallocate(jspin)    
      return
      end
      
!> @brief
!!    Construct the MEBFs from the fragment functions	
!!
!!  @author
!!    Remco Havenith
!!
!!  @todo
!!    To be replaced with gronor_make_basestate
!!
      subroutine gronor_old_makebasestate
     &     (lfnout,civm,vecsm,ioccm,nbasm,nactm,
     &     inactm,idetm,maxci,maxvec,maxnact,mstates,
     &     civb,vecsc,ioccb,inactb,nactb,idetb,
     &     maxcib,nbasis,mactb,
     &     nmol,ncomb,nspin,thresh_CI)

      implicit none
      integer :: maxact
      parameter (maxact=31)
      integer :: nspin,nmol,mstates,maxci,maxvec,maxnact
      integer :: inactb,nactb,idetb,maxcib,nbasis,mactb
      real (kind=8) :: civm(maxci,mstates),vecsm(maxvec,maxvec,mstates)
      integer :: ioccm(maxnact,maxci,mstates)
      integer :: nbasm(mstates),nactm(mstates),inactm(mstates)
      integer :: idetm(mstates)
      real (kind=8) :: civb(maxcib),vecsc(nbasis,nbasis)
      integer :: ioccb(maxcib),ncomb(nmol)
      real (kind=8),allocatable :: ciint(:,:),ciout(:,:)
      integer, allocatable :: intocc(:,:,:),intoccout(:,:,:),           &
     & jspin(:),idetout(:),kspin(:),jcomb(:),kcomb(:),idetoutout(:),    &
     & nactms(:)
      real (kind=8) :: dnorm,thresh_CI
      real (kind=8), external :: ddot

      integer :: nmstates,nci,mhighs,jjstate,jj,kk,is,iprod,ivecs
      integer :: nspintot,ioff,idetoutoutd,llspin,kocc,lfnout,ioff2

      logical :: oalls,onew,ocomparocc

      real (kind=8) :: dsi(3,2)
      integer :: i,j,k, kstart, jstart
      data dsi/ 1.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0/

      inactb=0
      nactb=0
      idetb=1
      do i=1,nmol
       inactb=inactb+inactm(ncomb(i))
       nactb=nactb+nactm(ncomb(i))
       idetb=idetb*idetm(ncomb(i))
      enddo
      if(nactb.gt.maxact)stop 'too many active orbitals'
      if(idetb.gt.maxcib)then
       print *,idetb,maxcib
       stop 'you fool!!!'
      endif
      do i=1,nbasis
       do j=1,nbasis
        vecsc(i,j)=0.0d0
       enddo
      enddo
      jstart=1
      kstart=1

      do i=1,nmol
       do j=1,inactm(ncomb(i))
        do k=1,nbasm(ncomb(i))
         vecsc(kstart+k-1,jstart+j-1)=vecsm(k,j,ncomb(i))
        enddo
       enddo
       kstart=kstart+nbasm(ncomb(i))
       jstart=jstart+inactm(ncomb(i))
      enddo

      jstart=inactb+1
      kstart=1

      do i=1,nmol
       do j=1,nactm(ncomb(i))
        do k=1,nbasm(ncomb(i))
         vecsc(kstart+k-1,jstart+j-1)=vecsm(k,j+inactm(ncomb(i)),
     &                                                        ncomb(i))
        enddo
       enddo
       kstart=kstart+nbasm(ncomb(i))
       jstart=jstart+nactm(ncomb(i))
      enddo

      allocate(jspin(nmol))

      do i=1,nmol
       jspin(i)=0
       do j=1,nactm(ncomb(i))
        if(abs((ioccm(j,1,ncomb(i)))).eq.1) jspin(i)=jspin(i)+
     &                                           ioccm(j,1,ncomb(i))
       enddo
      enddo

      oalls=.true.

      do i=2,nmol
       if(jspin(i).ne.jspin(1))oalls=.false.
      enddo

      if(jspin(1).ne.0)oalls=.false.

!     Singlet-singlet couple

      if(oalls) then
        call old_productf(
     &       civm,ioccm,nactm,idetm,nmol,ncomb,civb,ioccb,maxci,
     &               mstates,maxnact,mactb,maxcib,idetb,thresh_CI)
      else
       nmstates=0
       nci=0
       mhighs=0
       do i=1,nmol
        nmstates=nmstates+jspin(i)+1
        nci=max(nci,idetm(ncomb(i)))
        mhighs=max(mhighs,jspin(i))
       enddo

!     Need better way to calculate nci

       nci=(mhighs+2)*nci
       if(mhighs.gt.2)stop 'program the rest yourself!'
       allocate(ciint(nci,nmstates),intocc(maxnact,nci,nmstates),
     &          idetout(nmstates),kspin(nmstates),nactms(nmstates))
       jjstate=0
       do i=1,nmol
        jjstate=jjstate+1
        do jj=1,idetm(ncomb(i))
         ciint(jj,jjstate)=civm(jj,ncomb(i))
         do kk=1,nactm(ncomb(i))
          intocc(kk,jj,jjstate)=ioccm(kk,jj,ncomb(i))
         enddo
        enddo
        idetout(jjstate)=idetm(ncomb(i))
        kspin(jjstate)=jspin(i)
        nactms(jjstate)=nactm(ncomb(i))
        if(jspin(i).eq.0) goto 10
        do is=jspin(i)-2,-jspin(i)+2,-2
         jjstate=jjstate+1
         call sminop(ciint(1,jjstate-1),intocc(1,1,jjstate-1),
     &               ciint(1,jjstate),intocc(1,1,jjstate),
     &               idetout(jjstate-1),nactm(ncomb(i)),
     &               idetout(jjstate),maxnact)
         kspin(jjstate)=is
         nactms(jjstate)=nactm(ncomb(i))
         if (idetout(jjstate).gt.nci)then
           print *,'nci::',nci,idetout(jjstate),jjstate
           stop 'nci error'
         endif

         dnorm=ddot(idetout(jjstate),ciint(1,jjstate),1,
     &              ciint(1,jjstate),1)
         dnorm=1.0d0/dsqrt(dnorm)
         call dscal(idetout(jjstate),dnorm,ciint(1,jjstate),1)
        enddo
        jjstate=jjstate+1
        do jj=1,idetm(ncomb(i))
         ciint(jj,jjstate)=civm(jj,ncomb(i))
         do kk=1,nactm(ncomb(i))
          if(abs(ioccm(kk,jj,ncomb(i))).eq.1) then
           intocc(kk,jj,jjstate)=-ioccm(kk,jj,ncomb(i))
          else
           intocc(kk,jj,jjstate)=ioccm(kk,jj,ncomb(i))
          endif
         enddo
         ciint(jj,jjstate)=civm(jj,ncomb(i))
        enddo
        idetout(jjstate)=idetm(ncomb(i))
        kspin(jjstate)=-jspin(i)
        nactms(jjstate)=nactm(ncomb(i))
   10   continue
       enddo

       iprod=1
       do i=1,nmol
        jspin(i)=jspin(i)+1
        iprod=iprod*(jspin(i))
       enddo
       allocate(jcomb(nmol),kcomb(nmol))
       if(nci**2.gt.maxcib) print*,'WARNING maxcib too small',
     &    maxcib,nci**2
         allocate(idetoutout(iprod))
       ivecs=0
       ioff2=1
       do i=1,iprod
        call calciloop(i,jcomb,nmol,jspin)
        nspintot=0
        ioff=0
        idetoutoutd=1
        do j=1,nmol
         nspintot=nspintot+kspin(ioff+jcomb(j))
         kcomb(j)=ioff+jcomb(j)
         ioff=ioff+jspin(j)
         idetoutoutd=idetoutoutd*idetout(kcomb(j))
        enddo
        if(nspintot.eq.nspin) then
         ivecs=ivecs+1
         call old_productf(
     &        ciint,intocc,nactms,idetout,nmol,kcomb,
     &                 civb(ioff2),ioccb(ioff2),nci,
     &                 nmstates,maxnact,mactb,maxcib-ioff2+1,
     &                 idetoutoutd,thresh_CI)

         idetoutout(ivecs)=idetoutoutd
         ioff2=ioff2+idetoutoutd
        endif
       enddo
       dnorm=1.0d0/(dsqrt(ivecs*1.0d0))
       idetb=0

!      TODO: Need to take care of the sign when forming the combinations!
!      Currently only trip-trip-singlet...

       llspin=ivecs-1
       ioff2=1
       do i=1,ivecs
        call dscal(idetoutout(i),dnorm*dsi(i,llspin),civb(ioff2),1)
        ioff2=ioff2+idetoutout(i)
       enddo
       idetb=ioff2-1

       deallocate(ciint)
       deallocate(nactms)
       deallocate(intocc)
       deallocate(idetout)
       deallocate(kspin)
       deallocate(jcomb)
       deallocate(kcomb)
       deallocate(idetoutout)
      endif

      deallocate (jspin)

      return
      end
