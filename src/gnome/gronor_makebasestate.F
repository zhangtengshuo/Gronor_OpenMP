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
      use cidist
      use gnome_data
      use gnome_parameters
      implicit none
      

      integer              :: ibase,state1,state2
      integer              :: i,j,k,kk,iact,jact
      integer              :: jstart,kstart
      integer              :: ndets,istate
      integer              :: ms1,ms2,det1,det2
      integer, external    :: perm_ab,isetsign4
      
!     Because of the packing, maxact is limited to 32 (I guess)
      integer, parameter   :: maxact = 32
      
      integer, allocatable :: spin(:)
      integer, allocatable :: ndet1_ms(:),ndet2_ms(:)
      integer, allocatable :: occ1_ms(:,:,:),occ2_ms(:,:,:)
      integer, allocatable :: iocc_tmp(:,:)

      logical              :: all_singlet

      real (kind = 8)              :: prod,dnorm
      real (kind = 8), external    :: clebsch_gordon
      real (kind = 8), allocatable :: ci1_ms(:,:)
      real (kind = 8), allocatable :: ci2_ms(:,:)
      real (kind = 8), allocatable :: cicoef_tmp(:)
      
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
         call errquit(me,9,"Too many active orbitals")
      endif

* Combining the vectors of the fragments

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

* Just one molecule
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

* Two singlet fragments
      if ( all_singlet ) then
        k = 0
        do i = 1, idetm(state1)
          do j = 1, idetm(state2)
            prod = civm(i,state1) * civm(j,state2)
            if ( abs(prod) .gt. tau_CI) then
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
     &                                   perm_ab(iocch,nactb(ibase))
              civb(k,ibase) = civb(k,ibase) * 
     &                                   isetsign4(iocch,nactb(ibase))
            end if
          end do
        end do
        ndets = k
        goto 1581
      end if

* All other cases
      allocate(ndet1_ms(spin(1)+1))
      allocate(ci1_ms(nci,spin(1)+1))
      allocate(occ1_ms(maxnact,nci,spin(1)+1))
      allocate(ndet2_ms(spin(2)+1))
      allocate(ci2_ms(nci,spin(2)+1))
      allocate(occ2_ms(maxnact,nci,spin(2)+1))
      ndet1_ms = 0
      ci1_ms   = 0.0
      occ1_ms  = 0
      ndet2_ms = 0
      ci2_ms   = 0.0
      occ2_ms  = 0
*   M_S of fragment 1
      istate=ncombv(1,ibase)
      allocate(iocc_tmp(nactm(istate),idetm(istate)))
      allocate(cicoef_tmp(idetm(istate)))
      do j = 1, idetm(istate)
        cicoef_tmp(j) = civm(j,istate)
        do k = 1, nactm(istate)
          iocc_tmp(k,j) = ioccm(k,j,istate)
        end do
      end do
      call sminop(cicoef_tmp,ci1_ms,iocc_tmp,occ1_ms,idetm(istate),
     &             ndet1_ms,spin(1),nactm(istate),nci,maxnact) 
      deallocate(cicoef_tmp)
      deallocate(iocc_tmp)
*   M_S of fragment 2
      istate=ncombv(2,ibase)
      allocate(iocc_tmp(nactm(istate),idetm(istate)))
      allocate(cicoef_tmp(idetm(istate)))
      do j = 1, idetm(istate)
        cicoef_tmp(j) = civm(j,istate)
        do k = 1, nactm(istate)
          iocc_tmp(k,j) = ioccm(k,j,istate)
        end do
      end do
      call sminop(cicoef_tmp,ci2_ms,iocc_tmp,occ2_ms,idetm(istate),
     &             ndet2_ms,spin(2),nactm(istate),nci,maxnact)
      deallocate(cicoef_tmp)
      deallocate(iocc_tmp)

      
      k = 0
      do i = 1, spin(1) + 1
        ms1  = spin(1) - 2*(i-1)
        do j = 1, spin(2) + 1
          ms2 = spin(2) - 2*(j-1)
          if ( ms1 + ms2 .eq. nspin ) then
            do det1 = 1, ndet1_ms(i)
              do det2 = 1, ndet2_ms(j)
                dnorm = clebsch_gordon(spin(1)/2.d0,spin(2)/2.d0,
     &                      ms1/2.d0,ms2/2.d0,
     &                      nspin/2.d0,nspin/2.d0)
                prod = ci1_ms(det1,i)*ci2_ms(det2,j)
                if (abs(prod) .gt. tau_CI) then
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
     &                              perm_ab(iocch,nactb(ibase))
                  civb(k,ibase) = civb(k,ibase) * 
     &                              isetsign4(iocch,nactb(ibase))
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


 1581 dnorm = 0.0
      idetb(ibase) = ndets
      do i = 1, idetb(ibase)
        dnorm = dnorm + civb(i,ibase)**2
      end do
      dnorm = 1/dsqrt(dnorm)
      do i = 1, idetb(ibase)
        civb(i,ibase) = civb(i,ibase) * dnorm
      end do
      deallocate(spin)
      return
      end 

