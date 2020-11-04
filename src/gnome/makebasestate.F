      subroutine makebasestate(ibase)
      use cidef
      use gnome_data
      use gnome_parameters
      implicit none
      
      integer, parameter   :: maxact=40
      integer, allocatable :: spin(:)
      integer, allocatable :: mspin(:),mscomb(:),mul_comb(:) ! ms spin and combination for the products
      integer, allocatable :: ndet_mspin(:),occ_mspin(:,:,:) ! mspin states from sminop
      integer, allocatable :: multiplicity(:),ndetprod(:),nact_mspin(:)
      integer              :: i,j,k,ibase
      integer              :: jstart,kstart
!      integer              :: nci ! max ci coefficients in input and number of monomer state      
      integer              :: mspin_states,msum,nproducts
      integer              :: ndetout,ivecs
      integer              :: offsetm,offsetdet ! offset monomer and offset determinants in the product
      
      real (kind=8),allocatable :: ci_mspin(:,:) ! mspin ci coeff from sminop      
      real (kind=8)        :: dnorm
      
      logical              :: all_singlet

      integer,external      :: perm_ab
      real(kind=8),external :: cg     

      inactb(ibase) = 0
      nactb(ibase)  = 0
      idetb(ibase)  = 1
      do i = 1, nmol
         inactb(ibase) = inactb(ibase) + inactm(ncombv(i,ibase))
         nactb(ibase)  = nactb(ibase) + nactm(ncombv(i,ibase))
         idetb(ibase)  = idetb(ibase) * idetm(ncombv(i,ibase))
      end do
      
      if ( nactb(ibase) .gt. maxact)then
         write(lfnout,*)'Error nact',nactb(ibase),maxact
         stop 'too many active orbitals'
      endif
      if ( idetb(ibase) .gt. maxcib)then
         write(lfnout,*)'Error ndet',idetb(ibase), maxcib
         stop 'too many determinants'
      endif
*     combining the vectors of the fragments
      do i=1,nbasis
         do j=1,nbasis
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
      spin = 0
      do i = 1, nmol
         do j = 1, nactm(ncombv(i,ibase))
            if ( abs((ioccm(j,1,ncombv(i,ibase)))) .eq. 1 ) then
               spin(i) = spin(i) + ioccm(j,1,ncombv(i,ibase))
            endif
         end do
         if ( all_singlet ) then
            if ( spin(i) .ne. 0 ) all_singlet = .false.
         end if
      end do
      

      k=1
      do i=1,nmol
         nci=idetm(ncombv(i,ibase))*(2**spin(i)) 
         if(nci.lt.k)then
            nci=k
         else                                           
            k=nci
         endif
      enddo
      do i=1,maxcib
         civb(i,ibase)=0.0d0
         ioccb(i,ibase)=0
      enddo
      
      if ( all_singlet .or. ( nmol .eq. 1 ) ) then         
         call productf2(civm,ioccm,nactm,idetm,nmol,ncombv(:,ibase),
     &        civb(:,ibase),ioccb(:,ibase),maxci,mstates,maxnact,
     &        mactb,maxcib,idetb(ibase),thresh_CI,lfnout)                              
      else
         mspin_states = 0
         do i = 1, nmol
            mspin_states = mspin_states + spin(i) + 1
         enddo
         allocate(mscomb(nmol))
         allocate(mul_comb(nmol))
         allocate(mspin(mspin_states))
         allocate(nact_mspin(mspin_states))
         allocate(ndet_mspin(mspin_states))
         allocate(ci_mspin(nci,mspin_states))
         allocate(occ_mspin(maxnact,nci,mspin_states))
         allocate(multiplicity(mspin_states))
         
         ci_mspin=0.0d0
         occ_mspin=0
         multiplicity=0

         offsetm=0

         do i=1,nmol                        
            call sminop2
     &           (civm(1,ncombv(i,ibase)),ci_mspin(1,offsetm+1),
     &           ioccm(1,1,ncombv(i,ibase)),occ_mspin(1,1,offsetm+1),
     &           idetm(ncombv(i,ibase)),ndet_mspin(offsetm+1),
     &           spin(i),nactm(ncombv(i,ibase)),nci,lfnout)
            
            flush(lfnout)
            mspin(offsetm+1)=spin(i)
            do j=2,spin(i)+1
               mspin(offsetm+j)=mspin(offsetm+j-1)-2
            enddo
            offsetm=offsetm+spin(i)+1                        
         enddo

         k=0
         do i=1,nmol
            do j=1,spin(i)+1
               k=k+1
               nact_mspin(k)=nactm(ncombv(i,ibase))
            enddo
         enddo
*     Product of the monomers
         ndetout=0         
         nproducts=1         
         do i=1,nmol
            multiplicity(i)=spin(i)+1
            nproducts=nproducts*(spin(i)+1)            
         enddo         
         allocate(ndetprod(nproducts))
         ndetprod=0
         ivecs=0
         offsetdet=1
         do i=1,nproducts
            call calciloop(i,mscomb,nmol,multiplicity)        
            msum=0
            offsetm=0
            do j=1,nmol
               msum=msum+mspin(mscomb(j)+offsetm)
               mul_comb(j)=offsetm+mscomb(j)
               offsetm=offsetm+multiplicity(j)       
            enddo
            
            
            if(msum.eq.nspin)then
               ivecs=ivecs+1
               if(nmol.le.2)then                                    
                  call productf2(ci_mspin,occ_mspin,nact_mspin,
     &                 ndet_mspin,nmol,mul_comb,civb(offsetdet,ibase),
     &                 ioccb(offsetdet,ibase),nci,mspin_states,
     &                 maxnact,mactb,maxcib-offsetdet+1,ndetout,
     &                 thresh_CI,lfnout)
                  
                  ndetprod(ivecs)=ndetout
                  dnorm=cg(spin(1)/2.0d0,spin(2)/2.0d0,
     &                 mspin(mscomb(1))/2.0d0,
     &                 mspin(multiplicity(1)+mscomb(2))/2.0d0,
     &                 nspin/2.0d0,nspin/2.0d0)
                  call dscal(ndetprod(ivecs),dnorm,civb(offsetdet,
     &                 ibase),1)                                    
               else
                  stop 'Three or more molecules not supported' 
!     TODO :  select intermediate couplings by user input
               endif
               offsetdet=offsetdet+ndetout
            endif
         enddo
         dnorm=0.0d0
         idetb(ibase)=0
         do i=1,ivecs
            idetb(ibase)=idetb(ibase)+ndetprod(i)
         enddo
         do i=1,idetb(ibase)
            dnorm=dnorm+civb(i,ibase)**2
         enddo
         do i=1,idetb(ibase)
            civb(i,ibase)=civb(i,ibase)*1.0d0/dsqrt(dnorm)
         enddo
         
         deallocate(mspin)
         deallocate(nact_mspin)
         deallocate(ndet_mspin)
         deallocate(ci_mspin)
         deallocate(mscomb)
         deallocate(mul_comb)
         deallocate(occ_mspin)
         deallocate(multiplicity)
         deallocate(ndetprod)
      endif
      deallocate(spin)            
      return
      end
