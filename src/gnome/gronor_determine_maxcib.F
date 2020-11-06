      subroutine gronor_determine_maxcib
      use cidef
      use gnome_parameters
      implicit none

      integer              :: istate,iact,ibase,i,j,k
      integer              :: det1,det2,cib
      integer              :: ms1,ms2,offsetm,mspin_states,spin(2)
      integer,allocatable  :: mspin(:),ndet_mspin(:,:),occ_mspin(:,:,:)

      real (kind = 8)               :: cg,prod,dnorm
      real (kind = 8), allocatable  :: ci_mspin(:,:)

      allocate( spinm(mstates) )
      spinm = 0
      do istate = 1, mstates
        do iact = 1, nactm(istate)
          if ( abs(ioccm(iact,1,istate)) .eq. 1 ) then
            spinm(istate) = spinm(istate) + ioccm(iact,1,istate)
          end if
        end do
      end do

      maxcib = 0

      if ( nmol .eq. 1 ) then
        do i = 1, mstates
          maxcib = max(maxcib,idetm(i))
        end do
        return
      end if

      do ibase = 1, nbase
        cib = 0
        if ( ( spinm(ncombv(1,ibase)) .eq. 0 ) .and.
     &             ( spinm(ncombv(2,ibase)) .eq. 0 ) ) then
          do det1 = 1, idetm(ncombv(1,ibase))
            do det2 = 1, idetm(ncombv(2,ibase))
              prod = civm(det1,ncombv(1,ibase)) *
     &                            civm(det2,ncombv(2,ibase))
              if ( abs(prod) .gt. thresh_CI ) cib = cib + 1
            end do
          end do
        else
          mspin_states = spinm(ncombv(1,ibase))+1 +
     &                        spinm(ncombv(2,ibase))+1
          allocate(mspin(mspin_states))                       ! ms values
          allocate(ndet_mspin(2,mspin_states))
          allocate(ci_mspin(nci,mspin_states))
          allocate(occ_mspin(maxnact,nci,mspin_states))       ! not really needed here
          ndet_mspin = 0
          offsetm = 0
          do i = 1, nmol 
            spin(i) = spinm(ncombv(i,ibase))
            call sminop2(civm(1,ncombv(i,ibase)),ci_mspin(1,offsetm+1),
     &           ioccm(1,1,ncombv(i,ibase)),occ_mspin(1,1,offsetm+1),
     &           idetm(ncombv(i:,ibase)),ndet_mspin(i,:),
     &           spin(i),nactm(ncombv(i,ibase)),nci,lfnout)
            mspin(offsetm+1)=spin(i)
            do j=2,spin(i)+1
               mspin(offsetm+j)=mspin(offsetm+j-1)-2
            enddo
            offsetm=offsetm+spin(i)+1
          end do

          do j = 1, spin(1) + 1
            ms1  = spin(1) - 2*(j-1)
            do k = 1, spin(2) + 1
              ms2 = spin(2) - 2*(k-1)
              if ( ms1 + ms2 .eq. nspin ) then
                do det1 = 1, ndet_mspin(1,j)
                  do det2 = 1, ndet_mspin(2,k)
*                    dnorm = cg(spin(1)/2.d0,spin(2)/2.d0,            ! filtering after multiplying with the
*     &                         ms1/2.d0,ms2/2.d0,                    ! CG-coefficients significantly reduces
*     &                         nspin/2.d0,nspin/2.d0)                ! the number of dets in the MEBF, but might
*                    prod = dnorm*ci_mspin(det1,j)*ci_mspin(det2,k)   ! introduce some spin contamination
                    prod = ci_mspin(det1,j)*ci_mspin(det2,k)
                    if (abs(prod) .gt. thresh_CI) cib = cib + 1
                  end do
                end do
              end if
            end do
          end do
          deallocate(mspin,ndet_mspin,ci_mspin,occ_mspin)
        end if
        if ( ipr. gt. 10 ) then
          write(lfnout,*)
          write(lfnout,'(a,i4,a,i10)') 
     &             ' number of dets in MEBF ',ibase,' :',cib
        end if
        idetb(ibase) = cib
        maxcib = max(cib,maxcib)
      end do
      write(lfnout,*)
      write(lfnout,'(2a,i10)')' Maximum number of determinants in ',
     &                         'MEBFs (maxcib) : ',maxcib
      flush(lfnout)
      return
      end subroutine gronor_determine_maxcib
