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
      
!> @brief   Counts the number of determinants in each MEBF
!>
!> @details Going through the list of determinants of the fragment
!>          wave function for each MEBF and decide if the product
!>          (and the products of the different MS-components) is
!>          above the threshold for being included in the calculation
!>          of the Hamiltonian. Some additional checks on the consistency
!>          of the input are made.
!>          
!>          
!> @author  Coen de Graaf, URV
!> @date    October 2020
!>
!> @todo    Decide whether we filter before or after multiplying with
!>          the Clebsch-Gordon coefficient. This may substantially reduce
!>          the number of determinants in the MEBF at the cost of possibly
!>          introducing some spin contamination.
!> 
      subroutine gronor_determine_maxcib
      use cidef
      use cidist
      use gnome_parameters
      implicit none

      integer              :: istate,iact,ibase,i,j,k
      integer              :: det1,det2,cib
      integer              :: ms1,ms2,offsetm,mspin_states,spin(2)
      integer              :: minspin,maxspin
      integer,allocatable  :: mspin(:),ndet_mspin(:,:),occ_mspin(:,:,:)
      integer,allocatable  :: iocc_tmp(:,:)

      real (kind=8)               :: prod
      real (kind=8), allocatable  :: ci_mspin(:,:)
      real (kind=8), allocatable  :: cicoef_tmp(:)

      allocate( spinm(mstates) )
      spinm=0
      do istate=1,mstates
        do iact=1,nactm(istate)
          if(abs(ioccm(iact,1,istate)).eq.1) then
            spinm(istate)=spinm(istate)+ioccm(iact,1,istate)
          endif
        enddo
      enddo

      maxcib=0

      if(nmol.eq.1) then
        do i=1,mstates
          maxcib=max(maxcib,idetm(i))
          if(spinm(i).ne.nspin) then
            write(lfnout,601) i,nspin
 601        format(' State ',i,' is incompatible with spin ',
     &           'of molecule (',i,')')
            call errquit(me,9,"Incompatible spin state")
          endif
        enddo
        return
      endif

      do ibase=1,nbase
        if (nmol.eq.2) then
          minspin=abs(spinm(ncombv(1,ibase))-spinm(ncombv(2,ibase)))
          maxspin=    spinm(ncombv(1,ibase))+spinm(ncombv(2,ibase))
          if (( maxspin.lt.nspin ).or.( minspin.gt.nspin)) then
            write(lfnout,602) ibase,spinm(1),spinm(2),nspin
 602        format(' Spins of fragments incompatible with total spin',
     &           ' for MEBF',i,/,' Spin(1) = ',i4,/,' Spin(2) = ',i4,/,
     &           ' Cannot be coupled to ',i)
            call errquit(me,9,"Incompatible spin state")
          end if
        end if
        cib=0
        if((nspin.eq.0).and.
     &       (spinm(ncombv(1,ibase)).eq.0).and.
     &       (spinm(ncombv(2,ibase)).eq.0)) then
          do det1=1,idetm(ncombv(1,ibase))
            do det2=1,idetm(ncombv(2,ibase))
              prod=civm(det1,ncombv(1,ibase)) *
     &                            civm(det2,ncombv(2,ibase))
              if(abs(prod).gt.thresh) cib=cib+1
            enddo
          enddo
        else
          mspin_states=spinm(ncombv(1,ibase))+1 +
     &                        spinm(ncombv(2,ibase))+1
          allocate(mspin(mspin_states))                       ! ms values
          allocate(ndet_mspin(2,mspin_states))
          allocate(ci_mspin(nci,mspin_states))
          allocate(occ_mspin(maxnact,nci,mspin_states))       ! not really needed here
          ndet_mspin=0
          offsetm=0
          do i=1,nmol 
            istate=ncombv(i,ibase)
            spin(i)=spinm(istate)
            allocate(iocc_tmp(nactm(istate),idetm(istate)))
            allocate(cicoef_tmp(idetm(istate)))
            iocc_tmp=0
            cicoef_tmp=0
            do j=1,idetm(istate)
              cicoef_tmp(j)=civm(j,istate)
              do k=1,nactm(istate)
                iocc_tmp(k,j)=ioccm(k,j,istate)
              end do
             end do
            call sminop(cicoef_tmp,ci_mspin(1,offsetm+1),
     &           iocc_tmp,
     &           occ_mspin(1,1,offsetm+1),
     &           idetm(istate),ndet_mspin(i,:),
     &           spin(i),nactm(istate),nci,maxnact)
             deallocate(cicoef_tmp)
             deallocate(iocc_tmp)
            mspin(offsetm+1)=spin(i)
            do j=2,spin(i)+1
               mspin(offsetm+j)=mspin(offsetm+j-1)-2
            enddo
            offsetm=offsetm+spin(i)+1
          enddo

          do j=1,spin(1)+1
            ms1 =spin(1)-2*(j-1)
            do k=1,spin(2)+1
              ms2=spin(2)-2*(k-1)
              if(ms1+ms2.eq.nspin) then
                do det1=1,ndet_mspin(1,j)
                  do det2=1,ndet_mspin(2,k)
*                    dnorm=cg(spin(1)/2.d0,spin(2)/2.d0,            ! filtering after multiplying with the
*     &                         ms1/2.d0,ms2/2.d0,                    ! CG-coefficients significantly reduces
*     &                         nspin/2.d0,nspin/2.d0)                ! the number of dets in the MEBF, but might
*                    prod=dnorm*ci_mspin(det1,j)*ci_mspin(det2,k)   ! introduce some spin contamination
                    prod=ci_mspin(det1,j)*ci_mspin(det2,k+spin(1)+1)
                    if (abs(prod).gt.thresh) cib=cib+1
                  enddo
                enddo
              endif
            enddo
          enddo
          deallocate(mspin,ndet_mspin,ci_mspin,occ_mspin)
        endif
        if(cib.eq.0) then 
          write(lfnout,603) spin(1)+1,spin(2)+1,nspin+1
 603      format(' Spin(1) = ',i,/,' Spin(2) = ',i,/,
     &         ' Total spin = ',i, ' is not possible')
          flush(lfnout)
          call errquit(me,9,"Incompatible spin state")
        endif

        if ((me.eq.master).and.(ipr.ge.30)) write(lfnout,604) ibase,cib
 604      format(' Number of determinants in MEBF ',i4,' : ',i10)

        idetb(ibase)=cib
        maxcib=max(cib,maxcib)
      enddo
      return
      end subroutine gronor_determine_maxcib
