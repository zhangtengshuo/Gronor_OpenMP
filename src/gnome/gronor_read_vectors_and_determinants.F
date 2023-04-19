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

!> @brief   Read the CI coefficients and the vectors
!>
!>          Duplicates are removed from the list of determinants by
!>          remove_duplicates_fildet. Two linear combinations of determinants
!>          of two CSFs sharing the same orbital occupation (that is,
!>          arising from the same electronic configuration) may contain
!>          some determinants that are identical.
!> @author  Aitor Sanchez, URV
!> @author  Coen de Graaf, URV
!> @date    October 2020
!> @param   idet_raw number of determinants for each monomer function
!> @param   maxci_raw maximum number of raw determinants
!> @param   civ_raw  CI coefficients of the raw list of determinants
!> @param   occ_raw ocupation string as read from the root_nnn.det
!>          file

      subroutine gronor_read_vectors_and_determinants()

      use inp
      use cidef
      use cidist
      use gnome_parameters

      implicit none

      integer                          :: i,j,k,nume
      integer                          :: maxci_raw
      integer           , allocatable  :: idet_raw(:)
      real      (kind=8)               :: sdum,sume,sump,sumc
      real      (kind=8), allocatable  :: civ_raw(:,:)
      character (len=80)               :: dumstr
      character (len=255), allocatable :: occ_raw(:,:)
      logical                          :: ostate

      maxci=0
      maxci_raw=0
      maxnact=0
      maxvec=0
      idetm=0
      allocate ( idet_raw(mstates) )

!     Starting with the determinants

!     Filling idetm, inactm,nactm and determine the maximum dets and active orbitals.
!     For reading from vecdet a raw value is calculated as duplicates are still to be removed

      allocate(nFrozen(mstates))
      allocate(fragLabel(mstates))
      allocate(ecasscf(mstates))
      allocate(ecaspt2(mstates))
      allocate(ecorr(mstates))
      allocate(nElectrons(mstates))
      allocate(spinm(mstates))
      ncorr=1
      do i=1,mstates
        write(fildet,2203) trim(fragfile(i))
 2203   format(a,'.det')
        open(unit=lfndet,file=trim(fildet),form='formatted',                    &
     &       status='old',err=999)
        read(lfndet,1600) inactm(i),nFrozen(i),idet_raw(i),                     &
     &       fragLabel(i),ecasscf(i),nElectrons(i),tau_MO, ecaspt2(i)
        if(ecaspt2(i).ge.0.0d0) ncorr=0
        ecorr(i)=ecaspt2(i)-ecasscf(i)
        read(lfndet,*) sdum,dumstr
        nactm(i)=len_trim(dumstr)
        maxci_raw=max(maxci_raw,idet_raw(i))
        maxnact=max(maxnact,nactm(i))
! determine the spin of the fragment wave function
! spin multiplicity = (n_alpha - n_beta) + 1
        spinm(i) = 0
        do j = 1, nactm(i)
          if ( dumstr(j:j) .eq. 'a' ) spinm(i) = spinm(i) + 1
          if ( dumstr(j:j) .eq. 'b' ) spinm(i) = spinm(i) - 1
        end do
        spinm(i) = spinm(i) + 1
        inactm(i) = inactm(i) - nFrozen(i)
        close(lfndet)
      enddo
 1600 format(2i5,i12,4x,a6,f22.12,i5,e10.3,f22.12)
      allocate(mebfLabel(nbase))
      if ((fragLabel(1)(1:8).eq.'no_label').or.(nmol.gt.15)) then
         mebfLabels = .false.
         do j = 1, nbase
            write(mebfLabel(j),618) j
         end do
 618     format(i4)
      else
         mebfLabels = .true.
         if(ncorr.ne.0) then
            write(lfnout,607)                                                   &
     &           ' MEBF ','Electrons   Sum of fragment energies',               &
     &           'E(CASSCF)','E(CASPT2)','E_corr'
         else
            write(lfnout,627)                                                   &
     &           ' MEBF ','Electrons   Sum of fragment energies',               &
     &           'E(CASSCF)'
         endif
 607     format(/,a,t16,a,/,t34,a,t56,a,t78,a,/)
 627     format(/,a,t16,a,/,t34,a,/)
         do j = 1,nbase
            write(mebfLabel(j),616)                                             &
     &           (trim(fragLabel(ncombv(k,j))),k=1,nmol)
 616        format(10a)
            nume=0
            sume=0.0d0
            sump=0.0d0
            sumc=0.0d0
            do k=1,nmol
               nume=nume+nElectrons(ncombv(k,j))
               sume=sume+ecasscf(ncombv(k,j))
               sump=sump+ecaspt2(ncombv(k,j))
               sumc=sumc+ecorr(ncombv(k,j))
            enddo
            if(ncorr.ne.0) then
               write(lfnout,617) mebflabel(j),nume,sume,sump,sumc
            else
               write(lfnout,637) mebflabel(j),nume,sume
            endif
 617        format(1x,a,t16,i5,3x,3f22.12)
 637        format(1x,a,t16,i5,3x,f22.12)
         enddo
      endif
      write(lfnarx,701) nbase,nmol
      write(lfnxrx,701) nbase,nmol
 701  format('State',2i10)
      do j=1,nbase
            nume=0
            sume=0.0d0
            do k=1,nmol
               nume=nume+nElectrons(ncombv(k,j))
               sume=sume+ecasscf(ncombv(k,j))
            enddo
         write(lfnarx,617) mebflabel(j),nume,sume
         write(lfnxrx,617) mebflabel(j),nume,sume
      enddo

! Read the CI vectors from the det files. The first loop over mstates 
! counts the number of unique determinants, and the second loop
! constructs the list by summing the coefficients of the duplicates.
! This could be done in one shot at the cost of allocating civm and
! ioccm with maxci_raw instead of maxci.
! 
      allocate(civ_raw(maxci_raw,mstates))
      allocate(occ_raw(maxci_raw,mstates))          ! a string is printed in the vecdet
      civ_raw = 0.0d0
      occ_raw = ''
      do i=1,mstates
        write(fildet,2203) trim(fragfile(i))
        open(unit=lfndet,file=trim(fildet),form='formatted',                    &
     &         status='old',err=999)
        read(lfndet,*)
        do j=1,idet_raw(i)
          read(lfndet,*) civ_raw(j,i),occ_raw(j,i)
          occ_raw(j,i)=adjustl(occ_raw(j,i))
        enddo
        close(unit=lfndet)
        if (idet_raw(i) .gt. 1 ) then
          call quicksort_string(civ_raw(:,i),occ_raw(:,i),idet_raw(i))
        end if
        idetm(i) = 1
        do j = 2, idet_raw(i)
          if ( occ_raw(j,i) .ne. occ_raw(j-1,i) ) then
            idetm(i) = idetm(i) + 1
          endif
        end do
        maxci=max(maxci,idetm(i))
      enddo
! Now that the number of unique dets is known, civm and occm_string can
! be allocated and the unique list of determinants constructed.
! For the moment we keep ioccm, but this array is temporary
      allocate( civm(maxci,mstates) )
      allocate( occm_string(maxci,mstates) )
      allocate( ioccm(maxnact,maxci,mstates) )
      do i=1,mstates
        k = 1
        civm(1,i) = civ_raw(1,i)
        occm_string(1,i) = occ_raw(1,i)
        do j=2,idet_raw(i)
          if ( occ_raw(j,i) .ne. occ_raw(j-1,i) ) then
            k = k + 1
            civm(k,i) = civ_raw(j,i)
            occm_string(k,i) = occ_raw(j,i)
          else
            civm(k,i) = civm(k,i) + civ_raw(j,i)
          end if
        end do
        do j = 1, idetm(i)
          dumstr=trim(adjustl(occm_string(j,i)))
          do k=1,nactm(i)
            if(dumstr(k:k).eq.'2') then
              ioccm(k,j,i)=2
            elseif(dumstr(k:k).eq.'a') then
              ioccm(k,j,i)=1
            elseif(dumstr(k:k).eq.'b') then
              ioccm(k,j,i)=-1
            else 
              ioccm(k,j,i)=0
            endif
          enddo
        enddo
      enddo
      deallocate(civ_raw)
      deallocate(idet_raw)
      deallocate(occ_raw)

!     Next, reading the vectors

      do i=1,mstates
        write(filvec,2201) trim(vecfile(i))
 2201   format(a,'.vec')
        open(unit=lfnvec,file=trim(filvec),form='formatted',                    &
     &        status='old',err=997)
        read(lfnvec,*)nbasm(i)
        maxvec=max(maxvec,nbasm(i))
        close(unit=lfnvec)
      enddo
      allocate( vecsm(maxvec,maxvec,mstates) )
      do i=1,mstates
        write(filvec,2201) trim(vecfile(i))
        open(unit=lfnvec,file=trim(filvec),form='formatted',                    &
     &        status='old',err=997)
        read(lfnvec,*)
        do j=1,nbasm(i)
          read(lfnvec,1003) (vecsm(k,j,i),k=1,nbasm(i))
 1003     format(4F18.14)
        enddo
        close(unit=lfnvec)
      enddo

!     Finally, the output section

      if(ipr.ge.20) then
         if(ncorr.ne.0) then
            write(lfnout,603)
         else
            write(lfnout,623)
         endif
 603    format(/,' Dimensions',//,                                              &
     &       '       istate   ndetm   nactm  ninatm   nbasm',                   &
     &       '  nelecs         CASSCF energy  ',                                &
     &       '       CASPT2 energy    files',/)
 623    format(/,' Dimensions',//,                                              &
     &       '       istate   ndetm   nactm  ninatm   nbasm',                   &
     &       '  nelecs         CASSCF energy    files',/)
        do i=1,mstates
           ostate=.true.
           do j=1,i-1
              if(trim(fragfile(i)).eq.trim(fragfile(j))) ostate=.false.
           enddo
!          do j=1,nbase
!            do k=1,nmol
!              if(ncombv(k,j).eq.i) ostate=.true.
!            enddo
!          enddo
           if(ostate) then
              write(fildet,203) trim(fragfile(i))
              write(filvec,201) trim(vecfile(i))
              if(ncorr.ne.0) then
                 write(lfnout,604) fragLabel(i),                                &
     &                idetm(i),nactm(i),inactm(i),nbasm(i),                     &
     &                nElectrons(i),ecasscf(i),ecaspt2(i),                      &
     &                trim(fildet),                                             &
     &                trim(filvec)
              else
                 write(lfnout,634) fragLabel(i),                                &
     &                idetm(i),nactm(i),inactm(i),nbasm(i),                     &
     &                nElectrons(i),ecasscf(i),                                 &
     &                trim(fildet),                                             &
     &                trim(filvec)
              endif
 604        format(8x,2x,a2,5i8,3x,2f22.12,2x,a,', ',a)
 634        format(8x,2x,a2,5i8,3x,f22.12,2x,a,', ',a)
          endif
        enddo
        write(lfnout,605) maxci,maxnact,maxvec
 605    format(13x,' ------  ------',8x,'  ------',/,                           &
     &     '        max:',i8,i8,8x,i8)
      endif
 201  format(a,'.vec')
 203  format(a,'.det')

      if(ipr.gt.0) then
        write(lfnout,606) (i,i=1,nbase)
 606    format(/,' Molecular states included in this calculation',//,           &
     &       ' State        : ',20i4)
        do i=1,nmol
          if(i.eq.1) then
            write(lfnout,608) i,(fragLabel(ncombv(i,j)),j=1,nbase)
          else
            write(lfnout,614) i,(fragLabel(ncombv(i,j)),j=1,nbase)
          endif
 608      format(/,' Fragment',i4,' : ',20(2x,a2))
 614      format(' Fragment',i4,' : ',20(2x,a2))
        enddo
      endif

      if(ipr.ge.30) then
        write(lfnout,609)
 609    format(//,'  istate   idetm           ci coefficient')
        do i=1,mstates
          do j=1,idetm(i)
            if(j.eq.1) then
              write(lfnout,610) i,j,civm(j,i),occm_string(j,i)
 610          format(/,2i8,f25.14,6x,a32)
            else
              write(lfnout,611) j,civm(j,i),occm_string(j,i)
 611          format(8x,i8,f25.14,6x,a32)
            endif
          enddo
! format statements 610 and 611 assume a maximum of 32 active orbitals
          if(ipr.ge.40) then
            write(filvec,201) trim(vecfile(i))
            write(lfnout,612) trim(filvec),nbasm(i)
 612        format(/,' Vector File ',a,//,                                      &
     &          ' Number of basis functions ',i8,/)
            do j=1,nbasm(i)
              write(lfnout,613) j,(vecsm(k,j,i),k=1,nbasm(i))
 613          format(i5,5f25.14,/,(5x,5f25.14))
            enddo
          endif
        enddo
      endif
      flush(lfnout)

      write(filsys,202) trim(root)
 202  format(a,'.sys')
      open(unit=lfnsys,file=trim(filsys),form='formatted',                      &
     &      status='old',err=995)

      if(ipr.gt.0) write(lfnout,615) trim(filsys)
 615  format(/,' System information read from ',a)

      return

 995  write(lfnout,985) trim(filsys)
      call gronor_abort(261,trim(filsys))
 997  write(lfnout,987) trim(filvec)
      call gronor_abort(262,trim(filvec))
 999  write(lfnout,989) trim(fildet)
      call gronor_abort(263,trim(fildet))
 985  format('Unable to open system file ',a)
 987  format('Unable to open vects  file ',a)
 989  format('Unable to open vecdet file ',a)

      end subroutine gronor_read_vectors_and_determinants

