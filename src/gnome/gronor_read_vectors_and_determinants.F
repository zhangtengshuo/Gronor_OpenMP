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
!> @param   civ_raw  CI coefficients of the raw list of determinants
!> @param   maxci_raw maximum number of raw determinants
!> @param   occ_string ocupation string as read from the root_nnn.det
!>          file

      subroutine gronor_read_vectors_and_determinants() 

      use inp
      use cidef
      use cidist
      use gnome_parameters

      implicit none

      integer                          :: i,j,k,l,istat,length
      integer                          :: maxci_raw
      integer           , allocatable  :: idet_raw(:),iocc_tmp(:)
      real      (kind=8)               :: sdum
      real      (kind=8), allocatable  :: civ_raw(:,:)
      character (len=80)               :: dumstr
      character (len=255), allocatable :: occ_string(:,:)
      character (len=255)              :: string

      maxci=0
      maxci_raw=0
      maxnact=0
      maxvec=0
      idetm=0
      allocate ( idet_raw(mstates) )

!     Starting with the determinants
      
!     Filling idetm, inactm,nactm and determine the maximum dets and active orbitals.
!     For reading from vecdet a raw value is calculated as duplicates are still to be removed

      do i=1,mstates
        write(fildet,203) root(1:inp_strlen(root)),i
        length=inp_strlen(fildet)
        open(unit=lfndet,file=fildet(1:length),form='formatted',
     &       status='old',err=999)
        read(lfndet,*) inactm(i)
        read(lfndet,*) sdum,dumstr
        nactm(i)=len_trim(dumstr)
        idet_raw(i)=1
        do 
          read(lfndet,*,iostat=istat)
          if(istat .ne. 0) then
            rewind(lfndet)
            exit
          else
            idet_raw(i)=idet_raw(i) + 1
          endif
        enddo 
        close(unit=lfndet)
        maxci_raw=max(maxci_raw,idet_raw(i))
        maxnact=max(maxnact,nactm(i))
      enddo
      allocate(civ_raw(maxci_raw,mstates))
      allocate(occ_string(maxci_raw,mstates))          ! a string is printed in the vecdet
      do i=1,mstates
        write(fildet,203) root(1:inp_strlen(root)),i
        length=inp_strlen(fildet)
        open(unit=lfndet,file=fildet(1:length),form='formatted',
     &         status='old',err=999)            
        read(lfndet,*)
        do j=1,idet_raw(i)
          read(lfndet,*) civ_raw(j,i),occ_string(j,i)
          occ_string(j,i)=adjustl(occ_string(j,i))
        enddo                                    
        call remove_duplicates_fildet
     &         (idet_raw(i),idetm(i),nactm(i),
     &         civ_raw(:,i),occ_string(:,i),thresh_CI)
        close(unit=lfndet)                        
        maxci=max(maxci,idetm(i))     
      enddo
      allocate( civm(maxci,mstates) )               ! Now that the number of unique dets is known
      allocate( ioccm(maxnact,maxci,mstates) )      ! civm and ioccm can be allocated
      do i=1,mstates
        do j=1,idetm(i)
          civm(j,i)=civ_raw(j,i)
          dumstr=trim(adjustl(occ_string(j,i)))
          do k=1,nactm(i)
            if(dumstr(k:k).eq.'2') then
              ioccm(k,j,i)=2
            elseif(dumstr(k:k).eq.'a') then
              ioccm(k,j,i)=1
            elseif(dumstr(k:k).eq.'b') then
              ioccm(k,j,i)=-1
            elseif(dumstr(k:k).eq.'0') then
              ioccm(k,j,i)=0
            else
              write(lfnout,'(2A)') 'weird occupation in ',dumstr
              write(lfnout,'(A,I4,A,I4)')'determinant ',j, 'of state ',i
              call errquit(me,9,"Inconsistent occupation")
              stop
            endif
          enddo
        enddo
      enddo
      deallocate(civ_raw)
      deallocate(idet_raw)
      deallocate(occ_string)

!     Next, reading the vectors
      
      do i=1,mstates
        write(filvec,201) root(1:inp_strlen(root)),i
        length=inp_strlen(filvec)
        open(unit=lfnvec,file=filvec(1:length),form='formatted',
     &        status='old',err=997)
        read(lfnvec,*)nbasm(i)
        maxvec=max(maxvec,nbasm(i))
        close(unit=lfnvec)
      enddo
      allocate( vecsm(maxvec,maxvec,mstates) )
      do i=1,mstates
        write(filvec,201) root(1:inp_strlen(root)),i
        length=inp_strlen(filvec)
        open(unit=lfnvec,file=filvec(1:length),form='formatted',
     &        status='old',err=997)
        read(lfnvec,*)
        do j=1,nbasm(i)
          read(lfnvec,1003) (vecsm(k,j,i),k=1,nbasm(i))
 1003     format(4F18.14)
        enddo
        close(unit=lfnvec)
      enddo

!     Finally, the output section
      
      if(ipr.ge.1) then
        write(lfnout,604)
 604    format(/,' Dimensions',//,
     &     '  istate   ndetm   nactm  ninatm   nbasm   files',/)
        do i=1,mstates
          write(fildet,203) root(1:inp_strlen(root)),i
          write(filvec,201) root(1:inp_strlen(root)),i
          write(lfnout,605) i,idetm(i),nactm(i),inactm(i),nbasm(i),
     &              fildet(1:inp_strlen(fildet)),
     &              filvec(1:inp_strlen(filvec))
 605      format(5i8,3x,a,', ',a) 
        enddo
        write(lfnout,606) maxci,maxnact,maxvec
 606    format(8x,'  ------  ------',8x,'  ------',/,
     &     '    max:',i8,i8,8x,i8)
      endif
 201  format(a,'_',i3.3,'.vec')
 203  format(a,'_',i3.3,'.det')

      write(lfnout,613) (i,i=1,nbase)
 613  format(/,' Molecular states included in this calculation',//,
     &      ' State        : ',20i4)
      write(lfnout,624)
 624  format(/)
      do i=1,nmol
        write(lfnout,617) i,(ncombv(i,j),j=1,nbase)
 617    format(' Molecule',i4,' : ',20i4)
      enddo

      if(ipr.ge.3) then
        write(lfnout,607)
 607    format(//,'  istate   idetm           ci coefficient')
        do i=1,mstates
          do j=1,idetm(i)
            if(j.eq.1) then
              write(lfnout,608) i,j,civm(j,i),
     &             (ioccm(k,j,i),k=1,nactm(i))
 608          format(/,2i8,f25.14,6x,32(i2,1x))
            else
              write(lfnout,609) j,civm(j,i),
     &             (ioccm(k,j,i),k=1,nactm(i))
 609          format(8x,i8,f25.14,6x,32(i2,1x))
            endif
          enddo
          if(ipr.ge.4) then
            write(filvec,201) root(1:inp_strlen(root)),i
            length=inp_strlen(filvec)
            write(lfnout,2002) filvec(1:length),nbasm(i)
 2002       format(/,' Vector File ',a,//,
     &          ' Number of basis functions ',i8,/)
            do j=1,nbasm(i)
              write(lfnout,2003) j,(vecsm(k,j,i),k=1,nbasm(i))
 2003         format(i5,5f25.14,/,(5x,5f25.14))
            enddo
          endif
        enddo
      endif
      flush(lfnout)

      write(filsym,202) root(1:inp_strlen(root))
 202  format(a,'.sym')
      length=inp_strlen(filsym)
      open(unit=lfnsym,file=filsym(1:length),form='formatted',
     &      status='old',err=995)

      if(ipr.ge.3) write(lfnout,610) filsym(1:inp_strlen(filsym))
 610  format(/,' System information read from ',a,/)

      return

 995  write(lfnout,985) filsym(1:inp_strlen(filsym))
      call errquit(me,6,"Error reading determinants")
 997  write(lfnout,987) filvec(1:inp_strlen(filvec))
      call errquit(me,6,"Error reading determinants")
 999  write(lfnout,989) fildet(1:inp_strlen(fildet))
      call errquit(me,6,"Error reading determinants")
 985  format('Unable to open symmetry file ',a)
 987  format('Unable to open vects  file ',a)
 989  format('Unable to open vecdet file ',a)

      end subroutine gronor_read_vectors_and_determinants

