!> @brief   Read the CI coefficients and the vectors
!>
!> @details For the moment, both civ-files and det-files can be used.
!>          Duplicates are removed from the list of determinants by
!>          remove_duplicates_fildet. Two linear combinations of
!determinants
!>          of two CSFs sharing the same orbital occupation (that is,
!>          arising from the same electronic configuration) may contain
!>          some determinants that are identical.
!> @author  Aitor Sanchez, URV
!> @author  Coen de Graaf, URV
!> @date    October 2020
!> @todo    Remove the possibility to read from civ files
!> @param   idet_raw number of determinants for each monomer function
!> @param   civ_raw  CI coefficients of the raw list of determinants
!> @param   maxci_raw maximum number of raw determinants
!> @param   occ_string ocupation string as read from the root_nnn.det
!>          file

      subroutine gronor_check_vectors_and_determinants()
      use inp
      use cidef
      use gnome_parameters

      implicit none
      
      integer :: i,length
      
      do i = 1, mstates
        write(fildet,203) root(1:inp_strlen(root)),i
        length=inp_strlen(fildet)
        open(unit=lfndet,file=fildet(1:length),form='formatted',
     &       status='old',err=998)
        close(unit=lfndet,status='keep')
 203  format(a,'_',i3.3,'.det')
      enddo
      idevel=1
      return
 998  continue
      
      do i = 1, mstates
        write(filvec,204) root(1:inp_strlen(root)),i
        length=inp_strlen(filvec)
        open(unit=lfnvec,file=filvec(1:length),form='formatted',
     &       status='old',err=999)
        close(unit=lfnvec,status='keep')
 204    format(a,'_',i3.3,'.vec')
        write(filciv,205) root(1:inp_strlen(root)),i
        length=inp_strlen(filciv)
        open(unit=lfnciv,file=filciv(1:length),form='formatted',
     &       status='old',err=999)
        close(unit=lfnciv,status='keep')
 205    format(a,'_',i3.3,'.civ')
      enddo
      idevel=0
      return
 999  continue
      
      stop 'No vector files'
      return
      end subroutine gronor_check_vectors_and_determinants
      
      subroutine gronor_read_vectors_and_determinants() 

      use inp
      use cidef
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

      maxci     = 0
      maxci_raw = 0
      maxnact   = 0
      maxvec    = 0
      idetm     = 0
      allocate ( idet_raw(mstates) )

* Starting with the determinants
      
! Filling idetm, inactm,nactm and determine the maximum dets and active orbitals.
! For reading from vecdet a raw value is calculated as duplicates are still to be removed

      if ( vecdet ) then
        do i = 1, mstates
          write(fildet,203) root(1:inp_strlen(root)),i
          length=inp_strlen(fildet)
          open(unit=lfndet,file=fildet(1:length),form='formatted',
     &         status='old',err=999)
          read(lfndet,*) inactm(i)
          read(lfndet,*) sdum,dumstr
          nactm(i)=len_trim(dumstr)
          idet_raw(i) = 1
          do 
            read(lfndet,*,iostat=istat)
            if ( istat .ne. 0 ) then
              rewind(lfndet)
              exit
            else
              idet_raw(i) = idet_raw(i) + 1
            endif
          enddo 
          close(unit=lfndet)
          maxci_raw = max(maxci_raw,idet_raw(i))
          maxnact = max(maxnact,nactm(i))
        end do
        allocate(civ_raw(maxci_raw,mstates))
        allocate(occ_string(maxci_raw,mstates))          ! a string is printed in the vecdet
      else
        do i = 1, mstates
          write(filciv,200) root(1:inp_strlen(root)),i
          length=inp_strlen(filciv)
          open(unit=lfnciv,file=filciv(1:length),form='formatted',
     &           status='old',err=998)
          read(lfnciv,*) idetm(i),nactm(i),inactm(i)
          close(unit=lfnciv)
          maxci = max(maxci,idetm(i))
          maxnact = max(maxnact,nactm(i))
        end do
        allocate(civm(maxci,mstates))
        allocate(ioccm(maxnact,maxci,mstates))
      end if
      if ( vecdet ) then 
        do i = 1, mstates
          write(fildet,203) root(1:inp_strlen(root)),i
          length=inp_strlen(fildet)
          open(unit=lfndet,file=fildet(1:length),form='formatted',
     &           status='old',err=999)            
          read(lfndet,*)
          do j = 1, idet_raw(i)
            read(lfndet,*) civ_raw(j,i),occ_string(j,i)
            occ_string(j,i)=adjustl(occ_string(j,i))
          enddo                                    
          call remove_duplicates_fildet
     &           (idet_raw(i),idetm(i),nactm(i),
     &           civ_raw(:,i),occ_string(:,i),thresh_CI)
          close(unit=lfndet)                        
          maxci = max(maxci,idetm(i))     
        end do
        allocate( civm(maxci,mstates) )               ! Now that the number of unique dets is known
        allocate( ioccm(maxnact,maxci,mstates) )      ! civm and ioccm can be allocated
        do i = 1, mstates
          do j = 1, idetm(i)
            civm(j,i)=civ_raw(j,i)
            dumstr=trim(adjustl(occ_string(j,i)))
            do k = 1, nactm(i)
              if ( dumstr(k:k) .eq. '2' ) then
                ioccm(k,j,i) = 2
              elseif ( dumstr(k:k) .eq. 'a' ) then
                ioccm(k,j,i) = 1
              elseif ( dumstr(k:k) .eq. 'b' ) then
                ioccm(k,j,i) = -1
              elseif ( dumstr(k:k) .eq. '0' ) then
                ioccm(k,j,i) = 0
              else
                stop 'weird occupation'
              endif
            enddo
          enddo
        enddo
        deallocate(civ_raw)
        deallocate(idet_raw)
        deallocate(occ_string)
      else
        allocate( iocc_tmp(maxnact) )
        do i = 1, mstates
          iocc_tmp = 0
          write(filciv,200) root(1:inp_strlen(root)),i
          length=inp_strlen(filciv)
          open(unit=lfnciv,file=filciv(1:length),form='formatted',
     &           status='old',err=998)
          read(lfnciv,*)
          do j = 1, idetm(i)
            read(lfnciv,1002) l,civm(j,i),(iocc_tmp(k),k=1,nactm(i))
 1002       format(i7,f25.14,6x,20(i2,x))
            do k = 1, nactm(i)
              if ( iocc_tmp(k) .eq. 11 ) then
                ioccm(k,j,i)=2
              elseif( iocc_tmp(k) .eq. 1 ) then
                ioccm(k,j,i)=-1
              elseif( iocc_tmp(k) .eq. 10 ) then
                ioccm(k,j,i)=1
              elseif( iocc_tmp(k) .eq. 0 ) then
                ioccm(k,j,i)=0
              else
                print *,(iocc_tmp(l),l=1,nactm(i))
                stop 'weird occupation'
              endif
            enddo
          enddo
          close(unit=lfnciv)
        end do
        deallocate( iocc_tmp )
      end if

* Next, reading the vectors 
      do i = 1, mstates
        write(filvec,201) root(1:inp_strlen(root)),i
        length=inp_strlen(filvec)
        open(unit=lfnvec,file=filvec(1:length),form='formatted',
     &        status='old',err=997)
        read(lfnvec,*)nbasm(i)
        maxvec = max(maxvec,nbasm(i))
        close(unit=lfnvec)
      end do
      allocate( vecsm(maxvec,maxvec,mstates) )
      do i = 1, mstates
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

* Finally, the output section
      if ( ipr .ge. 1 ) then
        write(lfnout,604)
 604    format(/,' Dimensions',//,
     &     '  istate   ndetm   nactm  ninatm   nbasm   files',/)
        do i = 1, mstates
          if ( vecdet ) then
            write(fildet,203) root(1:inp_strlen(root)),i
          else
            write(fildet,200) root(1:inp_strlen(root)),i
          end if
          write(filvec,201) root(1:inp_strlen(root)),i
          write(lfnout,605) i,idetm(i),nactm(i),inactm(i),nbasm(i),
     &              fildet(1:inp_strlen(fildet)),
     &              filvec(1:inp_strlen(filvec))
 605      format(5i8,3x,a,', ',a) 
        end do
        write(lfnout,606) maxci,maxnact,maxvec
 606    format(8x,'  ------  ------',8x,'  ------',/,
     &     '    max:',i8,i8,8x,i8)
      end if
 200  format(a,'_',i3.3,'.civ')
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

      if ( ipr .ge. 3 ) then
        write(lfnout,607)
 607    format(//,'  istate   idetm           ci coefficient')
        do i = 1, mstates
          do j = 1, idetm(i)
            if ( j .eq. 1 ) then
              write(lfnout,608) i,j,civm(j,i),
     &             (ioccm(k,j,i),k=1,nactm(i))
 608          format(/,2i8,f25.14,6x,32(i2,1x))
            else
              write(lfnout,609) j,civm(j,i),
     &             (ioccm(k,j,i),k=1,nactm(i))
 609          format(8x,i8,f25.14,6x,32(i2,1x))
            endif
          end do
          if ( ipr .ge. 4) then
            write(filvec,201) root(1:inp_strlen(root)),i
            length=inp_strlen(filvec)
            write(lfnout,2002) filvec(1:length),nbasm(i)
 2002       format(/,' Vector File ',a,//,
     &          ' Number of basis functions ',i8,/)
            do j = 1, nbasm(i)
              write(lfnout,2003) j,(vecsm(k,j,i),k=1,nbasm(i))
 2003         format(i5,5f25.14,/,(5x,5f25.14))
            end do
          endif
        end do
      end if
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
      call errquit()
 997  write(lfnout,987) filvec(1:inp_strlen(filvec))
      call errquit()
 998  write(lfnout,988) filciv(1:inp_strlen(filciv))
      call errquit()
 999  write(lfnout,989) fildet(1:inp_strlen(fildet))
      call errquit()
 985  format('Unable to open symmetry file ',a)
 987  format('Unable to open vects  file ',a)
 988  format('Unable to open civec  file ',a)
 989  format('Unable to open vecdet file ',a)

      end subroutine gronor_read_vectors_and_determinants

