      subroutine print_results(hbase,sbase,nbase)
      
      use gnome_parameters, only : itest,ncols
      use cidef           , only : filf75,lfnf75,lfnout,lfntst

      implicit none

      integer,intent(in)            :: nbase
      integer                       :: i,j,k,ii,ik,il
      integer                       :: nk,info,lwork
      real(kind=8),intent(in)       :: hbase(nbase,nbase)
      real(kind=8),intent(in)       :: sbase(nbase,nbase)
      real(kind=8)                  :: hbt(nbase),tc(nbase,nbase)
      real(kind=8)                  :: hsave(nbase,nbase)
      real(kind=8)                  :: ssave(nbase,nbase)
      real(kind=8)                  :: hev(nbase)
      real(kind=8),allocatable      :: work(:)
      character(len=4)              :: token

      write(*,670)
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,670)
        write(lfntst,670)
 670    format(//,' Hamiltonian Matrix',/)
        write(*,671) (i,i=ii,il)
        write(lfnout,671) (i,i=ii,il)
        write(lfntst,671) (i,i=ii,il)
 671    format(6x,7(6x,i8,6x))
        do j=1,nbase
          write(*,672) j,(hbase(i,j),i=ii,il)
          write(lfnout,672) j,(hbase(i,j),i=ii,il)
 672      format(i5,1x,10f20.10)
          do i=ii,il
            hbt(i)=hbase(i,j)
            if(dabs(hbt(i)).lt.1.0d-08) hbt(i)=0.0d0
          enddo
          write(lfntst,1672) j,(hbt(i),i=ii,il)
 1672     format(i5,1x,10f20.8)
        enddo
        write(lfnout,*) ' '
        write(*,*) ' '
      enddo
      write(*,673)
      nk=nbase/ncols
      if(mod(nbase,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nbase,k*ncols)
        write(lfnout,673)
        write(lfntst,673)
 673    format(//,' Overlap Matrix',/)
        write(*,671) (i,i=ii,il)
        write(lfnout,671) (i,i=ii,il)
        write(lfntst,671) (i,i=ii,il)
        do j=1,nbase
          write(*,672) j,(sbase(i,j),i=ii,il)
          write(lfnout,672) j,(sbase(i,j),i=ii,il)
          do i=ii,il
            hbt(i)=sbase(i,j)
            if(dabs(hbt(i)).lt.1.0d-08) hbt(i)=0.0d0
          enddo
          write(lfntst,1672) j,(hbt(i),i=ii,il)
        enddo
        write(lfnout,*) ' '
        write(*,*) ' '
      enddo
      
      
      open(unit=lfnf75,file=filf75,form='unformatted')
      write(75)nbase
      write(75)((hbase(j,i),i=1,j),j=1,nbase)
      write(75)((sbase(j,i),i=1,j),j=1,nbase)
      close(unit=lfnf75,status='keep')
      
      if(nbase.gt.1) then
        do i=2,nbase
          do j=1,i-1
            tc(i,j)=(hbase(i,j)
     &           -0.5d0*(hbase(i,i)+hbase(j,j))*sbase(i,j))/
     &           (1.0d0-sbase(i,j)*sbase(i,j))
          enddo
        enddo

        write(lfnout,676)
        write(*,676)
 676    format(//,' Electronic Couplings',/)
        nk=nbase/ncols
        if(mod(nbase,ncols).ne.0) nk=nk+1
        do k=1,nk
          ii=(k-1)*ncols+1
          il=min(nbase,k*ncols)
          write(lfnout,671) (i,i=ii,il)
          write(*,671) (i,i=ii,il)
          do j=1,nbase
            if(j.eq.1) then
              write(lfnout,672) 1
              write(*,672) 1
            else
              ik=min(j-1,il)
              write(lfnout,672) j,(tc(j,i),i=ii,ik)
              write(*,672) j,(tc(j,i),i=ii,ik)
            endif
          enddo
          write(lfnout,*) ' '
          write(*,*) ' '
        enddo

!     Diagonalize the NOCI matrix and print the energies and wave functions
!     Don't touch hbase and sbase, needed again when shifting
          
        info=0
        lwork = 4 * nbase
        hsave = hbase
        ssave = sbase
        allocate(work(lwork))
        call dsygv(1,'V','L',nbase,hsave,nbase,
     &                   ssave,nbase,hev,work,lwork,info)
        deallocate(work)
        if ( info .ne. 0 ) write(*,*) 'something went wrong in dyegv'
        
! Coen 2020/07/18
! Check if the eigenvectors come in rows or columns!
! a transpose(hsave) might be needed
! printing takes care of it now: hsave(j,i) instead of hsave(i,j)
          
        write(*,677)
        write(lfnout,677)
 677    format(//,' NOCI energies and wave functions',/)

        nk=nbase/ncols
        if(mod(nbase,ncols).ne.0) nk=nk+1
        do k=1,nk
          ii=(k-1)*ncols+1
          il=min(nbase,k*ncols)
          
          write(*,'(a,12i20)')'    State:',(i,i=ii,il)
          write(*,678)'   Energy:    ', (hev(i),i=ii,il)
          write(*,*) ' '
          write(lfnout,'(a,12i20)')'     State:',(i,i=ii,il)
          write(lfnout,678)'    Energy:   ', (hev(i),i=ii,il)
          write(lfnout,*) ' '
          if(itest.eq.2) then
            write(lfntst,677)
            write(lfntst,'(a,12i20)')'     State:',(i,i=ii,il)
            write(lfntst,1678)'    Energy:   ', (hev(i),i=ii,il)
            write(lfntst,*) ' '
          endif
 678      format((a,12f20.10),//)
 1678     format((a,12f20.8),//)
          token="MEBF"
          do j=1,nbase
            write(*,679) token,j,(hsave(j,i),i=ii,il)
            write(lfnout,679) token,j,(hsave(j,i),i=ii,il)
 679        format(a4,1x,i5,4x,12f20.10)
            if(itest.eq.2) then
              write(lfntst,1679) token,j,(hsave(j,i),i=ii,il)
 1679         format(a4,1x,i5,4x,12f20.8)
            endif
            token="    "
          enddo
        enddo
      endif
      
      return
      
      end subroutine print_results
