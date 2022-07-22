!>  @brief
!>    Construct the MEBF vectors from the fragment vectors
!>
!>  @authors
!>    Aitor Sanchez-Mansilla, Coen de Graaf
!>    (URV, Tarragona)
!>    Tjerk Straatsma
!>    (ORNL)
!>
!>  @date
!>    July 2022
!>

      subroutine gronor_mebf_vectors(ibase)
      use cidef
      use cidist
      use gnome_data
      use gnome_parameters
      implicit none

      integer :: ibase
      integer :: i,j,k,jstart,kstart
      integer, parameter   :: maxact=32 ! because of the packing, maxact is limited to 32 (I guess)

      if(idbg.ge.1.and.me.eq.master) idbg=100

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i4)') date(1:8),time(1:8),
     &       ' makebasestate for nfrags ',nfrags
        flush(lfndbg)
      endif
      
      inactb(ibase)=0
      nactb(ibase)=0
      do i=1,nfrags
        inactb(ibase)=inactb(ibase)+inactm(ncombv(i,ibase))
        nactb(ibase)=nactb(ibase)+nactm(ncombv(i,ibase))
      enddo

      if(nactb(ibase).gt.maxact) then
        call errquit(me,3,"Too many active orbitals")
!     write(lfnout,*)' * * *  Error  * * *'
!     write(lfnout,*)'Max. number of active orbitals: ',maxact
!     write(lfnout,*)'Actual active orbitals        : ',nactb(ibase)
!     stop 'Too many active orbitals'
      endif

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a,i4,a,i6)') date(1:8),time(1:8),
     &       ' number of inactive orbitals for state',ibase,
     &       ' is ',inactb(ibase)
        write(lfndbg,'(a,1x,a,1x,a,i4,a,i6)') date(1:8),time(1:8),
     &       ' number of active orbitals for state  ',ibase,
     &       ' is ',nactb(ibase)
        flush(lfndbg)
      endif
      
!     Combining the vectors of the fragments

      do i=1,nbasis
        do j=1,nbasis
          vecsb(i,j,ibase)=0.0d0
        enddo
      enddo
      jstart=1
      kstart=1
      do i=1,nfrags
        do j=1,inactm(ncombv(i,ibase))
          do k=1,nbasm(ncombv(i,ibase))
            vecsb(kstart+k-1,jstart+j-1,ibase)=
     &           vecsm(k,j,ncombv(i,ibase))
          enddo
        enddo
        kstart=kstart+nbasm(ncombv(i,ibase))
        jstart=jstart+inactm(ncombv(i,ibase))
      enddo

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' vecsb step 1 completed'
        flush(lfndbg)
      endif
      
      jstart=inactb(ibase)+1
      kstart=1

      do i=1,nfrags
        do j=1,nactm(ncombv(i,ibase))
          do k=1,nbasm(ncombv(i,ibase))
            vecsb(kstart+k-1,jstart+j-1,ibase)=vecsm(k,j+
     &           inactm(ncombv(i,ibase)),ncombv(i,ibase))
          enddo
        enddo
        kstart=kstart+nbasm(ncombv(i,ibase))
        jstart=jstart+nactm(ncombv(i,ibase))
      enddo

      if(idbg.ge.50) then
        call swatch(date,time)
        write(lfndbg,'(a,1x,a,1x,a)') date(1:8),time(1:8),
     &       ' vecsb step 2 completed'
        flush(lfndbg)
      endif
      
      return
      end 
