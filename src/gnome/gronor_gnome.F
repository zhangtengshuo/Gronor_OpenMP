
!> @brief
!! General Non-Orthogonal Matrix Element calculation using GNOME structure
!!     
!! @author  R. Broer, RUG
!! @author  T. P. Straatsma, ORNL
!! @date    2016
!!
  
#include "tgt_guard.fh"
      
      subroutine gronor_gnome(lfndbg,ihc,nhc)

      use mpi
      use cidist
      use gnome_integrals
      use gnome_parameters
      use gnome_data

      implicit none

      integer :: lfndbg
      real (kind=8) :: vecs1(nbasis,nbasis),vecs2(nbasis,nbasis)

      integer :: nworv(2)
      integer :: idet

      character (len=25) :: idfmess, line

      integer (kind=4) :: ierr
      integer :: nwor, n, k, iv, ib, i, j, ihc, nhc

      integer :: ntvc,ivc,ibas

      if(idbg.ge.12) then
      write(lfndbg,600) nbasis
 600  format(/,' Number of basis functions is',t45,i8)
      write(lfndbg,603) nelec(1)
 603  format(' Number of electrons is',t45,i8)
      endif

      e1=0.0d0
      e2=0.0d0
      e2c=0.0d0
      etot=0.0d0
      ttest=0

      do idet=1,2

        if(idbg.ge.13) then
          write(lfndbg,604) idet
 604      format(/,' Transformation of MO set',i8,/)
          write(lfndbg,605) n,(ioccup(k,idet),k=1,nact(idet))
 605      format(' Active orbital occupation in symrep',i4,' : ',32i3)
        endif
        
        call timer_start(11)
        call gronor_transvc(lfndbg,idet)
        call timer_stop(11)

        call timer_start(12)
!       call gronor_order(lfndbg,idet)
        call timer_stop(12)

        call timer_start(13)
        call gronor_tranout(lfndbg,idet)
        call timer_stop(13)

        if(idbg.ge.14) write(lfndbg,606) idet
 606    format(/,' Construction of M.O.set',i2,' completed')

      enddo

      ntcla=ntcl(1)
      ntclb=ntcl(2)
      ntopa=ntop(1)
      ntopb=ntop(2)
      
      nveca=ntcla+ntopa
      nvecb=ntclb+ntopb
      ntesta=nveca+ntcla
      ntestb=nvecb+ntclb
      
      if(ntesta.ne.ntestb) stop'the number of electrons are not equal'
      
      nelecs=ntesta
      n1bas=nbas*(nbas+1)/2
      nstdim=max(1,nelecs*nelecs,n1bas)
      mbasel=max(nelecs,nbas)
      
      if(nveca.ne.ntcl(1)+ntop(1)) stop 'here1'
      if(nvecb.ne.ntcl(2)+ntop(2)) stop 'here2'
      
      if(idbg.gt.10) then
        do idet=1,2
          ntvc=ntcl(idet)+ntop(idet)
          write(lfndbg,1603) ntvc
 1603     format(/,' Closed shell M.O.''s',i5,/)
          do ivc=1,ntvc
            if(ivc.eq.ntcl(idet)+1) write(lfndbg,1604)
 1604       format(/,' Open shell M.O.'' s:')
            write(lfndbg,1605)  ' (',ivc,')',
     &           (vec(ivc,ibas,idet),ibas=1,nbas)
 1605       format(a2,i3,a1,(t9,10f12.8))
          enddo
        enddo
      endif

      do ib=1,nbas
        do iv=1,nveca
          va(iv,ib)=vec(iv,ib,1)
        enddo
      enddo
      do ib=1,nbas
        do iv=1,nvecb
          vb(iv,ib)=vec(iv,ib,2)
        enddo
      enddo

      if(iamacc.gt.0) then
        
!_ACCTGT_($acc data copyin(s,ta,va,vb,tat) create(tb,w1,w2,taa,u,w,ev))
!_ACCTGT_($acc& create(diag,bdiag,cdiag,bsdiag,csdiag,sdiag,aaa,tt,aat,sm))
!_ACCTGT_($acc& create(diagl,bdiagl,bsdiagl,csdiagl,sml,aaal,ttl,aatl,tatl))
        
!_OMPTGT_($omp target data map(to:s,ta,va,vb,tat) map(alloc:tb,w1,w2,taa,u,w,ev))
!_OMPTGT_($omp& map(alloc:diag,bdiag,cdiag,bsdiag,csdiag,sdiag,aaa,tt,aat,sm))
!_OMPTGT_($omp& map(alloc: diagl,bdiagl,bsdiagl,csdiagl,sml,aaal,ttl,aatl,tatl)) 

!     Calculations of the overlap matrices

        call timer_start(14)
        call gronor_moover(lfndbg)
        call timer_stop(14)

!     Calculation of the cofactor matrices and arrays corresponding
!     to the total overlap

        call timer_start(15)
        call gronor_cofac1(lfndbg)
        call timer_stop(15)
        
        if(idbg.ge.12) then
          if(ising.eq.0) write(lfndbg,609)
          if(ising.eq.1) write(lfndbg,610)
          if(ising.eq.2) write(lfndbg,611)
          if(ising.eq.3) write(lfndbg,612)
 609      format(/,' A has no singularities')
 610      format(/,' A has a single singularity')
 611      format(/,' A has two singularities: ',
     &         'one electron matrix elements are zero')
 612      format(/,' A has more than two singularities: ',
     &         'one and two electron matrix elements are zero')
        endif

        if(ising.lt.3) then
        
          if(corres) then
            call timer_start(16)
            call gronor_cororb(lfndbg)
            call timer_stop(16)
          endif      
          
          call timer_start(17)
          call gronor_tramat(lfndbg)
          call timer_stop(17)
          
          if(idipole.ne.0) then
            call timer_start(18)
            call gronor_dipole(lfndbg)
            call timer_stop(18)
          endif

          if(idbg.ge.13) then
            if(icalc.eq.1.or.icalc.eq.3) write(lfndbg,613) icalc
 613        format(' No calculation of two-electron matrix elements',i4)
          endif

!     Transformation of the  m.o.'s into the bassisset of the two
!     electron integrals

          call timer_start(19)
!         call gronor_trsym(lfndbg)
          call timer_stop(19)

!     transformation of the first order cofactor matrix
!     (x-matrix to f-matrix in terms of the basis set of the 2-el.integr

          call timer_start(20)
          call gronor_tramat2(lfndbg)
          call timer_stop(20)

!     Calculation of the one electron Hamiltonian matrix elements

          call timer_start(21)
          if(icalc.le.1.and.ising.le.1) call gronor_gnone(lfndbg)
          call timer_stop(21)
        endif

!     Calculation of the two-electron matrix elements

        if((icalc.eq.2.or.icalc.eq.0)) then
          call timer_start(22)
          if(lbatch.lt.0) then
            if(labels.eq.0) then
              call gronor_gntwo_batch(lfndbg,ihc,nhc)
            else
              call gronor_gntwo_batch_mlab(lfndbg,ihc,nhc)
            endif
          elseif(lbatch.gt.0) then
            if(labels.eq.0) then
              call gronor_gntwo_batch2(lfndbg,ihc,nhc)
            else
              call gronor_gntwo_batch2_mlab(lfndbg,ihc,nhc)
            endif
          else
            if(labels.eq.0) then
              if(ising.le.2) call gronor_gntwo(lfndbg)
            else
              if(ising.le.2) call gronor_gntwo_mlab(lfndbg)
            endif
          endif
          call timer_stop(22)
        endif

!_ACCTGT_($acc end data)
!_OMPTGT_($omp end target data)

      else

!     OpenMP routines

!     Calculations of the overlap matrices

        call timer_start(14)
        call gronor_moover_omp(lfndbg)
        call timer_stop(14)

!     Calculation of the cofactor matrices and arrays corresponding
!     to the total overlap

        call timer_start(15)
        call gronor_cofac1_omp(lfndbg)
        call timer_stop(15)
        
        if(idbg.ge.12) then
          if(ising.eq.0) write(lfndbg,609)
          if(ising.eq.1) write(lfndbg,610)
          if(ising.eq.2) write(lfndbg,611)
          if(ising.eq.3) write(lfndbg,612)
        endif

        if(ising.lt.3) then

          if(corres) then
            call timer_start(16)
            call gronor_cororb(lfndbg)
            call timer_stop(16)
          endif
          
          call timer_start(17)
          call gronor_tramat_omp(lfndbg)
          call timer_stop(17)
          
          if(idipole.ne.0) then
            call timer_start(18)
            call gronor_dipole(lfndbg)
            call timer_stop(18)
          endif
          
          if(idbg.ge.13) then
            if(icalc.eq.1.or.icalc.eq.3) write(lfndbg,613) icalc
          endif
          
!     Transformation of the  m.o.'s into the bassisset of the two
!     electron integrals

          call timer_start(19)
!     call gronor_trsym_omp(lfndbg)
          call timer_stop(19)

!     transformation of the first order cofactor matrix
!     (x-matrix to f-matrix in terms of the basis set of the 2-el.integr

          call timer_start(20)
          call gronor_tramat2_omp(lfndbg)
          call timer_stop(20)

!     Calculation of the one electron Hamiltonian matrix elements

          call timer_start(21)
          if(icalc.le.1.and.ising.le.1) call gronor_gnone_omp(lfndbg)
          call timer_stop(21)

          endif
!     Calculation of the two-electron matrix elements

          if((icalc.eq.2.or.icalc.eq.0)) then
            call timer_start(22)
            if(lbatch.lt.0) then
              if(labels.eq.0) then
                call gronor_gntwo_omp_batch(lfndbg,ihc,nhc)
              else
                call gronor_gntwo_omp_batch_mlab(lfndbg,ihc,nhc)
              endif
            elseif(lbatch.gt.0) then
              if(labels.eq.0) then
                call gronor_gntwo_omp_batch2(lfndbg,ihc,nhc)
              else
                call gronor_gntwo_omp_batch2_mlab(lfndbg,ihc,nhc)
              endif
            else
              if(labels.eq.0) then
                call gronor_gntwo_omp(lfndbg)
              else
                call gronor_gntwo_omp_mlab(lfndbg)
              endif
            endif
            call timer_stop(22)
          endif

      endif

      return
      end
