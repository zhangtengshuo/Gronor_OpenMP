      subroutine gronor_makebasestate                                     
     &     (lfnout,civm,vecsm,ioccm,nbasm,nactm,                            
     &     inactm,idetm,maxci,maxvec,maxnact,mstates,                     
     &     civb,vecsc,ioccb,inactb,nactb,idetb,                           
     &     maxcib,nbasis,mactb,                        
     &     nmol,ncomb,nspin,thresh_CI)

      implicit none
      integer :: maxact
      parameter (maxact=31)
      integer :: nspin,nmol,mstates,maxci,maxvec,maxnact
      integer :: inactb,nactb,idetb,maxcib,nbasis,mactb
      integer :: ioccm(maxnact,maxci,mstates)
      integer :: nbasm(mstates),nactm(mstates),inactm(mstates)
      integer :: idetm(mstates)
      integer :: ioccb(maxcib),ncomb(nmol)
      integer :: nmstates,nci,mhighs,jjstate,is,iprod,ivecs
      integer :: nspintot,ioff,idetoutoutd,llspin,kocc,lfnout,ioff2
      integer :: addspin,spin_sum,ns,opp_jjstate
      integer, allocatable :: intocc(:,:,:),intoccout(:,:,:),           
     &     jspin(:),idetout(:),kspin(:),jcomb(:),kcomb(:),
     &     nactms(:),jjspin(:),spin_counter(:),idetoutout(:),
     &     multiplicity(:) 
      integer :: i,j,k,l,m,jj,kk,ll,kstart,jstart
      real (kind=8) :: dnorm,thresh_CI
      real (kind=8) :: civb(maxcib),vecsc(nbasis,nbasis)
      real (kind=8),allocatable :: ciint(:,:),ciout(:,:)
      real (kind=8) :: civm(maxci,mstates),vecsm(maxvec,maxvec,mstates)
      real (kind=8), external :: ddot,cg
      logical :: oalls,onew,ocomparocc

      inactb=0
      nactb=0
      idetb=1
      do i=1,nmol
         inactb=inactb+inactm(ncomb(i))
         nactb=nactb+nactm(ncomb(i))
         idetb=idetb*idetm(ncomb(i))
      enddo

      if(nactb.gt.maxact) stop 'too many active orbitals > maxact'
      if(idetb.gt.maxcib) stop 'too many determinants > maxcib'

      do i=1,nbasis
         do j=1,nbasis
            vecsc(i,j)=0.0d0
         enddo
      enddo
      jstart=1
      kstart=1
      do i=1,nmol
         do j=1,inactm(ncomb(i))
            do k=1,nbasm(ncomb(i))
               vecsc(kstart+k-1,jstart+j-1)=vecsm(k,j,ncomb(i))
            enddo
         enddo
         kstart=kstart+nbasm(ncomb(i))
         jstart=jstart+inactm(ncomb(i))
      enddo
      
      jstart=inactb+1
      kstart=1
      
      do i=1,nmol
         do j=1,nactm(ncomb(i))
            do k=1,nbasm(ncomb(i))
               vecsc(kstart+k-1,jstart+j-1)=vecsm(k,j+inactm(ncomb(i)),
     &              ncomb(i))
            enddo
         enddo
         kstart=kstart+nbasm(ncomb(i))
         jstart=jstart+nactm(ncomb(i))
      enddo

      allocate(jspin(nmol))
      
      do i=1,nmol
        jspin(i)=0
        do j=1,nactm(ncomb(i))
          if(abs((ioccm(j,1,ncomb(i)))).eq.1)
     &         jspin(i)=jspin(i)+ioccm(j,1,ncomb(i))
        enddo
      enddo

      oalls=.true.
      
      do i=2,nmol
         if(jspin(i).ne.jspin(1))oalls=.false.
      enddo
      
      if(jspin(1).ne.0)oalls=.false.
      
!     Singlet-singlet couple
      
      if(oalls.or.nmol.eq.1) then
        call productf(civm,ioccm,nactm,idetm,nmol,ncomb,civb,ioccb,
     &       maxci,mstates,maxnact,mactb,maxcib,idetb,thresh_CI,lfnout)
      else
        nmstates=0
        nci=0
        mhighs=0
        spin_sum=0
        do i=1,nmol
          nmstates=nmstates+jspin(i)+1
          nci=max(nci,idetm(ncomb(i)))
          mhighs=max(mhighs,jspin(i))
          spin_sum=spin_sum+jspin(i)+1
        enddo
        
!     Construction of ms states

        nci=(spin_sum+2)*nci
!     safe also for > 2 fragments 

!     nci=(mhighs+2)*nci
!     nci=int(sqrt(maxcib)/(2*mhighs))

        allocate(ciint(nci,nmstates),intocc(maxnact,nci,nmstates),
     &       idetout(nmstates),kspin(nmstates),nactms(nmstates))
        jjstate=0
        do i=1,nmol
          jjstate=jjstate+1
          do jj=1,idetm(ncomb(i))
            ciint(jj,jjstate)=civm(jj,ncomb(i)) ! first ms state from molcas -> civec
            do kk=1,nactm(ncomb(i))
              intocc(kk,jj,jjstate)=ioccm(kk,jj,ncomb(i))
            enddo
          enddo
          
          idetout(jjstate)=idetm(ncomb(i))
          kspin(jjstate)=jspin(i)
          nactms(jjstate)=nactm(ncomb(i))

          if(jspin(i).eq.0) goto 10
          
          if (mod(jspin(i)+1,2).eq.0)  then
!     two cases, odd or even multiplicity
            ns=((jspin(i)+1)/2)-1
!     multiplicity/2 - 1
          else
            ns=jspin(i)/2
!     (multiplicity-1)/2
          endif
          
          is=jspin(i)           ! ms spin of the current state

          do j=1,ns
!     do is=jspin(i)-2,-jspin(i)+2,-2 ! build all intermediate ms states with sminop
            jjstate=jjstate+1
            is=is-2
               call sminop(ciint(1,jjstate-1),intocc(1,1,jjstate-1),
     &              ciint(1,jjstate),intocc(1,1,jjstate),
     &              idetout(jjstate-1),nactm(ncomb(i)),
     &              idetout(jjstate),maxnact,lfnout)
               kspin(jjstate)=is ! ms spin of jjstate
               nactms(jjstate)=nactm(ncomb(i))
               if (idetout(jjstate).gt.nci)  stop 'error: nci too small'
               dnorm=ddot(idetout(jjstate),ciint(1,jjstate),1,
     &              ciint(1,jjstate),1)
               dnorm=1.0d0/dsqrt(dnorm)
               call dscal(idetout(jjstate),dnorm,ciint(1,jjstate),1)
            enddo
            
            if(mod(jspin(i)+1,2).eq.0)then
               opp_jjstate=jjstate
               ns=(jspin(i)+1)/2
            else
               opp_jjstate=jjstate-1
               ns=jspin(i)/2
            endif

            do j=1,ns
               jjstate=jjstate+1
               do jj=1,idetout(opp_jjstate)
                  ciint(jj,jjstate)=ciint(jj,opp_jjstate)
                  do kk=1,nactm(ncomb(i))
                     if(abs(intocc(kk,jj,opp_jjstate)).eq.1) then
                        intocc(kk,jj,jjstate)=-intocc(kk,jj,opp_jjstate)
                     else
                        intocc(kk,jj,jjstate)=intocc(kk,jj,opp_jjstate)
                     endif
                  enddo
               enddo
               idetout(jjstate)=idetout(opp_jjstate)
               kspin(jjstate)=-kspin(opp_jjstate)
               nactms(jjstate)=nactm(ncomb(i))
               opp_jjstate=opp_jjstate-1
            enddo

 10         continue
         enddo
         
!     Product of the states

         iprod=1
         allocate(multiplicity(nmol))
         do i=1,nmol
            iprod=iprod*(jspin(i)+1)
            multiplicity(i)=jspin(i)+1
         enddo

         allocate(jcomb(nmol),kcomb(nmol),jjspin(nmol),
     &        spin_counter(nmol))
         if(nci**2.gt.maxcib) then
           write(*,*)' Warning maxcib too small',maxcib,nci**2
           stop
         endif
         allocate(idetoutout(iprod))
         ivecs=0
         ioff2=1

         do i=1,iprod
            call calciloop(i,jcomb,nmol,multiplicity)
            nspintot=0
            ioff=0
            idetoutoutd=1
            do j=1,nmol
               nspintot=nspintot+kspin(ioff+jcomb(j))
               spin_counter(j)=nspintot
               jjspin(j)=kspin(ioff+jcomb(j))
               kcomb(j)=ioff+jcomb(j)
               ioff=ioff+multiplicity(j)
               idetoutoutd=idetoutoutd*idetout(kcomb(j))
            enddo
            if(nspintot.eq.nspin) then
               ivecs=ivecs+1
               if(nmol.le.2)then ! TWO MOLECULES CASE
                  call productf(ciint,intocc,nactms,idetout,nmol,kcomb,
     &                 civb(ioff2),ioccb(ioff2),nci,
     &                 nmstates,maxnact,mactb,maxcib-ioff2+1,
     &                 idetoutoutd,thresh_CI,lfnout)   
                  idetoutout(ivecs)=idetoutoutd
                  dnorm=cg(jspin(1)/2.0d0,jspin(2)/2.0d0,
     &                 jjspin(1)/2.0d0,jjspin(2)/2.0d0,
     &                 nspin/2.0d0,nspin/2.0d0)
                  call dscal(idetoutout(ivecs),dnorm,civb(ioff2),1)
                else
                  
!     Three or more molecules
                  write(*,*)'work in progress'
                  
               endif
               ioff2=ioff2+idetoutoutd
            endif
         enddo

         idetb=0
         llspin=ivecs-1
         ioff2=1
         do i=1,ivecs
            ioff2=ioff2+idetoutout(i)
         enddo
         idetb=ioff2-1
         deallocate(ciint)
         deallocate(nactms)
         deallocate(intocc)
         deallocate(idetout)
         deallocate(kspin)
         deallocate(jjspin)
         deallocate(spin_counter)
         deallocate(jcomb)
         deallocate(kcomb)
         deallocate(idetoutout)
         deallocate(multiplicity)
      endif
      deallocate(jspin)    
      return
      end
      
      subroutine gronor_old_makebasestate
     &     (lfnout,civm,vecsm,ioccm,nbasm,nactm,
     &     inactm,idetm,maxci,maxvec,maxnact,mstates,
     &     civb,vecsc,ioccb,inactb,nactb,idetb,
     &     maxcib,nbasis,mactb,
     &     nmol,ncomb,nspin,thresh_CI)

      implicit none
      integer :: maxact
      parameter (maxact=31)
      integer :: nspin,nmol,mstates,maxci,maxvec,maxnact
      integer :: inactb,nactb,idetb,maxcib,nbasis,mactb
      real (kind=8) :: civm(maxci,mstates),vecsm(maxvec,maxvec,mstates)
      integer :: ioccm(maxnact,maxci,mstates)
      integer :: nbasm(mstates),nactm(mstates),inactm(mstates)
      integer :: idetm(mstates)
      real (kind=8) :: civb(maxcib),vecsc(nbasis,nbasis)
      integer :: ioccb(maxcib),ncomb(nmol)
      real (kind=8),allocatable :: ciint(:,:),ciout(:,:)
      integer, allocatable :: intocc(:,:,:),intoccout(:,:,:),           &
     & jspin(:),idetout(:),kspin(:),jcomb(:),kcomb(:),idetoutout(:),    &
     & nactms(:)
      real (kind=8) :: dnorm,thresh_CI
      real (kind=8), external :: ddot

      integer :: nmstates,nci,mhighs,jjstate,jj,kk,is,iprod,ivecs
      integer :: nspintot,ioff,idetoutoutd,llspin,kocc,lfnout,ioff2

      logical :: oalls,onew,ocomparocc

      real (kind=8) :: dsi(3,2)
      integer :: i,j,k, kstart, jstart
      data dsi/ 1.0d0, 1.0d0, 0.0d0, 1.0d0, 1.0d0, 1.0d0/

      inactb=0
      nactb=0
      idetb=1
      do i=1,nmol
       inactb=inactb+inactm(ncomb(i))
       nactb=nactb+nactm(ncomb(i))
       idetb=idetb*idetm(ncomb(i))
      enddo
      if(nactb.gt.maxact)stop 'too many active orbitals'
      if(idetb.gt.maxcib)then
       print *,idetb,maxcib
       stop 'you fool!!!'
      endif
      do i=1,nbasis
       do j=1,nbasis
        vecsc(i,j)=0.0d0
       enddo
      enddo
      jstart=1
      kstart=1

      do i=1,nmol
       do j=1,inactm(ncomb(i))
        do k=1,nbasm(ncomb(i))
         vecsc(kstart+k-1,jstart+j-1)=vecsm(k,j,ncomb(i))
        enddo
       enddo
       kstart=kstart+nbasm(ncomb(i))
       jstart=jstart+inactm(ncomb(i))
      enddo

      jstart=inactb+1
      kstart=1

      do i=1,nmol
       do j=1,nactm(ncomb(i))
        do k=1,nbasm(ncomb(i))
         vecsc(kstart+k-1,jstart+j-1)=vecsm(k,j+inactm(ncomb(i)),
     &                                                        ncomb(i))
        enddo
       enddo
       kstart=kstart+nbasm(ncomb(i))
       jstart=jstart+nactm(ncomb(i))
      enddo

      allocate(jspin(nmol))

      do i=1,nmol
       jspin(i)=0
       do j=1,nactm(ncomb(i))
        if(abs((ioccm(j,1,ncomb(i)))).eq.1) jspin(i)=jspin(i)+
     &                                           ioccm(j,1,ncomb(i))
       enddo
      enddo

      oalls=.true.

      do i=2,nmol
       if(jspin(i).ne.jspin(1))oalls=.false.
      enddo

      if(jspin(1).ne.0)oalls=.false.

!     Singlet-singlet couple

      if(oalls) then
        call old_productf(
     &       civm,ioccm,nactm,idetm,nmol,ncomb,civb,ioccb,maxci,
     &               mstates,maxnact,mactb,maxcib,idetb,thresh_CI)
      else
       nmstates=0
       nci=0
       mhighs=0
       do i=1,nmol
        nmstates=nmstates+jspin(i)+1
        nci=max(nci,idetm(ncomb(i)))
        mhighs=max(mhighs,jspin(i))
       enddo

!     Need better way to calculate nci

       nci=(mhighs+2)*nci
       if(mhighs.gt.2)stop 'program the rest yourself!'
       allocate(ciint(nci,nmstates),intocc(maxnact,nci,nmstates),
     &          idetout(nmstates),kspin(nmstates),nactms(nmstates))
       jjstate=0
       do i=1,nmol
        jjstate=jjstate+1
        do jj=1,idetm(ncomb(i))
         ciint(jj,jjstate)=civm(jj,ncomb(i))
         do kk=1,nactm(ncomb(i))
          intocc(kk,jj,jjstate)=ioccm(kk,jj,ncomb(i))
         enddo
        enddo
        idetout(jjstate)=idetm(ncomb(i))
        kspin(jjstate)=jspin(i)
        nactms(jjstate)=nactm(ncomb(i))
        if(jspin(i).eq.0) goto 10
        do is=jspin(i)-2,-jspin(i)+2,-2
         jjstate=jjstate+1
         call sminop(ciint(1,jjstate-1),intocc(1,1,jjstate-1),
     &               ciint(1,jjstate),intocc(1,1,jjstate),
     &               idetout(jjstate-1),nactm(ncomb(i)),
     &               idetout(jjstate),maxnact)
         kspin(jjstate)=is
         nactms(jjstate)=nactm(ncomb(i))
         if (idetout(jjstate).gt.nci)then
           print *,'nci::',nci,idetout(jjstate),jjstate
           stop 'nci error'
         endif

         dnorm=ddot(idetout(jjstate),ciint(1,jjstate),1,
     &              ciint(1,jjstate),1)
         dnorm=1.0d0/dsqrt(dnorm)
         call dscal(idetout(jjstate),dnorm,ciint(1,jjstate),1)
        enddo
        jjstate=jjstate+1
        do jj=1,idetm(ncomb(i))
         ciint(jj,jjstate)=civm(jj,ncomb(i))
         do kk=1,nactm(ncomb(i))
          if(abs(ioccm(kk,jj,ncomb(i))).eq.1) then
           intocc(kk,jj,jjstate)=-ioccm(kk,jj,ncomb(i))
          else
           intocc(kk,jj,jjstate)=ioccm(kk,jj,ncomb(i))
          endif
         enddo
         ciint(jj,jjstate)=civm(jj,ncomb(i))
        enddo
        idetout(jjstate)=idetm(ncomb(i))
        kspin(jjstate)=-jspin(i)
        nactms(jjstate)=nactm(ncomb(i))
   10   continue
       enddo

       iprod=1
       do i=1,nmol
        jspin(i)=jspin(i)+1
        iprod=iprod*(jspin(i))
       enddo
       allocate(jcomb(nmol),kcomb(nmol))
       if(nci**2.gt.maxcib) print*,'WARNING maxcib too small',
     &    maxcib,nci**2
         allocate(idetoutout(iprod))
       ivecs=0
       ioff2=1
       do i=1,iprod
        call calciloop(i,jcomb,nmol,jspin)
        nspintot=0
        ioff=0
        idetoutoutd=1
        do j=1,nmol
         nspintot=nspintot+kspin(ioff+jcomb(j))
         kcomb(j)=ioff+jcomb(j)
         ioff=ioff+jspin(j)
         idetoutoutd=idetoutoutd*idetout(kcomb(j))
        enddo
        if(nspintot.eq.nspin) then
         ivecs=ivecs+1
         call old_productf(
     &        ciint,intocc,nactms,idetout,nmol,kcomb,
     &                 civb(ioff2),ioccb(ioff2),nci,
     &                 nmstates,maxnact,mactb,maxcib-ioff2+1,
     &                 idetoutoutd,thresh_CI)

         idetoutout(ivecs)=idetoutoutd
         ioff2=ioff2+idetoutoutd
        endif
       enddo
       dnorm=1.0d0/(dsqrt(ivecs*1.0d0))
       idetb=0

!      TODO: Need to take care of the sign when forming the combinations!
!      Currently only trip-trip-singlet...

       llspin=ivecs-1
       ioff2=1
       do i=1,ivecs
        call dscal(idetoutout(i),dnorm*dsi(i,llspin),civb(ioff2),1)
        ioff2=ioff2+idetoutout(i)
       enddo
       idetb=ioff2-1

       deallocate(ciint)
       deallocate(nactms)
       deallocate(intocc)
       deallocate(idetout)
       deallocate(kspin)
       deallocate(jcomb)
       deallocate(kcomb)
       deallocate(idetoutout)
      endif

      deallocate (jspin)

      return
      end
