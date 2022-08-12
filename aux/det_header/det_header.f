      program det_header

      implicit none

      integer :: iarg,ifr,ito,i,k
      integer :: ninact,nfrozen,nraw,nelecs
      character (len=255) :: card,filout,fildet
      real(kind=8) :: ecasscf, ecaspt2
      character (len=6) :: label
      real (kind=8) :: energy,ecorr,tau
      real (kind=8), allocatable :: coef(:)
      character (len=50), allocatable :: occ(:)
      character (len=50) :: string3


      iarg=1
      call getarg(iarg,filout)
      iarg=2
      call getarg(iarg,fildet)

      ifr=index(filout,'_')+1
      ito=index(filout,'.')-1
      
      open(unit=10,file=trim(filout),form="formatted",status="old")

      ecasscf=0.0d0
      ecaspt2=0.0d0

      do while(.true.)
         read(10,'(a)',end=1,err=1) card
         if(card(1:29).eq."      RASSCF energy for state")
     &        read(card(49:64),'(f16.8)') ecasscf
         if(card(1:23).eq."      Total energy:    ")
     &        read(card(29:46),'(f18.10)') ecaspt2
      enddo
 1    continue

      close(unit=10)

      write(*,'(a,a)') "State=",filout(ifr:ito)
      write(*,'(a,f20.10)') "E(CASSCF)=",ecasscf
      write(*,'(a,f20.10)') "E(CASPT2)=",ecaspt2


      open(unit=11,file=trim(fildet),form="formatted",status="old")

      read(11,100) ninact,nfrozen,nraw,label,energy,nelecs,tau
 100  format(2i5,i12,4x,a6,f22.12,i5,e10.3,f22.12)

      allocate(coef(nraw),occ(nraw))

      do i=1,nraw
         read(11,101) coef(i),occ(i)
 101     format(e15.8,a)
      enddo

      string3=occ(1)
      do i=1,50
         if(string3(i:i).ne.' ') k=i
      enddo

      label="      "
      label(1:ito-ifr+1)=filout(ifr:ito)

      rewind(unit=11)
      write(11,102) ninact,nfrozen,nraw,label,ecasscf,nelecs,tau,ecaspt2
 102  format(2i5,i12,4x,a6,f22.12,i5,1pe10.3,0pf22.12)
      do i=1,nraw
         write(11,101) coef(i),occ(i)(1:k)
      enddo

      close(unit=11)

      end
