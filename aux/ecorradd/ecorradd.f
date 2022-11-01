      program ecorradd

      implicit none

      integer :: i,k,iarg,ninact,nfrozen,nraw,nelecs
      character (len=6) :: label
      character (len=255) :: string1, string2, fildet
      real (kind=8) :: energy,ecorr,tau
      real (kind=8), allocatable :: coef(:)
      character (len=50), allocatable :: occ(:)
      character (len=50) :: string3
      
      iarg=1
      call getarg(iarg,string1)
      iarg=2
      call getarg(iarg,string2)

      fildet=trim(string1)//".det"
      read(string2,*) ecorr

      open(unit=10,file=trim(fildet),form="formatted",status="old")
      read(10,100) ninact,nfrozen,nraw,label,energy,nelecs,tau
 100  format(2i5,i12,a6,f26.12,i5,e10.3)

      allocate(coef(nraw),occ(nraw))

      do i=1,nraw
         read(10,101) coef(i),occ(i)
 101     format(e15.8,a)
      enddo

      string3=occ(1)
      do i=1,50
         if(string3(i:i).ne.' ') k=i
      enddo

      rewind(unit=10)
      write(10,102) ninact,nfrozen,nraw,label,energy,nelecs,tau,ecorr
 102  format(2i5,i12,a6,f26.12,i5,1pe10.3,f26.12)
      do i=1,nraw
         write(10,101) coef(i),occ(i)(1:k)
      enddo

      close(unit=10)

      end
