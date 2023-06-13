      program merge_xyz

      implicit none

      integer :: i, k, iarg, nfil, n, natoms
      character (len=255) :: f,filxyz(10)
      character (len=2) :: atom
      real (kind=8) :: x,y,z

      nfil=0
      do i=1,10
         iarg=i
         call getarg(iarg,filxyz(i))
         if(len(trim(filxyz(i))).eq.0) exit
         nfil=i
      enddo

      natoms=0
      do i=2,nfil
         f=trim(filxyz(i))//".xyz"
         open(unit=10,file=trim(f),form="formatted",status="old")
         read(10,*) n
         close(unit=10,status="keep")
         natoms=natoms+n
      enddo

      f=trim(filxyz(1))//".xyz"
      open(unit=11,file=trim(f),form="formatted",status="unknown")
      write(11,'(i6)') natoms
      write(11,'(a,i6,a)') "Merged file from",nfil-1," xyz files"
      do i=2,nfil
         f=trim(filxyz(i))//".xyz"
         open(unit=10,file=trim(f),form="formatted",status="old")
         read(10,*) n
         read(10,*) f
         do k=1,n
            read(10,*) atom, x, y, z
            write(11,'(a,t3,3f16.8)') trim(atom),x,y,z
         enddo
         close(unit=10,status="keep")
      enddo
      close(unit=11,status="keep")
      

      end
