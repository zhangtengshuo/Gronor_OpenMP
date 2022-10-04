      program reorder

      implicit none

      integer :: iarg,i,j,k,m,ic,jc,it,jt,num,ir,jr,iq,jq,ip,jp
      integer :: lfn1,lfn2,lfn3
      integer :: natoms
      integer :: natoms2
      integer :: iswap,jswap,nswap
      real (kind=8) :: d,dmax,dnew

      character (len=255) :: file1,file2,file3
      character (len=255) :: card1,card2

      character (len=2), allocatable :: name1(:)
      character (len=2), allocatable :: name2(:)
      character (len=2), allocatable :: namet(:)
      
      real (kind=8), allocatable :: coord1(:,:)
      real (kind=8), allocatable :: coord2(:,:)
      real (kind=8), allocatable :: coordt(:,:)

      
      real (kind=8), allocatable :: dist1(:,:)
      real (kind=8), allocatable :: dist2(:,:)

      integer (kind=8), allocatable :: ndx(:)
      integer (kind=8), allocatable :: ndx2(:)
      integer (kind=8), allocatable :: ndone(:)

      integer (kind=8), allocatable :: nb(:)
      integer (kind=8), allocatable :: ib(:,:)
      integer (kind=8), allocatable :: na(:)
      integer (kind=8), allocatable :: nh(:)
      integer (kind=8), allocatable :: nc(:)
      integer (kind=8), allocatable :: nn(:)
      integer (kind=8), allocatable :: no(:)
      integer (kind=8), allocatable :: id(:)

      integer (kind=8), allocatable :: nb2(:)
      integer (kind=8), allocatable :: ib2(:,:)
      integer (kind=8), allocatable :: na2(:)
      integer (kind=8), allocatable :: nh2(:)
      integer (kind=8), allocatable :: nc2(:)
      integer (kind=8), allocatable :: nn2(:)
      integer (kind=8), allocatable :: no2(:)
      integer (kind=8), allocatable :: id2(:)
      
      logical :: done

      iarg=1
      call getarg(iarg,file1)
      iarg=2
      call getarg(iarg,file2)
      iarg=3
      call getarg(iarg,file3)

      lfn1=10
      lfn2=11
      lfn3=12
      
      open(unit=lfn1,file=trim(file1),form='formatted',status='old')
      read(lfn1,*) natoms
      read(lfn1,'(a)') card1
      
      open(unit=lfn2,file=trim(file2),form='formatted',status='old')
      read(lfn2,*) natoms2
      read(lfn2,'(a)') card2

      if(natoms.ne.natoms2) stop "Number of atoms not equal"
      
      allocate(name1(natoms))
      allocate(coord1(natoms,3))
      allocate(dist1(natoms,natoms))
      allocate(name2(natoms))
      allocate(coord2(natoms,3))
      allocate(dist2(natoms,natoms))
      allocate(namet(natoms))
      allocate(coordt(natoms,3))
      
      allocate(ndx(natoms))
      allocate(ndx2(natoms))
      allocate(ndone(natoms))
      
      allocate(nb(natoms))
      allocate(ib(natoms,5))
      allocate(na(natoms))
      allocate(nh(natoms))
      allocate(nc(natoms))
      allocate(nn(natoms))
      allocate(no(natoms))
      allocate(id(natoms))
      
      allocate(nb2(natoms))
      allocate(ib2(natoms,5))
      allocate(na2(natoms))
      allocate(nh2(natoms))
      allocate(nc2(natoms))
      allocate(nn2(natoms))
      allocate(no2(natoms))
      allocate(id2(natoms))
      
      do i=1,natoms
         ndx(i)=i
         ndx2(i)=i
         ndone(i)=0
      enddo

      do i=1,natoms
         read(lfn1,*) name1(i),(coord1(i,j),j=1,3)
         read(lfn2,*) name2(i),(coord2(i,j),j=1,3)
      enddo

      close(unit=lfn1)
      close(unit=lfn2)

!      write(*,'(a)') "Reference"
      call analyze(natoms,name1,coord1,nb,ib,id)

      do i=1,natoms
         if(trim(name1(i)).ne.trim(name2(ndx2(i)))) then
            do j=i+1,natoms
               if(trim(name1(i)).eq.trim(name2(ndx2(j)))) then
                  k=ndx2(i)
                  ndx2(i)=ndx2(j)
                  ndx2(j)=k
                  exit
               endif
            enddo
         endif
      enddo

      do i=1,natoms
         namet(i)=name2(i)
         coordt(i,1)=coord2(i,1)
         coordt(i,2)=coord2(i,2)
         coordt(i,3)=coord2(i,3)
      enddo

      do i=1,natoms
         name2(i)=namet(ndx2(i))
         coord2(i,1)=coordt(ndx2(i),1)
         coord2(i,2)=coordt(ndx2(i),2)
         coord2(i,3)=coordt(ndx2(i),3)
      enddo

!      write(*,'(a)') "Target"
      call analyze(natoms,name2,coord2,nb2,ib2,id2)
      
      m=0
      do i=1,natoms
         if(m.eq.0.and.nb(i).eq.1.and.id(i).eq.id2(i)) m=i
      enddo

      do i=1,natoms
         ndx(i)=0
         ndx2(i)=0
         ndone(i)=0
      enddo

      if(m.ne.0) then
         ndx(m)=m
         num=1
         do while(num.gt.0)
            num=0
            do i=1,natoms
               if(ndx(i).eq.0) num=num+1
               if(ndx(i).gt.0) then
                  do j=1,nb2(i)
!                     write(*,'(a,5i3,2i16)') "Test ",i,j,ib(i,j),ib2(ndx(i),j),ndx(ib(i,j)),id2(ib2(ndx(i),j)),id(ib(i,j))
                     if(ndx(ib(i,j)).eq.0.and.id2(ib2(ndx(i),j)).eq.id(ib(i,j))) then
                        ndx(ib(i,j))=ib2(ndx(i),j)
!                        write(*,'(a,40i3)') "From ",(m,m=1,natoms)
!                        write(*,'(a,40i3)') "To   ",(ndx(m),m=1,natoms)
                     endif
                  enddo
               endif
            enddo
         enddo

         do i=1,natoms
            namet(i)=name2(i)
            coordt(i,1)=coord2(i,1)
            coordt(i,2)=coord2(i,2)
            coordt(i,3)=coord2(i,3)
         enddo

         do i=1,natoms
            name2(i)=namet(ndx(i))
            coord2(i,1)=coordt(ndx(i),1)
            coord2(i,2)=coordt(ndx(i),2)
            coord2(i,3)=coordt(ndx(i),3)
         enddo

      endif



      open(unit=lfn3,file=trim(file3),form='formatted',status='unknown')
      write(lfn3,'(i4)') natoms
      write(lfn3,'(a)') card2(1:len_trim(card1))
      do i=1,natoms
         write(lfn3,'(a,t4,3f15.8)') trim(name2(i)),(coord2(i,j),j=1,3)
      enddo

      close(unit=lfn3)
            
      end
      
      subroutine analyze(natoms,name,coord,nb,ib,id)
        
        implicit none

        integer :: natoms
        
        character (len=2) :: name(natoms)

        real (kind=8) :: coord(natoms,3)

        integer (kind=8) :: nb(natoms)
        integer (kind=8) :: id(natoms)
        integer (kind=8) :: ib(natoms,5)
        
        integer (kind=8), allocatable :: ia(:)
        
        integer (kind=8), allocatable :: n2a(:)
        integer (kind=8), allocatable :: n2h(:)
        integer (kind=8), allocatable :: n2c(:)
        integer (kind=8), allocatable :: n2n(:)
        integer (kind=8), allocatable :: n2o(:)
        
        integer (kind=8), allocatable :: n3(:)
        integer (kind=8), allocatable :: n3a(:)
        integer (kind=8), allocatable :: n3h(:)
        integer (kind=8), allocatable :: n3c(:)
        integer (kind=8), allocatable :: n3n(:)
        integer (kind=8), allocatable :: n3o(:)
        
        integer (kind=8), allocatable :: n4(:)
        integer (kind=8), allocatable :: n4a(:)
        integer (kind=8), allocatable :: n4h(:)
        integer (kind=8), allocatable :: n4c(:)
        integer (kind=8), allocatable :: n4n(:)
        integer (kind=8), allocatable :: n4o(:)

        real (kind=8), allocatable :: dm(:)
        
        integer :: i,j,k,m,id2,id3,id4
        real (kind=8) :: d,dcut

        allocate(ia(natoms))
        
        allocate(n2a(natoms))
        allocate(n2h(natoms))
        allocate(n2c(natoms))
        allocate(n2n(natoms))
        allocate(n2o(natoms))

        allocate(n3(natoms))
        allocate(n3a(natoms))
        allocate(n3h(natoms))
        allocate(n3c(natoms))
        allocate(n3n(natoms))
        allocate(n3o(natoms))

        allocate(n4(natoms))
        allocate(n4a(natoms))
        allocate(n4h(natoms))
        allocate(n4c(natoms))
        allocate(n4n(natoms))
        allocate(n4o(natoms))

        allocate(dm(natoms))

        
        do i=1,natoms
           nb(i)=0
           n2a(i)=0
           n2h(i)=0
           n2c(i)=0
           n2n(i)=0
           n2o(i)=0
           n3a(i)=0
           n3h(i)=0
           n3c(i)=0
           n3n(i)=0
           n3o(i)=0
           n4a(i)=0
           n4h(i)=0
           n4c(i)=0
           n4n(i)=0
           n4o(i)=0
           dm(i)=999.999d9
        enddo

        do i=1,natoms
           if(trim(name(i)).eq."H") ia(i)=1
           if(trim(name(i)).eq."C") ia(i)=6
           if(trim(name(i)).eq."N") ia(i)=7
           if(trim(name(i)).eq."O") ia(i)=8
           if(trim(name(i)).eq."F") ia(i)=9
           if(trim(name(i)).eq."Cl") ia(i)=17
           if(trim(name(i)).eq."Br") ia(i)=35
           if(trim(name(i)).eq."I") ia(i)=53
        enddo
        
        do i=1,natoms-1
           do j=i+1,natoms
              d=0.0d0
              do k=1,3
                 d=d+(coord(i,k)-coord(j,k))**2
              enddo
              d=dsqrt(d)
              if(d.lt.dm(i)) dm(i)=d
              if(d.lt.dm(j)) dm(j)=d
              dcut=1.6d0
              if(ia(i).gt.15.or.ia(j).gt.15) dcut=1.99d0
              if(d.lt.dcut) then
                 nb(i)=nb(i)+1
                 nb(j)=nb(j)+1
                 ib(i,nb(i))=j
                 ib(j,nb(j))=i
                 if(ia(i).eq.1) n2h(j)=n2h(j)+1
                 if(ia(j).eq.1) n2h(i)=n2h(i)+1
                 if(ia(i).eq.6) n2c(j)=n2c(j)+1
                 if(ia(j).eq.6) n2c(i)=n2c(i)+1
                 if(ia(i).eq.7) n2n(j)=n2n(j)+1
                 if(ia(j).eq.7) n2n(i)=n2n(i)+1
                 if(ia(i).eq.8) n2o(j)=n2o(j)+1
                 if(ia(j).eq.8) n2o(i)=n2o(i)+1
              endif
           enddo
        enddo

        do i=1,natoms
           do j=1,nb(i) ! loop over number of neighbors of i
              do k=1,nb(ib(i,j)) ! loop over number of second-neighbors of i
                 n3(i)=n3(i)+1
                 if(ia(ib(ib(i,j),k)).eq.1) n3h(i)=n3h(i)+1
                 if(ia(ib(ib(i,j),k)).eq.6) n3c(i)=n3c(i)+1
                 if(ia(ib(ib(i,j),k)).eq.7) n3n(i)=n3n(i)+1
                 if(ia(ib(ib(i,j),k)).eq.8) n3o(i)=n3o(i)+1
                 do m=1,nb(ib(ib(i,j),k)) ! loop over number of third-neighbors of i
                    n4(i)=n4(i)+1
                    if(ia(ib(ib(ib(i,j),k),m)).eq.1) n4h(i)=n4h(i)+1
                    if(ia(ib(ib(ib(i,j),k),m)).eq.6) n4c(i)=n4c(i)+1
                    if(ia(ib(ib(ib(i,j),k),m)).eq.7) n4n(i)=n4n(i)+1
                    if(ia(ib(ib(ib(i,j),k),m)).eq.8) n4o(i)=n4o(i)+1
                 enddo
              enddo
           enddo
        enddo
        
        do i=1,natoms
           id2=4*(4*(4*(4*n2a(i)+n2h(i))+n2c(i))+n2n(i))+n2o(i)
           id3=16*(16*(16*(16*n3a(i)+n3h(i))+n3c(i))+n3n(i))+n3o(i)
           id4=64*(64*(64*(64*n4a(i)+n4h(i))+n4c(i))+n4n(i))+n4o(i)
           id(i)=id4
        enddo
        
        do i=1,natoms
           if(nb(i).gt.1) then
              do j=1,nb(i)-1
                 do k=j+1,nb(i)
                    if(id(ib(i,j)).gt.id(ib(i,k))) then
                       m=ib(i,j)
                       ib(i,j)=ib(i,k)
                       ib(i,k)=m
                    endif
                 enddo
              enddo
           endif
        enddo

!        do i=1,natoms
!           write(*,'(i3,1x,a,t6,i10,3f12.6,6i4," : ",10i4)') i,trim(name(i)),id(i),coord(i,1),coord(i,2),coord(i,3), &
!                nb(i),n2a(i),n2h(i),n2c(i),n2n(i),n2o(i),(ib(i,m),m=1,nb(i))
!        enddo

!        do i=1,natoms
!           write(*,'(i5," : ",10i10)') i,(ib(i,m),m=1,nb(i)),(id(ib(i,m)),m=1,nb(i))
!        enddo
        
        deallocate(n2a,n2h,n2c,n2n,n2o)
        deallocate(n3,n3a,n3h,n3c,n3n,n3o)
        deallocate(n4,n4a,n4h,n4c,n4n,n4o)
        deallocate(dm)
        
        return

      end subroutine analyze

