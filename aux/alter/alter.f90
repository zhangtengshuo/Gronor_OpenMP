program alter

  character (len=255) :: card,filout,filalt
  real (kind=8) :: occ(1000)
  character (len=5) :: orb(1000)
  integer :: num, iarg, ndx, nocc, nalt, ndxz,ndxn, ncas
  integer :: iz(1000),inz(1000)
  integer :: ialt(12,2)

  iarg=1
  call getarg(iarg,filout)
  iarg=2
  call getarg(iarg,filalt)

  ncas=4
  
  do i=1,1000
     occ(i)=0.0
     orb(i)="     "
  enddo
  
  open(unit=10,file=trim(filout),form='formatted')

  card(1:2)="       "
  do while(card(1:7).ne."  Index")
     read(10,'(a)') card
  enddo

  do while(card(1:2).ne."--")
     read(10,'(a)') card
     if(card(1:2).ne."--") then
        read(card,'(i5)') ndx
        if(ndx.gt.0) then
           read(card,'(15x,f10.4,19x,a5)') occ(ndx),orb(ndx)
           num=ndx
        endif
     endif
  enddo

  close(unit=10)

  do i=1,num
     if(occ(i).gt.0.0) nocc=i
  enddo

  ndxz=0
  do i=nocc-ncas,1,-1
     if(orb(i)(2:3).eq.'pz') then
        ndxz=ndxz+1
        iz(ndxz)=i
     endif
  enddo
  ndx=0
  nalt=0
  do i=nocc,nocc+1-ncas,-1
     if(orb(i)(2:3).ne.'pz') then
        ndx=ndx+1
        nalt=nalt+1
        ialt(nalt,1)=iz(ndx)
        ialt(nalt,2)=i
     endif
  enddo
  
  ndxz=0
  do i=nocc+1+ncas,num
     if(orb(i)(2:3).eq.'pz') then
        ndxz=ndxz+1
        iz(ndxz)=i
     endif
  enddo
  ndx=0
  do i=nocc+1,nocc+ncas
     if(orb(i)(2:3).ne.'pz') then
        ndx=ndx+1
        nalt=nalt+1
        ialt(nalt,1)=i
        ialt(nalt,2)=iz(ndx)
     endif
  enddo

  open(unit=10,file=trim(filalt),form='formatted')
  
  write(10,'(i6)') nalt
  do i=1,nalt
     write(10,'(2i6)') ialt(i,1),ialt(i,2)
  enddo
  close(unit=10)
  
end program alter
