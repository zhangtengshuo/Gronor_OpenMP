

      subroutine sminop(ciin,ioccin,ciout,ioccout,idetin,nact,idetout,
     & maxn)
      implicit none
      integer :: idetin,nact,idetout,maxn
      real (kind=8) :: ciin(idetin)
      integer :: ioccin(maxn,idetin)
      real (kind=8) :: ciout(*)
      integer :: ioccout(maxn,*)
      integer, allocatable :: iopen(:),nocc(:)
      logical :: onew,ocomparocc
      integer :: k,i,nopen,j,jocc
      integer :: ioccb(32),ioccb2(32),iocca(32),na,nb,ib,ia
      integer :: nperm, jperm, ifind, index, ipar

      idetout=0
      allocate (iopen(nact),nocc(nact))

      do i=1,idetin
       nopen=0
       na=0
       nb=0
       ib=0
       ia=0
       do j=1,nact
        if(ioccin(j,i).eq.2) then
         na=na+1
         nb=nb+1
         ia=ia+1
         iocca(ia)=j
         ib=ib+1
         ioccb(ib)=j
        else if(ioccin(j,i).eq.1) then
         nopen=nopen+1
         iopen(nopen)=j
         na=na+1
         ia=ia+1
         iocca(ia)=j
        else if(ioccin(j,i).eq.-1) then
         nb=nb+1
         ib=ib+1
         ioccb(ib)=j
        endif
       enddo
       do j=1,nopen
        do k=1,nact
         nocc(k)=ioccin(k,i)
        enddo
        do k=1,nb
         ioccb2(k+1)=ioccb(k)
        enddo
        nocc(iopen(j))=-1
        index=ifind(iopen(j),iocca,na)
        nperm=index-1+na-1
        ioccb2(1)=iopen(j)
        call bubble_sort(ioccb2,nb+1,jperm)
        nperm=nperm+jperm
        ipar=(-1)**nperm
         idetout=idetout+1
         ciout(idetout)=ipar*ciin(i)
         do k=1,nact
          ioccout(k,idetout)=nocc(k)
         enddo
       enddo
      enddo

      deallocate(iopen,nocc)

      return
      end

      logical function ocomparocc(ia,ib,nact)
      implicit none
      integer :: ia(nact),ib(nact),nact
      integer :: i
      ocomparocc=.true.

      do i=1,nact
       if(ia(i).ne.ib(i))ocomparocc=.false.
      enddo

      return
      end
