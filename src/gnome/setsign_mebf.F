
      integer function setsign_mebf(iocc,nact1,nact2,nact)
!     Sets the sign the MEBF
      implicit none
      integer :: nact1,nact2,nact
      integer :: iocc1(nact1),iocc2(nact2)
      integer :: iocc1aux(nact1),iocc2aux(nact2)
      integer :: iocc(nact),ioccaux(nact)
      integer :: i,j,nperm
      
      setsign_mebf=0
      j=0
      nperm=0
      do i=1,nact1
         j=j+1
         iocc1(i)=iocc(j)
      enddo
      do i=1,nact2
         j=j+1
         iocc2(i)=iocc(j)
      enddo      
      j=0
      do i=1,nact1
         if(iocc1(i).eq.2)then
            j=j+1
            iocc1aux(j)=1
            j=j+1
            iocc1aux(j)=-1
         elseif(iocc1(i).eq.1)then
            j=j+1
            iocc1aux(j)=1
         elseif(iocc1(i).eq.-1)then
            j=j+1
            iocc1aux(j)=-1
         endif
      enddo
      j=0
      do i=1,nact2
         if(iocc2(i).eq.2)then
            j=j+1
            iocc2aux(j)=1
            j=j+1
            iocc2aux(j)=-1
         elseif(iocc2(i).eq.1)then
            j=j+1
            iocc2aux(j)=1
         elseif(iocc2(i).eq.-1)then
            j=j+1
            iocc2aux(j)=-1
         endif
      enddo
!     order iocc1 and iocc2 separately
      call bubble_sort(iocc1aux,nact1,nperm)
      call bubble_sort(iocc2aux,nact2,nperm)                        
      j=0
      do i=1,nact1
         j=j+1
         ioccaux(j)=iocc1aux(i)
      enddo
      do i=1,nact2
         j=j+1
         ioccaux(j)=iocc2aux(i)
      enddo     
      call bubble_sort(ioccaux,nact,nperm)
      if(mod(nperm,2).eq.1)then
         setsign_mebf=-1
      else
         setsign_mebf=1
      endif
      return
      end
