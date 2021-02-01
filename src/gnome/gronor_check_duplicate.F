      subroutine check_duplicate(list,det,korb,nact,ndetlist,dup,i_dup,
     &     lfnout)

      implicit none
      
      integer :: nact,ndetlist,i,j,korb
      integer :: list(nact,ndetlist),det(nact),newdet(nact)
      integer :: i_dup,lfnout
      logical :: dup
      do i=1,nact
         newdet(i)=det(i)
      enddo
      newdet(korb)=-1
      dup=.false.
      do i=1,ndetlist         
         if(all(newdet.eq.list(:,i)))then
            dup=.true.
            i_dup=i
            goto 99
         else
            dup=.false.
         endif
      enddo
 99   continue
      return
      end
