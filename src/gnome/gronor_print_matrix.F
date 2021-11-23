      subroutine gronor_print_matrix(lfn,r,nrdim)
      
      use gnome_parameters, only : ncols
      
      implicit none
      integer, intent(in) :: lfn,nrdim
      real(kind=8), intent(in) :: r(nrdim,nrdim)
      integer :: i,j,k,nk,ii,il

      write(lfn,*) ' '
      nk=nrdim/ncols
      if(mod(nrdim,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nrdim,k*ncols)
        write(lfn,600) (i,i=ii,il)
 600    format(6x,7(6x,i8,6x))
 601    format(i5,1x,10f20.10)
        do j=1,nrdim
          write(lfn,601) j,(r(i,j),i=ii,il)
        enddo
        write(lfn,*) ' '
      enddo
      flush(lfn)

      return

      end subroutine gronor_print_matrix
