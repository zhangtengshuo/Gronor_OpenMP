      subroutine gronor_print_matrix                                      &
     &     (lfn1,lfn2,lfn3,header,key,olabel,labels,olower,scale,         &
     &     rmat,nrdim,ncols)

      implicit none
      
      integer, intent(in) :: lfn1,lfn2,lfn3,nrdim,ncols
      logical, intent(in) :: olabel,olower
      real(kind=8), intent(in) :: scale,rmat(nrdim,nrdim)
      character(len=128) :: header,key,labels(nrdim)

      real(kind=8), allocatable :: rt(:)
      integer :: i,j,k,nk,ii,il,ik

      allocate(rt(nrdim))

      write(lfn1,600) trim(header)
 600  format(//,1x,a)
      if(lfn2.gt.0) write(lfn2,601) trim(key(1:5)),nrdim,ncols
 601  format(a5,2i10)
      if(lfn3.gt.0) write(lfn3,600) trim(header)
      nk=nrdim/ncols
      if(mod(nrdim,ncols).ne.0) nk=nk+1
      do k=1,nk
        ii=(k-1)*ncols+1
        il=min(nrdim,k*ncols)
        if(.not.olabel) then
          write(lfn1,602) (i,i=ii,il)
 602      format(/,6x,'|',7(6x,i8,6x))
          write(lfn1,603) ' -----|',                                      &
     &         ('--------------------',i=ii,il)
 603      format(10a)
        else
          write(lfn1,604) '|',(trim(labels(i)),i=ii,il)
 604      format(/,12x,a,t22,a,t42,a,t62,a,t82,a,t102,a,t122,a,t142,a)
          write(lfn1,603) ' -----------|',                                &
     &         ('--------------------',i=ii,il)
        endif
        if(lfn2.gt.0) write(lfn2,603) (i,i=ii,il)
        if(lfn3.gt.0) write(lfn3,603) (i,i=ii,il)
        do j=1,nrdim
          ik=il
          if(olower) ik=min(j-1,il)
          if (.not.olabel) then
            write(lfn1,605) j,(scale*rmat(i,j),i=ii,ik)
 605        format(1x,i4,1x,'|',10f20.10)
          else
            write(lfn1,606) trim(labels(j)),(scale*rmat(i,j),i=ii,ik)
 606        format(1x,a,t13,'|',                                          &
     &           t14,f20.10,t34,f20.10,t54,f20.10,t74,f20.10,             &
     &           t94,f20.10,t114,f20.10,t134,f20.10)
          endif
          if(lfn2.gt.0) write(lfn2,607) j,(scale*rmat(i,j),i=ii,ik)
 607      format(i5,1x,10e20.13)
        enddo
        do i=ii,il
          rt(i)=scale*rmat(i,j)
          if(dabs(rt(i)).lt.1.0d-08) rt(i)=0.0d0
        enddo
          if(lfn3.gt.0) write(lfn3,608) j,(rt(i),i=ii,il)
 608      format(i5,1x,10f20.6)
      enddo
      flush(lfn1)
      if(lfn2.gt.0) flush(lfn2)
      if(lfn3.gt.0) flush(lfn3)

      deallocate(rt)
      
      return

      end subroutine gronor_print_matrix
