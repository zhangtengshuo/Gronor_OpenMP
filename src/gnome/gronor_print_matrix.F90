subroutine gronor_print_matrix(lfno,lfna,lfnr,lfnt,header,key,olabel,labels,olower,scale, &
    rmat,nrdim,ncols,nds)

  implicit none

  integer, intent(in) :: lfno,lfna,lfnr,lfnt,nrdim,ncols,nds
  logical, intent(in) :: olabel,olower
  real(kind=8), intent(in) :: scale,rmat(nrdim,nrdim)
  character(len=128) :: header,key,labels(nrdim)

  real(kind=8), allocatable :: rt(:),rt0(:)
  integer :: i,j,k,nk,ii,il,ik,lenlab

  allocate(rt(nrdim))
  allocate(rt0(nrdim))

  write(lfno,600) trim(header)
600 format(//,1x,a)
  if(lfna.gt.0) write(lfna,601) trim(key(1:5)),nrdim,ncols
  if(lfnr.gt.0) write(lfnr,601) trim(key(1:5)),nrdim,ncols
601 format(a5,2i10)
  if(lfnt.gt.0) write(lfnt,600) trim(header)
  nk=nrdim/ncols
  if(mod(nrdim,ncols).ne.0) nk=nk+1
  lenlab=0
  do i=1,nrdim
    lenlab=max(lenlab,len(trim(labels(i))))
  enddo
  do k=1,nk
    ii=(k-1)*ncols+1
    il=min(nrdim,k*ncols)
    if(.not.olabel) then
      write(lfno,602) (i,i=ii,il)
602   format(/,6x,'|',7(6x,i8,6x))
      write(lfno,603) ' -----|',('--------------------',i=ii,il)
603   format(10a)
    else
      if(lenlab.le.10) then
        write(lfno,604) '|',(trim(labels(i)),i=ii,il)
        write(lfno,603) ' -----------|',('--------------------',i=ii,il)
      elseif(lenlab.le.20) then
        write(lfno,1604) '|',(trim(labels(i)),i=ii,il)
        write(lfno,603) ' ---------------------|',('--------------------',i=ii,il)
      else
        write(lfno,2604) '|',(trim(labels(i)),i=ii,il)
        write(lfno,603) ' -------------------------------|',('--------------------',i=ii,il)
      endif
604   format(/,12x,a,t22,a,t42,a,t62,a,t82,a,t102,a,t122,a,t142,a)
1604  format(/,22x,a,t32,a,t52,a,t72,a,t92,a,t112,a,t132,a,t152,a)
2604  format(/,32x,a,t42,a,t62,a,t82,a,t102,a,t122,a,t142,a,t162,a)
    endif
    if(lfna.gt.0) write(lfna,'(17x,10(i4,16x))') (i,i=ii,il)
    if(lfnr.gt.0) write(lfnr,'(17x,10(i4,16x))') (i,i=ii,il)
    if(lfnt.gt.0) write(lfnt,'(17x,10(i4,16x))') (i,i=ii,il)
    do j=1,nrdim
      ik=il
      if(olower) ik=min(j-1,il)
      if (.not.olabel) then
        write(lfno,605) j,(scale*rmat(i,j),i=ii,ik)
605     format(1x,i4,1x,'|',10f20.10)
      else
        if(lenlab.le.10) then
          write(lfno,606) trim(labels(j)),(scale*rmat(i,j),i=ii,ik)
        elseif(lenlab.le.20) then
          write(lfno,1606) trim(labels(j)),(scale*rmat(i,j),i=ii,ik)
        else
          write(lfno,2606) trim(labels(j)),(scale*rmat(i,j),i=ii,ik)
        endif
606     format(1x,a,t13,'|',t14,f20.10,t34,f20.10,t54,f20.10,t74,f20.10, &
            t94,f20.10,t114,f20.10,t134,f20.10)
1606    format(1x,a,t23,'|',t24,f20.10,t44,f20.10,t64,f20.10,t84,f20.10, &
            t104,f20.10,t124,f20.10,t144,f20.10)
2606    format(1x,a,t33,'|',t34,f20.10,t54,f20.10,t74,f20.10,t94,f20.10, &
            t114,f20.10,t134,f20.10,t154,f20.10)
      endif
      if(lfna.gt.0) write(lfna,607) j,(scale*rmat(i,j),i=ii,ik)
      if(lfnr.gt.0) write(lfnr,607) j,(scale*rmat(i,j),i=ii,ik)
607   format(i5,1x,10e20.13)
      do i=ii,il
        rt(i)=scale*rmat(i,j)
        if(dabs(rt(i)).lt.1.0d-08) rt(i)=0.0d0
      enddo
      if(lfnt.gt.0) then
        do i=ii,il
          rt0(i)=rt(i)
          if(abs(rt0(i)).lt.1.0d-3) rt0(i)=0.0d0
        enddo
        if(nds.eq.3) then
          write(lfnt,608) j,(rt0(i),i=ii,il)
        else
          write(lfnt,609) j,(rt0(i),i=ii,il)
        endif
      endif
608   format(i5,1x,10f20.3)
609   format(i5,1x,10f20.6)
    enddo
  enddo
  flush(lfno)
  if(lfna.gt.0) flush(lfna)
  if(lfnr.gt.0) flush(lfnr)
  if(lfnt.gt.0) flush(lfnt)

  deallocate(rt)
  deallocate(rt0)

  return

end subroutine gronor_print_matrix
