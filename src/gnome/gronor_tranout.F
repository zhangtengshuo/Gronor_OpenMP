

      subroutine gronor_tranout(lfndbg,idet)
      use cidist
      use gnome_parameters
      use gnome_data
      implicit none
      integer :: idet, lfndbg
      integer :: ivc,ntvc,ibas,i

      ntvc=ntcl(idet)+ntop(idet)

      if(idbg.ge.15)
     &  write(lfndbg,601)  ntcl(idet),ntop(idet),
     &                     (ioccn(i,idet),i=1,ntop(idet))
 601  format(/,' Tranout: nclose:',i3,',nopen:',i3,',ocopen:',20i4)
      if(idbg.ge.13) write(lfndbg,602)
 602  format(/,' M.O.''s transformed and ordered')

      if(idbg.gt.13) then
       write(lfndbg,603)
 603   format(/,' closed shell M.O.''s',/)
       do ivc=1,ntvc
        if(ivc.eq.ntcl(idet)+1) write(lfndbg,604)
 604    format(/,' open shell M.O.'' s:')
        write(lfndbg,605)  ' (',ivc,')',(vec(ivc,ibas,idet),ibas=1,nbas)
 605    format(a2,i3,a1,(t9,10f12.8))
       enddo
      endif

      return
      end
