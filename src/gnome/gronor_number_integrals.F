
      subroutine gronor_number_integrals(numone,numtwo)

      use inp
      use cidef
      use gnome_data
      use gnome_integrals
      use gnome_parameters

      implicit none

      integer :: mdum, numone, numtwo
      real (kind=8) :: sdum

      open(unit=lfnint,file=filone,form='unformatted',status='old',
     & err=1994)
      goto 1995
 1994 continue
      open(unit=lfnint,file=filint,form='unformatted',status='old',
     & err=994)
 1995 continue
      rewind(lfnint)
      if(ione.eq.0) then
       read(lfnint,err=993) namint,intone,potnuc,nbas,mbuf,mdum
      else
       read(lfnint,err=993) namint,nbas,intone,potnuc
      endif
      if(ione.eq.0) then
       read(lfnint,err=993) sdum
       read(lfnint,err=993) sdum
       read(lfnint,err=993) sdum
      else
       call gronor_read0(lfnint,intone)
       call gronor_read0(lfnint,intone)
       call gronor_read0(lfnint,intone)
       call gronor_read0(lfnint,intone)
      endif
      int2=0
      read(lfnint,end=993) nt

      numone=nt(1)
      numtwo=nt(2)

      rewind(unit=lfnint)
      close(unit=lfnint)

      return
 993  write(lfnout,983) filint(1:inp_strlen(filint))
      call errquit()
 994  write(lfnout,984) filint(1:inp_strlen(filint))
      call errquit()
 983  format('Error reading integral file ',a)
 984  format('Unable to open integral file ',a)
      end
