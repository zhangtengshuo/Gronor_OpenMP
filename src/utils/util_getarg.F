      subroutine util_getarg(iarg_in,arg)
* $Id: util_getarg.F 19707 2010-10-29 17:59:36Z d3y133 $
*
* The statement call util_getarg( k , arg ) gets the  kth  command-
* line argument and puts it into the string arg.
*
      implicit none
c#include "errquit.fh"
      integer iarg_in
      character*(*) arg
#if defined(EXT_INT) && !defined(GCC4)
      integer*4 iarg
#else
      integer iarg
#if defined(WIN32)
      integer idum
#endif
#endif
#if defined(CRAY)
      integer*4 length,ierror
#endif
#if defined(HPUX)
#ifdef EXT_INT
!$HP$ ALIAS getarg='getarg__'
#else
!$HP$ ALIAS getarg='getarg_'
#endif
#endif
*
      iarg=iarg_in
#if defined(CRAY)
      call pxfgetarg(iarg,arg,length,ierror)
c
c      if (ierror.ne.0)
c     + call errquit('util_getarg:pxfgetarg fatal error',ierror,
c     &       UNKNOWN_ERR)
#elif defined(WIN32)
      call getarg(iarg,arg,idum)
#else
      call getarg(iarg,arg)
#endif
*
      return
      end
