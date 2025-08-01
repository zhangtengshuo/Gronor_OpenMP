      double precision function util_cpusec()
      implicit none
! $Id: util_cpusec.F 19707 2010-10-29 17:59:36Z d3y133 $
!
! routine to give cpu seconds since start of execution
!  delta times are handled by all applications
!
#if defined(IBM) || defined(SP1) ||defined(HPUX)
#define __UTIL_CPUSEC__
      double precision ibm_cputime
      external ibm_cputime
      util_cpusec = ibm_cputime()
#endif
#if defined(KSR)
#define __UTIL_CPUSEC__
!--KSR-start
      real ksr_cputime
      external ksr_cputime
!
      util_cpusec = ksr_cputime()
!--KSR-end
#endif
#if defined(SUN) || defined(__PARAGON__) || defined(DECOSF) || defined(SOLARIS) || defined(FUJITSU_VPP)|| defined(LINUXALPHA) || defined(MACX)
#define __UTIL_CPUSEC__
      real*4 tarray(2)
      real*4 etime
#ifndef GCC4
      external etime
#endif
!
      util_cpusec = etime(tarray)
!
#endif

#if defined(SGI)
#define __UTIL_CPUSEC__
#if defined(SGITFP)
      real*4 tarray(2)        ! nb: real == real*8 on TFP
      real*4 etime
      external etime
#else
      real tarray(2)
      real etime
      external etime
#endif
      util_cpusec = etime(tarray)
#endif

#if defined(DELTA)
#define __UTIL_CPUSEC__
      double precision dclock
      external dclock
      double precision first
      logical ofirst
      data ofirst /.true./
!
      if (ofirst) then
!first = dclock()
!ofirst = .false.
      endif
!
      util_cpusec = dclock() - first
#endif

#if (defined(LINUX) || defined(CYGNUS)) && !defined(LINUXIA64) && !defined(CATAMOUNT) && !defined(IBM)
#define __UTIL_CPUSEC__
      double precision linux_cputime
      external linux_cputime
!
      util_cpusec = linux_cputime()
#endif

#if defined(WIN32) || defined(LINUXIA64)
#define __UTIL_CPUSEC__
      double precision win32_cputime
      external win32_cputime
!
      util_cpusec = win32_cputime()
#endif

!#if defined(CRAY)
!#define __UTIL_CPUSEC__
!      double precision secondr
!      external secondr
!      double precision first
!      logical ofirst
!      save first, ofirst
!      data ofirst /.true./
!c
!      if (ofirst) then
!	first = secondr()
!	ofirst = .false.
!      endif
!c
!      util_cpusec = secondr() - first ! Actually wall time
!#endif

#if defined(FUJITSU_SOLARIS) || defined(CATAMOUNT)
#define __UTIL_CPUSEC__
#ifdef FUJITSU_SOLARIS
      intrinsic cpu_time
#endif
!
      call cpu_time (util_cpusec)
#endif

!
#if !defined(__UTIL_CPUSEC__)
      util_cpusec = 0.0d00
#endif
      end


