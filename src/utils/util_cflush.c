/*
 $Id: util_cflush.c 19707 2010-10-29 17:59:36Z d3y133 $
 */
#include <stdio.h>
/* this flushes stderr and stdout
   temporary hack for system where fortran flsuh is broken*/
int util_cflush_()
{
  fflush(stdout);
  fflush(stderr);
}
