#include <stdlib.h>
#include <stdio.h>
#include <string.h>

typedef long logical;		/* Equivalent C type to FORTRAN logical */
#ifdef EXT_INT
typedef long integer;		/* Equivalent C type to FORTRAN integer */
#else
typedef int integer;		/* Equivalent C type to FORTRAN integer */
#endif
#define FORTRAN_TRUE  ((logical) 1)
#define FORTRAN_FALSE ((logical) 0)

#ifndef WIN32
#define FATR 
#endif

#define MAX_CLEN 4096
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
static int fortchar_to_string(_fcd f, int flen, char *buf, 
			      const int buflen)
#else
static int fortchar_to_string(const char *f, int flen, char *buf, 
			      const int buflen)
#endif
{
#if (defined(CRAY) || defined(USE_FCD))&& !defined(__crayx1)
  char *fstring = _fcdtocp(f);
  flen = _fcdlen(f);

  while (flen-- && fstring[flen] == ' ')
    ;

  if ((flen+1) >= buflen)
    return 0;			/* Won't fit */

  flen++;
  buf[flen] = 0;
  while(flen--)
    buf[flen] = fstring[flen];

#else

  while (flen-- && f[flen] == ' ')
    ;

  if ((flen+1) >= buflen)
    return 0;			/* Won't fit */

  flen++;
  buf[flen] = 0;
  while(flen--)
    buf[flen] = f[flen];
#endif
  return 1;
}
static int string_to_fortchar( f, flen, buf)
  int flen;
  char *buf;
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  _fcd f;
#else
  char *f;
#endif
{
  int len = strlen(buf), i;
#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  flen = _fcdlen(f);
#endif

  if (len > flen) 
    return 0;			/* Won't fit */

#if (defined(CRAY) || defined(USE_FCD)) && !defined(__crayx1)
  for (i=0; i<len; i++)
    _fcdtocp(f)[i] = buf[i];
  for (i=len; i<flen; i++)
    _fcdtocp(f)[i] = ' ';
#else
  for (i=0; i<len; i++)
    f[i] = buf[i];
  for (i=len; i<flen; i++)
    f[i] = ' ';
#endif
  return 1;
}
