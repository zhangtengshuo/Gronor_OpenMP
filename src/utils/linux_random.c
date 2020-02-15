/*-----------------------------------------------------*\
 $Id: linux_random.c 19707 2010-10-29 17:59:36Z d3y133 $
\*-----------------------------------------------------*/
#include <stdlib.h>
void linux_sran_(int *input_seed)
{
  unsigned int seed;

  seed = (unsigned) *input_seed;
  (void) srandom(seed);
}
double linux_rand_(void)
{
  return (double) (((double) random())/(double) RAND_MAX);
}
