#include <time.h>
#include <unistd.h>

#include "xorgens.h"
#include "xorgens.c"

double random_double (void) {
  return xor4096r(0);
}

unsigned long seed_rng (void) {  
  /* pick "random" seed based on time */
  unsigned long rseed;
  long pid, tim;
  pid = getpid();
  pid = pid*65538; /* so it fills in bits 2 through 31 */
  tim = time(0);
  rseed = pid^tim;
  xor4096r(rseed);
  return rseed;
}

void set_rng_seed(unsigned long newseed) {
  xor4096r(newseed);
}
