#include <stdlib.h>

double random_(void) { 
  static long int i=0;
  if(i == 0) {
    i = (long int) time( NULL );
    srand48(i); 
  }
  return drand48();
}

