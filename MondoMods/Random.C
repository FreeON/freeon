#include <stdlib.h>
/*
double random_(void){return rand(); }
*/
double random_(void){ 
srand( (unsigned int) time( NULL ));
return drand48(); }
