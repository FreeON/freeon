#ifdef IPM_TIMER

#include "IPM_timer.h"
double mondo_timer() { return IPM_timer_get_time();}
double mondo_timer_() { return IPM_timer_get_time();}
double mondo_timer__() { return IPM_timer_get_time();}

#else

#ifdef PARALLEL
#include "mpi.h"

double mondo_timer() {return MPI_Wtime();}
double mondo_timer_() {return MPI_Wtime();}
double mondo_timer__() {return MPI_Wtime();}

#else 
double mondo_timer() {return wall_seconds_();}
double mondo_timer_() {return wall_seconds_();}
double mondo_timer__() {return wall_seconds_();}

#endif
#endif
