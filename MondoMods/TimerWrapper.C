#ifdef IPM
#include "IPM_timer.h"
double mondotimer() { return IPM_timer_get_time();}
double mondotimer_() { return IPM_timer_get_time();}
double mondotimer__() { return IPM_timer_get_time();}
#else
#ifdef PARALLEL
#include "mpi.h"
double mondotimer() {return MPI_Wtime();}
double mondotimer_() {return MPI_Wtime();}
double mondotimer__() {return MPI_Wtime();}
#else
double mondotimer() {return wall_seconds_();}
double mondotimer_() {return wall_seconds_();}
double mondotimer__() {return wall_seconds_();}
#endif 
#endif
