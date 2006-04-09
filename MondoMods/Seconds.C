/*
  Author: Matt Challacombe
  CPU AND WALL TIMERS BASED ON C INTRINSICS 
*/
#include <time.h>
#include <sys/times.h>
#include <unistd.h>
double cpu_seconds_(void)
{
 double CPS=1.0/CLOCKS_PER_SEC;
 double CLOCKS=clock();
 return CLOCKS*CPS;
}
double cpu_seconds(void){return cpu_seconds_();}
double cpu_seconds__(void){return cpu_seconds_();}

double wall_seconds_(void)
{
 struct tms tbuff;
 double CTK=1.0/sysconf(_SC_CLK_TCK);
 /* double CTK=1.0/CLK_TCK; This is obsolete. */
 double WALL=times(&tbuff);
 return WALL*CTK;
}
double wall_seconds(void){return wall_seconds_();}
double wall_seconds__(void){return wall_seconds_();}
