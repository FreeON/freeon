/*
  Author: Matt Challacombe
  CPU AND WALL TIMERS BASED ON C INTRINSICS 
*/
#include <time.h>
#include <sys/times.h>
double cpu_seconds_(void)
{
 double CPS=1.0/CLOCKS_PER_SEC;
 double CLOCKS=clock();
 return CLOCKS*CPS;
}
double cpu_seconds(void){return cpu_seconds_();}

double wall_seconds_(void)
{
 struct tms tbuff;
 double CTK=1.0/CLK_TCK;
 double WALL=times(&tbuff);
 return WALL*CTK;
}
double wall_seconds(void){return wall_seconds_();}
