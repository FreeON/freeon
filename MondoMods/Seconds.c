/*
!
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!    CPU AND WALL TIMERS BASED ON C INTRINSICS 
!
*/
#include <time.h>
#include <sys/times.h>
double cpu_seconds_(void)
{
 double CPS=1.0/CLOCKS_PER_SEC;
 double CLOCKS=clock();
 return CLOCKS*CPS;
}

double cpu_seconds__(void)
{
 double CPS=1.0/CLOCKS_PER_SEC;
 double CLOCKS=clock();
 return CLOCKS*CPS;
}
double wall_seconds_(void)
{
 struct tms tbuff;
 double CTK=1.0/CLK_TCK;
 double WALL=times(&tbuff);
 return WALL*CTK;
}

double wall_seconds__(void)
{
 struct tms tbuff;
 double CTK=1.0/CLK_TCK;
 double WALL=times(&tbuff);
 return WALL*CTK;
}
