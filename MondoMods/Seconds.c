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
 double CLOCKS=clock();
 double CPS=CLOCKS_PER_SEC;
 return CLOCKS/CPS;
}
double cpu_seconds(void){return cpu_seconds_();}
double cpu_seconds__(void){return cpu_seconds_();}

double wall_seconds_(void)
{
 struct tms tbuff;
 double WALL=times(&tbuff);
 double CTK=CLK_TCK;
 return WALL/CTK;
}
double wall_seconds(void){return wall_seconds_();}
double wall_seconds__(void){return wall_seconds_();}
