/*
!--  This source code is part of the MondoSCF suite of 
!--  linear scaling electronic structure codes.  
!
!--  Matt Challacombe
!--  Los Alamos National Laboratory
!--  Copyright 1999, The University of California
!
!    ROUTINE TO GENERATE A TRACE BACK
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
void trap_(void){abort();}
void trap(void){trap_();}
