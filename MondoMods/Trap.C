/*
    ROUTINE TO GENERATE A TRACE BACK
    Author:  Matt Challacombe
*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
void trap_(void){abort();}
void trap(void){trap_();}
