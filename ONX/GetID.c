#include <stdio.h>
#include <stdlib.h>
#include <string.h>



void getid(int *mondoid)
{
    char* IDName;

    if((IDName = getenv("MONDO_ID")) == NULL) {
       printf(" The MONDO_ID environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }

    *mondoid = atoi(IDName);
}


void getid_(int *mondoid)
{
    getid(mondoid);
}
    
void GetID(int *mondoid)
{
    getid(mondoid);
}

