#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void cconvert(int *CDrv, int *LngDrv)
{
    char* DriverHome;
    char FileName[100];
    char* Mode = "r";
    FILE* InFile;
    int i;

    if ((DriverHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in cconvert.\n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }
    sprintf(FileName,"%s/MakeK/Drivers/CDriver.ascii",DriverHome);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find the ONX contraction driver data file.\n");
       exit( EXIT_FAILURE );
    }
    fscanf(InFile,"%d",LngDrv);
    for (i=0; i < (*LngDrv); ++i) {
       fscanf(InFile,"%d",CDrv+i);
    }
    fclose(InFile);
    sprintf(FileName,"%s/MakeK/Drivers/CDriver.binary",DriverHome);
    Mode = "wb";
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create the ONX contraction driver data file\n");
       exit( EXIT_FAILURE );
    }
    fwrite(LngDrv,sizeof(int),1,InFile);
    fwrite(CDrv,sizeof(int),(*LngDrv),InFile);
    fclose(InFile);
}


void ccdriver(int *CDrv, int *LngDrv)
{
    char* DriverHome;
    char* Mode = "rb";
    FILE* InFile;
    char FileName[100];

    if ((DriverHome = getenv("MONDO_HOME")) == NULL)
       exit( EXIT_FAILURE );
    sprintf(FileName,"%s/MakeK/Drivers/CDriver.binary",DriverHome);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf("ONX binary contraction drivers do not exist.\n");
       printf("Creating the binary contraction driver.\n");
       cconvert(CDrv,LngDrv); 
    }
    else {
       fread(LngDrv,sizeof(int),1,InFile);
       fread(CDrv,sizeof(int),(*LngDrv),InFile);
       fclose(InFile);
    }
}

void ccdriver_(int *CDrv, int *LngDrv)
{
   ccdriver(CDrv, LngDrv);
}

 
