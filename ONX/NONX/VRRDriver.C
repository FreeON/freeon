#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void vconvert(int *VDrv, int *LDrv, int *LngDrv, int *LngLoc)
{
    char* DriverHome;
    char FileName[100];
    char* Mode = "r";
    FILE* InFileD;
    FILE* InFileL;
    int i;

    if ((DriverHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in cconvert.\n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }
    sprintf(FileName,"%s/ONX/Drivers/VRRDriver.ascii",DriverHome);
    if ((InFileD = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find the ONX VRR driver data file.\n");
       exit( EXIT_FAILURE );
    }
    sprintf(FileName,"%s/ONX/Drivers/VRRLoc.ascii",DriverHome);
    if ((InFileL = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find the ONX VRR address data file.\n");
       exit( EXIT_FAILURE );
    }

    fscanf(InFileL,"%d",LngDrv);
    fscanf(InFileL,"%d",LngLoc);

    for (i=0; i < (*LngDrv); ++i) {
       fscanf(InFileD,"%d",VDrv+i);
    }
    fclose(InFileD);

    for (i=0; i < (*LngLoc); ++i) {
       fscanf(InFileL,"%d",LDrv+i);
    }
    fclose(InFileL);

    sprintf(FileName,"%s/ONX/Drivers/VRRDriver.binary",DriverHome);
    Mode = "wb";
    if ((InFileD = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create the ONX VRR driver data file\n");
       exit( EXIT_FAILURE );
    }
    fwrite(LngDrv,sizeof(int),1,InFileD);
    fwrite(VDrv,sizeof(int),(*LngDrv),InFileD);
    fclose(InFileD);

    sprintf(FileName,"%s/ONX/Drivers/VRRLoc.binary",DriverHome);
    Mode = "wb";
    if ((InFileL = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create the ONX address data file\n");
       exit( EXIT_FAILURE );
    }
    fwrite(LngLoc,sizeof(int),1,InFileL);
    fwrite(LDrv,sizeof(int),(*LngLoc),InFileL);
    fclose(InFileL);
}


void vrrdriver(int *VDrv, int *LDrv, int *LngDrv, int *LngLoc)
{
    char* DriverHome;
    char* Mode = "rb";
    FILE* InFile;
    char FileName[100];

    if ((DriverHome = getenv("MONDO_HOME")) == NULL)
       exit( EXIT_FAILURE );
    sprintf(FileName,"%s/ONX/Drivers/VRRDriver.binary",DriverHome);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf("ONX binary VRR drivers do not exist.\n");
       printf("Creating the binary contraction driver.\n");
       vconvert(VDrv,LDrv,LngDrv,LngLoc); 
    }
    else {
       fread(LngDrv,sizeof(int),1,InFile);
       fread(VDrv,sizeof(int),(*LngDrv),InFile);
       fclose(InFile);
       sprintf(FileName,"%s/ONX/Drivers/VRRLoc.binary",DriverHome);
       InFile = fopen(FileName,Mode);
       fread(LngLoc,sizeof(int),1,InFile);
       fread(LDrv,sizeof(int),(*LngLoc),InFile);
       fclose(InFile);
    }
}

void vrrlng(int *LngDrv, int *LngLoc)
{
    char* DriverHome;
    char FileName[100];
    char* Mode = "r";
    FILE* InFile;

    if ((DriverHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in cconvert.\n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }
    sprintf(FileName,"%s/ONX/Drivers/VRRLoc.ascii",DriverHome);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find the ONX VRR address data file.\n");
       exit( EXIT_FAILURE );
    }
    fscanf(InFile,"%d",LngDrv);
    fscanf(InFile,"%d",LngLoc);
    fclose(InFile);
}

void vrrdriver_(int *VDrv, int *LDrv, int *LngDrv, int *LngLoc)
{
   vrrdriver(VDrv, LDrv, LngDrv, LngLoc);
}

void vrrlng_(int *LngDrv, int *LngLoc)
{
   vrrlng(LngDrv,LngLoc);
} 
