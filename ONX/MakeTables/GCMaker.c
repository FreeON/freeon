#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void gcmaker(int *GCDrv, int *BC, int *KC, int *L1, int *L2,
             int *L3, int *L4, int *L5)
{
    char* DriverHome;
    char FileName[100];
    char SB[4];
    char SK[4];
    char* Mode = "r";
    FILE* InFile;
    int i;

    if((DriverHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in gcdriver. \n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }
    sprintf(SB,"%d",*BC);
    sprintf(SK,"%d",*KC);

/* G1 */

    sprintf(FileName,"%s/MondoMods/MMA/Functions/G1.%sx%s.ascii",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find the gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }

    fscanf(InFile,"%d",L1);
    fscanf(InFile,"%d",L2);
    fscanf(InFile,"%d",L3);
    fscanf(InFile,"%d",L4);
    fscanf(InFile,"%d",L5);
    for (i=0; i < (*L1); ++i) {
       fscanf(InFile,"%d",GCDrv+i);
    }
    fclose(InFile);
    sprintf(FileName,"%s/MondoMods/MMA/Functions/G1.%sx%s.bin",DriverHome,SB,SK);
    Mode="wb";
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create the gradient contraction driver\n");
       exit( EXIT_FAILURE );
    }
    fwrite(L1,sizeof(int),1,InFile);
    fwrite(L2,sizeof(int),1,InFile);
    fwrite(L3,sizeof(int),1,InFile);
    fwrite(L4,sizeof(int),1,InFile);
    fwrite(L5,sizeof(int),1,InFile);
    fwrite(GCDrv,sizeof(int),(*L1),InFile);
    fclose(InFile);

/* G2 */

    Mode="r";
    sprintf(FileName,"%s/MondoMods/MMA/Functions/G2.%sx%s.ascii",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find the gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }

    fscanf(InFile,"%d",L1);
    for (i=0; i < (*L1)*4; ++i) {
       fscanf(InFile,"%d",GCDrv+i);
    }
    fclose(InFile);
    sprintf(FileName,"%s/MondoMods/MMA/Functions/G2.%sx%s.bin",DriverHome,SB,SK);
    Mode="wb";
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create the gradient contraction driver\n");
       exit( EXIT_FAILURE );
    }
    fwrite(L1,sizeof(int),1,InFile);
    fwrite(GCDrv,sizeof(int),(*L1)*4,InFile);
    fclose(InFile);

/* G3 */

    Mode="r";
    sprintf(FileName,"%s/MondoMods/MMA/Functions/G3.%sx%s.ascii",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find the gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }

    fscanf(InFile,"%d",L1);
    for (i=0; i < (*L1)*4; ++i) {
       fscanf(InFile,"%d",GCDrv+i);
    }
    fclose(InFile);
    sprintf(FileName,"%s/MondoMods/MMA/Functions/G3.%sx%s.bin",DriverHome,SB,SK);
    Mode="wb";
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create the gradient contraction driver\n");
       exit( EXIT_FAILURE );
    }
    fwrite(L1,sizeof(int),1,InFile);
    fwrite(GCDrv,sizeof(int),(*L1)*4,InFile);
    fclose(InFile);

/* G4 */

    Mode="r";
    sprintf(FileName,"%s/MondoMods/MMA/Functions/G4.%sx%s.ascii",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find the gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }

    fscanf(InFile,"%d",L1);
    for (i=0; i < (*L1)*4; ++i) {
       fscanf(InFile,"%d",GCDrv+i);
    }
    fclose(InFile);
    sprintf(FileName,"%s/MondoMods/MMA/Functions/G4.%sx%s.bin",DriverHome,SB,SK);
    Mode="wb";
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create the gradient contraction driver\n");
       exit( EXIT_FAILURE );
    }
    fwrite(L1,sizeof(int),1,InFile);
    fwrite(GCDrv,sizeof(int),(*L1)*4,InFile);
    fclose(InFile);

/* G5 */

    Mode="r";
    sprintf(FileName,"%s/MondoMods/MMA/Functions/G5.%sx%s.ascii",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find the gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }

    fscanf(InFile,"%d",L1);
    for (i=0; i < (*L1)*6; ++i) {
       fscanf(InFile,"%d",GCDrv+i);
    }
    fclose(InFile);
    sprintf(FileName,"%s/MondoMods/MMA/Functions/G5.%sx%s.bin",DriverHome,SB,SK);
    Mode="wb";
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create the gradient contraction driver\n");
       exit( EXIT_FAILURE );
    }
    fwrite(L1,sizeof(int),1,InFile);
    fwrite(GCDrv,sizeof(int),(*L1)*6,InFile);
    fclose(InFile);






}


void gcmaker_(int *GCDrv,int *BC, int *KC, int *L1,int *L2,
               int *L3, int *L4, int *L5)
{
    gcmaker(GCDrv,BC,KC,L1,L2,L3,L4,L5);
}
    

