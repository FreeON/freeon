#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void gdloader(int *GDrv1, int *GDrv2, int *GDrv3, int *GDrv4, int *GDrv5,
              int *BC, int *KC, int *L1, int *L2, int *L3, int *L4, int *L5,
              int *La, int *Lc, int *Lb, int *Ln)
{
    char* MondoHome;
    char DriverHome[100];
    char FileName[100];
    char SB[4];
    char SK[4];
    char* Mode = "rb";
    FILE* InFile;
    int i;

    if((MondoHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in gdloader. \n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }
    sprintf(DriverHome,"%s/ONX/GDrivers",MondoHome);
    sprintf(SB,"%d",*BC);
    sprintf(SK,"%d",*KC);

/* G1 */

    sprintf(FileName,"%s/G1.%sx%s.binary",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find this gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }
    fread(L1,sizeof(int),1,InFile);
    fread(La,sizeof(int),1,InFile);
    fread(Lc,sizeof(int),1,InFile);
    fread(Lb,sizeof(int),1,InFile);
    fread(Ln,sizeof(int),1,InFile);
    fread(GDrv1,sizeof(int),(*L1),InFile);
    fclose(InFile);

/* G2 */

    sprintf(FileName,"%s/G2.%sx%s.binary",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find this gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }
    fread(L2,sizeof(int),1,InFile);
    fread(GDrv2,sizeof(int),(*L2)*4,InFile);

/* G3 */

    sprintf(FileName,"%s/G3.%sx%s.binary",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find this gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }
    fread(L3,sizeof(int),1,InFile);
    fread(GDrv3,sizeof(int),(*L3)*4,InFile);
 
/* G4 */

    sprintf(FileName,"%s/G4.%sx%s.binary",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find this gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }
    fread(L4,sizeof(int),1,InFile);
    fread(GDrv4,sizeof(int),(*L4)*4,InFile);

/* G5 */

    sprintf(FileName,"%s/G5.%sx%s.binary",DriverHome,SB,SK);
    if ((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Could not find this gradient contraction driver. \n");
       exit( EXIT_FAILURE );
    }
    fread(L5,sizeof(int),1,InFile);
    fread(GDrv5,sizeof(int),(*L5)*6,InFile);
}


void gdloader_(int *GDrv1, int *GDrv2, int *GDrv3, int *GDrv4, int *GDrv5,
              int *BC, int *KC, int *L1, int *L2, int *L3, int *L4, int *L5,
              int *La, int *Lc, int *Lb, int *Ln)
{
    gdloader(GDrv1, GDrv2, GDrv3, GDrv4, GDrv5,
             BC, KC, L1, L2, L3, L4, L5,
             La, Lc, Lb, Ln);
}
    

