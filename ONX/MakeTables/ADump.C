#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void ADump(int *Ltot, int *GammaA)
{
    char* GammaHome;
    char FileName[100];
    char* Mode = "wb";
    FILE* InFile;

    if ((GammaHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in ADump.\n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }
    
    sprintf(FileName,"%s/MondoMods/MMA/Functions/Gamma_A.bin",GammaHome);

    if((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create a binary Gamma asymptotics table.\n");
       exit( EXIT_FAILURE );
    }

     fwrite(Ltot,sizeof(int),1,InFile);
     fwrite(GammaA,sizeof(double),(*Ltot)+1,InFile);
     fclose(InFile);
}

void adump_(int *Ltot, int *GammaA)
{
    ADump(Ltot,GammaA);
}

void adump(int *Ltot, int *GammaA)
{
    ADump(Ltot,GammaA);
}

