#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void GammaHeader(int *Mesh, int *Switch, int *Grid)
{
    char* GammaHome;
    char FileName[100];
    char* Mode = "rb";
    FILE* InFile;

    if ((GammaHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in GDump.\n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }

    sprintf(FileName,"%s/MondoMods/MMA/Functions/Gamma_Header.bin",GammaHome);

    if((InFile = fopen(FileName,Mode)) == NULL){
       printf("Binary gamma header table does not exist.\n");
       exit( EXIT_FAILURE );
    }

    fread(Mesh,sizeof(int),1,InFile);
    fread(Switch,sizeof(double),1,InFile);
    fread(Grid,sizeof(double),1,InFile);
}

void gammaheader_(int *Mesh, int *Switch, int *Grid)
{
   GammaHeader(Mesh,Switch,Grid);
}

void gammaheader(int *Mesh, int *Switch, int *Grid)
{
   GammaHeader(Mesh,Switch,Grid);
}

