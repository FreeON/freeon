#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void HDump(int *Mesh, int *Switch, int *Grid)
{
    char* GammaHome;
    char FileName[100];
    char* Mode = "wb";
    FILE* InFile;

    if ((GammaHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in GDump.\n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }
    
    sprintf(FileName,"%s/MondoMods/MMA/Functions/Gamma_Header.bin",GammaHome);

    if((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create a binary Gamma header table.\n");
       exit( EXIT_FAILURE );
    }

     fwrite(Mesh,sizeof(int),1,InFile);
     fwrite(Switch,sizeof(double),1,InFile);
     fwrite(Grid,sizeof(double),1,InFile);
     fclose(InFile);
}

void hdump_(int *Mesh, int *Switch, int *Grid)
{
    HDump(Mesh,Switch,Grid);
}

void hdump(int *Mesh, int *Switch, int *Grid)
{
    HDump(Mesh,Switch,Grid);
}

