#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void GDump(int *L, int *Length, int *Asymp, int *Table)
{
    char* GammaHome;
    char SL[3];
    char FileName[100];
    char* Mode = "wb";
    FILE* InFile;

    if ((GammaHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in GDump.\n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }
    
    sprintf(SL,"%d",*L);
    sprintf(FileName,"%s/MondoMods/MMA/Functions/Gamma_%s.bin",GammaHome,SL);

    if((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create a binary Gamma table.\n");
       exit( EXIT_FAILURE );
    }

     fwrite(Asymp,sizeof(double),1,InFile);
     fwrite(Table,sizeof(double),(*Length),InFile);
     fclose(InFile);
}

void gdump_(int *L, int *Length, int *Asymp, int *Table)
{
    GDump(L, Length, Asymp, Table);
}

void gdump(int *L, int *Length, int *Asymp, int *Table)
{
    GDump(L, Length, Asymp, Table);
}

