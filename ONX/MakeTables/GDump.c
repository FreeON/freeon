#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void GDump(int *L, int *Length, int *Table)
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

     fwrite(Table,sizeof(double),(*Length),InFile);
     fclose(InFile);
}

void gdump_(int *L, int *Length, int *Table)
{
    GDump(L, Length, Table);
}

void gdump(int *L, int *Length, int *Table)
{
    GDump(L, Length, Table);
}

