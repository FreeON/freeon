#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void EDump(int *Length, int *Table)
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
    
    sprintf(FileName,"%s/MondoMods/MMA/Functions/Exponential.bin",GammaHome);

    if((InFile = fopen(FileName,Mode)) == NULL) {
       printf(" %s\n",FileName);
       printf(" Error trying to create a binary Exp table.\n");
       exit( EXIT_FAILURE );
    }

     fwrite(Table,sizeof(double),(*Length),InFile);
     fclose(InFile);
}

void edump_(int *Length, int *Table)
{
    EDump(Length, Table);
}

void edump(int *Length, int *Table)
{
    EDump(Length, Table);
}

