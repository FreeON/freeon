#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void GetGammaTable(int *L, int *Length, int *Table)
{
    char* GammaHome;
    char SL[3];
    char FileName[100];
    char* Mode = "rb";
    FILE* InFile;

    if ((GammaHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in GDump.\n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }

    sprintf(SL,"%d",*L);
    sprintf(FileName,"%s/MondoMods/MMA/Functions/Gamma_%s.bin",GammaHome,SL);

    if((InFile = fopen(FileName,Mode)) == NULL){
       printf("Binary Gamma tables do not exist.\n");
       exit( EXIT_FAILURE );
    }
    fread(Table,sizeof(double),(*Length),InFile);
    fclose(InFile);
}

void getgammatable_(int *L, int *Length, int *Table)
{
   GetGammaTable(L, Length, Table);
}

void getgammatable(int *L, int *Length, int *Table)
{
   GetGammaTable(L, Length, Table);
}

