#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void GammaTable(int *L, int *Length, int *Table)
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
    fread(Table,sizeof(double),4*(*Length+1),InFile);
    fclose(InFile);
}

void gammatable_(int *L, int *Length, int *Table)
{
   GammaTable(L, Length, Table);
}

void gammatable(int *L, int *Length, int *Table)
{
   GammaTable(L, Length, Table);
}

/*-----------------------------------------------------------------------*/

void ExpTable(int *Length, int *Table)
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

    sprintf(FileName,"%s/MondoMods/MMA/Functions/Exponential.bin",GammaHome);

    if((InFile = fopen(FileName,Mode)) == NULL){
       printf("Binary Gamma tables do not exist.\n");
       exit( EXIT_FAILURE );
    }
    fread(Table,sizeof(double),4*(*Length+1),InFile);
    fclose(InFile);
}

void exptable_(int *Length, int *Table)
{
   ExpTable(Length, Table);
}

void exptable(int *Length, int *Table)
{
   ExpTable(Length, Table);
}

/*-----------------------------------------------------------------------*/

void GammaHeader(int *Mesh, int *Switch, int *Grid)
{
    char* GammaHome;
    char FileName[100];
    char* Mode = "rb";
    FILE* InFile;

    if ((GammaHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in GammaHeader.\n");
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

/*-----------------------------------------------------------------------*/


void GammaAsymptotics(int *Lr, int *GammaA)
{
    char* GammaHome;
    char FileName[100];
    char* Mode = "rb";
    FILE* InFile;

    if ((GammaHome = getenv("MONDO_HOME")) == NULL) {
       printf(" Error in GammaAsymptotics.\n");
       printf(" The MONDO_HOME environment variable is not set.\n");
       exit( EXIT_FAILURE );
    }

    sprintf(FileName,"%s/MondoMods/MMA/Functions/Gamma_A.bin",GammaHome);

    if((InFile = fopen(FileName,Mode)) == NULL){
       printf("Binary gamma asymptotics table does not exist.\n");
       exit( EXIT_FAILURE );
    }

    fread(Lr,sizeof(int),1,InFile);
    fread(GammaA,sizeof(double),(*Lr)+1,InFile);
}

void gammaasymptotics_(int *Lr, int *GammaA)
{
   GammaAsymptotics(Lr,GammaA);
}

void gammaasymptotics(int *Lr, int *GammaA)
{
   GammaAsymptotics(Lr,GammaA);
}

