#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void onxio(int *MyID, int *Action, int *FileN, int *Type, int *Leng, void *Buff)
{
    int PID = 0;
    char* MondoScratch;
    char FileName[100];
    char* Mode;
    int size;
    FILE* FileP;

    GetID(&PID);

    if(( MondoScratch = getenv("MONDO_SCRATCH")) == NULL) 
       exit( EXIT_FAILURE );

    switch (*FileN)
    {
       case 0:
          {
          sprintf(FileName,"%s/M%d_%d_ONX0",MondoScratch,PID,*MyID);
          break;
          }
       case 1:
          {
          sprintf(FileName,"%s/M%d_%d_ONX1",MondoScratch,PID,*MyID);
          break;
          }
       case 2:
          {
          sprintf(FileName,"%s/M%d_%d_ONX2",MondoScratch,PID,*MyID);
          break;
          }
       case 3:
          {
          sprintf(FileName,"%s/M%d_%d_ONX3",MondoScratch,PID,*MyID);
          break;
          }
       case 4:
          {
          sprintf(FileName,"%s/M%d_%d_ONX4",MondoScratch,PID,*MyID);
          break;
          }
       case 5:
          {
          sprintf(FileName,"%s/M%d_%d_ONX5",MondoScratch,PID,*MyID);
          break;
          }
       case 6:
          {
          sprintf(FileName,"%s/M%d_%d_ONX6",MondoScratch,PID,*MyID);
          break;
          }
       case 7:
          {
          sprintf(FileName,"%s/M%d_%d_ONX7",MondoScratch,PID,*MyID);
          break;
          }
       case 8:
          {
          sprintf(FileName,"%s/M%d_%d_ONX8",MondoScratch,PID,*MyID);
          break;
          }
       default:
          {
          printf("file number = %d\n",*FileN);
          printf("bad file number in onxio\n");
          exit( EXIT_FAILURE );
          }
    }
    switch (*Type)
    {
       case 0:
          {
          size = sizeof(int);
          break;
          }
       case 1:
          {
          size = sizeof(double);
          break;
          }
       default:
          {printf("bad Type in onxio\n");
          exit( EXIT_FAILURE );
          }
    }

    switch (*Action)
    {
       case 0:
          {
          Mode="wb";
          if ((FileP = fopen(FileName,Mode)) == NULL) {
             printf(" Error writing to file %s\n",FileName);
             exit( EXIT_FAILURE );
             }
          fwrite(Buff,size,(*Leng),FileP);
          fclose(FileP);
          break;
          }
       case 1:
          {
          Mode="rb";
          if ((FileP = fopen(FileName,Mode)) == NULL) {
             printf(" Error reading from file %s\n",FileName);
             exit( EXIT_FAILURE );
             }
          fread(Buff,size,(*Leng),FileP);
          fclose(FileP);
          break;
          }
       case 2:
          {
          if (remove(FileName) == -1) {
             printf(" Error erasing file %s\n",FileName);
             exit( EXIT_FAILURE );
             }
          break;
          }
       default:
          {printf("bad action in onxio\n");
          exit( EXIT_FAILURE );
          }
    }
}


void onxio_(int *MyID, int *Action, int *FileN, int *Type, int *Leng, void *Buff)
{
onxio(MyID,Action,FileN,Type,Leng,Buff);
}


