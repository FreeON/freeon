/*   FORKS A CHILD PROCESS    */
/*   Author: Matt Challacombe */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/wait.h>
#include <sys/types.h>
#include <unistd.h>
/*
       pid_t fork(void);

       #include <sys/types.h>
       #include <sys/wait.h>

       pid_t wait(int *status)
       pid_t waitpid(pid_t pid, int *status, int options);
*/


int spawn(int* nc,int* maxlen,int* ichr){return spawn_(nc,maxlen,ichr); }

int spawn_(int* nc,int* maxlen,int* ichr)
{
  int i,j,k,setnull,status;
  int clen=*maxlen-1;
  char* argv[100];
  pid_t pid,wpid;
  int ZERO_ERROR=0;
  int EXIT_ERROR=-120384;
  int FORK_ERROR=-320498;
  int DUMP_ERROR=-580234;
  int SGNL_ERROR=-674034;
  int MISC_ERROR=-843503;
  int ierr=MISC_ERROR;
  k=0; 
  for(i=0; i<*nc; i++){
     argv[i] = (char*)calloc(clen+1,sizeof(char*));
     for(j=0; j<=clen; j++){
         argv[i][j]=(char) ichr[k]; k=k+1;
     }
     setnull=clen+1;
     for(j=0; j<=clen; j++){
       if(argv[i][j]==' '){
          setnull=j; 
          break;
       }
     }
     argv[i][setnull]='\0';
  }
  argv[*nc]=NULL;
  /*
  printf(".... ing Spawning  Spawning  Spawning  Spawning Spawn ....\n");
  for(i=0; i<=*nc; i++){printf("argv[%i]=<%s>\n",i,argv[i]);} 
  */
  pid=fork();
  if(pid==0)
    {
      execvp(argv[0],argv);
      printf("Unable to EXECVP <%s>.\n",argv[0]);
      for(i=0; i<=*nc; i++){printf("argv[%i]=<%s>\n",i,argv[i]);}
      _exit(EXIT_FAILURE);
    }
  else if(pid<0)
    {
      printf("Unable to FORK <%s>.\n",argv[0]);
      for(i=0; i<=*nc; i++){printf("argv[%i]=<%s>\n",i,argv[i]);}
      ierr=FORK_ERROR;
    }
  else
    {
      wpid=waitpid(pid,&status,0);
      if(wpid==pid){
         ierr=ZERO_ERROR;
         return(ierr);
      }
      if(WSTOPSIG   (status)!=0)ierr=SGNL_ERROR;
#ifdef AIX
#else     
      if(WCOREDUMP  (status)!=0)ierr=DUMP_ERROR; 
#endif
      if(WEXITSTATUS(status)!=0)ierr=EXIT_ERROR;
    }  
  for(i=0; i<*nc; i++){free(argv[i]);}
  return(ierr);
}
