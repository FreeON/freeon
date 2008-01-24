/*
     This code is part of the MondoSCF suite of programs for linear scaling
     electronic structure theory and ab initio molecular dynamics.

     Copyright (2004). The Regents of the University of California. This
     material was produced under U.S. Government contract W-7405-ENG-36
     for Los Alamos National Laboratory, which is operated by the University
     of California for the U.S. Department of Energy. The U.S. Government has
     rights to use, reproduce, and distribute this software.  NEITHER THE
     GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
     OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.

     This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by the
     Free Software Foundation; either version 2 of the License, or (at your
     option) any later version. Accordingly, this program is distributed in
     the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE. See the GNU General Public License at www.gnu.org for details.

     While you may do as you like with this software, the GNU license requires
     that you clearly mark derivative software.  In addition, you are encouraged
     to return derivative works to the MondoSCF group for review, and possible
     disemination in future releases.
*/
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

  //printf(".... ing Spawning  Spawning  Spawning  Spawning Spawn ....\n");
  //for(i=0; i<=*nc; i++){printf("argv[%i]=<%s>\n",i,argv[i]);} 

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

int spawn(int* nc,int* maxlen,int* ichr){return spawn_(nc,maxlen,ichr); }
