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
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*    F90 interface to the C HDF5 API                                        */
/*    Author: Matt Challacombe and CK Gan                                    */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
#include <stdio.h>
#include <stdlib.h>

#include "config.h"

#if defined (HAVE_INTERNAL_HDF5)
#include "hdf5.h"
#else
#include <hdf5.h>
#endif

#include <string.h>

#include <sys/types.h>
#include <regex.h>

/*
  THESE TURN ON VARIOUS LEVELS OF DEBUG
  #define debug_interface
  #define debug_all
*/

/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
char* IntToChar(int* NC, int* IntArray)
{
  int   j;
  char* VarName;

  VarName = (char*) malloc((*NC+1)*sizeof(char));
  for(j=0; j < *NC; j++) { VarName[j] = (char) IntArray[j]; }
  VarName[*NC] = '\0';
  return VarName;
}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5createfile_(int* NC, int* IChr)
{
  char* filename;
  hid_t fid;
  int FileID;

  filename = IntToChar(NC, IChr);
#ifdef debug_interface
  printf("IN CREATE_HDF5_FILE: Creating %s \n",filename);
#endif
  fid=H5Fcreate(filename,H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#ifdef debug_interface
#ifdef debug_all
  printf("IN CREATE_HDF5_FILE: FileID = %d \n",fid);
#endif
#endif
  FileID=fid;
  free(filename);
  return FileID;
}
int hdf5createfile(int* NC, int* IChr){return hdf5createfile_(NC,IChr);}
int hdf5createfile__(int* NC, int* IChr){return hdf5createfile_(NC,IChr);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5openfile_(int* NC, int* IChr)
{
  char* filename;
  hid_t fid;
  int   FileID;

  filename = IntToChar(NC,IChr);
#ifdef debug_interface
  printf("IN OPEN_HDF5_FILE: Opening <%s> \n", filename);
#endif
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
  fid=H5Fopen(filename,H5F_ACC_RDWR,H5P_DEFAULT);
#ifdef debug_interface
#ifdef debug_all
  printf("IN OPEN_HDF5_FILE: FileID = %d \n",fid);
#endif
#endif
  FileID=fid;
  free(filename);
  return FileID;
}
int hdf5openfile(int* NC, int* IChr){return hdf5openfile_(NC,IChr);}
int hdf5openfile__(int* NC, int* IChr){return hdf5openfile_(NC,IChr);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5closefile_(int* FileID)
{
  herr_t stat;
  hid_t fid;
  int   STATUS;

  fid=*FileID;
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_FILE: Flushing FileID= %d \n",FileID);
#endif
#endif
  stat=H5Fflush(fid,H5F_SCOPE_GLOBAL);
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_FILE: Closing  FileID= %d \n",FileID);
#endif
#endif
  stat=H5Fclose(fid);
  if(stat==-1){STATUS=stat; return STATUS;}
  /*printf("HDF5CloseFile: Could not close FileID= %d \n",*FileID);}*/
  stat=H5close();
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_FILE: STATUS = %d \n",stat);
#endif
#endif
  STATUS=stat;
  return STATUS;
}
int hdf5closefile(int* FileID){return hdf5closefile_(FileID);}
int hdf5closefile__(int* FileID){return hdf5closefile_(FileID);}
/*=================================================================================*/


/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5creategroup_(int* FileID, int* NC, int* IChr)
{
  char* filename;
  hid_t fid;
  hid_t gid;
  int GroupID;

  filename = IntToChar(NC,IChr);
#ifdef debug_interface
  printf("IN CREATE_HDF5_GROUP: Creating %s \n",filename);
#endif
  fid=*FileID;
  gid=H5Gcreate1(fid,filename,0);
#ifdef debug_interface
#ifdef debug_all
  printf("IN CREATE_HDF5_GROUP: FileID = %d, GroupID= %d \n",fid,gid);
#endif
#endif
  GroupID=gid;
  free(filename);
  return GroupID;
}
int hdf5creategroup(int* FileID, int* NC, int* IChr){return hdf5creategroup_(FileID,NC,IChr);}
int hdf5creategroup__(int* FileID, int* NC, int* IChr){return hdf5creategroup_(FileID,NC,IChr);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5opengroup_(int* FileID, int* NC, int* IChr)
{
  char* filename;
  hid_t fid;
  hid_t gid;
  int   GroupID;

  filename = IntToChar(NC,IChr);
#ifdef debug_interface
  printf("IN OPEN_HDF5_GROUP: Opening <%s> \n",filename);
#endif
  fid=*FileID;
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
  gid=H5Gopen1(fid,filename);
#ifdef debug_interface
#ifdef debug_all
  printf("IN OPEN_HDF5_GROUP: FileID = %d, GroupID = %d \n",fid,gid);
#endif
#endif
  GroupID=gid;
  free(filename);
  return GroupID;
}
int hdf5opengroup(int* FileID, int* NC, int* IChr){return hdf5opengroup_(FileID,NC,IChr);}
int hdf5opengroup__(int* FileID, int* NC, int* IChr){return hdf5opengroup_(FileID,NC,IChr);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5closegroup_(int* GroupID)
{
  herr_t stat;
  int   STATUS;
  hid_t gid;
  gid=*GroupID;
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_GROUP: Flushing GroupID= %d \n",GroupID);
#endif
#endif
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_GROUP: Closing  GroupID= %d \n",GroupID);
#endif
#endif
  stat=H5Gclose(gid);
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_GROUP: STATUS = %d \n",stat);
#endif
#endif
  STATUS=stat;
  return STATUS;
}
int hdf5closegroup(int* GroupID){return hdf5closegroup_(GroupID);}
int hdf5closegroup__(int* GroupID){return hdf5closegroup_(GroupID);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
void hdf5opendata_(int* FileID, int* NC, int* IChr, int* DataId, int* DataSpc)
{
  char* filename;
  hid_t fid,did,dspc;

  fid=*FileID;
  filename = IntToChar(NC,IChr);
#ifdef debug_interface
#ifdef debug_all
  printf("OPEN_HDF5_DATA:  FileID = %d \n",fid);
#endif
  printf("OPEN_HDF5_DATA:  VarName= <%s> \n",filename);
#endif
  did=H5Dopen1(fid,filename);
  dspc=H5Dget_space(did);
#ifdef debug_interface
#ifdef debug_all
  printf("OPEN_HDF5_DATA:  DataId = %d \n",did);
  printf("OPEN_HDF5_DATA:  DataSpc= %d \n",dspc);
#endif
#endif
  free(filename);
  *DataId=did;
  *DataSpc=dspc;
}
void hdf5opendata(int* FileID, int* NC, int* IChr, int* DataId, int* DataSpc)
{hdf5opendata_(FileID,NC,IChr,DataId,DataSpc);}
void hdf5opendata__(int* FileID, int* NC, int* IChr, int* DataId, int* DataSpc)
{hdf5opendata_(FileID,NC,IChr,DataId,DataSpc);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
void hdf5extenddata_(int* DataId, int* DataSpc, int* N)
{
  hid_t did,dspc,stat;
  hsize_t siz[1];
  int   STATUS;
  siz[0]=*N;
  did=*DataId;
  dspc=*DataSpc;
#ifdef debug_interface
  printf("EXTEND_HDF5_DATA:  DataId = %d \n",did);
  printf("EXTEND_HDF5_DATA:  NewSiz = %d \n",siz[0]);
#endif
  H5Dextend(did,siz);
  H5Sclose(dspc);
  dspc=H5Dget_space(did);
  *DataSpc=dspc;
#ifdef debug_interface
  printf("EXTEND_HDF5_DATA:  Old DataSpc = %d \n",dspc);
  printf("EXTEND_HDF5_DATA:  New DataSpc = %d \n",dspc);
#endif
}
void hdf5extenddata(int* DataId, int* DataSpc, int* N)
{hdf5extenddata_(DataId,DataSpc,N);}
void hdf5extenddata__(int* DataId, int* DataSpc, int* N)
{hdf5extenddata_(DataId,DataSpc,N);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
void hdf5selectdata_(int* DataId, int* DataSpc, int* NewSize)
{
  hid_t did,dspc;
  hsize_t      cnt[1];
  hsize_t      off[1];
  did=*DataId;
  dspc=*DataSpc;
  off[0]=0;
  cnt[0]=*NewSize;
  H5Sselect_hyperslab(dspc,H5S_SELECT_SET,off,NULL,cnt,NULL);
}
void hdf5selectdata(int* DataId, int* DataSpc, int* NewSize)
{hdf5selectdata_(DataId,DataSpc,NewSize);}
void hdf5selectdata__(int* DataId, int* DataSpc, int* NewSize)
{hdf5selectdata_(DataId,DataSpc,NewSize);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5closedata_(int* DataId,int* DataSpc)
{
  hid_t did,dspc;
  herr_t stat1,stat2;
  int   STATUS;
  did=*DataId;
  dspc=*DataSpc;
#ifdef debug_interface
#ifdef debug_all
  printf("CLOSE_HDF5_DATA:  DataId = %d \n",did);
  printf("CLOSE_HDF5_DATA:  DataSpc= %d \n",dspc);
#endif
#endif
  stat1=H5Dclose(did);
  stat2=H5Sclose(dspc);
  STATUS=stat1||stat2;
  return STATUS;
}
int hdf5closedata(int* DataId,int* DataSpc){return hdf5closedata_(DataId,DataSpc);}
int hdf5closedata__(int* DataId,int* DataSpc){return hdf5closedata_(DataId,DataSpc);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
void hdf5createdata_(int* FileID,int* Type,int* N,int* NC,int* IChr,int* ULimit, int* DataId, int* DataSpc)
{
  char* filename;
  hsize_t dms[1],mxdms[1],chnk[1];
  hid_t fid,did,dtyp,dspc,dprp,stat;

  filename = IntToChar(NC,IChr);
  fid=*FileID;
  dms[0]=*N;
#ifdef debug_interface
#ifdef debug_all
  printf("IN CREATE_HDF5_DATA: FileID = %d \n",fid);
  printf("IN CREATE_HDF5_DATA: N      = %d \n",dms[0]);
#endif
#endif
  if(*Type==24){dtyp=H5T_NATIVE_INT;}
  else if(*Type==6){dtyp=H5T_NATIVE_DOUBLE;}
  else{printf("Waring, unknown type = %d in HDF5DataId!",*Type); *DataId=-1; return;}
  dprp=H5Pcreate(H5P_DATASET_CREATE);
  if(*ULimit==1)
    {
      mxdms[0]=H5S_UNLIMITED;
      dspc=H5Screate_simple(1,dms,mxdms);
      chnk[0]=8192; /* 2048; */ /* 512; */
      stat=H5Pset_chunk(dprp,1,chnk);
#ifdef debug_interface
      printf("IN CREATE_HDF5_DATA: UNLIMITED DIMENSION = %d\n",dms[0]);
      printf("IN CREATE_HDF5_DATA: CHUNKING SIZE       = %d\n",chnk[0]);
#endif
      if(stat==-1){*DataId=-1; return;}
    }
  else
    {dspc=H5Screate_simple(1,dms,NULL);
#ifdef debug_interface
#ifdef debug_all
    printf("IN CREATE_HDF5_DATA: DIMENSIONS LIMITED TO %d\n",dms[0]);
#endif
#endif
    }
  did=H5Dcreate1(fid,filename,dtyp,dspc,dprp);
#ifdef debug_interface
  printf("IN CREATE_HDF5_DATA: VarName =<%s> \n",filename);
  printf("IN CREATE_HDF5_DATA: DataId  = %d \n",did);
#ifdef debug_all
  printf("IN CREATE_HDF5_DATA: DataTyp = %d \n",dtyp);
  printf("IN CREATE_HDF5_DATA: DataSpc = %d \n",dspc);
  printf("IN CREATE_HDF5_DATA: DataProp= %d \n",dprp);
#endif
#endif
  *DataId=did;
  *DataSpc=dspc;
  free(filename);
}
void hdf5createdata(int* FileID,int* Type,int* N,int* NC,int* IChr,int* ULimit, int* DataId, int* DataSpc)
{hdf5createdata_(FileID,Type,N,NC,IChr,ULimit,DataId,DataSpc);}
void hdf5createdata__(int* FileID,int* Type,int* N,int* NC,int* IChr,int* ULimit, int* DataId, int* DataSpc)
{hdf5createdata_(FileID,Type,N,NC,IChr,ULimit,DataId,DataSpc);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5sizeofdata_(int* DataSpc)
{
  hid_t dspc;
  hsize_t dsiz;
  int DataSize;
  dspc=*DataSpc;
#ifdef debug_interface
#ifdef debug_all
  printf("IN SIZE_OF_DATA: DataSpc = %d \n",*DataSpc);
#endif
#endif
  dsiz=H5Sget_simple_extent_npoints(dspc);
#ifdef debug_interface
#ifdef debug_all
  printf("IN SIZE_OF_DATA: DataSiz = %d \n",dsiz);
#endif
#endif
  DataSize=dsiz;
  return DataSize;
}
int hdf5sizeofdata(int* DataSpc)
{return hdf5sizeofdata_(DataSpc);}
int hdf5sizeofdata__(int* DataSpc)
{return hdf5sizeofdata_(DataSpc);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5writeintegervector_(int* DataId, int* DataSpc, int* Data)
{
  hid_t did,dspc;
  herr_t stat;
  int STATUS;
  did=*DataId;
  dspc=*DataSpc;
  stat=H5Dwrite(did,H5T_NATIVE_INT,H5S_ALL,dspc,H5P_DEFAULT,Data);
  STATUS=stat;
  return STATUS;
}
int hdf5writeintegervector(int* DataId, int* DataSpc, int* Data)
{return hdf5writeintegervector_(DataId,DataSpc,Data);}
int hdf5writeintegervector__(int* DataId, int* DataSpc, int* Data)
{return hdf5writeintegervector_(DataId,DataSpc,Data);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5writedoublevector_(int* DataId, int* DataSpc, double* Data)
{
  hid_t did,dspc;
  herr_t stat;
  int STATUS;
  did=*DataId;
  dspc=*DataSpc;
  stat=H5Dwrite(did,H5T_NATIVE_DOUBLE,H5S_ALL,dspc,H5P_DEFAULT,Data);
  STATUS=stat;
  return STATUS;
}
int hdf5writedoublevector(int* DataId, int* DataSpc, double* Data)
{return hdf5writedoublevector_(DataId,DataSpc,Data);}
int hdf5writedoublevector__(int* DataId, int* DataSpc, double* Data)
{return hdf5writedoublevector_(DataId,DataSpc,Data);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5readintegervector_(int* DataId, int* DataSpc, int* Data)
{
  hid_t did,dspc;
  herr_t stat;
  int STATUS;
  did=*DataId;
  dspc=*DataSpc;
  stat=H5Dread(did,H5T_NATIVE_INT,H5S_ALL,dspc,H5P_DEFAULT,Data);
  STATUS=stat;
  return STATUS;
}
int hdf5readintegervector(int* DataId, int* DataSpc, int* Data)
{return hdf5readintegervector_(DataId,DataSpc,Data);}
int hdf5readintegervector__(int* DataId, int* DataSpc, int* Data)
{return hdf5readintegervector_(DataId,DataSpc,Data);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5readdoublevector_(int* DataId, int* DataSpc, double* Data)
{
  hid_t did,dspc;
  herr_t stat;
  int   STATUS;
  did=*DataId;
  dspc=*DataSpc;
  stat=H5Dread(did,H5T_NATIVE_DOUBLE,H5S_ALL,dspc,H5P_DEFAULT,Data);
  STATUS=stat;
  return STATUS;
}
int hdf5readdoublevector(int* DataId, int* DataSpc, double* Data)
{return hdf5readdoublevector_(DataId,DataSpc,Data);}
int hdf5readdoublevector__(int* DataId, int* DataSpc, double* Data)
{return hdf5readdoublevector_(DataId,DataSpc,Data);}

/* Get the library version. */
int
hdf5version (unsigned *majnum, unsigned *minnum, unsigned *relnum)
{
  return H5get_libversion(majnum, minnum, relnum);
}
int
hdf5version_ (unsigned *majnum, unsigned *minnum, unsigned *relnum)
{
  return hdf5version(majnum, minnum, relnum);
}

typedef struct op_data_t
{
  hid_t hdfID;
  regex_t pattern;
}
op_data_t;

/* Deletes object from group given a name. The name can contain wildcards. */
herr_t
hdf5delete_op (hid_t groupID, const char *name, const H5L_info_t *info, void *op_data_arg)
{
  struct op_data_t *op_data = (struct op_data_t*) op_data_arg;
  int result;

  /* Debugging output. */
  fprintf(stderr, "[hdf5delete_op] group ID %i, hdf ID %i, object name = %s... ", groupID, op_data->hdfID, name);

  /* Match link name with objectName pattern. */
  result = regexec(&(op_data->pattern), name, 0, NULL, 0);
  if (result == 0)
  {
    fprintf(stderr, "found match... deleting object... ");
    fprintf(stderr, "link type %i... ", info->type);

    if ((result = H5Ldelete(op_data->hdfID, name, H5P_DEFAULT)) < 0)
    {
      fprintf(stderr, "%s:%i error\n", __FILE__, __LINE__);
      H5Eprint(H5Eget_current_stack(), stderr);
    }

    else
    {
      fprintf(stderr, "ok\n");
      return result;
    }
  }

  else if (result == REG_ESPACE)
  {
    fprintf(stderr, "ran out of memory\n");
  }

  else
  {
    fprintf(stderr, "no match\n");
  }

  /* Return an error. */
  return -1;
}

int
hdf5delete (int* id, int* groupNameN, int* groupNameC, int* objectNameN, int* objectNameC)
{
  char *groupName = IntToChar(groupNameN, groupNameC);
  char *objectName = IntToChar(objectNameN, objectNameC);
  hid_t *hdfID = (hid_t*) id;
  herr_t result;
  hsize_t index = 0;
  struct op_data_t op_data;

  /* Set some data. */
  op_data.hdfID = *hdfID;

  fprintf(stderr, "[hdf5delete] hdfID %i\n", *hdfID);
  fprintf(stderr, "[hdf5delete] searching for objects in group \"%s\"\n", groupName);
  fprintf(stderr, "[hdf5delete] object name filter \"%s\"\n", objectName);

  /* Compile regular expression pattern. */
  if (regcomp(&(op_data.pattern), objectName, 0) != 0)
  {
    fprintf(stderr, "[hdf5delete] error compiling pattern %s\n", objectName);
    exit(1);
  }

  /* More debugging. */
  //fprintf(stderr, "[hdf5delete] pattern.allocated = %u\n", op_data.pattern.allocated);

  /* Start to search. */
  while ((result = H5Literate_by_name(*hdfID, groupName, H5_INDEX_NAME,
          H5_ITER_NATIVE, &index, hdf5delete_op, (void*) &op_data, H5P_DEFAULT)) >= 0)
  {}

  exit(1);

  /* Clean up memory. */
  free(groupName);
  free(objectName);

  /* Return. */
  return 0;
}

int
hdf5delete_ (int* id, int* groupNameN, int* groupNameC, int* objectNameN, int* objectNameC)
{
  return hdf5delete(id, groupNameN, groupNameC, objectNameN, objectNameC);
}
