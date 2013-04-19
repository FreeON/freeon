/** @file
 */

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

#include <execinfo.h>
#include <stdarg.h>
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

/* THESE TURN ON VARIOUS LEVELS OF DEBUG */
//#define debug_interface
//#define debug_all

#define HDF5ERROR(...) hdf5error(__FILE__, __LINE__, __VA_ARGS__)

typedef struct op_data_t
{
  hid_t hdfID;
  regex_t pattern;
}
op_data_t;

/** Print an error message and exit.
 *
 * @param filename The filename the error occurred.
 * @param line The linenumber the error occurred.
 * @param format The format of the error message, see the man-page for
 * printf() for details.
 */
void
hdf5error (const char *const filename, const int line, ...)
{
  va_list va;
  int format_length;

  void *array[50];
  size_t size;
  char **strings;
  size_t i;

  char *old_format;
  char *new_format;

  /* Initialize variadic argument list. */
  va_start(va, line);

  /* Get format string. */
  old_format = va_arg(va, char*);

  format_length = 6+strlen(filename)+14+strlen(old_format);
  new_format = calloc(format_length, sizeof(char));
  snprintf(new_format, format_length, "[%s:%i] ", filename, line);
  strncat(new_format, old_format, format_length);

  /* Print error. */
  vprintf(new_format, va);
  va_end(va);

  /* Cleanup. */
  free(new_format);

  /* Print backtrace. */
  printf("\n");

  size = backtrace(array, 50);
  strings = backtrace_symbols(array, size);

  printf("Obtained %zd stack frames.\n", size-2);

  /* Omit the first 2 frames since we know those are the error handlers. */
  for(i = 2; i < size; i++)
  {
    printf("%s\n", strings[i]);
  }

  free(strings);

  /* Exit with signal so debuggers can produce a backtrace. */
  abort();
}

/** Convert an integer array into a string. Necessary for Fortran/C interface.
 *
 * @param NC The number of characters.
 * @param IntArray The integer array.
 *
 * @return A C-String. The string has to be free'ed by the caller.
 */
char *
IntToChar (int* NC, int* IntArray)
{
  int   j;
  char* VarName;

  VarName = (char*) malloc((*NC+1)*sizeof(char));
  for(j=0; j < *NC; j++) { VarName[j] = (char) IntArray[j]; }
  VarName[*NC] = '\0';
  return VarName;
}

/* Wrapper for H5Fcreate().
 *
 * @param NC The number of characters in the filename string.
 * @param IChar The filename string.
 *
 * @return The file descriptor of the created file.
 */
int
F77_FUNC (hdf5createfile, HDF5CREATEFILE) (int* NC, int* IChr)
{
  char* filename;
  hid_t fid;

  filename = IntToChar(NC, IChr);
#ifdef debug_interface
  printf("IN CREATE_HDF5_FILE: Creating %s \n",filename);
#endif
  if ((fid = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) < 0)
  {
    HDF5ERROR("failed to create hdf5 file\n");
  }
#ifdef debug_interface
#ifdef debug_all
  printf("IN CREATE_HDF5_FILE: FileID = %d \n",fid);
#endif
#endif
  free(filename);
  return fid;
}

/** Wrapper for H5FOpen().
 *
 * @param NC The number of characters in the filename string.
 * @param IChar The filename string.
 *
 * @return The file descriptor of the created file.
 */
int
F77_FUNC (hdf5openfile, HDF5OPENFILE) (int* NC, int* IChr)
{
  char* filename;
  hid_t fid;

  filename = IntToChar(NC,IChr);
#ifdef debug_interface
  printf("IN OPEN_HDF5_FILE: Opening <%s> \n", filename);
#endif
  H5Eset_auto(H5E_DEFAULT, NULL, NULL);
  if((fid = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT)) < 0)
  {
    HDF5ERROR("failed to open hdf5 file\n");
  }
#ifdef debug_interface
#ifdef debug_all
  printf("IN OPEN_HDF5_FILE: FileID = %d \n",fid);
#endif
#endif
  free(filename);
  return fid;
}

/** Wrapper for H5Fclose().
 *
 * @param NC The number of characters in the filename string.
 * @param IChar The filename string.
 *
 * @return Returns a positive number on success, and a negative one on
 * failure.
 */
int
F77_FUNC (hdf5closefile, HDF5CLOSEFILE) (int* FileID)
{
  herr_t result;
  hid_t fid;

  fid=*FileID;
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_FILE: Flushing FileID= %d \n",*FileID);
#endif
#endif
  if ((result = H5Fflush(fid,H5F_SCOPE_GLOBAL)) < 0)
  {
    HDF5ERROR("error flushing hdf5 file\n");
  }
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_FILE: Closing  FileID= %d \n",*FileID);
#endif
#endif
  if((result = H5Fclose(fid)) < 0)
  {
    return result;
  }
  result = H5close();
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_FILE: STATUS = %d \n", result);
#endif
#endif
  return result;
}

/** Wrapper for H5Gcreate().
 *
 * @param NC The number of characters in the filename string.
 * @param IChar The filename string.
 *
 * @return The group id.
 */
int
F77_FUNC (hdf5creategroup, HDF5CREATEGROUP) (int* FileID, int* NC, int* IChr)
{
  char* filename;
  hid_t fid;
  hid_t gid;

  filename = IntToChar(NC,IChr);
#ifdef debug_interface
  printf("IN CREATE_HDF5_GROUP: Creating %s \n",filename);
#endif
  fid=*FileID;
  if((gid = H5Gcreate(fid, filename, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) < 0)
  {
    HDF5ERROR("error created group\n");
  }
#ifdef debug_interface
#ifdef debug_all
  printf("IN CREATE_HDF5_GROUP: FileID = %d, GroupID= %d \n",fid,gid);
#endif
#endif
  free(filename);
  return gid;
}

/** Wrapper for H5Gopen().
 *
 * @param NC The number of characters in the filename string.
 * @param IChar The filename string.
 *
 * @return The group id.
 */
int
F77_FUNC (hdf5opengroup, HDF5OPENGROUP) (int* FileID, int* NC, int* IChr)
{
  char* filename;
  hid_t fid;
  hid_t gid;

  filename = IntToChar(NC,IChr);
#ifdef debug_interface
  printf("IN OPEN_HDF5_GROUP: Opening <%s> \n",filename);
#endif
  fid=*FileID;
  if(H5Eset_auto(H5E_DEFAULT, NULL, NULL) < 0)
  {
    HDF5ERROR("failed to turn off error printing\n");
  }
  if((gid = H5Gopen(fid, filename, H5P_DEFAULT)) < 0)
  {
    HDF5ERROR("error opening group\n");
  }
#ifdef debug_interface
#ifdef debug_all
  printf("IN OPEN_HDF5_GROUP: FileID = %d, GroupID = %d \n",fid,gid);
#endif
#endif
  free(filename);
  return gid;
}

/** Wrapper for H5Gclose().
 *
 * @param GroupID The group ID.
 *
 * @return The return status of the close operation.
 */
int
F77_FUNC (hdf5closegroup, HDF5CLOSEGROUP) (int* GroupID)
{
  herr_t result;
  hid_t gid;

  gid=*GroupID;
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_GROUP: Flushing GroupID= %d \n",*GroupID);
#endif
#endif
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_GROUP: Closing  GroupID= %d \n",*GroupID);
#endif
#endif
  result = H5Gclose(gid);
#ifdef debug_interface
#ifdef debug_all
  printf("IN CLOSE_HDF5_GROUP: STATUS = %d \n", result);
#endif
#endif
  return result;
}

/** Wrapper for H5Dopen().
 */
void
F77_FUNC (hdf5opendata, HDF5OPENDATA) (int* FileID, int* NC, int* IChr, int* DataId, int* DataSpc)
{
  char* filename;
  hid_t fid, did, dspc;

  fid=*FileID;
  filename = IntToChar(NC,IChr);
#ifdef debug_interface
#ifdef debug_all
  printf("OPEN_HDF5_DATA:  FileID = %d \n",fid);
#endif
  printf("OPEN_HDF5_DATA:  VarName = <%s> \n",filename);
#endif
  if((did = H5Dopen(fid, filename, H5P_DEFAULT)) < 0)
  {
    /* We will ignore this error. In InOut.F90 in the Put interface we will
     * always go through the open -> write -> close pattern, where the open
     * will sometimes fail since we haven't actually created the data yet. The
     * write however, will create the entry in the hdf5 file so that
     * subsequent reads will succeed. */

    //HDF5ERROR("error opening dataset\n");
  }
  if((dspc = H5Dget_space(did)) < 0)
  {
    //HDF5ERROR("error getting space\n");
  }
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

/** Wrapper for H5Dextend().
 */
void
F77_FUNC (hdf5extenddata, HDF5EXTENDDATA) (int* DataId, int* DataSpc, int* N)
{
  hid_t did, dspc;
  hsize_t siz[1];

  siz[0]=*N;
  did=*DataId;
  dspc=*DataSpc;
#ifdef debug_interface
  printf("EXTEND_HDF5_DATA:  DataId = %d \n",did);
  printf("EXTEND_HDF5_DATA:  NewSiz = %lld \n",siz[0]);
#endif
  if(H5Dset_extent(did, siz) < 0)
  {
    HDF5ERROR("error setting data extent\n");
  }
  if(H5Sclose(dspc) < 0)
  {
    HDF5ERROR("error closing data set\n");
  }
#ifdef debug_interface
  printf("EXTEND_HDF5_DATA:  Old DataSpc = %d \n",dspc);
#endif
  if((dspc = H5Dget_space(did)) < 0)
  {
    HDF5ERROR("error getting space\n");
  }
  *DataSpc=dspc;
#ifdef debug_interface
  printf("EXTEND_HDF5_DATA:  New DataSpc = %d \n",dspc);
#endif
}

/** Wrapper for H5Sselect().
 */
void
F77_FUNC (hdf5selectdata, HDF5SELECTDATA) (int* DataId, int* DataSpc, int* NewSize)
{
  hid_t dspc;
  hsize_t cnt[1];
  hsize_t off[1];

  dspc=*DataSpc;
  off[0]=0;
  cnt[0]=*NewSize;
  if(H5Sselect_hyperslab(dspc, H5S_SELECT_SET, off, NULL, cnt, NULL) < 0)
  {
    HDF5ERROR("error selecting hyperslab\n");
  }
}

/** Wrapper for H5Dclose().
 */
int
F77_FUNC (hdf5closedata, HDF5CLOSEDATA) (int* DataId,int* DataSpc)
{
  hid_t did, dspc;
  herr_t result;

  did=*DataId;
  dspc=*DataSpc;
#ifdef debug_interface
#ifdef debug_all
  printf("CLOSE_HDF5_DATA:  DataId = %d \n",did);
  printf("CLOSE_HDF5_DATA:  DataSpc= %d \n",dspc);
#endif
#endif
  if((result = H5Dclose(did)) < 0)
  {
    return result;
  }
  result = H5Sclose(dspc);
  return result;
}

/** Wrapper for H5Pcreate().
 */
void
F77_FUNC (hdf5createdata, HDF5CREATEDATA) (int* FileID,int* Type,int* N,int* NC,int* IChr,int* ULimit, int* DataId, int* DataSpc)
{
  char* filename;
  hsize_t dms[1],mxdms[1],chnk[1];
  hid_t fid,did,dtyp,dspc,dprp,result;

  filename = IntToChar(NC,IChr);
  fid=*FileID;
  dms[0]=*N;
#ifdef debug_interface
#ifdef debug_all
  printf("IN CREATE_HDF5_DATA: FileID = %d \n",fid);
  printf("IN CREATE_HDF5_DATA: N      = %lld \n",dms[0]);
#endif
#endif
  if(*Type == 24) dtyp = H5T_NATIVE_INT;
  else if(*Type == 6) dtyp = H5T_NATIVE_DOUBLE;
  else
  {
    printf("Warning, unknown type = %d in HDF5DataId!",*Type);
    *DataId=-1;
    return;
  }
  if((dprp = H5Pcreate(H5P_DATASET_CREATE)) < 0)
  {
    HDF5ERROR("error creating property list\n");
  }
  if(*ULimit==1)
  {
    mxdms[0]=H5S_UNLIMITED;
    if((dspc = H5Screate_simple(1, dms, mxdms)) < 0)
    {
      HDF5ERROR("error creating simple dataspace\n");
    }
    chnk[0]=8192; /* 2048; */ /* 512; */
    result = H5Pset_chunk(dprp, 1, chnk);
#ifdef debug_interface
    printf("IN CREATE_HDF5_DATA: UNLIMITED DIMENSION = %lld\n",dms[0]);
    printf("IN CREATE_HDF5_DATA: CHUNKING SIZE       = %lld\n",chnk[0]);
#endif
    if(result < 0)
    {
      *DataId=-1;
      return;
    }
  }
  else
  {
    if((dspc = H5Screate_simple(1, dms, NULL)) < 0)
    {
      HDF5ERROR("error creating simple dataspace\n");
    }
#ifdef debug_interface
#ifdef debug_all
    printf("IN CREATE_HDF5_DATA: DIMENSIONS LIMITED TO %lld\n",dms[0]);
#endif
#endif
  }
  if((did = H5Dcreate(fid, filename, dtyp, dspc, H5P_DEFAULT, dprp, H5P_DEFAULT)) < 0)
  {
    HDF5ERROR("error creating new dataset\n");
  }
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

/** Wrapper for H5Sget_simple_extent_npoints().
 */
int
F77_FUNC (hdf5sizeofdata, HDF5SIZEOFDATA) (int* DataSpc)
{
  hid_t dspc;
  hsize_t dsiz;

  dspc=*DataSpc;
#ifdef debug_interface
#ifdef debug_all
  printf("IN SIZE_OF_DATA: DataSpc = %d \n",*DataSpc);
#endif
#endif
  if((dsiz = H5Sget_simple_extent_npoints(dspc)) < 0)
  {
    HDF5ERROR("error getting simple extent\n");
  }
#ifdef debug_interface
#ifdef debug_all
  printf("IN SIZE_OF_DATA: DataSiz = %lld \n",dsiz);
#endif
#endif
  return dsiz;
}

/** Wrapper for H5Dwrite().
 */
int
F77_FUNC (hdf5writeintegervector, HDF5WRITEINTEGERVECTOR) (int* DataId, int* DataSpc, int* Data)
{
  hid_t did,dspc;
  herr_t result;

  did=*DataId;
  dspc=*DataSpc;
  if((result = H5Dwrite(did, H5T_NATIVE_INT, H5S_ALL, dspc, H5P_DEFAULT,Data)) < 0)
  {
    HDF5ERROR("error writing data\n");
  }
  return result;
}

/** Wrapper for H5Dwrite().
 */
int
F77_FUNC (hdf5writedoublevector, HDF5WRITEDOUBLEVECTOR) (int* DataId, int* DataSpc, double* Data)
{
  hid_t did,dspc;
  herr_t result;

  did=*DataId;
  dspc=*DataSpc;
  if((result = H5Dwrite(did, H5T_NATIVE_DOUBLE, H5S_ALL, dspc, H5P_DEFAULT, Data)) < 0)
  {
    HDF5ERROR("error writing data\n");
  }
  return result;
}

/** Wrapper for H5Dread().
 */
int
F77_FUNC (hdf5readintegervector, HDF5READINTEGERVECTOR) (int* DataId, int* DataSpc, int* Data)
{
  hid_t did,dspc;
  herr_t result;

  did=*DataId;
  dspc=*DataSpc;
  if((result = H5Dread(did, H5T_NATIVE_INT, H5S_ALL, dspc, H5P_DEFAULT, Data)) < 0)
  {
    HDF5ERROR("error reading data\n");
  }
  return result;
}

/** Wrapper for H5Dread().
 */
int
F77_FUNC (hdf5readdoublevector, HDF5READDOUBLEVECTOR) (int* DataId, int* DataSpc, double* Data)
{
  hid_t did,dspc;
  herr_t result;

  did=*DataId;
  dspc=*DataSpc;
  if((result = H5Dread(did, H5T_NATIVE_DOUBLE, H5S_ALL, dspc, H5P_DEFAULT, Data)) < 0)
  {
    HDF5ERROR("error reading data\n");
  }
  return result;
}

/** Get the library version. */
int
F77_FUNC (hdf5version, HDF5VERSION) (unsigned *majnum, unsigned *minnum, unsigned *relnum)
{
  return H5get_libversion(majnum, minnum, relnum);
}

/** Delete object from group given a name. The name can contain wildcards.
 */
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

/** Wrapper for H5LDelete().
 */
int
F77_FUNC (hdf5delete, HDF5DELETE) (int* id, int* groupNameN, int* groupNameC, int* objectNameN, int* objectNameC)
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
