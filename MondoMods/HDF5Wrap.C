/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*--  This source code is part of the MondoSCF suite of                      */ 
/*--  linear scaling electronic structure codes.                             */
/*--  Matt Challacombe                                                       */
/*--  Los Alamos National Laboratory                                         */
/*--  Copyright 200, The University of California                            */
/*                                                                           */
/*    F90 interface to the C HDF5 API                                        */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
#include <stdio.h>
#include <stdlib.h>  
#include <hdf5.h>
#include <string.h>
/*
THESE TURN ON VARIOUS LEVELS OF DEBUG
#define debug_interface
#define debug_all

*/
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
const char* IntToChar(int* NC, int* IntArray)
{
   int   j;
   char* VarName;
   char  CharTemp[132];
   VarName = (char*)calloc(*NC,sizeof(char*));
   for(j=0; j<=*NC-1; j++){CharTemp[j]=(char) IntArray[j];}
   CharTemp[*NC]=NULL;
   memcpy(VarName,CharTemp,*NC); 
   return VarName;
}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5createfile_(int* NC, int* IChr)
{
   hid_t fid;
   int FileID;
#ifdef debug_interface
   printf("IN CREATE_HDF5_FILE: Creating %s \n",IntToChar(NC,IChr));
#endif
   fid=H5Fcreate(IntToChar(NC,IChr),H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
#ifdef debug_interface
#ifdef debug_all
   printf("IN CREATE_HDF5_FILE: FileID = %d \n",fid);
#endif
#endif
   FileID=fid;
   return FileID;
}
int hdf5createfile(int* NC, int* IChr)
{return hdf5createfile_(NC,IChr);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
int hdf5openfile_(int* NC, int* IChr) 
{
   hid_t fid;
   int   FileID;
#ifdef debug_interface
   printf("IN OPEN_HDF5_FILE: Opening <%s> \n",IntToChar(NC,IChr));
#endif
   H5Eset_auto(NULL, NULL);
   fid=H5Fopen(IntToChar(NC,IChr),H5F_ACC_RDWR,H5P_DEFAULT);
#ifdef debug_interface
#ifdef debug_all
   printf("IN OPEN_HDF5_FILE: FileID = %d \n",fid);
#endif
#endif
   FileID=fid;
   return FileID;
}
int hdf5openfile(int* NC, int* IChr)
{return hdf5openfile_(NC,IChr);}
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
   if(stat==-1){printf("HDF5CloseFile: Could not close FileID= %d \n",*FileID);}
   stat=H5close();
#ifdef debug_interface
#ifdef debug_all
   printf("IN CLOSE_HDF5_FILE: STATUS = %d \n",stat);
#endif
#endif
   STATUS=stat;
   return STATUS;
}
int hdf5closefile(int* FileID)
{return hdf5closefile_(FileID);}
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
void hdf5opendata_(int* FileID, int* NC, int* IChr, int* DataId, int* DataSpc)
{
   hid_t fid,did,dspc;
   fid=*FileID;
#ifdef debug_interface
#ifdef debug_all
   printf("OPEN_HDF5_DATA:  FileID = %d \n",fid);
#endif
   printf("OPEN_HDF5_DATA:  VarName= <%s> \n",IntToChar(NC,IChr));
#endif
   did=H5Dopen(fid,IntToChar(NC,IChr));
   dspc=H5Dget_space(did);
#ifdef debug_interface
#ifdef debug_all
   printf("OPEN_HDF5_DATA:  DataId = %d \n",did);
   printf("OPEN_HDF5_DATA:  DataSpc= %d \n",dspc);
#endif
#endif
   *DataId=did;
   *DataSpc=dspc;
}
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
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
void hdf5selectdata_(int* DataId, int* DataSpc, int* NewSize)
{
   hid_t did,dspc;
   hsize_t      cnt[1]; 
   hssize_t     off[1]; 
   did=*DataId;
   dspc=*DataSpc;
   off[0]=0;
   cnt[0]=*NewSize;
   H5Sselect_hyperslab(dspc,H5S_SELECT_SET,off,NULL,cnt,NULL);
}
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
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
/*                                                                           */
/*+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+*/
void hdf5createdata_(int* FileID,int* Type,int* N,int* NC,int* IChr,int* ULimit, int* DataId, int* DataSpc)
{ 
   hsize_t dms[1],mxdms[1],chnk[1];
   hid_t fid,did,dtyp,dspc,dprp,stat;
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
    chnk[0]=512;
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
   did=H5Dcreate(fid,IntToChar(NC,IChr),dtyp,dspc,dprp);
#ifdef debug_interface
   printf("IN CREATE_HDF5_DATA: VarName =<%s> \n",IntToChar(NC,IChr));
   printf("IN CREATE_HDF5_DATA: DataId  = %d \n",did);
#ifdef debug_all
   printf("IN CREATE_HDF5_DATA: DataTyp = %d \n",dtyp);
   printf("IN CREATE_HDF5_DATA: DataSpc = %d \n",dspc);
   printf("IN CREATE_HDF5_DATA: DataProp= %d \n",dprp);
#endif
#endif
   *DataId=did;
   *DataSpc=dspc;
}
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

