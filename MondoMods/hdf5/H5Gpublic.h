/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdf.ncsa.uiuc.edu/HDF5/doc/Copyright.html.  If you do not have     *
 * access to either file, you may request a copy from hdfhelp@ncsa.uiuc.edu. *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*-------------------------------------------------------------------------
 *
 * Created:             H5Gproto.h
 *                      Jul 11 1997
 *                      Robb Matzke <matzke@llnl.gov>
 *
 * Purpose:             Public declarations for the H5G package (symbol
 *                      tables).
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
#ifndef _H5Gpublic_H
#define _H5Gpublic_H

/* Public headers needed by this file */
#include <sys/types.h>

#include "H5public.h"
#include "H5Ipublic.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Types of links */
typedef enum H5G_link_t {
    H5G_LINK_ERROR	= -1,
    H5G_LINK_HARD	= 0,
    H5G_LINK_SOFT	= 1
} H5G_link_t;

/*
 * An object has a certain type. The first few numbers are reserved for use
 * internally by HDF5. Users may add their own types with higher values.  The
 * values are never stored in the file -- they only exist while an
 * application is running.  An object may satisfy the `isa' function for more
 * than one type.
 */
#define H5G_UNKNOWN	(-1)		/* Unknown object type		*/
#define H5G_LINK	0		/* Object is a symbolic link	*/
#define H5G_GROUP	1		/* Object is a group		*/
#define H5G_DATASET	2		/* Object is a dataset		*/
#define H5G_TYPE	3		/* Object is a named data type	*/
#define H5G_RESERVED_4	4		/* Reserved for future use	*/
#define H5G_RESERVED_5	5		/* Reserved for future use	*/
#define H5G_RESERVED_6	6		/* Reserved for future use	*/
#define H5G_RESERVED_7	7		/* Reserved for future use	*/

#define H5G_NTYPES	256		/* Max possible number of types	*/
#define H5G_NLIBTYPES	8		/* Number of internal types	*/
#define H5G_NUSERTYPES	(H5G_NTYPES-H5G_NLIBTYPES)
#define H5G_USERTYPE(X)	(8+(X))		/* User defined types		*/

/* Information about an object */
typedef struct H5G_stat_t {
    unsigned long 	fileno[2];	/*file number			*/
    unsigned long 	objno[2];	/*object number			*/
    unsigned 		nlink;		/*number of hard links to object*/
    int 		type;		/*basic object type		*/
    time_t		mtime;		/*modification time		*/
    size_t		linklen;	/*symbolic link value length	*/
} H5G_stat_t;
    

typedef herr_t (*H5G_iterate_t)(hid_t group, const char *name,
				void *op_data);

H5_DLL hid_t H5Gcreate(hid_t loc_id, const char *name, size_t size_hint);
H5_DLL hid_t H5Gopen(hid_t loc_id, const char *name);
H5_DLL herr_t H5Gclose(hid_t group_id);
H5_DLL herr_t H5Giterate(hid_t loc_id, const char *name, int *idx,
			  H5G_iterate_t op, void *op_data);
H5_DLL herr_t H5Gmove(hid_t loc_id, const char *src, const char *dst);
H5_DLL herr_t H5Glink(hid_t loc_id, H5G_link_t type, const char *cur_name,
		       const char *new_name);
H5_DLL herr_t H5Gunlink(hid_t loc_id, const char *name);
H5_DLL herr_t H5Gget_objinfo(hid_t loc_id, const char *name,
			      hbool_t follow_link, H5G_stat_t *statbuf/*out*/);
H5_DLL herr_t H5Gget_linkval(hid_t loc_id, const char *name, size_t size,
			      char *buf/*out*/);
H5_DLL herr_t H5Gset_comment(hid_t loc_id, const char *name,
			      const char *comment);
H5_DLL int H5Gget_comment(hid_t loc_id, const char *name, size_t bufsize,
			   char *buf);

#ifdef __cplusplus
}
#endif
#endif
