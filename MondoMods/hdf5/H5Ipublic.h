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

/*
 * This file contains function prototypes for each exported function in
 * the H5I module.
 */
#ifndef _H5Ipublic_H
#define _H5Ipublic_H

/* Public headers needed by this file */
#include "H5public.h"

/*
 * Group values allowed.  Start with `1' instead of `0' because it makes the
 * tracing output look better when hid_t values are large numbers.  Change the
 * GROUP_BITS in H5I.c if the MAXID gets larger than 32 (an assertion will
 * fail otherwise).
 */
typedef enum {
    H5I_BADID		= (-1),	/*invalid Group				    */
    H5I_FILE		= 1,	/*group ID for File objects		    */
    H5I_FILE_CLOSING,		/*files pending close due to open objhdrs   */
    H5I_TEMPLATE_0,		    /*group ID for Template objects		    */
    H5I_TEMPLATE_1,	       	/*group ID for Template objects		    */
    H5I_TEMPLATE_2,	        /*group ID for Template objects		    */
    H5I_TEMPLATE_3,	        /*group ID for Template objects		    */
    H5I_TEMPLATE_4,	        /*group ID for Template objects		    */
    H5I_TEMPLATE_5,	        /*group ID for Template objects		    */
    H5I_TEMPLATE_6,	        /*group ID for Template objects		    */
    H5I_TEMPLATE_7,	        /*group ID for Template objects		    */
    H5I_TEMPLATE_MAX,	    /*not really a group ID			        */
    H5I_GROUP,		        /*group ID for Group objects		    */
    H5I_DATATYPE,	        /*group ID for Datatype objects		    */
    H5I_DATASPACE,	        /*group ID for Dataspace objects	    */
    H5I_DATASET,	        /*group ID for Dataset objects		    */
    H5I_ATTR,		        /*group ID for Attribute objects	    */
    H5I_TEMPBUF,	        /*group ID for Temporary buffer objects	*/
    H5I_REFERENCE,	        /*group ID for Reference objects	    */
    H5I_VFL,			    /*group ID for virtual file layer	    */
    H5I_GENPROP_CLS,        /*group ID for generic property list classes */
    H5I_GENPROP_LST,        /*group ID for generic property lists   */
    
    H5I_NGROUPS		        /*number of valid groups, MUST BE LAST!	    */
} H5I_type_t;

/* Type of atoms to return to users */
typedef int hid_t;

/* An invalid object ID. This is also negative for error return. */
#define H5I_INVALID_HID         (-1)

#ifdef __cplusplus
extern "C" {
#endif

H5_DLL H5I_type_t H5Iget_type(hid_t id);


#ifdef __cplusplus
}
#endif
#endif
