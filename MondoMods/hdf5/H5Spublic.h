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
 * This file contains public declarations for the H5S module.
 */
#ifndef _H5Spublic_H
#define _H5Spublic_H

/* Public headers needed by this file */
#include "H5public.h"
#include "H5Ipublic.h"

/* Define atomic datatypes */
#define H5S_ALL         0
#define H5S_UNLIMITED	((hsize_t)(hssize_t)(-1))

/* Define user-level maximum number of dimensions */
#define H5S_MAX_RANK    32

/* Different types of dataspaces */
typedef enum H5S_class_t {
    H5S_NO_CLASS         = -1,  /*error                                      */
    H5S_SCALAR           = 0,   /*scalar variable                            */
    H5S_SIMPLE           = 1,   /*simple data space                          */
    H5S_COMPLEX          = 2    /*complex data space                         */
} H5S_class_t;

/* Different ways of combining selections */
typedef enum H5S_seloper_t {
    H5S_SELECT_NOOP      = -1,  /* error                                     */
    H5S_SELECT_SET       = 0,   /* Select "set" operation 		     */
    H5S_SELECT_OR,              /* Binary "or" operation for hyperslabs
                                 * (add new selection to existing selection)		
                                 */
    H5S_SELECT_APPEND,          /* Append elements to end of point selection */
    H5S_SELECT_PREPEND,         /* Prepend elements to beginning of point selection */
    H5S_SELECT_INVALID          /* Invalid upper bound on selection operations */
} H5S_seloper_t;

#ifdef __cplusplus
extern "C" {
#endif

/* Functions in H5S.c */
H5_DLL hid_t H5Screate(H5S_class_t type);
H5_DLL hid_t H5Screate_simple(int rank, const hsize_t dims[],
			       const hsize_t maxdims[]);
H5_DLL herr_t H5Sset_extent_simple(hid_t space_id, int rank,
				    const hsize_t dims[],
				    const hsize_t max[]);
H5_DLL hid_t H5Scopy(hid_t space_id);
H5_DLL herr_t H5Sclose(hid_t space_id);
H5_DLL hssize_t H5Sget_simple_extent_npoints(hid_t space_id);
H5_DLL int H5Sget_simple_extent_ndims(hid_t space_id);
H5_DLL int H5Sget_simple_extent_dims(hid_t space_id, hsize_t dims[],
				      hsize_t maxdims[]);
H5_DLL htri_t H5Sis_simple(hid_t space_id);
H5_DLL herr_t H5Sset_space(hid_t space_id, int rank, const hsize_t *dims);
H5_DLL hssize_t H5Sget_select_npoints(hid_t spaceid);
H5_DLL herr_t H5Sselect_hyperslab(hid_t space_id, H5S_seloper_t op,
				   const hssize_t start[],
				   const hsize_t _stride[],
				   const hsize_t count[],
				   const hsize_t _block[]);
H5_DLL herr_t H5Sselect_elements(hid_t space_id, H5S_seloper_t op,
				  size_t num_elemn,
				  const hssize_t **coord);
H5_DLL H5S_class_t H5Sget_simple_extent_type(hid_t space_id);
H5_DLL herr_t H5Sset_extent_none(hid_t space_id);
H5_DLL herr_t H5Sextent_copy(hid_t dst_id,hid_t src_id);
H5_DLL herr_t H5Sselect_all(hid_t spaceid);
H5_DLL herr_t H5Sselect_none(hid_t spaceid);
H5_DLL herr_t H5Soffset_simple(hid_t space_id, const hssize_t *offset);
H5_DLL htri_t H5Sselect_valid(hid_t spaceid);
H5_DLL hssize_t H5Sget_select_hyper_nblocks(hid_t spaceid);
H5_DLL hssize_t H5Sget_select_elem_npoints(hid_t spaceid);
H5_DLL herr_t H5Sget_select_hyper_blocklist(hid_t spaceid, hsize_t startblock, hsize_t numblocks, hsize_t *buf);
H5_DLL herr_t H5Sget_select_elem_pointlist(hid_t spaceid, hsize_t startpoint, hsize_t numpoints, hsize_t *buf);
H5_DLL herr_t H5Sget_select_bounds(hid_t spaceid, hsize_t *start, hsize_t *end);

#ifdef __cplusplus
}
#endif
#endif
