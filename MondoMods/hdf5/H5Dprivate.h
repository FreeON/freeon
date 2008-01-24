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
 * This file contains private information about the H5D module
 */
#ifndef _H5Dprivate_H
#define _H5Dprivate_H

#include "H5Dpublic.h"

/* Private headers needed by this file */
#include "H5private.h"
#include "H5Fprivate.h"		/*for the H5F_t type			     */
#include "H5Gprivate.h"		/*symbol tables				     */
#include "H5MMpublic.h"		/*for H5MM_allocate_t and H5MM_free_t types  */
#ifdef H5_HAVE_PARALLEL
#include "H5FDmpio.h"           /* MPI-I/O file driver info                  */
#endif /* H5_HAVE_PARALLEL */
#include "H5Oprivate.h"		/*object Headers			     */
#include "H5Sprivate.h"		/*for the H5S_t type			     */
#include "H5Tprivate.h"		/*for the H5T_t type			     */

/*
 * Feature: Define H5D_DEBUG on the compiler command line if you want to
 *	    debug dataset I/O. NDEBUG must not be defined in order for this
 *	    to have any effect.
 */
#ifdef NDEBUG
#  undef H5D_DEBUG
#endif

#define H5D_RESERVED_ATOMS  0

/* Set the minimum object header size to create objects with */
#define H5D_MINHDR_SIZE 256 

/* Dataset creation property list */
typedef struct H5D_create_t {
    H5D_layout_t	layout;		/*storage layout		     */
    int		chunk_ndims;	/*chunk dimensionality		     */
    hsize_t		chunk_size[32];	/*chunk size if chunked storage	     */
    H5O_fill_t		fill;		/*fill value			     */
    H5O_efl_t		efl;		/*external file list		     */
    H5O_pline_t		pline;		/*data filter pipeline		     */
} H5D_create_t;

/* Data transfer property list */
typedef struct H5D_xfer_t {
    hsize_t		buf_size;	/*max temp buffer size		     */
    void		*tconv_buf;	/*type conversion buffer or null     */
    void		*bkg_buf;	/*background buffer or null	     */
    H5T_bkg_t	need_bkg;	/*type of background buffer needed   */
    double		split_ratios[3];/*B-tree node splitting ratios	     */
    unsigned       cache_hyper;    /*cache hyperslab blocks during I/O? */
    unsigned       block_limit;    /*largest hyperslab block to cache   */
    H5MM_allocate_t 	vlen_alloc; 	/*VL datatype allocation function    */
    void		*alloc_info;    /*VL datatype allocation information */
    H5MM_free_t	vlen_free;      /*VL datatype free function	     */
    void		*free_info;	/*VL datatype free information	     */
    hid_t		driver_id;	/*File driver ID		     */
    void		*driver_info;	/*File driver specific information   */
#ifdef H5_HAVE_PARALLEL
    unsigned            use_view;       /* Whether to use a view for this I/O */
    MPI_Datatype        btype;          /* MPI type for buffer (memory) */
    MPI_Datatype        ftype;          /* MPI type for file */
    unsigned            block_before_meta_write;    /* Whether to block before metadata writes */
    H5FD_mpio_xfer_t	xfer_mode;	/*collective or independent I/O	*/
    unsigned            library_internal;/* Whether this property list is used internally to the library or not */
#endif /* H5_HAVE_PARALLEL */
} H5D_xfer_t;

typedef struct H5D_t H5D_t;

/* library variables */
H5_DLLVAR const H5D_create_t H5D_create_dflt;
H5_DLLVAR H5D_xfer_t H5D_xfer_dflt;

/* Functions defined in H5D.c */
H5_DLL herr_t H5D_init(void);
H5_DLL H5D_t *H5D_create(H5G_entry_t *loc, const char *name,
			  const H5T_t *type, const H5S_t *space,
			  hid_t dxpl_id, const H5D_create_t *create_parms);
H5_DLL H5D_t *H5D_open(H5G_entry_t *loc, const char *name, hid_t dxpl_id);
H5_DLL herr_t H5D_close(H5D_t *dataset);
H5_DLL htri_t H5D_isa(H5G_entry_t *ent, hid_t dxpl_id);
H5_DLL herr_t H5D_read(H5D_t *dataset, const H5T_t *mem_type,
			const H5S_t *mem_space, const H5S_t *file_space,
			hid_t dset_xfer_plist, void *buf/*out*/);
H5_DLL herr_t H5D_write(H5D_t *dataset, const H5T_t *mem_type,
			 const H5S_t *mem_space, const H5S_t *file_space,
			 hid_t dset_xfer_plist, const void *buf);
H5_DLL herr_t H5D_extend(H5D_t *dataset, const hsize_t *size, hid_t dxpl_id);
H5_DLL H5G_entry_t *H5D_entof(H5D_t *dataset);
H5_DLL H5T_t *H5D_typeof(H5D_t *dset);
H5_DLL H5S_t *H5D_get_space(H5D_t *dset, hid_t dxpl_id);
H5_DLL H5D_t * H5D_open_oid(H5G_entry_t *ent, hid_t dxpl_id);
H5_DLL H5F_t * H5D_get_file(const H5D_t *dset);
H5_DLL hsize_t H5D_get_storage_size(H5D_t *dset, hid_t dxpl_id);
H5_DLL void *H5D_vlen_get_buf_size_alloc(size_t size, void *info);
H5_DLL herr_t H5D_vlen_get_buf_size(void *elem, hid_t type_id, hsize_t ndim, hssize_t *point, void *op_data);

#endif
