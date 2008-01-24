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
 * This file contains function prototypes for each exported function in the
 * H5P module.
 */
#ifndef _H5Ppublic_H
#define _H5Ppublic_H

/* Default Template for creation, access, etc. templates */
#define H5P_DEFAULT     0

/* Public headers needed by this file */
#include "H5public.h"
#include "H5Ipublic.h"
#include "H5Dpublic.h"
#include "H5Fpublic.h"
#include "H5MMpublic.h"
#include "H5Zpublic.h"


/* Metroworks <sys/types.h> doesn't define off_t. */
#ifdef __MWERKS__
typedef long off_t;
/* Metroworks does not define EINTR in <errno.h> */
# define EINTR 4
#endif
/*__MWERKS__*/

/* Property list classes */
typedef enum H5P_class_t {
    H5P_NO_CLASS         = -1,  /*error return value                 */
    H5P_FILE_CREATE      = 0,   /*file creation properties           */
    H5P_FILE_ACCESS      = 1,   /*file access properties             */
    H5P_DATASET_CREATE   = 2,   /*dataset creation properties        */
    H5P_DATASET_XFER     = 3,   /*data transfer properties	     */
    H5P_MOUNT		 = 4,	/*file mounting properties	     */
    H5P_NCLASSES         = 5    /*this must be last!                 */
} H5P_class_t;

/* H5P_DATASET_XFER was the name from the beginning through 1.2.  It was
 * changed to H5P_DATA_XFER on v1.3.0.  Then it was changed back to
 * H5P_DATASET_XFER right before the release of v1.4.0-beta2.
 * Define an alias here to help applications that had ported to v1.3.
 * Should be removed in later version.
 */
#define H5P_DATA_XFER H5P_DATASET_XFER

#ifdef __cplusplus
extern "C" {
#endif

/* Define property list class callback function pointer types */
typedef herr_t (*H5P_cls_create_func_t)(hid_t prop_id, void *create_data);
typedef herr_t (*H5P_cls_close_func_t)(hid_t prop_id, void *close_data);

/* Define property list callback function pointer types */
typedef herr_t (*H5P_prp_create_func_t)(const char *name, void **def_value);
typedef herr_t (*H5P_prp_set_func_t)(hid_t prop_id, const char *name, void **value);
typedef herr_t (*H5P_prp_get_func_t)(hid_t prop_id, const char *name, void *value);
typedef herr_t (*H5P_prp_close_func_t)(const char *name, void *value);

/* Define property list iteration function type */
typedef herr_t (*H5P_iterate_t)(hid_t id, const char *name, void *iter_data); 

/*
 * The library created property list classes
 */
#define H5P_NO_CLASS_NEW   	(H5open(), H5P_NO_CLASS_g)
#define H5P_NO_CLASS_HASH_SIZE   1  /* 1, not 0, otherwise allocations get weird */
#define H5P_FILE_CREATE_NEW	(H5open(), H5P_FILE_CREATE_g)
#define H5P_FILE_CREATE_HASH_SIZE   17
#define H5P_FILE_ACCESS_NEW (H5open(), H5P_FILE_ACCESS_g)
#define H5P_FILE_ACCESS_HASH_SIZE   17
#define H5P_DATASET_CREATE_NEW  (H5open(), H5P_DATASET_CREATE_g)
#define H5P_DATASET_CREATE_HASH_SIZE   17
#define H5P_DATASET_XFER_NEW   (H5open(), H5P_DATASET_XFER_g)
#define H5P_DATASET_XFER_HASH_SIZE   17
#define H5P_MOUNT_NEW       (H5open(), H5P_MOUNT_g)
#define H5P_MOUNT_HASH_SIZE   17
H5_DLLVAR hid_t H5P_NO_CLASS_g;
H5_DLLVAR hid_t H5P_FILE_CREATE_g;
H5_DLLVAR hid_t H5P_FILE_ACCESS_g;
H5_DLLVAR hid_t H5P_DATASET_CREATE_g;
H5_DLLVAR hid_t H5P_DATASET_XFER_g;
H5_DLLVAR hid_t H5P_MOUNT_g;

/* Public functions */
H5_DLL hid_t H5Pcreate_class(hid_t parent, const char *name, unsigned hashsize,
            H5P_cls_create_func_t cls_create, void *create_data,
            H5P_cls_close_func_t cls_close, void *close_data);
H5_DLL char *H5Pget_class_name(hid_t pclass_id);
H5_DLL hid_t H5Pcreate_list(hid_t cls_id);
H5_DLL herr_t H5Pregister(hid_t cls_id, const char *name, size_t size,
            void *def_value, H5P_prp_create_func_t prp_create,
            H5P_prp_set_func_t prp_set, H5P_prp_get_func_t prp_get,
            H5P_prp_close_func_t prp_close);
H5_DLL herr_t H5Pinsert(hid_t plist_id, const char *name, size_t size,
            void *value, H5P_prp_set_func_t prp_set, H5P_prp_get_func_t prp_get,
            H5P_prp_close_func_t prp_close);
H5_DLL herr_t H5Pset(hid_t plist_id, const char *name, void *value);
H5_DLL htri_t H5Pexist(hid_t plist_id, const char *name);
H5_DLL herr_t H5Pget_size(hid_t id, const char *name, size_t *size);
H5_DLL herr_t H5Pget_nprops(hid_t id, size_t *nprops);
H5_DLL hid_t H5Pget_class_new(hid_t plist_id);
H5_DLL hid_t H5Pget_class_parent(hid_t pclass_id);
H5_DLL herr_t H5Pget(hid_t plist_id, const char *name, void * value);
H5_DLL htri_t H5Pequal(hid_t id1, hid_t id2);
H5_DLL int H5Piterate(hid_t id, int *idx, H5P_iterate_t iter_func,
            void *iter_data);
H5_DLL herr_t H5Premove(hid_t plist_id, const char *name);
H5_DLL herr_t H5Punregister(hid_t pclass_id, const char *name);
H5_DLL herr_t H5Pclose_list(hid_t plist_id);
H5_DLL herr_t H5Pclose_class(hid_t plist_id);
H5_DLL hid_t H5Pcreate(H5P_class_t type);
H5_DLL herr_t H5Pclose(hid_t plist_id);
H5_DLL hid_t H5Pcopy(hid_t plist_id);
H5_DLL H5P_class_t H5Pget_class(hid_t plist_id);
H5_DLL herr_t H5Pget_version(hid_t plist_id, int *boot/*out*/,
			      int *freelist/*out*/, int *stab/*out*/,
			      int *shhdr/*out*/);
H5_DLL herr_t H5Pset_userblock(hid_t plist_id, hsize_t size);
H5_DLL herr_t H5Pget_userblock(hid_t plist_id, hsize_t *size);
H5_DLL herr_t H5Pset_alignment(hid_t fapl_id, hsize_t threshold,
				hsize_t alignment);
H5_DLL herr_t H5Pget_alignment(hid_t fapl_id, hsize_t *threshold/*out*/,
				hsize_t *alignment/*out*/);
H5_DLL herr_t H5Pset_sizes(hid_t plist_id, size_t sizeof_addr,
			    size_t sizeof_size);
H5_DLL herr_t H5Pget_sizes(hid_t plist_id, size_t *sizeof_addr/*out*/,
			    size_t *sizeof_size/*out*/);
H5_DLL herr_t H5Pset_sym_k(hid_t plist_id, int ik, int lk);
H5_DLL herr_t H5Pget_sym_k(hid_t plist_id, int *ik/*out*/, int *lk/*out*/);
H5_DLL herr_t H5Pset_istore_k(hid_t plist_id, int ik);
H5_DLL herr_t H5Pget_istore_k(hid_t plist_id, int *ik/*out*/);
H5_DLL herr_t H5Pset_layout(hid_t plist_id, H5D_layout_t layout);
H5_DLL H5D_layout_t H5Pget_layout(hid_t plist_id);
H5_DLL herr_t H5Pset_chunk(hid_t plist_id, int ndims, const hsize_t dim[]);
H5_DLL int H5Pget_chunk(hid_t plist_id, int max_ndims, hsize_t dim[]/*out*/);
H5_DLL herr_t H5Pset_external(hid_t plist_id, const char *name, off_t offset,
			       hsize_t size);
H5_DLL int H5Pget_external_count(hid_t plist_id);
H5_DLL herr_t H5Pget_external(hid_t plist_id, int idx, size_t name_size,
			       char *name/*out*/, off_t *offset/*out*/,
			       hsize_t *size/*out*/);
H5_DLL herr_t H5Pset_driver(hid_t plist_id, hid_t driver_id,
			     const void *driver_info);
#if defined(WANT_H5_V1_2_COMPAT) || defined(H5_WANT_H5_V1_2_COMPAT)
H5_DLL H5F_driver_t H5Pget_driver(hid_t plist_id);
H5_DLL herr_t H5Pset_stdio(hid_t plist_id);
H5_DLL herr_t H5Pget_stdio(hid_t plist_id);
H5_DLL herr_t H5Pset_sec2(hid_t plist_id);
H5_DLL herr_t H5Pget_sec2(hid_t plist_id);
H5_DLL herr_t H5Pset_core(hid_t plist_id, size_t increment);
H5_DLL herr_t H5Pget_core(hid_t plist_id, size_t *increment/*out*/);
H5_DLL herr_t H5Pset_split(hid_t plist_id, const char *meta_ext, hid_t meta_plist_id,
	      const char *raw_ext, hid_t raw_plist_id);
H5_DLL herr_t H5Pget_split(hid_t plist_id, size_t meta_ext_size, char *meta_ext/*out*/,
	      hid_t *meta_properties/*out*/, size_t raw_ext_size,
	      char *raw_ext/*out*/, hid_t *raw_properties/*out*/);
H5_DLL herr_t H5Pset_family(hid_t plist_id, hsize_t memb_size, hid_t memb_plist_id);
H5_DLL herr_t H5Pget_family(hid_t plist_id, hsize_t *memb_size/*out*/,
	       hid_t *memb_plist_id/*out*/);
#if defined(H5_HAVE_PARALLEL)
H5_DLL herr_t H5Pset_mpi(hid_t plist_id, MPI_Comm comm, MPI_Info info);
H5_DLL herr_t H5Pget_mpi(hid_t plist_id, MPI_Comm *comm, MPI_Info *info);
H5_DLL herr_t H5Pset_xfer(hid_t plist_id, H5D_transfer_t data_xfer_mode);
H5_DLL herr_t H5Pget_xfer(hid_t plist_id, H5D_transfer_t *data_xfer_mode);
#endif /*H5_HAVE_PARALLEL*/
#else /* WANT_H5_V1_2_COMPAT */
H5_DLL hid_t H5Pget_driver(hid_t plist_id);
#endif /* WANT_H5_V1_2_COMPAT */
H5_DLL void *H5Pget_driver_info(hid_t plist_id);
H5_DLL herr_t H5Pset_buffer(hid_t plist_id, hsize_t size, void *tconv,
			     void *bkg);
H5_DLL hsize_t H5Pget_buffer(hid_t plist_id, void **tconv/*out*/,
			     void **bkg/*out*/);
H5_DLL herr_t H5Pset_preserve(hid_t plist_id, hbool_t status);
H5_DLL int H5Pget_preserve(hid_t plist_id);
H5_DLL herr_t H5Pset_filter(hid_t plist_id, H5Z_filter_t filter,
			     unsigned int flags, size_t cd_nelmts,
			     const unsigned int c_values[]);
H5_DLL int H5Pget_nfilters(hid_t plist_id);
H5_DLL H5Z_filter_t H5Pget_filter(hid_t plist_id, int filter,
				   unsigned int *flags/*out*/,
				   size_t *cd_nelmts/*out*/,
				   unsigned cd_values[]/*out*/,
				   size_t namelen, char name[]);
H5_DLL herr_t H5Pset_deflate(hid_t plist_id, unsigned aggression);
H5_DLL herr_t H5Pset_cache(hid_t plist_id, int mdc_nelmts, int rdcc_nelmts,
			    size_t rdcc_nbytes, double rdcc_w0);
H5_DLL herr_t H5Pget_cache(hid_t plist_id, int *mdc_nelmts/*out*/,
			    int *rdcc_nelmts/*out*/,
			    size_t *rdcc_nbytes/*out*/, double *rdcc_w0);
H5_DLL herr_t H5Pset_hyper_cache(hid_t plist_id, unsigned cache,
				  unsigned limit);
H5_DLL herr_t H5Pget_hyper_cache(hid_t plist_id, unsigned *cache,
				  unsigned *limit);
H5_DLL herr_t H5Pset_btree_ratios(hid_t plist_id, double left, double middle,
				   double right);
H5_DLL herr_t H5Pget_btree_ratios(hid_t plist_id, double *left/*out*/,
				   double *middle/*out*/,
				   double *right/*out*/);
H5_DLL herr_t H5Pset_fill_value(hid_t plist_id, hid_t type_id,
				 const void *value);
H5_DLL herr_t H5Pget_fill_value(hid_t plist_id, hid_t type_id,
				 void *value/*out*/);
H5_DLL herr_t H5Pset_gc_references(hid_t fapl_id, unsigned gc_ref);
H5_DLL herr_t H5Pget_gc_references(hid_t fapl_id, unsigned *gc_ref/*out*/);
H5_DLL herr_t H5Pset_vlen_mem_manager(hid_t plist_id,
                                       H5MM_allocate_t alloc_func,
                                       void *alloc_info, H5MM_free_t free_func,
                                       void *free_info);
H5_DLL herr_t H5Pget_vlen_mem_manager(hid_t plist_id,
                                       H5MM_allocate_t *alloc_func,
                                       void **alloc_info,
                                       H5MM_free_t *free_func,
                                       void **free_info);
H5_DLL herr_t H5Pset_meta_block_size(hid_t fapl_id, hsize_t size);
H5_DLL herr_t H5Pget_meta_block_size(hid_t fapl_id, hsize_t *size/*out*/);
H5_DLL herr_t H5Pset_sieve_buf_size(hid_t fapl_id, hsize_t size);
H5_DLL herr_t H5Pget_sieve_buf_size(hid_t fapl_id, hsize_t *size/*out*/);
H5_DLL herr_t H5Pset_small_data_block_size(hid_t fapl_id, hsize_t size);
H5_DLL herr_t H5Pget_small_data_block_size(hid_t fapl_id, hsize_t *size/*out*/);

#ifdef __cplusplus
}
#endif
#endif


