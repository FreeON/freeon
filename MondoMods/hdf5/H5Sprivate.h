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
 * This file contains private information about the H5S module
 */
#ifndef _H5Sprivate_H
#define _H5Sprivate_H

#include "H5Spublic.h"

/* Private headers needed by this file */
#include "H5private.h"
#include "H5Dpublic.h"
#include "H5Fprivate.h"
#include "H5Oprivate.h"

#define H5S_RESERVED_ATOMS  2

/* Flags to indicate special dataspace features are active */
#define H5S_VALID_MAX	0x01
#define H5S_VALID_PERM	0x02

/* Flags for H5S_find */
#define H5S_CONV_PAR_IO_POSSIBLE        0x0001
/* The storage options are mutually exclusive */
/* (2-bits reserved for storage type currently) */
#define H5S_CONV_STORAGE_COMPACT        0x0000  /* i.e. '0' */
#define H5S_CONV_STORAGE_CONTIGUOUS     0x0002  /* i.e. '1' */
#define H5S_CONV_STORAGE_CHUNKED        0x0004  /* i.e. '2' */
#define H5S_CONV_STORAGE_MASK           0x0006

/* Forward references of common typedefs */
typedef struct H5S_t H5S_t;
typedef struct H5S_pnt_node_t H5S_pnt_node_t;

/* Enumerated type for the type of selection */
typedef enum {
    H5S_SEL_ERROR	= -1, 		/* Error			*/
    H5S_SEL_NONE	= 0,           	/* Nothing selected 		*/
    H5S_SEL_POINTS	= 1,         	/* Sequence of points selected	*/
    H5S_SEL_HYPERSLABS	= 2,     	/* Hyperslab selection defined	*/
    H5S_SEL_ALL		= 3,            /* Entire extent selected	*/
    H5S_SEL_N		= 4		/*THIS MUST BE LAST		*/
}H5S_sel_type;

/* Point selection iteration container */
typedef struct {
    hsize_t elmt_left;      /* Number of elements left to iterate over */
    H5S_pnt_node_t *curr;   /* Pointer to next node to output */
} H5S_point_iter_t;

/* Hyperslab selection iteration container */
typedef struct {
    hsize_t elmt_left;      /* Number of elements left to iterate over */
    hssize_t *pos;          /* Position to start iterating at */
} H5S_hyper_iter_t;

/* "All" selection iteration container */
typedef struct {
    hsize_t elmt_left;      /* Number of elements left to iterate over */
    hsize_t offset;         /* Next element to output */
} H5S_all_iter_t;

/* Selection iteration container */
typedef union {
    H5S_point_iter_t pnt;   /* Point selection iteration information */
    H5S_hyper_iter_t hyp;   /* Hyperslab selection iteration information */
    H5S_all_iter_t all;     /* "All" selection iteration information */
} H5S_sel_iter_t;

/*
 * Data space conversions usually take place in two halves.  One half
 * transfers data points between memory and a data type conversion array
 * where the points are contiguous, and the other half transfers points
 * between the type conversion array and the file.
 */
typedef struct H5S_fconv_t {
    /* Identification */
    const char 		*name;
    H5S_sel_type	type;
    
    /* Initialize file element numbering information */
    herr_t (*init)(const struct H5O_layout_t *layout, const H5S_t *space,
		   H5S_sel_iter_t *iter);

    /* Determine optimal number of elements to transfer */
    hsize_t (*avail)(const H5S_t *file_space, const H5S_sel_iter_t *file_iter,
		    hsize_t max);

    /* Gather elements from disk to type conversion buffer */
    hsize_t (*gath)(H5F_t *f, const struct H5O_layout_t *layout,
		   const struct H5O_pline_t *pline,
		   const struct H5O_fill_t *fill,
		   const struct H5O_efl_t *efl, size_t elmt_size,
		   const H5S_t *file_space, H5S_sel_iter_t *file_iter,
		   hsize_t nelmts, hid_t dxpl_id, void *tconv_buf/*out*/);

    /* Scatter elements from type conversion buffer to disk */
    herr_t (*scat)(H5F_t *f, const struct H5O_layout_t *layout,
		   const struct H5O_pline_t *pline,
		   const struct H5O_fill_t *fill,
		   const struct H5O_efl_t *efl, size_t elmt_size,
		   const H5S_t *file_space, H5S_sel_iter_t *file_iter,
		   hsize_t nelmts, hid_t dxpl_id, const void *tconv_buf);
} H5S_fconv_t;

typedef struct H5S_mconv_t {
    /* Identification */
    const char		*name;
    H5S_sel_type	type;
    
    /* Initialize memory element numbering information */
    herr_t (*init)(const struct H5O_layout_t *layout, const H5S_t *space,
		   H5S_sel_iter_t *iter);

    /* Gather elements from app buffer to type conversion buffer */
    hsize_t (*gath)(const void *buf, size_t elmt_size,
		   const H5S_t *mem_space, H5S_sel_iter_t *mem_iter,
		   hsize_t nelmts, void *tconv_buf/*out*/);

    /* Scatter elements from type conversion buffer to application buffer */
    herr_t (*scat)(const void *tconv_buf, size_t elmt_size,
		   const H5S_t *mem_space, H5S_sel_iter_t *mem_iter,
		   hsize_t nelmts, void *buf/*out*/);
} H5S_mconv_t;

typedef struct H5S_conv_t {
    const H5S_fconv_t	*f;
    const H5S_mconv_t	*m;

    /*
     * If there is no data type conversion then it might be possible to
     * transfer data points between application memory and the file in one
     * step without going through the data type conversion buffer.
     *
     * rky 980918
     * If the direct read or write function determines that the transfer
     * must be done indirectly, i.e., through the conversion buffer or
     * (in the case of parallel MPI-IO) in block-by-block transfers
     * then the function returns with the value of must_convert!=0,
     * the function's return value is SUCCEED,
     * and no transfer of data is attempted.
     * Otherwise the direct read or write function returns must_convert 0,
     * with the function's return value being SUCCEED or FAIL
     * depending on whether or not the transfer of data was successful.
     */
    
    /* Read from file to application w/o intermediate scratch buffer */
    herr_t (*read)(H5F_t *f, const struct H5O_layout_t *layout,
		   const struct H5O_pline_t *pline,
		   const struct H5O_fill_t *fill,
		   const struct H5O_efl_t *efl, size_t elmt_size,
		   const H5S_t *file_space, const H5S_t *mem_space,
		   hid_t dxpl_id, void *buf/*out*/);


    /* Write directly from app buffer to file */
    herr_t (*write)(H5F_t *f, const struct H5O_layout_t *layout,
		   const struct H5O_pline_t *pline,
		   const struct H5O_fill_t *fill,
		   const struct H5O_efl_t *efl, size_t elmt_size,
		   const H5S_t *file_space, const H5S_t *mem_space,
		   hid_t dxpl_id, const void *buf);
    
#ifdef H5S_DEBUG
    struct {
	H5_timer_t	scat_timer;		/*time spent scattering	*/
	hsize_t		scat_nbytes;		/*scatter throughput	*/
	hsize_t		scat_ncalls;		/*number of calls	*/
	H5_timer_t	gath_timer;		/*time spent gathering	*/
	hsize_t		gath_nbytes;		/*gather throughput	*/
	hsize_t		gath_ncalls;		/*number of calls	*/
	H5_timer_t	bkg_timer;		/*time for background	*/
	hsize_t		bkg_nbytes;		/*background throughput	*/
	hsize_t		bkg_ncalls;		/*number of calls	*/
	H5_timer_t	read_timer;		/*time for read calls	*/
	hsize_t		read_nbytes;		/*total bytes read	*/
	hsize_t		read_ncalls;		/*number of calls	*/
	H5_timer_t	write_timer;		/*time for write calls	*/
	hsize_t		write_nbytes;		/*total bytes written	*/
	hsize_t		write_ncalls;		/*number of calls	*/
    } stats[2];		/* 0=output, 1=input */
#endif
} H5S_conv_t;

/* Conversion information for the various data space selection types */
H5_DLLVAR const H5S_fconv_t	H5S_POINT_FCONV[];
H5_DLLVAR const H5S_mconv_t	H5S_POINT_MCONV[];
H5_DLLVAR const H5S_fconv_t	H5S_ALL_FCONV[];
H5_DLLVAR const H5S_mconv_t	H5S_ALL_MCONV[];
H5_DLLVAR const H5S_fconv_t	H5S_HYPER_FCONV[];
H5_DLLVAR const H5S_mconv_t	H5S_HYPER_MCONV[];
H5_DLLVAR const H5S_fconv_t	H5S_NONE_FCONV[];
H5_DLLVAR const H5S_mconv_t	H5S_NONE_MCONV[];

/* We get the declaration of H5G_entry_t from the H5Oprivate.h file */

H5_DLL H5S_t *H5S_create(H5S_class_t type);
H5_DLL H5S_t *H5S_copy(const H5S_t *src);
H5_DLL herr_t H5S_close(H5S_t *ds);
H5_DLL H5S_conv_t *H5S_find(const H5S_t *mem_space, const H5S_t *file_space,
                unsigned flags);
H5_DLL H5S_class_t H5S_get_simple_extent_type(const H5S_t *ds);
H5_DLL hssize_t H5S_get_simple_extent_npoints(const H5S_t *ds);
H5_DLL hsize_t H5S_get_npoints_max(const H5S_t *ds);
H5_DLL int H5S_get_simple_extent_ndims(const H5S_t *ds);
H5_DLL int H5S_get_simple_extent_dims(const H5S_t *ds, hsize_t dims[]/*out*/,
					hsize_t max_dims[]/*out*/);
H5_DLL herr_t H5S_modify(struct H5G_entry_t *ent, const H5S_t *space, hid_t dxpl_id);
H5_DLL H5S_t *H5S_read(struct H5G_entry_t *ent, hid_t dxpl_id);
H5_DLL int H5S_cmp(const H5S_t *ds1, const H5S_t *ds2);
H5_DLL htri_t H5S_is_simple(const H5S_t *sdim);
H5_DLL unsigned H5S_nelem(const H5S_t *space);
H5_DLL herr_t H5S_select_hyperslab(H5S_t *space, H5S_seloper_t op,
				    const hssize_t start[],
				    const hsize_t *_stride,
				    const hsize_t count[],
				    const hsize_t *_block);
H5_DLL int H5S_get_hyperslab(const H5S_t *ds, hssize_t offset[]/*out*/,
			       hsize_t size[]/*out*/, hsize_t stride[]/*out*/);
H5_DLL herr_t H5S_select_copy(H5S_t *dst, const H5S_t *src);
H5_DLL herr_t H5S_extent_release(H5S_t *space);
H5_DLL herr_t H5S_select_release(H5S_t *space);
H5_DLL hssize_t H5S_get_select_npoints(const H5S_t *space);
H5_DLL int H5S_extend(H5S_t *space, const hsize_t *size);
H5_DLL herr_t H5S_set_extent_simple(H5S_t *space, int rank,
				     const hsize_t *dims, const hsize_t *max);
H5_DLL htri_t H5S_select_valid(const H5S_t *space);
H5_DLL herr_t H5S_debug(H5F_t *f, hid_t dxpl_id, const void *_mesg, FILE *stream,
			 int indent, int fwidth);
H5_DLL hssize_t H5S_select_serial_size(const H5S_t *space);
H5_DLL herr_t H5S_select_serialize(const H5S_t *space, uint8_t *buf);
H5_DLL herr_t H5S_select_deserialize(H5S_t *space, const uint8_t *buf);
H5_DLL htri_t H5S_select_contiguous(const H5S_t *space);
H5_DLL htri_t H5S_select_single(const H5S_t *space);
H5_DLL htri_t H5S_select_regular(const H5S_t *space);
H5_DLL htri_t H5S_select_shape_same(const H5S_t *space1, const H5S_t *space2);
H5_DLL herr_t H5S_select_iterate(void *buf, hid_t type_id, H5S_t *space,
				  H5D_operator_t op, void *operator_data);
H5_DLL herr_t H5S_sel_iter_release(const H5S_t *space,
				    H5S_sel_iter_t *sel_iter);

#ifdef H5_HAVE_PARALLEL
/* MPI-IO function to read directly from app buffer to file rky980813 */
H5_DLL herr_t H5S_mpio_spaces_read(H5F_t *f,
				    const struct H5O_layout_t *layout,
				    const struct H5O_pline_t *pline,
		                    const struct H5O_fill_t *fill,
				    const struct H5O_efl_t *efl,
				    size_t elmt_size, const H5S_t *file_space,
				    const H5S_t *mem_space, hid_t dxpl_id,
				    void *buf/*out*/);

/* MPI-IO function to write directly from app buffer to file rky980813 */
H5_DLL herr_t H5S_mpio_spaces_write(H5F_t *f,
				     const struct H5O_layout_t *layout,
				     const struct H5O_pline_t *pline,
		                     const struct H5O_fill_t *fill,
				     const struct H5O_efl_t *efl,
				     size_t elmt_size, const H5S_t *file_space,
				     const H5S_t *mem_space, hid_t dxpl_id,
				     const void *buf);

/* MPI-IO function to check if a direct I/O transfer is possible between
 * memory and the file */
H5_DLL htri_t H5S_mpio_opt_possible(const H5S_t *mem_space,
                                     const H5S_t *file_space, const unsigned flags);

#ifndef _H5S_IN_H5S_C
/* Global vars whose value comes from environment variable */
/* (Defined in H5S.c) */
H5_DLLVAR hbool_t		H5S_mpi_opt_types_g;
H5_DLLVAR hbool_t		H5S_mpi_prefer_derived_types_g;
#endif /* _H5S_IN_H5S_C */

#endif /* H5_HAVE_PARALLEL */

#endif /* _H5Sprivate_H */
