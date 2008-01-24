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

/* $Id: H5Pdxpl.c,v 1.1.2.1 2002/08/12 18:14:10 koziol Exp $ */

/* Private header files */
#include "H5private.h"		/* Generic Functions			*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5Iprivate.h"		/* IDs			  		*/
#include "H5Pprivate.h"		/* Property lists		  	*/

#define PABLO_MASK	H5Pdxpl_mask

/* Is the interface initialized? */
static int		interface_initialize_g = 0;
#define INTERFACE_INIT NULL

/* Local types */

/* Local static functions */

#ifdef WANT_H5_V1_2_COMPAT
#ifdef H5_HAVE_PARALLEL

/*-------------------------------------------------------------------------
 * Function:	H5Pset_xfer
 *
 * Signature:	herr_t H5Pset_xfer(hid_t plist_id,
 *		                   H5D_transfer_t data_xfer_mode)
 *
 * Purpose:	Set the transfer mode of the dataset transfer property list.
 *		The list can then be used to control the I/O transfer mode
 *		during dataset accesses.  This function is available only
 *		in the parallel HDF5 library and is not a collective function.
 *
 * Parameters:
 *		hid_t plist_id 
 *		    ID of a dataset transfer property list
 *		H5D_transfer_t data_xfer_mode
 *		    Data transfer modes: 
 *			H5D_XFER_INDEPENDENT 
 *			    Use independent I/O access. 
 *			H5D_XFER_COLLECTIVE 
 *			    Use MPI collective I/O access. 
 *			H5D_XFER_DFLT 
 *			    Use default I/O access.  Currently,
 *			    independent is the default mode.
 *			
 *	     
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Albert Cheng
 *		April 2, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_xfer(hid_t plist_id, H5D_transfer_t data_xfer_mode)
{
    H5FD_mpio_xfer_t	dx_xfer_mode;

    FUNC_ENTER(H5Pset_xfer, FAIL);
    H5TRACE2("e","iDt",plist_id,data_xfer_mode);

    /* Check arguments */
    if (H5P_DATASET_XFER != H5P_get_class(plist_id) ||
            NULL == H5I_object(plist_id)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset transfer property list");
    }
    
    switch (data_xfer_mode) {
        case H5D_XFER_COLLECTIVE:
            dx_xfer_mode = H5FD_MPIO_COLLECTIVE;
            break;

        case H5D_XFER_INDEPENDENT:
        case H5D_XFER_DFLT:
            dx_xfer_mode = H5FD_MPIO_INDEPENDENT;
            break;

        default:
            HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "invalid dataset transfer mode");
    }

    FUNC_LEAVE(H5Pset_dxpl_mpio(plist_id,dx_xfer_mode));
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_xfer
 *
 * Purpose:	Reads the transfer mode current set in the property list.
 *		This function is available only in the parallel HDF5 library
 *		and is not a collective function.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Albert Cheng
 *		April 2, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_xfer(hid_t plist_id, H5D_transfer_t *data_xfer_mode)
{
    H5FD_mpio_xfer_t	dx_xfer_mode;

    FUNC_ENTER (H5Pget_xfer, FAIL);
    H5TRACE2("e","i*Dt",plist_id,data_xfer_mode);

    /* Check arguments */
    if (H5P_DATASET_XFER != H5P_get_class(plist_id) ||
            NULL == H5I_object(plist_id)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset transfer property list");
    }

    if(data_xfer_mode) {
        if(H5Pget_dxpl_mpio(plist_id,&dx_xfer_mode)<0) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
                  "not a dataset transfer property list");
        }
        switch(dx_xfer_mode) {
            H5FD_MPIO_INDEPENDENT:
                *data_xfer_mode = H5D_XFER_INDEPENDENT;
                break;

            H5FD_MPIO_COLLECTIVE:
                *data_xfer_mode = H5D_XFER_COLLECTIVE;
                break;
	    default:
		HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		    "unknown transfer mode");
        } /* end switch */
    } /* end if */

    FUNC_LEAVE (SUCCEED);
}
#endif /* H5_HAVE_PARALLEL */
#endif /* WANT_H5_V1_2_COMPAT */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_buffer
 *
 * Purpose:	Given a dataset transfer property list, set the maximum size
 *		for the type conversion buffer and background buffer and
 *		optionally supply pointers to application-allocated buffers.
 *		If the buffer size is smaller than the entire amount of data
 *		being transfered between application and file, and a type
 *		conversion buffer or background buffer is required then
 *		strip mining will be used.
 *
 *		If TCONV and/or BKG are null pointers then buffers will be
 *		allocated and freed during the data transfer.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Monday, March 16, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_buffer(hid_t plist_id, hsize_t size, void *tconv, void *bkg)
{
    H5D_xfer_t		*plist = NULL;
    
    FUNC_ENTER (H5Pset_buffer, FAIL);
    H5TRACE4("e","ihxx",plist_id,size,tconv,bkg);

    /* Check arguments */
    if (H5P_DATASET_XFER != H5P_get_class (plist_id) ||
            NULL == (plist = H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a dataset transfer property list");
    }
    if (size<=0) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "buffer size must not be zero");
    }

    /* Update property list */
    plist->buf_size = size;
    plist->tconv_buf = tconv;
    plist->bkg_buf = bkg;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_buffer
 *
 * Purpose:	Reads values previously set with H5Pset_buffer().
 *
 * Return:	Success:	Buffer size.
 *
 *		Failure:	0
 *
 * Programmer:	Robb Matzke
 *              Monday, March 16, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
hsize_t
H5Pget_buffer(hid_t plist_id, void **tconv/*out*/, void **bkg/*out*/)
{
    H5D_xfer_t		*plist = NULL;
    
    FUNC_ENTER (H5Pget_buffer, 0);
    H5TRACE3("h","ixx",plist_id,tconv,bkg);

    /* Check arguments */
    if (H5P_DATASET_XFER != H5P_get_class (plist_id) ||
            NULL == (plist = H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, 0,
		       "not a dataset transfer property list");
    }

    /* Return values */
    if (tconv)
        *tconv = plist->tconv_buf;
    if (bkg)
        *bkg = plist->bkg_buf;

    FUNC_LEAVE (plist->buf_size);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_hyper_cache
 *
 * Purpose:	Given a dataset transfer property list, indicate whether to
 *		cache the hyperslab blocks during the I/O (which speeds
 *		things up) and the maximum size of the hyperslab block to
 *		cache.  If a block is smaller than to limit, it may still not
 *		be cached if no memory is available. Setting the limit to 0
 *		indicates no limitation on the size of block to attempt to
 *		cache.
 *
 *		The default is to cache blocks with no limit on block size
 *		for serial I/O and to not cache blocks for parallel I/O
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Monday, September 21, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_hyper_cache(hid_t plist_id, unsigned cache, unsigned limit)
{
    H5D_xfer_t		*plist = NULL;
    
    FUNC_ENTER (H5Pset_hyper_cache, FAIL);
    H5TRACE3("e","iIuIu",plist_id,cache,limit);

    /* Check arguments */
    if (H5P_DATASET_XFER != H5P_get_class (plist_id) ||
            NULL == (plist = H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a dataset transfer property list");
    }

    /* Update property list */
    plist->cache_hyper = (cache>0) ? 1 : 0;
    plist->block_limit = limit;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_hyper_cache
 *
 * Purpose:	Reads values previously set with H5Pset_hyper_cache().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Monday, September 21, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_hyper_cache(hid_t plist_id, unsigned *cache/*out*/,
		   unsigned *limit/*out*/)
{
    H5D_xfer_t		*plist = NULL;
    
    FUNC_ENTER (H5Pget_hyper_cache, FAIL);
    H5TRACE3("e","ixx",plist_id,cache,limit);

    /* Check arguments */
    if (H5P_DATASET_XFER != H5P_get_class (plist_id) ||
            NULL == (plist = H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a dataset transfer property list");
    }

    /* Return values */
    if (cache)
        *cache = plist->cache_hyper;
    if (limit)
        *limit = plist->block_limit;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_preserve
 *
 * Purpose:	When reading or writing compound data types and the
 *		destination is partially initialized and the read/write is
 *		intended to initialize the other members, one must set this
 *		property to TRUE.  Otherwise the I/O pipeline treats the
 *		destination datapoints as completely uninitialized.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Tuesday, March 17, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_preserve(hid_t plist_id, hbool_t status)
{
    H5D_xfer_t		*plist = NULL;
    
    FUNC_ENTER (H5Pset_preserve, FAIL);
    H5TRACE2("e","ib",plist_id,status);

    /* Check arguments */
    if (H5P_DATASET_XFER != H5P_get_class (plist_id) ||
            NULL == (plist = H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a dataset transfer property list");
    }

    /* Update property list */
    plist->need_bkg = status ? H5T_BKG_YES : H5T_BKG_NO;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_preserve
 *
 * Purpose:	The inverse of H5Pset_preserve()
 *
 * Return:	Success:	TRUE or FALSE
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *              Tuesday, March 17, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
H5Pget_preserve(hid_t plist_id)
{
    H5D_xfer_t		*plist = NULL;
    
    FUNC_ENTER (H5Pget_preserve, FAIL);
    H5TRACE1("Is","i",plist_id);

    /* Check arguments */
    if (H5P_DATASET_XFER != H5P_get_class (plist_id) ||
            NULL == (plist = H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a dataset transfer property list");
    }

    FUNC_LEAVE (plist->need_bkg?TRUE:FALSE);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_btree_ratios
 *
 * Purpose:	Queries B-tree split ratios.  See H5Pset_btree_ratios().
 *
 * Return:	Success:	Non-negative with split ratios returned through
 *				the non-null arguments.
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *              Monday, September 28, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_btree_ratios(hid_t plist_id, double *left/*out*/, double *middle/*out*/,
		    double *right/*out*/)
{
    H5D_xfer_t		*plist = NULL;

    FUNC_ENTER(H5Pget_btree_ratios, FAIL);
    H5TRACE4("e","ixxx",plist_id,left,middle,right);

    /* Check arguments */
    if (H5P_DATASET_XFER!=H5P_get_class(plist_id) ||
            NULL==(plist=H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset transfer property list");
    }

    /* Get values */
    if (left)
        *left = plist->split_ratios[0];
    if (middle)
        *middle = plist->split_ratios[1];
    if (right)
        *right = plist->split_ratios[2];

    FUNC_LEAVE(SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_btree_ratios
 *
 * Purpose:	Sets B-tree split ratios for a dataset transfer property
 *		list. The split ratios determine what percent of children go
 *		in the first node when a node splits.  The LEFT ratio is
 *		used when the splitting node is the left-most node at its
 *		level in the tree; the RIGHT ratio is when the splitting node
 *		is the right-most node at its level; and the MIDDLE ratio for
 *		all other cases.  A node which is the only node at its level
 *		in the tree uses the RIGHT ratio when it splits.  All ratios
 *		are real numbers between 0 and 1, inclusive.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Monday, September 28, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_btree_ratios(hid_t plist_id, double left, double middle,
		    double right)
{
    H5D_xfer_t		*plist = NULL;

    FUNC_ENTER(H5Pget_btree_ratios, FAIL);
    H5TRACE4("e","iddd",plist_id,left,middle,right);

    /* Check arguments */
    if (H5P_DATASET_XFER!=H5P_get_class(plist_id) ||
            NULL==(plist=H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset transfer property list");
    }
    if (left<0.0 || left>1.0 || middle<0.0 || middle>1.0 ||
            right<0.0 || right>1.0) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL,
		      "split ratio must satisfy 0.0<=X<=1.0");
    }
    
    /* Set values */
    plist->split_ratios[0] = left;
    plist->split_ratios[1] = middle;
    plist->split_ratios[2] = right;

    FUNC_LEAVE(SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_vlen_mem_manager
 *
 * Purpose:	Sets the memory allocate/free pair for VL datatypes.  The
 *		allocation routine is called when data is read into a new
 *		array and the free routine is called when H5Dvlen_reclaim is
 *		called.  The alloc_info and free_info are user parameters
 *		which are passed to the allocation and freeing functions
 *		respectively.  To reset the allocate/free functions to the
 *		default setting of using the system's malloc/free functions,
 *		call this routine with alloc_func and free_func set to NULL.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 1, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_vlen_mem_manager(hid_t plist_id, H5MM_allocate_t alloc_func,
        void *alloc_info, H5MM_free_t free_func, void *free_info)
{
    H5D_xfer_t		*plist = NULL;
    
    FUNC_ENTER(H5Pset_vlen_mem_manager, FAIL);
    H5TRACE5("e","ixxxx",plist_id,alloc_func,alloc_info,free_func,free_info);

    /* Check arguments */
    if (H5P_DATASET_XFER!=H5P_get_class(plist_id) ||
            NULL==(plist=H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset transfer property list");
    }

    /* Update property list */
    plist->vlen_alloc = alloc_func;
    plist->alloc_info = alloc_info;
    plist->vlen_free = free_func;
    plist->free_info = free_info;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_vlen_mem_manager
 *
 * Purpose:	The inverse of H5Pset_vlen_mem_manager()
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 1, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_vlen_mem_manager(hid_t plist_id, H5MM_allocate_t *alloc_func/*out*/,
			void **alloc_info/*out*/,
			H5MM_free_t *free_func/*out*/,
			void **free_info/*out*/)
{
    H5D_xfer_t		*plist = NULL;
    
    FUNC_ENTER(H5Pget_vlen_mem_manager, FAIL);
    H5TRACE5("e","ixxxx",plist_id,alloc_func,alloc_info,free_func,free_info);

    /* Check arguments */
    if (H5P_DATASET_XFER!=H5P_get_class(plist_id) ||
            NULL==(plist=H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset transfer property list");
    }

    if(alloc_func!=NULL)
        *alloc_func= plist->vlen_alloc;
    if(alloc_info!=NULL)
        *alloc_info= plist->alloc_info;
    if(free_func!=NULL)
        *free_func= plist->vlen_free;
    if(free_info!=NULL)
        *free_info= plist->free_info;

    FUNC_LEAVE (SUCCEED);
}

