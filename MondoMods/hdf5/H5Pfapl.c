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

/* $Id: H5Pfapl.c,v 1.1.2.2 2003/03/04 14:26:15 koziol Exp $ */

/* Private header files */
#include "H5private.h"		/* Generic Functions			*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5FDprivate.h"	/* File drivers				*/
#include "H5Iprivate.h"		/* IDs			  		*/
#include "H5Pprivate.h"		/* Property lists		  	*/

/* Default file driver - see H5Pget_driver() */
#include "H5FDsec2.h"		/* Posix unbuffered I/O	file driver	*/

#ifdef WANT_H5_V1_2_COMPAT
/* Other predefined file drivers */
#include "H5FDcore.h"		/* Files stored entirely in memory	*/
#include "H5FDfamily.h"		/* File families 			*/
#include "H5FDmpio.h"		/* Parallel files using MPI-2 I/O	*/
#include "H5FDstdio.h"		/* Standard C buffered I/O		*/
#include "H5FDsrb.h"            /* Remote access using SRB              */
#include "H5FDgass.h"		/* Remote files using GASS I/O		*/
#include "H5FDstream.h"         /* in-memory files streamed via sockets */
#include "H5FDmulti.h"		/* Usage-partitioned file family	*/
#include "H5FDlog.h"            /* sec2 driver with I/O logging (for debugging) */
#endif /* WANT_H5_V1_2_COMPAT */

#define PABLO_MASK	H5Pfapl_mask

/* Is the interface initialized? */
static int		interface_initialize_g = 0;
#define INTERFACE_INIT NULL

/* Local types */

/* Local static functions */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_alignment
 *
 * Purpose:	Sets the alignment properties of a file access property list
 *		so that any file object >= THRESHOLD bytes will be aligned on
 *		an address which is a multiple of ALIGNMENT.  The addresses
 *		are relative to the end of the user block; the alignment is
 *		calculated by subtracting the user block size from the
 *		absolute file address and then adjusting the address to be a
 *		multiple of ALIGNMENT.
 *
 *		Default values for THRESHOLD and ALIGNMENT are one, implying
 *		no alignment.  Generally the default values will result in
 *		the best performance for single-process access to the file.
 *		For MPI-IO and other parallel systems, choose an alignment
 *		which is a multiple of the disk block size.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Tuesday, June  9, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_alignment(hid_t fapl_id, hsize_t threshold, hsize_t alignment)
{
    H5F_access_t	*fapl = NULL;
    
    FUNC_ENTER (H5Pset_alignment, FAIL);
    H5TRACE3("e","ihh",fapl_id,threshold,alignment);

    /* Check args */
    if (H5P_FILE_ACCESS != H5P_get_class (fapl_id) ||
            NULL == (fapl = H5I_object (fapl_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }
    if (alignment<1) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "alignment must be positive");
    }

    /* Set values */
    fapl->threshold = threshold;
    fapl->alignment = alignment;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_alignment
 *
 * Purpose:	Returns the current settings for alignment properties from a
 *		file access property list.  The THRESHOLD and/or ALIGNMENT
 *		pointers may be null pointers.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Tuesday, June  9, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_alignment(hid_t fapl_id, hsize_t *threshold/*out*/,
		  hsize_t *alignment/*out*/)
{
    H5F_access_t	*fapl = NULL;

    FUNC_ENTER (H5Pget_alignment, FAIL);
    H5TRACE3("e","ixx",fapl_id,threshold,alignment);

    /* Check args */
    if (H5P_FILE_ACCESS != H5P_get_class (fapl_id) ||
            NULL == (fapl = H5I_object (fapl_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }

    /* Get values */
    if (threshold)
        *threshold = fapl->threshold;
    if (alignment)
        *alignment = fapl->alignment;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_driver
 *
 * Purpose:	Set the file driver (DRIVER_ID) for a file access or data
 *		transfer property list (PLIST_ID) and supply an optional
 *		struct containing the driver-specific properites
 *		(DRIVER_INFO).  The driver properties will be copied into the
 *		property list and the reference count on the driver will be
 *		incremented, allowing the caller to close the driver ID but
 *		still use the property list.
 *
 * Return:	Success:	Non-negative
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *              Tuesday, August  3, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_driver(hid_t plist_id, hid_t driver_id, const void *driver_info)
{
    H5F_access_t	*fapl=NULL;
    H5D_xfer_t		*dxpl=NULL;
    
    FUNC_ENTER(H5Pset_driver, FAIL);
    H5TRACE3("e","iix",plist_id,driver_id,driver_info);

    if (H5I_VFL!=H5I_get_type(driver_id) ||
            NULL==H5I_object(driver_id)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a file driver ID");
    }

    if (H5P_FILE_ACCESS==H5P_get_class(plist_id)) {
        if (NULL==(fapl=H5I_object(plist_id))) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
                  "not a file access property list");
        }
	
        /* Remove old driver */
        assert(fapl->driver_id>=0);
        H5FD_fapl_free(fapl->driver_id, fapl->driver_info);
        H5I_dec_ref(fapl->driver_id);

        /* Add new driver */
        H5I_inc_ref(driver_id);
        fapl->driver_id = driver_id;
        fapl->driver_info = H5FD_fapl_copy(driver_id, driver_info);
        
    } else if (H5P_DATASET_XFER==H5P_get_class(plist_id)) {
        if (NULL==(dxpl=H5I_object(plist_id))) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
                  "not a file access property list");
        }

        /* Remove old driver */
        assert(dxpl->driver_id>=0);
        H5FD_dxpl_free(dxpl->driver_id, dxpl->driver_info);
        H5I_dec_ref(dxpl->driver_id);

        /* Add new driver */
        H5I_inc_ref(driver_id);
        dxpl->driver_id = driver_id;
        dxpl->driver_info = H5FD_dxpl_copy(driver_id, driver_info);
        
    } else {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access or data transfer property list");
    }

    FUNC_LEAVE(SUCCEED);
}

#ifdef WANT_H5_V1_2_COMPAT

/*-------------------------------------------------------------------------
 * Function:	H5Pget_driver
 *
 * Purpose:	Return the ID of the low-level file driver.  PLIST_ID should
 *		be a file access property list.
 *
 * Return:	Success:	A low-level driver ID
 *
 *		Failure:	H5F_LOW_ERROR (a negative value)
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
H5F_driver_t
H5Pget_driver(hid_t plist_id)
{
    H5F_access_t	*plist = NULL;
    H5F_driver_t    ret_value=H5F_LOW_ERROR;

    FUNC_ENTER (H5Pget_driver, H5F_LOW_ERROR);
    H5TRACE1("Fd","i",plist_id);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class (plist_id) ||
            NULL == (plist=H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, H5F_LOW_ERROR,
		       "not a file access property list");
    }
    
    if(plist->driver_id==H5FD_SEC2 || plist->driver_id==H5P_DEFAULT)
        ret_value=H5F_LOW_SEC2;
    else if(plist->driver_id==H5FD_STDIO)
        ret_value=H5F_LOW_STDIO;
    else if(plist->driver_id==H5FD_MPIO)
        ret_value=H5F_LOW_MPIO;
    else if(plist->driver_id==H5FD_CORE)
        ret_value=H5F_LOW_CORE;
    else if(plist->driver_id==H5FD_FAMILY)
        ret_value=H5F_LOW_FAMILY;
    else if(plist->driver_id==H5FD_MULTI) { /* Need to check if it's a split or multi file */
        H5FD_mem_t		mt;
        H5FD_mem_t		memb_map[H5FD_MEM_NTYPES];
        haddr_t         memb_addr[H5FD_MEM_NTYPES];
        unsigned           multi=0;

        /* Get the information from the multi file driver */
        if (H5Pget_fapl_multi(plist_id,memb_map,NULL,NULL,memb_addr,NULL)<0) {
            HRETURN_ERROR (H5E_PLIST, H5E_NOTFOUND, H5F_LOW_ERROR,
                   "can't get multi file information");
        }

        /* Check whether all of the meta data is in one file & the raw data in another */
        for (mt=H5FD_MEM_DEFAULT; mt<H5FD_MEM_NTYPES; mt++) {
            if(mt==H5FD_MEM_DRAW) {
                if(memb_map[mt]!=H5FD_MEM_DRAW) {
                    multi=1;
                    break;
                } /* end if */
            } /* end if */
            else {
                if(memb_map[mt]!=H5FD_MEM_SUPER) {
                    multi=1;
                    break;
                } /* end if */
            } /* end else */
        } /* end for */

        /* Check further if things look like a split file currently */
        if(!multi) {
            if(memb_addr[H5FD_MEM_SUPER]!=0 || memb_addr[H5FD_MEM_DRAW] != HADDR_MAX/2)
                multi=1;
        } /* end if */

        if(multi)
            ret_value=H5F_LOW_ERROR;    /* v1.2 didn't have multi-file driver */
        else
            ret_value=H5F_LOW_SPLIT;
    } /* end if */
    else
        ret_value=H5F_LOW_ERROR;    /* error, or driver unknown to v1.2 */

    FUNC_LEAVE (ret_value);
}

#else /* WANT_H5_V1_2_COMPAT */

/*-------------------------------------------------------------------------
 * Function:	H5Pget_driver
 *
 * Purpose:	Return the ID of the low-level file driver.  PLIST_ID should
 *		be a file access property list or data transfer propert list.
 *
 * Return:	Success:	A low-level driver ID which is the same ID
 *				used when the driver was set for the property
 *				list. The driver ID is only valid as long as
 *				the file driver remains registered.
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 * Modifications:
 *		Robb Matzke, 1999-08-03
 *		Rewritten to use the virtual file layer.
 *
 * 		Robb Matzke, 1999-08-05
 *		If the driver ID is H5FD_VFD_DEFAULT then substitute the current value of
 *		H5FD_SEC2.
 *
 * 		Quincey Koziol 2000-11-28
 *		Added internal function..
 *-------------------------------------------------------------------------
 */
hid_t
H5Pget_driver(hid_t plist_id)
{
    hid_t		ret_value=-1;

    FUNC_ENTER (H5Pget_driver, FAIL);
    H5TRACE1("i","i",plist_id);

    ret_value = H5P_get_driver(plist_id);

    FUNC_LEAVE(ret_value);
}
#endif /* WANT_H5_V1_2_COMPAT */


/*-------------------------------------------------------------------------
 * Function:	H5P_get_driver
 *
 * Purpose:	Return the ID of the low-level file driver.  PLIST_ID should
 *		be a file access property list or data transfer propert list.
 *
 * Return:	Success:	A low-level driver ID which is the same ID
 *				used when the driver was set for the property
 *				list. The driver ID is only valid as long as
 *				the file driver remains registered.
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 * Modifications:
 *		Robb Matzke, 1999-08-03
 *		Rewritten to use the virtual file layer.
 *
 * 		Robb Matzke, 1999-08-05
 *		If the driver ID is H5FD_VFD_DEFAULT then substitute the current value
 *      of H5FD_SEC2.
 *
 * 		Quincey Koziol 2000-11-28
 *		Added internal function..
 *-------------------------------------------------------------------------
 */
hid_t
H5P_get_driver(hid_t plist_id)
{
    H5F_access_t	*fapl=NULL;
    H5D_xfer_t		*dxpl=NULL;
    hid_t		ret_value=-1;

    FUNC_ENTER (H5P_get_driver, FAIL);

    if (H5P_FILE_ACCESS==H5P_get_class(plist_id) &&
            (fapl=H5I_object(plist_id))) {
        ret_value = fapl->driver_id;
	
    } else if (H5P_DATASET_XFER==H5P_get_class(plist_id) &&
	       (dxpl=H5I_object(plist_id))) {
        ret_value = dxpl->driver_id;
	
    } else {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access or data transfer property list");
    }

    if (H5FD_VFD_DEFAULT==ret_value)
        ret_value = H5FD_SEC2;

    FUNC_LEAVE(ret_value);
}

/*-------------------------------------------------------------------------
 * Function:	H5Pget_driver_info
 *
 * Purpose:	Returns a pointer directly to the file driver-specific
 *		information of a file access or data transfer property list.
 *
 * Return:	Success:	Ptr to *uncopied* driver specific data
 *				structure if any.
 *
 *		Failure:	NULL. Null is also returned if the driver has
 *				not registered any driver-specific properties
 *				although no error is pushed on the stack in
 *				this case.
 *
 * Programmer:	Robb Matzke
 *              Wednesday, August  4, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void *
H5Pget_driver_info(hid_t plist_id)
{
    H5F_access_t	*fapl=NULL;
    H5D_xfer_t		*dxpl=NULL;
    void		*ret_value=NULL;

    FUNC_ENTER(H5Pget_driver_info, NULL);

    if (H5P_FILE_ACCESS==H5P_get_class(plist_id) &&
            (fapl=H5I_object(plist_id))) {
        ret_value = fapl->driver_info;
	
    } else if (H5P_DATASET_XFER==H5P_get_class(plist_id) &&
	       (dxpl=H5I_object(plist_id))) {
        ret_value = dxpl->driver_info;
	
    } else {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, NULL,
		      "not a file access or data transfer property list");
    }
    
    FUNC_LEAVE(ret_value);
}

#ifdef WANT_H5_V1_2_COMPAT

/*-------------------------------------------------------------------------
 * Function:	H5Pset_stdio
 *
 * Purpose:	Set the low level file driver to use the functions declared
 *		in the stdio.h file: fopen(), fseek() or fseek64(), fread(),
 *		fwrite(), and fclose().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 19, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_stdio(hid_t plist_id)
{
    FUNC_ENTER (H5Pset_stdio, FAIL);
    H5TRACE1("e","i",plist_id);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class(plist_id) ||
            NULL == H5I_object(plist_id)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }

    FUNC_LEAVE(H5Pset_fapl_stdio(plist_id));
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_stdio
 *
 * Purpose:	If the file access property list is set to the stdio driver
 *		then this function returns zero; otherwise it returns a
 *		negative value.	 In the future, additional arguments may be
 *		added to this function to match those added to H5Pset_stdio().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_stdio(hid_t plist_id)
{
    herr_t	ret_value=FAIL;

    FUNC_ENTER (H5Pget_stdio, FAIL);
    H5TRACE1("e","i",plist_id);

    /* Check arguments and test driver */
    if (H5P_FILE_ACCESS==H5P_get_class(plist_id) &&
            (H5FD_STDIO == H5P_get_driver(plist_id))) {
        ret_value=SUCCEED;
    }
    else {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }

    FUNC_LEAVE (ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_sec2
 *
 * Purpose:	Set the low-level file driver to use the functions declared
 *		in the unistd.h file: open(), lseek() or lseek64(), read(),
 *		write(), and close().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 19, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_sec2(hid_t plist_id)
{
    FUNC_ENTER (H5Pset_sec2, FAIL);
    H5TRACE1("e","i",plist_id);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class(plist_id) ||
            NULL == H5I_object(plist_id)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }

    FUNC_LEAVE(H5Pset_fapl_sec2(plist_id));
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_sec2
 *
 * Purpose:	If the file access property list is set to the sec2 driver
 *		then this function returns zero; otherwise it returns a
 *		negative value.	 In the future, additional arguments may be
 *		added to this function to match those added to H5Pset_sec2().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_sec2(hid_t plist_id)
{
    herr_t	ret_value=FAIL;

    FUNC_ENTER (H5Pget_sec2, FAIL);
    H5TRACE1("e","i",plist_id);

    /* Check arguments and test driver */
    if (H5P_FILE_ACCESS==H5P_get_class(plist_id) &&
            (H5FD_SEC2 == H5P_get_driver(plist_id))) {
        ret_value=SUCCEED;
    }
    else {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }

    FUNC_LEAVE (ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_core
 *
 * Purpose:	Set the low-level file driver to use malloc() and free().
 *		This driver is restricted to temporary files which are not
 *		larger than the amount of virtual memory available. The
 *		INCREMENT argument determines the file block size and memory
 *		will be allocated in multiples of INCREMENT bytes. A liberal
 *		INCREMENT results in fewer calls to realloc() and probably
 *		less memory fragmentation.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 19, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_core(hid_t plist_id, size_t increment)
{
    FUNC_ENTER (H5Pset_core, FAIL);
    H5TRACE2("e","iz",plist_id,increment);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class(plist_id) ||
            NULL == H5I_object(plist_id)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }
    if (increment<1) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "increment must be positive");
    }

    FUNC_LEAVE(H5Pset_fapl_core(plist_id,increment,0));
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_core
 *
 * Purpose:	If the file access property list is set to the core driver
 *		then this function returns zero; otherwise it returns a
 *		negative value.	 On success, the block size is returned
 *		through the INCREMENT argument if it isn't the null pointer.
 *		In the future, additional arguments may be added to this
 *		function to match those added to H5Pset_core().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_core(hid_t plist_id, size_t *increment/*out*/)
{
    herr_t	ret_value=FAIL;

    FUNC_ENTER (H5Pget_core, FAIL);
    H5TRACE2("e","ix",plist_id,increment);

    /* Check arguments */
    if (H5P_FILE_ACCESS==H5P_get_class(plist_id) &&
            (H5FD_CORE == H5P_get_driver(plist_id)) &&
            H5Pget_fapl_core(plist_id,increment,NULL)>=0) {
        ret_value=SUCCEED;
    }
    else {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }

    FUNC_LEAVE (ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_split
 *
 * Purpose:	Set the low-level driver to split meta data from raw data,
 *		storing meta data in one file and raw data in another file.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 19, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_split(hid_t plist_id, const char *meta_ext, hid_t meta_plist_id,
	      const char *raw_ext, hid_t raw_plist_id)
{
    FUNC_ENTER (H5Pset_split, FAIL);
    H5TRACE5("e","isisi",plist_id,meta_ext,meta_plist_id,raw_ext,raw_plist_id);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class(plist_id) ||
            NULL == H5I_object(plist_id)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }
    if (H5P_DEFAULT!=meta_plist_id &&
            (H5P_FILE_ACCESS != H5P_get_class(meta_plist_id) ||
             NULL == H5I_object(meta_plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }
    if (H5P_DEFAULT!=raw_plist_id &&
            (H5P_FILE_ACCESS != H5P_get_class(raw_plist_id) ||
             NULL == H5I_object(raw_plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }

    /* Set driver */
    FUNC_LEAVE (H5Pset_fapl_split(plist_id,meta_ext,meta_plist_id,raw_ext,raw_plist_id));
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_split
 *
 * Purpose:	If the file access property list is set to the sec2 driver
 *		then this function returns zero; otherwise it returns a
 *		negative value.	 On success, at most META_EXT_SIZE characters
 *		are copied to the META_EXT buffer if non-null and at most
 *		RAW_EXT_SIZE characters are copied to the RAW_EXT buffer if
 *		non-null.  If the actual extension is larger than the number
 *		of characters requested then the buffer will not be null
 *		terminated (that is, behavior like strncpy()).	In addition,
 *		if META_PROPERTIES and/or RAW_PROPERTIES are non-null then
 *		the file access property list of the meta file and/or raw
 *		file is copied and its OID returned through these arguments.
 *		In the future, additional arguments may be added to this
 *		function to match those added to H5Pset_sec2().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_split(hid_t plist_id, size_t meta_ext_size, char *meta_ext/*out*/,
	      hid_t *meta_properties/*out*/, size_t raw_ext_size,
	      char *raw_ext/*out*/, hid_t *raw_properties/*out*/)
{
    H5FD_mem_t		mt;
    H5FD_mem_t		_memb_map[H5FD_MEM_NTYPES];
    hid_t		_memb_fapl[H5FD_MEM_NTYPES];
    char		_memb_name[H5FD_MEM_NTYPES][16];
    char		*_memb_name_ptrs[H5FD_MEM_NTYPES];
    haddr_t		_memb_addr[H5FD_MEM_NTYPES];

    FUNC_ENTER (H5Pget_split, FAIL);
    H5TRACE7("e","izxxzxx",plist_id,meta_ext_size,meta_ext,meta_properties,
             raw_ext_size,raw_ext,raw_properties);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class (plist_id) ||
            NULL == H5I_object (plist_id)) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }
    if (H5FD_MULTI != H5P_get_driver(plist_id)) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "the split driver is not set");
    }

    /* Reset output args for error handling */
    if (meta_ext && meta_ext_size>0)
        *meta_ext = '\0';
    if (raw_ext && raw_ext_size>0)
        *raw_ext = '\0';
    if (meta_properties)
        *meta_properties = FAIL;
    if (raw_properties)
        *raw_properties = FAIL;

    /* Set up the member extention pointers */
	for (mt=H5FD_MEM_DEFAULT; mt<H5FD_MEM_NTYPES; mt++)
	    _memb_name_ptrs[mt] = &_memb_name[mt][0];

    /* Get the information from the multi file driver */
    if (H5Pget_fapl_multi(plist_id,_memb_map,_memb_fapl,_memb_name_ptrs,_memb_addr,NULL)<0) {
        HRETURN_ERROR (H5E_PLIST, H5E_NOTFOUND, FAIL,
		       "can't get split file information");
    }

    /* Output arguments */
    if (meta_ext && meta_ext_size>0) {
        if (_memb_name[H5FD_MEM_SUPER]) {
            HDstrncpy (meta_ext, _memb_name[H5FD_MEM_SUPER], meta_ext_size);
        } else {
            HDstrncpy (meta_ext, ".meta", meta_ext_size);
        }
    }
    if (raw_ext && raw_ext_size>0) {
        if (_memb_name[H5FD_MEM_DRAW]) {
            HDstrncpy (raw_ext, _memb_name[H5FD_MEM_DRAW], raw_ext_size);
        } else {
            HDstrncpy (raw_ext, ".raw", raw_ext_size);
        }
    }
    if (meta_properties) {
        *meta_properties = _memb_fapl[H5FD_MEM_SUPER];
    }
    if (raw_properties) {
        *raw_properties = _memb_fapl[H5FD_MEM_DRAW];
    }
    
    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_family
 *
 * Purpose:	Sets the low-level driver to stripe the hdf5 address space
 *		across a family of files.  The MEMB_SIZE argument indicates
 *		the size in bytes of each family member and is only
 *		meaningful when creating new files or opening families that
 *		have only one member.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 19, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_family(hid_t plist_id, hsize_t memb_size, hid_t memb_plist_id)
{
    FUNC_ENTER (H5Pset_family, FAIL);
    H5TRACE3("e","ihi",plist_id,memb_size,memb_plist_id);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class(plist_id) ||
            NULL == H5I_object(plist_id)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }
    if (memb_size && memb_size<1024) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADRANGE, FAIL,
		       "family member size is too small");
    }
    if (H5P_DEFAULT!=memb_plist_id &&
            (H5P_FILE_ACCESS != H5P_get_class(memb_plist_id) ||
             NULL == H5I_object(memb_plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }

    /* Set driver */
    FUNC_LEAVE(H5Pset_fapl_family(plist_id,memb_size,memb_plist_id));
}
    

/*-------------------------------------------------------------------------
 * Function:	H5Pget_family
 *
 * Purpose:	If the file access property list is set to the family driver
 *		then this function returns zero; otherwise it returns a
 *		negative value.	 On success, if MEMB_PLIST_ID is a non-null
 *		pointer it will be initialized with the id of an open
 *		property list: the file access property list for the family
 *		members.  In the future, additional arguments may be added to
 *		this function to match those added to H5Pset_family().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_family(hid_t plist_id, hsize_t *memb_size/*out*/,
	       hid_t *memb_plist_id/*out*/)
{
    FUNC_ENTER (H5Pget_family, FAIL);
    H5TRACE3("e","ixx",plist_id,memb_size,memb_plist_id);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class (plist_id) ||
            NULL == H5I_object (plist_id)) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }
    if (H5FD_FAMILY == H5P_get_driver(plist_id)) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "the family driver is not set");
    }

    /* Retrieve args */
    FUNC_LEAVE (H5Pget_fapl_family(plist_id,memb_size,memb_plist_id));
}

#ifdef H5_HAVE_PARALLEL

/*-------------------------------------------------------------------------
 * Function:	H5Pset_mpi
 *
 * Signature:	herr_t H5Pset_mpi(hid_t plist_id, MPI_Comm comm, MPI_Info info)
 *
 * Purpose:	Store the access mode for MPIO call and the user supplied
 *		communicator and info in the access property list which can
 *		then be used to open file.  This function is available only
 *		in the parallel HDF5 library and is not a collective
 *		function.
 *
 * Parameters:
 *		hid_t plist_id 
 *		    ID of property list to modify 
 *		MPI_Comm comm 
 *		    MPI communicator to be used for file open as defined in
 *		    MPI_FILE_OPEN of MPI-2.  This function  does not make a
 *		    duplicated communicator. Any modification to comm after
 *		    this function call returns may have undetermined effect
 *		    to the access property list.  Users should call this
 *		    function again to setup the property list.
 *		MPI_Info info 
 *		    MPI info object to be used for file open as defined in
 *		    MPI_FILE_OPEN of MPI-2.  This function  does not make a
 *		    duplicated info. Any modification to info after
 *		    this function call returns may have undetermined effect
 *		    to the access property list.  Users should call this
 *		    function again to setup the property list.
 *	     
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Albert Cheng
 *		Feb 3, 1998
 *
 * Modifications:
 *
 *	Robb Matzke, 18 Feb 1998
 *	Check all arguments before the property list is updated so we don't
 *	leave the property list in a bad state if something goes wrong.  Also,
 *	the property list data type changed to allow more generality so all
 *	the mpi-related stuff is in the `u.mpi' member.  The `access_mode'
 *	will contain only mpi-related flags defined in H5Fpublic.h.
 *
 *	Albert Cheng, Apr 16, 1998
 *	Removed the access_mode argument.  The access_mode is changed
 *	to be controlled by data transfer property list during data
 *	read/write calls.
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_mpi(hid_t plist_id, MPI_Comm comm, MPI_Info info)
{
    FUNC_ENTER(H5Pset_mpi, FAIL);
    H5TRACE3("e","iMcMi",plist_id,comm,info);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class(plist_id) ||
            NULL == H5I_object(plist_id)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }
    
    /* Set driver */
    FUNC_LEAVE(H5Pset_fapl_mpio(plist_id,comm,info));
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_mpi
 *
 * Purpose:	If the file access property list is set to the mpi driver
 *		then this function returns zero; otherwise it returns a
 *		negative value.	 In the future, additional arguments may be
 *		added to this function to match those added to H5Pset_mpi().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Thursday, February 26, 1998
 *
 * Modifications:
 *
 *	Albert Cheng, Apr 16, 1998
 *	Removed the access_mode argument.  The access_mode is changed
 *	to be controlled by data transfer property list during data
 *	read/write calls.
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_mpi(hid_t plist_id, MPI_Comm *comm, MPI_Info *info)
{
    FUNC_ENTER (H5Pget_mpi, FAIL);
    H5TRACE3("e","i*Mc*Mi",plist_id,comm,info);

    /* Check arguments */
    if (H5P_FILE_ACCESS != H5P_get_class (plist_id) ||
            NULL == H5I_object (plist_id)) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }
    if (H5FD_MPIO == H5P_get_driver(plist_id)) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "the mpi driver is not set");
    }

    /* Retrieve args */
    FUNC_LEAVE (H5Pget_fapl_mpio(plist_id,comm,info));
}

#endif /* H5_HAVE_PARALLEL */
#endif /* WANT_H5_V1_2_COMPAT */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_cache
 *
 * Purpose:	Set the number of objects in the meta data cache and the
 *		maximum number of chunks and bytes in the raw data chunk
 *		cache.
 *
 * 		The RDCC_W0 value should be between 0 and 1 inclusive and
 *		indicates how much chunks that have been fully read or fully
 *		written are favored for preemption.  A value of zero means
 *		fully read or written chunks are treated no differently than
 *		other chunks (the preemption is strictly LRU) while a value
 *		of one means fully read chunks are always preempted before
 *		other chunks.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Tuesday, May 19, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_cache(hid_t plist_id, int mdc_nelmts,
	     int rdcc_nelmts, size_t rdcc_nbytes, double rdcc_w0)
{
    H5F_access_t	*fapl = NULL;
    
    FUNC_ENTER (H5Pset_cache, FAIL);
    H5TRACE5("e","iIsIszd",plist_id,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,
             rdcc_w0);

    /* Check arguments */
    if (H5P_FILE_ACCESS!=H5P_get_class (plist_id) ||
            NULL==(fapl=H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }
    if (mdc_nelmts<0) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "meta data cache size must be non-negative");
    }
    if (rdcc_nelmts<0) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL,
		      "raw data chunk cache nelmts must be non-negative");
    }
    if (rdcc_w0<0.0 || rdcc_w0>1.0) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "raw data cache w0 value must be between 0.0 and 1.0 "
		       "inclusive");
    }

    /* Set sizes */
    fapl->mdc_nelmts = mdc_nelmts;
    fapl->rdcc_nelmts = rdcc_nelmts;
    fapl->rdcc_nbytes = rdcc_nbytes;
    fapl->rdcc_w0 = rdcc_w0;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_cache
 *
 * Purpose:	Retrieves the maximum possible number of elements in the meta
 *		data cache and the maximum possible number of elements and
 *		bytes and the RDCC_W0 value in the raw data chunk cache.  Any
 *		(or all) arguments may be null pointers in which case the
 *		corresponding datum is not returned.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Tuesday, May 19, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_cache(hid_t plist_id, int *mdc_nelmts,
	     int *rdcc_nelmts, size_t *rdcc_nbytes, double *rdcc_w0)
{
    H5F_access_t	*fapl = NULL;
    
    FUNC_ENTER (H5Pget_cache, FAIL);
    H5TRACE5("e","i*Is*Is*z*d",plist_id,mdc_nelmts,rdcc_nelmts,rdcc_nbytes,
             rdcc_w0);

    /* Check arguments */
    if (H5P_FILE_ACCESS!=H5P_get_class (plist_id) ||
            NULL==(fapl=H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }

    /* Get sizes */
    if (mdc_nelmts)
        *mdc_nelmts = fapl->mdc_nelmts;
    if (rdcc_nelmts)
        *rdcc_nelmts = fapl->rdcc_nelmts;
    if (rdcc_nbytes)
        *rdcc_nbytes = fapl->rdcc_nbytes;
    if (rdcc_w0)
        *rdcc_w0 = fapl->rdcc_w0;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_gc_references
 *
 * Purpose:	Sets the flag for garbage collecting references for the file.
 *		Dataset region references (and other reference types
 *		probably) use space in the file heap.  If garbage collection
 *		is on and the user passes in an uninitialized value in a
 *		reference structure, the heap might get corrupted.  When
 *		garbage collection is off however and the user re-uses a
 *		reference, the previous heap block will be orphaned and not
 *		returned to the free heap space.  When garbage collection is
 *		on, the user must initialize the reference structures to 0 or
 *		risk heap corruption.
 *
 *		Default value for garbage collecting references is off, just
 *		to be on the safe side.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *		June, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_gc_references(hid_t fapl_id, unsigned gc_ref)
{
    H5F_access_t	*fapl = NULL;
    
    FUNC_ENTER(H5Pset_gc_references, FAIL);
    H5TRACE2("e","iIu",fapl_id,gc_ref);

    /* Check args */
    if (H5P_FILE_ACCESS!=H5P_get_class(fapl_id) ||
        NULL==(fapl=H5I_object(fapl_id))) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }

    /* Set values */
    fapl->gc_ref = (gc_ref!=0);

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_gc_references
 *
 * Purpose:	Returns the current setting for the garbage collection
 *		references property from a file access property list.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              June, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_gc_references(hid_t fapl_id, unsigned *gc_ref/*out*/)
{
    H5F_access_t	*fapl = NULL;

    FUNC_ENTER(H5Pget_gc_references, FAIL);
    H5TRACE2("e","ix",fapl_id,gc_ref);

    /* Check args */
    if (H5P_FILE_ACCESS!=H5P_get_class(fapl_id) ||
            NULL==(fapl=H5I_object(fapl_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a file access property list");
    }

    /* Get values */
    if (gc_ref)
        *gc_ref = fapl->gc_ref;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_meta_block_size
 *
 * Purpose:	Sets the minimum size of metadata block allocations when
 *      the H5FD_FEAT_AGGREGATE_METADATA is set by a VFL driver.
 *      Each "raw" metadata block is allocated to be this size and then
 *      specific pieces of metadata (object headers, local heaps, B-trees, etc)
 *      are sub-allocated from this block.
 *      
 *		The default value is set to 2048 (bytes), indicating that metadata
 *      will be attempted to be bunched together in (at least) 2K blocks in
 *      the file.  Setting the value to 0 with this API function will
 *      turn off the metadata aggregation, even if the VFL driver attempts to
 *      use that strategy.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Friday, August 25, 2000
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_meta_block_size(hid_t fapl_id, hsize_t size)
{
    H5F_access_t	*fapl = NULL;
    
    FUNC_ENTER (H5Pset_meta_block_size, FAIL);
    H5TRACE2("e","ih",fapl_id,size);

    /* Check args */
    if (H5P_FILE_ACCESS != H5P_get_class (fapl_id) ||
            NULL == (fapl = H5I_object (fapl_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }

    /* Set values */
    H5_ASSIGN_OVERFLOW(fapl->meta_block_size,size,hsize_t,size_t);

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_meta_block_size
 *
 * Purpose:	Returns the current settings for the metadata block allocation
 *      property from a file access property list.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Friday, August 29, 2000
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_meta_block_size(hid_t fapl_id, hsize_t *size/*out*/)
{
    H5F_access_t	*fapl = NULL;

    FUNC_ENTER (H5Pget_meta_block_size, FAIL);
    H5TRACE2("e","ix",fapl_id,size);

    /* Check args */
    if (H5P_FILE_ACCESS != H5P_get_class (fapl_id) ||
            NULL == (fapl = H5I_object (fapl_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }

    /* Get values */
    if (size)
        *size = fapl->meta_block_size;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_sieve_buf_size
 *
 * Purpose:	Sets the maximum size of the data seive buffer used for file
 *      drivers which are capable of using data sieving.  The data sieve
 *      buffer is used when performing I/O on datasets in the file.  Using a
 *      buffer which is large anough to hold several pieces of the dataset
 *      being read in for hyperslab selections boosts performance by quite a
 *      bit.
 *      
 *		The default value is set to 64KB, indicating that file I/O for raw data
 *      reads and writes will occur in at least 64KB blocks.
 *      Setting the value to 0 with this API function will turn off the
 *      data sieving, even if the VFL driver attempts to use that strategy.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, September 21, 2000
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_sieve_buf_size(hid_t fapl_id, hsize_t size)
{
    H5F_access_t	*fapl = NULL;
    
    FUNC_ENTER (H5Pset_sieve_buf_size, FAIL);
    H5TRACE2("e","ih",fapl_id,size);

    /* Check args */
    if (H5P_FILE_ACCESS != H5P_get_class (fapl_id) ||
            NULL == (fapl = H5I_object (fapl_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }

    /* Set values */
    fapl->sieve_buf_size = size;

    FUNC_LEAVE (SUCCEED);
} /* end H5Pset_sieve_buf_size() */


/*-------------------------------------------------------------------------
 * Function:	H5Pget_sieve_buf_size
 *
 * Purpose:	Returns the current settings for the data sieve buffer size
 *      property from a file access property list.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, September 21, 2000
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_sieve_buf_size(hid_t fapl_id, hsize_t *size/*out*/)
{
    H5F_access_t	*fapl = NULL;

    FUNC_ENTER (H5Pget_sieve_buf_size, FAIL);
    H5TRACE2("e","ix",fapl_id,size);

    /* Check args */
    if (H5P_FILE_ACCESS != H5P_get_class (fapl_id) ||
            NULL == (fapl = H5I_object (fapl_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }

    /* Get values */
    if (size)
        *size = fapl->sieve_buf_size;

    FUNC_LEAVE (SUCCEED);
} /* end H5Pget_sieve_buf_size() */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_small_data_block_size
 *
 * Purpose:	Sets the minimum size of "small" raw data block allocations
 *      when the H5FD_FEAT_AGGREGATE_SMALLDATA is set by a VFL driver.
 *      Each "small" raw data block is allocated to be this size and then
 *      specific pieces of raw data are sub-allocated from this block.
 *      
 *		The default value is set to 2048 (bytes), indicating that
 *      "small" raw data will be attempted to be bunched together in (at least)
 *      2K blocks in the file.  Setting the value to 0 with this API function
 *      will turn off the "small" raw data aggregation, even if the VFL driver
 *      attempts to use that strategy.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Monday, June 10, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_small_data_block_size(hid_t fapl_id, hsize_t size)
{
    H5F_access_t	*fapl = NULL;
    
    FUNC_ENTER (H5Pset_small_data_block_size, FAIL);
    H5TRACE2("e","ih",fapl_id,size);

    /* Check args */
    if (H5P_FILE_ACCESS != H5P_get_class (fapl_id) ||
            NULL == (fapl = H5I_object (fapl_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }

    /* Set values */
    H5_ASSIGN_OVERFLOW(fapl->sdata_block_size,size,hsize_t,size_t);

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_small_data_block_size
 *
 * Purpose:	Returns the current settings for the "small" raw data block
 *      allocation property from a file access property list.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Monday, June 10, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_small_data_block_size(hid_t fapl_id, hsize_t *size/*out*/)
{
    H5F_access_t	*fapl = NULL;

    FUNC_ENTER (H5Pget_small_data_block_size, FAIL);
    H5TRACE2("e","ix",fapl_id,size);

    /* Check args */
    if (H5P_FILE_ACCESS != H5P_get_class (fapl_id) ||
            NULL == (fapl = H5I_object (fapl_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a file access property list");
    }

    /* Get values */
    if (size)
        *size = fapl->sdata_block_size;

    FUNC_LEAVE (SUCCEED);
}

