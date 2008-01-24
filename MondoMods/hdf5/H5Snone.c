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
 * Programmer:  Quincey Koziol <koziol@ncsa.uiuc.edu>
 *              Tuesday, November 10, 1998
 *
 * Purpose:	"None" selection data space I/O functions.
 */

#define H5S_PACKAGE		/*suppress error about including H5Spkg	  */

#include "H5private.h"
#include "H5Eprivate.h"
#include "H5Iprivate.h"
#include "H5Spkg.h"
#include "H5Vprivate.h"
#include "H5Dprivate.h"

/* Interface initialization */
#define PABLO_MASK      H5Snone_mask
#define INTERFACE_INIT  NULL
static int             interface_initialize_g = 0;

static herr_t H5S_none_init (const struct H5O_layout_t *layout,
			    const H5S_t *space, H5S_sel_iter_t *iter);
static hsize_t H5S_none_favail (const H5S_t *space, const H5S_sel_iter_t *iter,
			      hsize_t max);
static hsize_t H5S_none_fgath (H5F_t *f, const struct H5O_layout_t *layout,
			     const struct H5O_pline_t *pline,
			     const struct H5O_fill_t *fill,
			     const struct H5O_efl_t *efl, size_t elmt_size,
			     const H5S_t *file_space,
			     H5S_sel_iter_t *file_iter, hsize_t nelmts,
			     hid_t dxpl_id, void *buf/*out*/);
static herr_t H5S_none_fscat (H5F_t *f, const struct H5O_layout_t *layout,
			     const struct H5O_pline_t *pline,
			     const struct H5O_fill_t *fill,
			     const struct H5O_efl_t *efl, size_t elmt_size,
			     const H5S_t *file_space,
			     H5S_sel_iter_t *file_iter, hsize_t nelmts,
			     hid_t dxpl_id, const void *buf);
static hsize_t H5S_none_mgath (const void *_buf, size_t elmt_size,
			     const H5S_t *mem_space, H5S_sel_iter_t *mem_iter,
			     hsize_t nelmts, void *_tconv_buf/*out*/);
static herr_t H5S_none_mscat (const void *_tconv_buf, size_t elmt_size,
			     const H5S_t *mem_space, H5S_sel_iter_t *mem_iter,
			     hsize_t nelmts, void *_buf/*out*/);
static herr_t H5S_select_none(H5S_t *space);

const H5S_fconv_t	H5S_NONE_FCONV[1] = {{
    "none", 					/*name			*/
    H5S_SEL_NONE,				/*selection type	*/
    H5S_none_init,				/*initialize		*/
    H5S_none_favail,				/*available		*/
    H5S_none_fgath,				/*gather		*/
    H5S_none_fscat,				/*scatter		*/
}};

const H5S_mconv_t	H5S_NONE_MCONV[1] = {{
    "none", 					/*name			*/
    H5S_SEL_NONE,				/*selection type	*/
    H5S_none_init,				/*initialize		*/
    H5S_none_mgath,				/*gather		*/
    H5S_none_mscat, 				/*scatter		*/
}};


/*-------------------------------------------------------------------------
 * Function:	H5S_none_init
 *
 * Purpose:	Initializes iteration information for none selection.
 *
 * Return:	non-negative on success, negative on failure.
 *
 * Programmer:	Quincey Koziol
 *              Tuesday, October 29, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5S_none_init(const struct H5O_layout_t UNUSED *layout,
	       const H5S_t UNUSED *space, H5S_sel_iter_t UNUSED *sel_iter)
{
    FUNC_ENTER(H5S_none_init, FAIL);

    /* Check args */
    assert(layout);
    assert(space && H5S_SEL_NONE==space->select.type);
    assert(sel_iter);

    FUNC_LEAVE (SUCCEED);
} /* H5S_none_init() */


/*-------------------------------------------------------------------------
 * Function:	H5S_none_favail
 *
 * Purpose:	Figure out the optimal number of elements to transfer to/from
 *		the file.
 *
 * Return:	zero always.
 *
 * Programmer:	Quincey Koziol
 *              Tuesday, October 29, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static hsize_t
H5S_none_favail(const H5S_t UNUSED *space, const H5S_sel_iter_t UNUSED *sel_iter, hsize_t UNUSED max)
{
    FUNC_ENTER(H5S_none_favail, 0);

    /* Check args */
    assert(space && H5S_SEL_NONE==space->select.type);
    assert(sel_iter);

    FUNC_LEAVE(0);
}   /* H5S_none_favail() */


/*-------------------------------------------------------------------------
 * Function:	H5S_none_fgath
 *
 * Purpose:	Gathers data points from file F and accumulates them in the
 *		type conversion buffer BUF.  The LAYOUT argument describes
 *		how the data is stored on disk and EFL describes how the data
 *		is organized in external files.  ELMT_SIZE is the size in
 *		bytes of a datum which this function treats as opaque.
 *		FILE_SPACE describes the data space of the dataset on disk
 *		and the elements that have been selected for reading (via
 *		hyperslab, etc).  This function will copy at most NELMTS
 *		elements.
 *
 * Return:	0 always
 *
 * Programmer:	Quincey Koziol
 *              Tuesday, October 29, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static hsize_t
H5S_none_fgath (H5F_t UNUSED *f, const struct H5O_layout_t UNUSED *layout,
	       const struct H5O_pline_t UNUSED *pline,
	       const struct H5O_fill_t UNUSED *fill, const struct H5O_efl_t UNUSED *efl,
	       size_t UNUSED elmt_size, const H5S_t UNUSED *file_space,
	       H5S_sel_iter_t UNUSED *file_iter, hsize_t UNUSED nelmts, hid_t UNUSED dxpl_id,
	       void UNUSED *buf/*out*/)
{
    FUNC_ENTER(H5S_none_fgath, 0);

    /* Check args */
    assert(f);
    assert(layout);
    assert(elmt_size>0);
    assert(file_space && H5S_SEL_NONE==file_space->select.type);
    assert(file_iter);
    assert(nelmts>0);
    assert(buf);

    FUNC_LEAVE(0);
} /* H5S_none_fgath() */


/*-------------------------------------------------------------------------
 * Function:	H5S_none_fscat
 *
 * Purpose:	Scatters dataset elements from the type conversion buffer BUF
 *		to the file F where the data points are arranged according to
 *		the file data space FILE_SPACE and stored according to
 *		LAYOUT and EFL. Each element is ELMT_SIZE bytes.
 *		The caller is requesting that NELMTS elements are copied.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Tuesday, October 29, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5S_none_fscat (H5F_t UNUSED *f, const struct H5O_layout_t UNUSED *layout,
	       const struct H5O_pline_t UNUSED *pline, const struct H5O_fill_t UNUSED *fill,
	       const struct H5O_efl_t UNUSED *efl, size_t UNUSED elmt_size,
	       const H5S_t UNUSED *file_space, H5S_sel_iter_t UNUSED *file_iter,
	       hsize_t UNUSED nelmts, hid_t UNUSED dxpl_id, const void UNUSED *buf)
{
    FUNC_ENTER(H5S_none_fscat, FAIL);

    /* Check args */
    assert(f);
    assert(layout);
    assert(elmt_size>0);
    assert(file_space && H5S_SEL_NONE==file_space->select.type);
    assert(file_iter);
    assert(nelmts>0);
    assert(buf);

    FUNC_LEAVE(SUCCEED);
}   /* H5S_none_fscat() */


/*-------------------------------------------------------------------------
 * Function:	H5S_none_mgath
 *
 * Purpose:	Gathers dataset elements from application memory BUF and
 *		copies them into the data type conversion buffer TCONV_BUF.
 *		Each element is ELMT_SIZE bytes and arranged in application
 *		memory according to MEM_SPACE.  
 *		The caller is requesting that at most NELMTS be gathered.
 *
 * Return:	0 always
 *
 * Programmer:	Quincey Koziol
 *              Tuesday, October 29, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static hsize_t
H5S_none_mgath (const void UNUSED *_buf, size_t UNUSED elmt_size,
	       const H5S_t UNUSED *mem_space, H5S_sel_iter_t UNUSED *mem_iter,
	       hsize_t UNUSED nelmts, void UNUSED *tconv_buf/*out*/)
{
    FUNC_ENTER(H5S_none_mgath, 0);

    /* Check args */
    assert(_buf);
    assert(elmt_size>0);
    assert(mem_space && H5S_SEL_NONE==mem_space->select.type);
    assert(mem_iter);
    assert(nelmts>0);
    assert(tconv_buf);

    FUNC_LEAVE(0);
}   /* H5S_none_mgath() */


/*-------------------------------------------------------------------------
 * Function:	H5S_none_mscat
 *
 * Purpose:	Scatters NELMTS data points from the type conversion buffer
 *		TCONV_BUF to the application buffer BUF.  Each element is
 *		ELMT_SIZE bytes and they are organized in application memory
 *		according to MEM_SPACE.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Tuesday, October 29, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5S_none_mscat (const void UNUSED *tconv_buf, size_t UNUSED elmt_size,
	       const H5S_t UNUSED *mem_space, H5S_sel_iter_t UNUSED *mem_iter,
	       hsize_t UNUSED nelmts, void UNUSED *_buf/*out*/)
{
    FUNC_ENTER (H5S_none_mscat, FAIL);

    /* Check args */
    assert(tconv_buf);
    assert(elmt_size>0);
    assert(mem_space && H5S_SEL_NONE==mem_space->select.type);
    assert(mem_iter);
    assert(nelmts>0);
    assert(_buf);

    FUNC_LEAVE(SUCCEED);
}   /* H5S_none_mscat() */


/*--------------------------------------------------------------------------
 NAME
    H5S_none_select_serialize
 PURPOSE
    Serialize the current selection into a user-provided buffer.
 USAGE
    herr_t H5S_none_select_serialize(space, buf)
        H5S_t *space;           IN: Dataspace pointer of selection to serialize
        uint8 *buf;             OUT: Buffer to put serialized selection into
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Serializes the current element selection into a buffer.  (Primarily for
    storing on disk).
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5S_none_select_serialize (const H5S_t *space, uint8_t *buf)
{
    herr_t ret_value=FAIL;  /* return value */

    FUNC_ENTER (H5S_none_select_serialize, FAIL);

    assert(space);

    /* Store the preamble information */
    UINT32ENCODE(buf, (uint32_t)space->select.type);  /* Store the type of selection */
    UINT32ENCODE(buf, (uint32_t)1);  /* Store the version number */
    UINT32ENCODE(buf, (uint32_t)0);  /* Store the un-used padding */
    UINT32ENCODE(buf, (uint32_t)0);  /* Store the additional information length */

    /* Set success */
    ret_value=SUCCEED;

    FUNC_LEAVE (ret_value);
}   /* H5S_none_select_serialize() */

/*--------------------------------------------------------------------------
 NAME
    H5S_none_select_deserialize
 PURPOSE
    Deserialize the current selection from a user-provided buffer.
 USAGE
    herr_t H5S_none_select_deserialize(space, buf)
        H5S_t *space;           IN/OUT: Dataspace pointer to place selection into
        uint8 *buf;             IN: Buffer to retrieve serialized selection from
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Deserializes the current selection into a buffer.  (Primarily for retrieving
    from disk).
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5S_none_select_deserialize (H5S_t *space, const uint8_t UNUSED *buf)
{
    herr_t ret_value=FAIL;  /* return value */

    FUNC_ENTER (H5S_none_select_deserialize, FAIL);

    assert(space);

    /* Change to "none" selection */
    if((ret_value=H5S_select_none(space))<0) {
        HGOTO_ERROR(H5E_DATASPACE, H5E_CANTDELETE, FAIL, "can't change selection");
    } /* end if */

done:
    FUNC_LEAVE (ret_value);
}   /* H5S_none_select_deserialize() */


/*--------------------------------------------------------------------------
 NAME
    H5S_select_none
 PURPOSE
    Specify that nothing is selected in the extent
 USAGE
    herr_t H5S_select_none(dsid)
        hid_t dsid;             IN: Dataspace ID of selection to modify
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    This function de-selects the entire extent for a dataspace.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5S_select_none (H5S_t *space)
{
    herr_t ret_value=SUCCEED;  /* return value */

    FUNC_ENTER (H5S_select_none, FAIL);

    /* Check args */
    assert(space);

    /* Remove current selection first */
    if(H5S_select_release(space)<0) {
        HGOTO_ERROR(H5E_DATASPACE, H5E_CANTDELETE, FAIL,
            "can't release hyperslab");
    } /* end if */

    /* Set selection type */
    space->select.type=H5S_SEL_NONE;

done:
    FUNC_LEAVE (ret_value);
}   /* H5S_select_none() */


/*--------------------------------------------------------------------------
 NAME
    H5Sselect_none
 PURPOSE
    Specify that nothing is selected in the extent
 USAGE
    herr_t H5Sselect_none(dsid)
        hid_t dsid;             IN: Dataspace ID of selection to modify
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    This function de-selects the entire extent for a dataspace.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Sselect_none (hid_t spaceid)
{
    H5S_t	*space = NULL;  /* Dataspace to modify selection of */
    herr_t ret_value=FAIL;  /* return value */

    FUNC_ENTER (H5Sselect_none, FAIL);

    /* Check args */
    if (H5I_DATASPACE != H5I_get_type(spaceid) ||
            NULL == (space=H5I_object(spaceid))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a data space");
    }

    /* Change to "none" selection */
    if((ret_value=H5S_select_none(space))<0) {
        HGOTO_ERROR(H5E_DATASPACE, H5E_CANTDELETE, FAIL, "can't change selection");
    } /* end if */

done:
    FUNC_LEAVE (ret_value);
}   /* H5Sselect_none() */


/*--------------------------------------------------------------------------
 NAME
    H5S_none_select_iterate
 PURPOSE
    Iterate over a none selection, calling a user's function for each
        element. (i.e. the user's function is not called because there are
        zero elements selected)
 USAGE
    herr_t H5S_none_select_iterate(buf, type_id, space, op, operator_data)
        void *buf;      IN/OUT: Buffer containing elements to iterate over
        hid_t type_id;  IN: Datatype ID of BUF array.
        H5S_t *space;   IN: Dataspace object containing selection to iterate over
        H5D_operator_t op; IN: Function pointer to the routine to be
                                called for each element in BUF iterated over.
        void *operator_data;    IN/OUT: Pointer to any user-defined data
                                associated with the operation.
 RETURNS
    Returns success (0).
 DESCRIPTION
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5S_none_select_iterate(void UNUSED *buf, hid_t UNUSED type_id, H5S_t UNUSED *space, H5D_operator_t UNUSED op,
        void UNUSED *operator_data)
{
    herr_t ret_value=SUCCEED;      /* return value */

    FUNC_ENTER (H5S_none_select_iterate, FAIL);

    assert(buf);
    assert(space);
    assert(op);
    assert(H5I_DATATYPE == H5I_get_type(type_id));

    FUNC_LEAVE (ret_value);
}   /* H5S_hyper_select_iterate() */
