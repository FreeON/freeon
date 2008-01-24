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

/* $Id: H5Pdcpl.c,v 1.1.2.1 2002/08/12 18:14:09 koziol Exp $ */

/* Private header files */
#include "H5private.h"		/* Generic Functions			*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5Iprivate.h"		/* IDs			  		*/
#include "H5MMprivate.h"	/* Memory management			*/
#include "H5Pprivate.h"		/* Property lists		  	*/

#define PABLO_MASK	H5Pdcpl_mask

/* Is the interface initialized? */
static int		interface_initialize_g = 0;
#define INTERFACE_INIT NULL

/* Local types */

/* Local static functions */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_layout
 *
 * Purpose:	Sets the layout of raw data in the file.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Tuesday, January  6, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_layout(hid_t plist_id, H5D_layout_t layout)
{
    H5D_create_t	   *plist = NULL;

    FUNC_ENTER(H5Pset_layout, FAIL);
    H5TRACE2("e","iDl",plist_id,layout);

    /* Check arguments */
    if (H5P_DATASET_CREATE != H5P_get_class(plist_id) ||
            NULL == (plist = H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset creation property list");
    }
    if (layout < 0 || layout >= H5D_NLAYOUTS) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADRANGE, FAIL,
		      "raw data layout method is not valid");
    }

    /* Set value */
    plist->layout = layout;

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:	H5Pget_layout
 *
 * Purpose:	Retrieves layout type of a dataset creation property list.
 *
 * Return:	Success:	The layout type
 *
 *		Failure:	H5D_LAYOUT_ERROR (negative)
 *
 * Programmer:	Robb Matzke
 *		Wednesday, January  7, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
H5D_layout_t
H5Pget_layout(hid_t plist_id)
{
    H5D_create_t	   *plist = NULL;

    FUNC_ENTER(H5Pget_layout, H5D_LAYOUT_ERROR);
    H5TRACE1("Dl","i",plist_id);

    /* Check arguments */
    if (H5P_DATASET_CREATE != H5P_get_class(plist_id) ||
            NULL == (plist = H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, H5D_LAYOUT_ERROR,
		      "not a dataset creation property list");
    }

    FUNC_LEAVE(plist->layout);
}

/*-------------------------------------------------------------------------
 * Function:	H5Pset_chunk
 *
 * Purpose:	Sets the number of dimensions and the size of each chunk to
 *		the values specified.  The dimensionality of the chunk should
 *		match the dimensionality of the data space.
 *
 *		As a side effect, the layout method is changed to
 *		H5D_CHUNKED.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Tuesday, January  6, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_chunk(hid_t plist_id, int ndims, const hsize_t dim[/*ndims*/])
{
    int			    i;
    H5D_create_t	   *plist = NULL;

    FUNC_ENTER(H5Pset_chunk, FAIL);
    H5TRACE3("e","iIs*[a1]h",plist_id,ndims,dim);

    /* Check arguments */
    if (H5P_DATASET_CREATE != H5P_get_class(plist_id) ||
            NULL == (plist = H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset creation property list");
    }
    if (ndims <= 0) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADRANGE, FAIL,
		      "chunk dimensionality must be positive");
    }
    if (ndims > H5S_MAX_RANK) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADRANGE, FAIL,
		      "chunk dimensionality is too large");
    }
    if (!dim) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL,
		      "no chunk dimensions specified");
    }
    for (i=0; i<ndims; i++) {
        if (dim[i] <= 0) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADRANGE, FAIL,
                  "all chunk dimensions must be positive");
        }
    }

    /* Set value */
    plist->layout = H5D_CHUNKED;
    plist->chunk_ndims = ndims;
    for (i = 0; i < ndims; i++)
        plist->chunk_size[i] = dim[i];

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:	H5Pget_chunk
 *
 * Purpose:	Retrieves the chunk size of chunked layout.  The chunk
 *		dimensionality is returned and the chunk size in each
 *		dimension is returned through the DIM argument.	 At most
 *		MAX_NDIMS elements of DIM will be initialized.
 *
 * Return:	Success:	Positive Chunk dimensionality.
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *		Wednesday, January  7, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
H5Pget_chunk(hid_t plist_id, int max_ndims, hsize_t dim[]/*out*/)
{
    int			i;
    H5D_create_t	*plist = NULL;

    FUNC_ENTER(H5Pget_chunk, FAIL);
    H5TRACE3("Is","iIsx",plist_id,max_ndims,dim);

    /* Check arguments */
    if (H5P_DATASET_CREATE != H5P_get_class(plist_id) ||
            NULL == (plist = H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset creation property list");
    }
    if (H5D_CHUNKED != plist->layout) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL,
		      "not a chunked storage layout");
    }

    for (i=0; i<plist->chunk_ndims && i<max_ndims && dim; i++)
        dim[i] = plist->chunk_size[i];

    FUNC_LEAVE(plist->chunk_ndims);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_external
 *
 * Purpose:	Adds an external file to the list of external files. PLIST_ID
 *		should be an object ID for a dataset creation property list.
 *		NAME is the name of an external file, OFFSET is the location
 *		where the data starts in that file, and SIZE is the number of
 *		bytes reserved in the file for the data.
 *
 *		If a dataset is split across multiple files then the files
 *		should be defined in order. The total size of the dataset is
 *		the sum of the SIZE arguments for all the external files.  If
 *		the total size is larger than the size of a dataset then the
 *		dataset can be extended (provided the data space also allows
 *		the extending).
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Tuesday, March	3, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_external(hid_t plist_id, const char *name, off_t offset, hsize_t size)
{
    int			idx;
    hsize_t		total, tmp;
    H5D_create_t	*plist = NULL;

    FUNC_ENTER(H5Pset_external, FAIL);
    H5TRACE4("e","isoh",plist_id,name,offset,size);

    /* Check arguments */
    if (H5P_DATASET_CREATE != H5P_get_class(plist_id) ||
            NULL == (plist = H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset creation property list");
    }
    if (!name || !*name) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL, "no name given");
    }
    if (offset<0) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "negative external file offset");
    }
    if (size<=0) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL, "zero size");
    }
    if (plist->efl.nused>0 &&
            H5O_EFL_UNLIMITED==plist->efl.slot[plist->efl.nused-1].size) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "previous file size is unlimited");
    }
    if (H5O_EFL_UNLIMITED!=size) {
        for (idx=0, total=size; idx<plist->efl.nused; idx++, total=tmp) {
            tmp = total + plist->efl.slot[idx].size;
            if (tmp <= total) {
                HRETURN_ERROR (H5E_EFL, H5E_OVERFLOW, FAIL,
                       "total external data size overflowed");
            }
        }
    }
    
    /* Add to the list */
    if (plist->efl.nused>=plist->efl.nalloc) {
        int na = plist->efl.nalloc + H5O_EFL_ALLOC;
        H5O_efl_entry_t *x = H5MM_realloc (plist->efl.slot,
                           na*sizeof(H5O_efl_entry_t));

        if (!x) {
            HRETURN_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL,
                   "memory allocation failed");
        }
        plist->efl.nalloc = na;
        plist->efl.slot = x;
    }
    idx = plist->efl.nused;
    plist->efl.slot[idx].name_offset = 0; /*not entered into heap yet*/
    plist->efl.slot[idx].name = H5MM_xstrdup (name);
    plist->efl.slot[idx].offset = offset;
    plist->efl.slot[idx].size = size;
    plist->efl.nused++;

    FUNC_LEAVE(SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_external_count
 *
 * Purpose:	Returns the number of external files for this dataset.
 *
 * Return:	Success:	Number of external files
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *              Tuesday, March  3, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
H5Pget_external_count(hid_t plist_id)
{
    H5D_create_t	*plist = NULL;
    
    FUNC_ENTER (H5Pget_external_count, FAIL);
    H5TRACE1("Is","i",plist_id);
    
    /* Check arguments */
    if (H5P_DATASET_CREATE != H5P_get_class(plist_id) ||
            NULL == (plist = H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset creation property list");
    }

    /* Return */
    FUNC_LEAVE (plist->efl.nused);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_external
 *
 * Purpose:	Returns information about an external file.  External files
 *		are numbered from zero to N-1 where N is the value returned
 *		by H5Pget_external_count().  At most NAME_SIZE characters are
 *		copied into the NAME array.  If the external file name is
 *		longer than NAME_SIZE with the null terminator, then the
 *		return value is not null terminated (similar to strncpy()).
 *
 *		If NAME_SIZE is zero or NAME is the null pointer then the
 *		external file name is not returned.  If OFFSET or SIZE are
 *		null pointers then the corresponding information is not
 *		returned.
 *
 * See Also:	H5Pset_external()
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Tuesday, March  3, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_external(hid_t plist_id, int idx, size_t name_size, char *name/*out*/,
		 off_t *offset/*out*/, hsize_t *size/*out*/)
{
    H5D_create_t	*plist = NULL;
    
    FUNC_ENTER (H5Pget_external, FAIL);
    H5TRACE6("e","iIszxxx",plist_id,idx,name_size,name,offset,size);
    
    /* Check arguments */
    if (H5P_DATASET_CREATE != H5P_get_class(plist_id) ||
            NULL == (plist = H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset creation property list");
    }
    if (idx<0 || idx>=plist->efl.nused) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADRANGE, FAIL,
		       "external file index is out of range");
    }

    /* Return values */
    if (name_size>0 && name)
        HDstrncpy (name, plist->efl.slot[idx].name, name_size);
    if (offset)
        *offset = plist->efl.slot[idx].offset;
    if (size)
        *size = plist->efl.slot[idx].size;

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_filter
 *
 * Purpose:	Adds the specified FILTER and corresponding properties to the
 *		end of the transient or permanent output filter pipeline
 *		depending on whether PLIST is a dataset creation or dataset
 *		transfer property list.  The FLAGS argument specifies certain
 *		general properties of the filter and is documented below.
 *		The CD_VALUES is an array of CD_NELMTS integers which are
 *		auxiliary data for the filter.  The integer vlues will be
 *		stored in the dataset object header as part of the filter
 *		information.
 *
 * 		The FLAGS argument is a bit vector of the following fields:
 *
 * 		H5Z_FLAG_OPTIONAL(0x0001)
 *		If this bit is set then the filter is optional.  If the
 *		filter fails during an H5Dwrite() operation then the filter
 *		is just excluded from the pipeline for the chunk for which it
 *		failed; the filter will not participate in the pipeline
 *		during an H5Dread() of the chunk.  If this bit is clear and
 *		the filter fails then the entire I/O operation fails.
 *
 * Note:	This function currently supports only the permanent filter
 *		pipeline.  That is, PLIST_ID must be a dataset creation
 *		property list.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Wednesday, April 15, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_filter(hid_t plist_id, H5Z_filter_t filter, unsigned int flags,
	       size_t cd_nelmts, const unsigned int cd_values[/*cd_nelmts*/])
{
    H5D_create_t	*plist = NULL;
    
    FUNC_ENTER (H5Pset_filter, FAIL);
    H5TRACE5("e","iZfIuz*[a3]Iu",plist_id,filter,flags,cd_nelmts,cd_values);

    /* Check arguments */
    if (H5P_DATASET_XFER==H5P_get_class(plist_id)) {
        HRETURN_ERROR(H5E_PLINE, H5E_UNSUPPORTED, FAIL,
		      "transient pipelines are not supported yet");
    }
    if (H5P_DATASET_CREATE!=H5P_get_class (plist_id) ||
            NULL==(plist=H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a dataset creation property list");
    }
    if (filter<0 || filter>H5Z_FILTER_MAX) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "invalid filter identifier");
    }
    if (flags & ~((unsigned)H5Z_FLAG_DEFMASK)) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL,
		      "invalid flags");
    }
    if (cd_nelmts>0 && !cd_values) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL,
		      "no client data values supplied");
    }

    /* Do it */
    if (H5Z_append(&(plist->pline), filter, flags, cd_nelmts, cd_values)<0) {
        HRETURN_ERROR(H5E_PLINE, H5E_CANTINIT, FAIL,
		      "unable to add filter to pipeline");
    }

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_nfilters
 *
 * Purpose:	Returns the number of filters in the permanent or transient
 *		pipeline depending on whether PLIST_ID is a dataset creation
 *		or dataset transfer property list.  In each pipeline the
 *		filters are numbered from zero through N-1 where N is the
 *		value returned by this function.  During output to the file
 *		the filters of a pipeline are applied in increasing order
 *		(the inverse is true for input).
 *
 * Note:	Only permanent filters are supported at this time.
 *
 * Return:	Success:	Number of filters or zero if there are none.
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *              Tuesday, August  4, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
H5Pget_nfilters(hid_t plist_id)
{
    H5D_create_t	*plist = NULL;
    
    FUNC_ENTER(H5Pget_nfilters, FAIL);
    H5TRACE1("Is","i",plist_id);

    if (H5P_DATASET_XFER==H5P_get_class(plist_id)) {
        HRETURN_ERROR(H5E_PLINE, H5E_UNSUPPORTED, FAIL,
		      "transient pipelines are not supported yet");
    }
    if (H5P_DATASET_CREATE!=H5P_get_class(plist_id) ||
            NULL==(plist=H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset creation property list");
    }

    FUNC_LEAVE((int)(plist->pline.nfilters));
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_filter
 *
 * Purpose:	This is the query counterpart of H5Pset_filter() and returns
 *		information about a particular filter number in a permanent
 *		or transient pipeline depending on whether PLIST_ID is a
 *		dataset creation or transfer property list.  On input,
 *		CD_NELMTS indicates the number of entries in the CD_VALUES
 *		array allocated by the caller while on exit it contains the
 *		number of values defined by the filter.  The IDX should be a
 *		value between zero and N-1 as described for H5Pget_nfilters()
 *		and the function will return failure if the filter number is
 *		out or range.
 * 
 * Return:	Success:	Filter identification number.
 *
 *		Failure:	H5Z_FILTER_ERROR (Negative)
 *
 * Programmer:	Robb Matzke
 *              Wednesday, April 15, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
H5Z_filter_t
H5Pget_filter(hid_t plist_id, int idx, unsigned int *flags/*out*/,
	       size_t *cd_nelmts/*in_out*/, unsigned cd_values[]/*out*/,
	       size_t namelen, char name[]/*out*/)
{
    H5D_create_t	*plist = NULL;
    size_t		i;
    
    FUNC_ENTER (H5Pget_filter, H5Z_FILTER_ERROR);
    H5TRACE7("Zf","iIsx*zxzx",plist_id,idx,flags,cd_nelmts,cd_values,namelen,
             name);
    
    /* Check arguments */
    if (H5P_DATASET_XFER==H5P_get_class(plist_id)) {
        HRETURN_ERROR(H5E_PLINE, H5E_UNSUPPORTED, H5Z_FILTER_ERROR,
		      "transient filters are not supported yet");
    }
    if (H5P_DATASET_CREATE!=H5P_get_class (plist_id) ||
            NULL==(plist=H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, H5Z_FILTER_ERROR,
		       "not a dataset creation property list");
    }
    if (idx<0 || (size_t)idx>=plist->pline.nfilters) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, H5Z_FILTER_ERROR,
		      "filter number is invalid");
    }
    if (cd_nelmts || cd_values) {
        if (cd_nelmts && *cd_nelmts>256) {
            /*
             * It's likely that users forget to initialize this on input, so
             * we'll check that it has a reasonable value.  The actual number
             * is unimportant because the H5O layer will detect when a message
             * is too large.
             */
            HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, H5Z_FILTER_ERROR,
                  "probable uninitialized *cd_nelmts argument");
        }
        if (cd_nelmts && *cd_nelmts>0 && !cd_values) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, H5Z_FILTER_ERROR,
                  "client data values not supplied");
        }

        /*
         * If cd_nelmts is null but cd_values is non-null then just ignore
         * cd_values
         */
        if (!cd_nelmts)
            cd_values = NULL;
    }

    /* Output values */
    if (flags)
        *flags = plist->pline.filter[idx].flags;
    if (cd_values) {
        for (i=0; i<plist->pline.filter[idx].cd_nelmts && i<*cd_nelmts; i++) {
            cd_values[i] = plist->pline.filter[idx].cd_values[i];
        }
    }
    if (cd_nelmts)
        *cd_nelmts = plist->pline.filter[idx].cd_nelmts;

    if (namelen>0 && name) {
        const char *s = plist->pline.filter[idx].name;
        if (!s) {
            H5Z_class_t *cls = H5Z_find(plist->pline.filter[idx].id);

            if (cls)
                s = cls->name;
        }
        if (s)
            HDstrncpy(name, s, namelen);
        else
            name[0] = '\0';
    }
    
    FUNC_LEAVE (plist->pline.filter[idx].id);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_deflate
 *
 * Purpose:	Sets the compression method for a permanent or transient
 *		filter pipeline (depending on whether PLIST_ID is a dataset
 *		creation or transfer property list) to H5Z_FILTER_DEFLATE
 *		and the compression level to LEVEL which should be a value
 *		between zero and nine, inclusive.  Lower compression levels
 *		are faster but result in less compression.  This is the same
 *		algorithm as used by the GNU gzip program.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Wednesday, April 15, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_deflate(hid_t plist_id, unsigned level)
{
    H5D_create_t	*plist = NULL;
    
    FUNC_ENTER (H5Pset_deflate, FAIL);
    H5TRACE2("e","iIu",plist_id,level);

    /* Check arguments */
    if (H5P_DATASET_XFER==H5P_get_class(plist_id)) {
        HRETURN_ERROR(H5E_PLINE, H5E_UNSUPPORTED, FAIL,
		      "transient filter pipelines are not supported yet");
    }
    if (H5P_DATASET_CREATE!=H5P_get_class (plist_id) ||
            NULL==(plist=H5I_object (plist_id))) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
		       "not a dataset creation property list");
    }
    if (level>9) {
        HRETURN_ERROR (H5E_ARGS, H5E_BADVALUE, FAIL,
		       "invalid deflate level");
    }

    /* Add the filter */
    if (H5Z_append(&(plist->pline), H5Z_FILTER_DEFLATE, H5Z_FLAG_OPTIONAL,
		   1, &level)<0) {
        HRETURN_ERROR(H5E_PLINE, H5E_CANTINIT, FAIL,
		      "unable to add deflate filter to pipeline");
    }

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pset_fill_value
 *
 * Purpose:	Set the fill value for a dataset creation property list. The
 *		VALUE is interpretted as being of type TYPE, which need not
 *		be the same type as the dataset but the library must be able
 *		to convert VALUE to the dataset type when the dataset is
 *		created.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, October  1, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_fill_value(hid_t plist_id, hid_t type_id, const void *value)
{
    H5D_create_t	*plist = NULL;
    H5T_t		*type = NULL;
    
    FUNC_ENTER(H5Pset_fill_value, FAIL);
    H5TRACE3("e","iix",plist_id,type_id,value);

    /* Check arguments */
    if (H5P_DATASET_CREATE!=H5P_get_class(plist_id) ||
            NULL==(plist=H5I_object(plist_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		      "not a dataset creation property list");
    }
    if (H5I_DATATYPE!=H5I_get_type(type_id) ||
            NULL==(type=H5I_object(type_id))) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a data type");
    }
    if (!value) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "no fill value specified");
    }

    /* Set the fill value */
    H5O_reset(H5O_FILL, &(plist->fill));
    if (NULL==(plist->fill.type=H5T_copy(type, H5T_COPY_TRANSIENT))) {
        HRETURN_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL,
		      "unable to copy data type");
    }
    plist->fill.size = H5T_get_size(type);
    if (NULL==(plist->fill.buf=H5MM_malloc(plist->fill.size))) {
        HRETURN_ERROR(H5E_RESOURCE, H5E_CANTINIT, FAIL,
		      "memory allocation failed for fill value");
    }
    HDmemcpy(plist->fill.buf, value, plist->fill.size);

    FUNC_LEAVE(SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_fill_value
 *
 * Purpose:	Queries the fill value property of a dataset creation
 *		property list.  The fill value is returned through the VALUE
 *		pointer and the memory is allocated by the caller.  The fill
 *		value will be converted from its current data type to the
 *		specified TYPE.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, October  1, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_fill_value(hid_t plist_id, hid_t type_id, void *value/*out*/)
{
    H5D_create_t	*plist = NULL;		/*property list		*/
    H5T_t		*type = NULL;		/*data type		*/
    H5T_path_t		*tpath = NULL;		/*type conversion info	*/
    void		*buf = NULL;		/*conversion buffer	*/
    void		*bkg = NULL;		/*conversion buffer	*/
    hid_t		src_id = -1;		/*source data type id	*/
    herr_t		ret_value = FAIL;	/*return value		*/
    
    FUNC_ENTER(H5Pget_fill_value, FAIL);
    H5TRACE3("e","iix",plist_id,type_id,value);

    /* Check arguments */
    if (H5P_DATASET_CREATE!=H5P_get_class(plist_id) ||
            NULL==(plist=H5I_object(plist_id))) {
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL,
		    "not a dataset creation proprety list");
    }
    if (H5I_DATATYPE!=H5I_get_type(type_id) ||
            NULL==(type=H5I_object(type_id))) {
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a data type");
    }
    if (!value) {
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL,
		    "no fill value output buffer");
    }

    /*
     * If no fill value is defined then return an error.  We can't even
     * return zero because we don't know the data type of the dataset and
     * data type conversion might not have resulted in zero.
     */
    if (NULL==plist->fill.buf) {
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "no fill value defined");
    }

    /*
     * Can we convert between the source and destination data types?
     */
    if (NULL==(tpath=H5T_path_find(plist->fill.type, type, NULL, NULL))) {
        HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL,
		    "unable to convert between src and dst data types");
    }
    src_id = H5I_register(H5I_DATATYPE,
			  H5T_copy (plist->fill.type, H5T_COPY_TRANSIENT));
    if (src_id<0) {
        HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL,
		    "unable to copy/register data type");
    }

    /*
     * Data type conversions are always done in place, so we need a buffer
     * other than the fill value buffer that is large enough for both source
     * and destination.  The app-supplied buffer might do okay.
     */
    if (H5T_get_size(type)>=H5T_get_size(plist->fill.type)) {
        buf = value;
        if (tpath->cdata.need_bkg>=H5T_BKG_TEMP &&
                NULL==(bkg=H5MM_malloc(H5T_get_size(type)))) {
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL,
                "memory allocation failed for type conversion");
        }
    } else {
        if (NULL==(buf=H5MM_malloc(H5T_get_size(plist->fill.type)))) {
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL,
                "memory allocation failed for type conversion");
        }
        if (tpath->cdata.need_bkg>=H5T_BKG_TEMP)
            bkg = value;
    }
    HDmemcpy(buf, plist->fill.buf, H5T_get_size(plist->fill.type));
        
    /* Do the conversion */
    if (H5T_convert(tpath, src_id, type_id, (hsize_t)1, 0, 0, buf, bkg, H5P_DEFAULT)<0) {
        HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL,
            "data type conversion failed");
    }
    if (buf!=value)
        HDmemcpy(value, buf, H5T_get_size(type));
    ret_value = SUCCEED;

done:
    if (buf!=value)
        H5MM_xfree(buf);
    if (bkg!=value)
        H5MM_xfree(bkg);
    if (src_id>=0)
        H5I_dec_ref(src_id);
    FUNC_LEAVE(ret_value);
}

