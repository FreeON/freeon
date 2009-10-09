/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Copyright by The HDF Group.                                               *
 * Copyright by the Board of Trustees of the University of Illinois.         *
 * All rights reserved.                                                      *
 *                                                                           *
 * This file is part of HDF5.  The full HDF5 copyright notice, including     *
 * terms governing use, modification, and redistribution, is contained in    *
 * the files COPYING and Copyright.html.  COPYING can be found at the root   *
 * of the source code distribution tree; Copyright.html can be found at the  *
 * root level of an installed copy of the electronic HDF5 document set and   *
 * is linked from the top-level documents page.  It can also be found at     *
 * http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
 * access to either file, you may request a copy from help@hdfgroup.org.     *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/*
 * Programmer: Quincey Koziol <koziol@ncsa.uiuc.edu>
 *	       Thursday, September 30, 2004
 */

/****************/
/* Module Setup */
/****************/

#define H5D_PACKAGE		/*suppress error about including H5Dpkg	  */


/***********/
/* Headers */
/***********/
#include "H5private.h"		/* Generic Functions			*/
#include "H5Dpkg.h"		/* Datasets				*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5Fprivate.h"		/* Files				*/


/****************/
/* Local Macros */
/****************/


/******************/
/* Local Typedefs */
/******************/


/********************/
/* Local Prototypes */
/********************/

/* Layout operation callbacks */
static herr_t H5D_efl_new(H5F_t *f, hid_t dapl_id, hid_t dxpl_id, H5D_t *dset,
    const H5P_genplist_t *dc_plist);
static hbool_t H5D_efl_is_space_alloc(const H5O_layout_t *layout);
static herr_t H5D_efl_io_init(const H5D_io_info_t *io_info, const H5D_type_info_t *type_info,
    hsize_t nelmts, const H5S_t *file_space, const H5S_t *mem_space,
    H5D_chunk_map_t *cm);
static ssize_t H5D_efl_readvv(const H5D_io_info_t *io_info,
    size_t dset_max_nseq, size_t *dset_curr_seq, size_t dset_len_arr[], hsize_t dset_offset_arr[],
    size_t mem_max_nseq, size_t *mem_curr_seq, size_t mem_len_arr[], hsize_t mem_offset_arr[]);
static ssize_t H5D_efl_writevv(const H5D_io_info_t *io_info,
    size_t dset_max_nseq, size_t *dset_curr_seq, size_t dset_len_arr[], hsize_t dset_offset_arr[],
    size_t mem_max_nseq, size_t *mem_curr_seq, size_t mem_len_arr[], hsize_t mem_offset_arr[]);

/* Helper routines */
static herr_t H5D_efl_read(const H5O_efl_t *efl, haddr_t addr, size_t size,
    uint8_t *buf);
static herr_t H5D_efl_write(const H5O_efl_t *efl, haddr_t addr, size_t size,
    const uint8_t *buf);


/*********************/
/* Package Variables */
/*********************/

/* External File List (EFL) storage layout I/O ops */
const H5D_layout_ops_t H5D_LOPS_EFL[1] = {{
    H5D_efl_new,
    H5D_efl_is_space_alloc,
    H5D_efl_io_init,
    H5D_contig_read,
    H5D_contig_write,
#ifdef H5_HAVE_PARALLEL
    NULL,
    NULL,
#endif /* H5_HAVE_PARALLEL */
    H5D_efl_readvv,
    H5D_efl_writevv,
    NULL
}};


/*******************/
/* Local Variables */
/*******************/



/*-------------------------------------------------------------------------
 * Function:	H5D_efl_new
 *
 * Purpose:	Constructs new EFL layout information for dataset
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, May 22, 2008
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_efl_new(H5F_t *f, hid_t UNUSED dapl_id, hid_t UNUSED dxpl_id, H5D_t *dset,
    const H5P_genplist_t *dc_plist)
{
    size_t dt_size;                     /* Size of datatype */
    hsize_t dim[H5O_LAYOUT_NDIMS];	/* Current size of data in elements */
    hsize_t max_dim[H5O_LAYOUT_NDIMS];  /* Maximum size of data in elements */
    hssize_t tmp_size;                  /* Temporary holder for raw data size */
    hsize_t max_points;                 /* Maximum elements */
    hsize_t max_storage;                /* Maximum storage size */
    int ndims;                          /* Rank of dataspace */
    int i;                              /* Local index variable */
    herr_t ret_value = SUCCEED;         /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_efl_new)

    /* Sanity checks */
    HDassert(f);
    HDassert(dset);
    HDassert(dc_plist);

    /*
     * The maximum size of the dataset cannot exceed the storage size.
     * Also, only the slowest varying dimension of a simple data space
     * can be extendible (currently only for external data storage).
     */

    /* Check for invalid dataset dimensions */
    if((ndims = H5S_get_simple_extent_dims(dset->shared->space, dim, max_dim)) < 0)
        HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "unable to initialize contiguous storage")
    for(i = 1; i < ndims; i++)
        if(max_dim[i] > dim[i])
            HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "only the first dimension can be extendible")

    /* Retrieve the size of the dataset's datatype */
    if(0 == (dt_size = H5T_get_size(dset->shared->type)))
        HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "unable to determine datatype size")

    /* Check for storage overflows */
    max_points = H5S_get_npoints_max(dset->shared->space);
    max_storage = H5O_efl_total_size(&dset->shared->dcpl_cache.efl);
    if(H5S_UNLIMITED == max_points) {
        if(H5O_EFL_UNLIMITED != max_storage)
            HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "unlimited dataspace but finite storage")
    } /* end if */
    else if((max_points * dt_size) < max_points)
        HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "dataspace * type size overflowed")
    else if((max_points * dt_size) > max_storage)
        HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "dataspace size exceeds external storage size")

    /* Compute the total size of dataset */
    tmp_size = H5S_GET_EXTENT_NPOINTS(dset->shared->space) * dt_size;
    H5_ASSIGN_OVERFLOW(dset->shared->layout.u.contig.size, tmp_size, hssize_t, hsize_t);

    /* Get the sieve buffer size for this dataset */
    dset->shared->cache.contig.sieve_buf_size = H5F_SIEVE_BUF_SIZE(f);

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_efl_new() */


/*-------------------------------------------------------------------------
 * Function:	H5D_efl_is_space_alloc
 *
 * Purpose:	Query if space is allocated for layout
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, January 15, 2009
 *
 *-------------------------------------------------------------------------
 */
static hbool_t
H5D_efl_is_space_alloc(const H5O_layout_t UNUSED *layout)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_efl_is_space_alloc)

    /* Sanity checks */
    HDassert(layout);

    /* EFL storage is currently treated as allocated */
    FUNC_LEAVE_NOAPI(TRUE)
} /* end H5D_efl_is_space_alloc() */


/*-------------------------------------------------------------------------
 * Function:	H5D_efl_io_init
 *
 * Purpose:	Performs initialization before any sort of I/O on the raw data
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, March 20, 2008
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_efl_io_init(const H5D_io_info_t *io_info, const H5D_type_info_t UNUSED *type_info,
    hsize_t UNUSED nelmts, const H5S_t UNUSED *file_space, const H5S_t UNUSED *mem_space,
    H5D_chunk_map_t UNUSED *cm)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_efl_io_init)

    HDmemcpy(&io_info->store->efl, &(io_info->dset->shared->dcpl_cache.efl), sizeof(H5O_efl_t));

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* end H5D_efl_io_init() */


/*-------------------------------------------------------------------------
 * Function:	H5D_efl_read
 *
 * Purpose:	Reads data from an external file list.  It is an error to
 *		read past the logical end of file, but reading past the end
 *		of any particular member of the external file list results in
 *		zeros.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Wednesday, March  4, 1998
 *
 * Modifications:
 *		Robb Matzke, 1999-07-28
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_efl_read(const H5O_efl_t *efl, haddr_t addr, size_t size, uint8_t *buf)
{
    int		fd = -1;
    size_t	to_read;
#ifndef NDEBUG
    hsize_t     tempto_read;
#endif /* NDEBUG */
    hsize_t     skip = 0;
    haddr_t     cur;
    ssize_t	n;
    size_t      u;                      /* Local index variable */
    herr_t      ret_value = SUCCEED;    /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_efl_read)

    /* Check args */
    HDassert(efl && efl->nused>0);
    HDassert(H5F_addr_defined(addr));
    HDassert(size < SIZET_MAX);
    HDassert(buf || 0 == size);

    /* Find the first efl member from which to read */
    for (u=0, cur=0; u<efl->nused; u++) {
	if (H5O_EFL_UNLIMITED==efl->slot[u].size || addr < cur+efl->slot[u].size) {
	    skip = addr - cur;
	    break;
	}
  	cur += efl->slot[u].size;
    }

    /* Read the data */
    while (size) {
        HDassert(buf);
	if (u>=efl->nused)
	    HGOTO_ERROR (H5E_EFL, H5E_OVERFLOW, FAIL, "read past logical end of file")
	if (H5F_OVERFLOW_HSIZET2OFFT (efl->slot[u].offset+skip))
	    HGOTO_ERROR (H5E_EFL, H5E_OVERFLOW, FAIL, "external file address overflowed")
	if ((fd=HDopen (efl->slot[u].name, O_RDONLY, 0))<0)
	    HGOTO_ERROR (H5E_EFL, H5E_CANTOPENFILE, FAIL, "unable to open external raw data file")
	if (HDlseek (fd, (off_t)(efl->slot[u].offset+skip), SEEK_SET)<0)
	    HGOTO_ERROR (H5E_EFL, H5E_SEEKERROR, FAIL, "unable to seek in external raw data file")
#ifndef NDEBUG
	tempto_read = MIN(efl->slot[u].size-skip,(hsize_t)size);
        H5_CHECK_OVERFLOW(tempto_read,hsize_t,size_t);
	to_read = (size_t)tempto_read;
#else /* NDEBUG */
	to_read = MIN((size_t)(efl->slot[u].size-skip), size);
#endif /* NDEBUG */
	if((n = HDread(fd, buf, to_read)) < 0)
	    HGOTO_ERROR(H5E_EFL, H5E_READERROR, FAIL, "read error in external raw data file")
	else if((size_t)n < to_read)
	    HDmemset(buf + n, 0, to_read - (size_t)n);
	HDclose(fd);
	fd = -1;
	size -= to_read;
	buf += to_read;
	skip = 0;
	u++;
    }

done:
    if (fd>=0)
        HDclose (fd);

    FUNC_LEAVE_NOAPI(ret_value)
}


/*-------------------------------------------------------------------------
 * Function:	H5D_efl_write
 *
 * Purpose:	Writes data to an external file list.  It is an error to
 *		write past the logical end of file, but writing past the end
 *		of any particular member of the external file list just
 *		extends that file.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Wednesday, March  4, 1998
 *
 * Modifications:
 *		Robb Matzke, 1999-07-28
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_efl_write(const H5O_efl_t *efl, haddr_t addr, size_t size, const uint8_t *buf)
{
    int		fd = -1;
    size_t	to_write;
#ifndef NDEBUG
    hsize_t	tempto_write;
#endif /* NDEBUG */
    haddr_t     cur;
    hsize_t     skip = 0;
    size_t	u;                      /* Local index variable */
    herr_t      ret_value = SUCCEED;    /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_efl_write)

    /* Check args */
    HDassert(efl && efl->nused>0);
    HDassert(H5F_addr_defined(addr));
    HDassert(size < SIZET_MAX);
    HDassert(buf || 0 == size);

    /* Find the first efl member in which to write */
    for (u=0, cur=0; u<efl->nused; u++) {
	if (H5O_EFL_UNLIMITED==efl->slot[u].size || addr < cur+efl->slot[u].size) {
	    skip = addr - cur;
	    break;
	}
	cur += efl->slot[u].size;
    }

    /* Write the data */
    while (size) {
        assert(buf);
	if (u>=efl->nused)
	    HGOTO_ERROR (H5E_EFL, H5E_OVERFLOW, FAIL, "write past logical end of file")
	if (H5F_OVERFLOW_HSIZET2OFFT (efl->slot[u].offset+skip))
	    HGOTO_ERROR (H5E_EFL, H5E_OVERFLOW, FAIL, "external file address overflowed")
	if ((fd=HDopen (efl->slot[u].name, O_CREAT|O_RDWR, 0666))<0) {
	    if (HDaccess (efl->slot[u].name, F_OK)<0) {
		HGOTO_ERROR (H5E_EFL, H5E_CANTOPENFILE, FAIL, "external raw data file does not exist")
	    } else {
		HGOTO_ERROR (H5E_EFL, H5E_CANTOPENFILE, FAIL, "unable to open external raw data file")
	    }
	}
	if (HDlseek (fd, (off_t)(efl->slot[u].offset+skip), SEEK_SET)<0)
	    HGOTO_ERROR (H5E_EFL, H5E_SEEKERROR, FAIL, "unable to seek in external raw data file")
#ifndef NDEBUG
	tempto_write = MIN(efl->slot[u].size-skip,(hsize_t)size);
        H5_CHECK_OVERFLOW(tempto_write,hsize_t,size_t);
        to_write = (size_t)tempto_write;
#else /* NDEBUG */
	to_write = MIN((size_t)(efl->slot[u].size-skip), size);
#endif /* NDEBUG */
	if ((size_t)HDwrite (fd, buf, to_write)!=to_write)
	    HGOTO_ERROR (H5E_EFL, H5E_READERROR, FAIL, "write error in external raw data file")
	HDclose (fd);
	fd = -1;
	size -= to_write;
	buf += to_write;
	skip = 0;
	u++;
    }

done:
    if (fd>=0)
        HDclose (fd);

    FUNC_LEAVE_NOAPI(ret_value)
}


/*-------------------------------------------------------------------------
 * Function:	H5D_efl_readvv
 *
 * Purpose:	Reads data from an external file list.  It is an error to
 *		read past the logical end of file, but reading past the end
 *		of any particular member of the external file list results in
 *		zeros.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Wednesday, May  7, 2003
 *
 *-------------------------------------------------------------------------
 */
static ssize_t
H5D_efl_readvv(const H5D_io_info_t *io_info,
    size_t dset_max_nseq, size_t *dset_curr_seq, size_t dset_len_arr[], hsize_t dset_offset_arr[],
    size_t mem_max_nseq, size_t *mem_curr_seq, size_t mem_len_arr[], hsize_t mem_offset_arr[])
{
    const H5O_efl_t *efl = &(io_info->store->efl); /* Pointer to efl info */
    unsigned char *buf;         /* Pointer to buffer to write */
    haddr_t addr;               /* Actual address to read */
    size_t total_size = 0;      /* Total size of sequence in bytes */
    size_t size;                /* Size of sequence in bytes */
    size_t u;                   /* Counting variable */
    size_t v;                   /* Counting variable */
    ssize_t ret_value;          /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_efl_readvv)

    /* Check args */
    HDassert(efl && efl->nused > 0);
    HDassert(io_info->u.rbuf);

    /* Work through all the sequences */
    for(u = *dset_curr_seq, v = *mem_curr_seq; u < dset_max_nseq && v < mem_max_nseq; ) {
        /* Choose smallest buffer to write */
        if(mem_len_arr[v] < dset_len_arr[u])
            size = mem_len_arr[v];
        else
            size = dset_len_arr[u];

        /* Compute offset on disk */
        addr = dset_offset_arr[u];

        /* Compute offset in memory */
        buf = (unsigned char *)io_info->u.rbuf + mem_offset_arr[v];

        /* Read data */
        if(H5D_efl_read(efl, addr, size, buf) < 0)
            HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "block write failed")

        /* Update memory information */
        mem_len_arr[v] -= size;
        mem_offset_arr[v] += size;
        if(mem_len_arr[v] == 0)
            v++;

        /* Update file information */
        dset_len_arr[u] -= size;
        dset_offset_arr[u] += size;
        if(dset_len_arr[u] == 0)
            u++;

        /* Increment number of bytes copied */
        total_size += size;
    } /* end for */

    /* Update current sequence vectors */
    *dset_curr_seq = u;
    *mem_curr_seq = v;

    /* Set return value */
    H5_ASSIGN_OVERFLOW(ret_value, total_size, size_t, ssize_t);

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_efl_readvv() */


/*-------------------------------------------------------------------------
 * Function:	H5D_efl_writevv
 *
 * Purpose:	Writes data to an external file list.  It is an error to
 *		write past the logical end of file, but writing past the end
 *		of any particular member of the external file list just
 *		extends that file.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Friday, May  2, 2003
 *
 *-------------------------------------------------------------------------
 */
static ssize_t
H5D_efl_writevv(const H5D_io_info_t *io_info,
    size_t dset_max_nseq, size_t *dset_curr_seq, size_t dset_len_arr[], hsize_t dset_offset_arr[],
    size_t mem_max_nseq, size_t *mem_curr_seq, size_t mem_len_arr[], hsize_t mem_offset_arr[])
{
    const H5O_efl_t *efl = &(io_info->store->efl); /* Pointer to efl info */
    const unsigned char *buf;   /* Pointer to buffer to write */
    haddr_t addr;               /* Actual address to read */
    size_t total_size = 0;      /* Total size of sequence in bytes */
    size_t size;                /* Size of sequence in bytes */
    size_t u;                   /* Counting variable */
    size_t v;                   /* Counting variable */
    ssize_t ret_value;          /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_efl_writevv)

    /* Check args */
    HDassert(efl && efl->nused > 0);
    HDassert(io_info->u.wbuf);

    /* Work through all the sequences */
    for(u = *dset_curr_seq, v = *mem_curr_seq; u < dset_max_nseq && v < mem_max_nseq; ) {
        /* Choose smallest buffer to write */
        if(mem_len_arr[v] < dset_len_arr[u])
            size = mem_len_arr[v];
        else
            size = dset_len_arr[u];

        /* Compute offset on disk */
        addr = dset_offset_arr[u];

        /* Compute offset in memory */
        buf = (const unsigned char *)io_info->u.wbuf + mem_offset_arr[v];

        /* Write data */
        if(H5D_efl_write(efl, addr, size, buf) < 0)
            HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "block write failed")

        /* Update memory information */
        mem_len_arr[v] -= size;
        mem_offset_arr[v] += size;
        if(mem_len_arr[v] == 0)
            v++;

        /* Update file information */
        dset_len_arr[u] -= size;
        dset_offset_arr[u] += size;
        if(dset_len_arr[u] == 0)
            u++;

        /* Increment number of bytes copied */
        total_size += size;
    } /* end for */

    /* Update current sequence vectors */
    *dset_curr_seq = u;
    *mem_curr_seq = v;

    /* Set return value */
    H5_ASSIGN_OVERFLOW(ret_value, total_size, size_t, ssize_t);

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_efl_writevv() */

