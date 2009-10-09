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

/*-------------------------------------------------------------------------
 *
 * Created:		H5HGcache.c
 *			Feb  5 2008
 *			Quincey Koziol <koziol@hdfgroup.org>
 *
 * Purpose:		Implement global heap metadata cache methods.
 *
 *-------------------------------------------------------------------------
 */

/****************/
/* Module Setup */
/****************/

#define H5F_PACKAGE		/*suppress error about including H5Fpkg	  */
#define H5HG_PACKAGE		/*suppress error about including H5HGpkg  */


/***********/
/* Headers */
/***********/
#include "H5private.h"		/* Generic Functions			*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5Fpkg.h"             /* File access				*/
#include "H5HGpkg.h"		/* Global heaps				*/
#include "H5MFprivate.h"	/* File memory management		*/
#include "H5MMprivate.h"	/* Memory management			*/


/****************/
/* Local Macros */
/****************/


/******************/
/* Local Typedefs */
/******************/


/********************/
/* Package Typedefs */
/********************/


/********************/
/* Local Prototypes */
/********************/

/* Metadata cache callbacks */
static H5HG_heap_t *H5HG_load(H5F_t *f, hid_t dxpl_id, haddr_t addr, const void *udata1,
			      void *udata2);
static herr_t H5HG_flush(H5F_t *f, hid_t dxpl_id, hbool_t dest, haddr_t addr,
			 H5HG_heap_t *heap, unsigned UNUSED * flags_ptr);
static herr_t H5HG_clear(H5F_t *f, H5HG_heap_t *heap, hbool_t destroy);
static herr_t H5HG_size(const H5F_t *f, const H5HG_heap_t *heap, size_t *size_ptr);


/*********************/
/* Package Variables */
/*********************/


/*****************************/
/* Library Private Variables */
/*****************************/


/*******************/
/* Local Variables */
/*******************/

/* H5HG inherits cache-like properties from H5AC */
const H5AC_class_t H5AC_GHEAP[1] = {{
    H5AC_GHEAP_ID,
    (H5AC_load_func_t)H5HG_load,
    (H5AC_flush_func_t)H5HG_flush,
    (H5AC_dest_func_t)H5HG_dest,
    (H5AC_clear_func_t)H5HG_clear,
    (H5AC_size_func_t)H5HG_size,
}};



/*-------------------------------------------------------------------------
 * Function:	H5HG_load
 *
 * Purpose:	Loads a global heap collection from disk.
 *
 * Return:	Success:	Ptr to a global heap collection.
 *
 *		Failure:	NULL
 *
 * Programmer:	Robb Matzke
 *              Friday, March 27, 1998
 *
 *-------------------------------------------------------------------------
 */
static H5HG_heap_t *
H5HG_load(H5F_t *f, hid_t dxpl_id, haddr_t addr, const void UNUSED * udata1,
	   void UNUSED * udata2)
{
    H5HG_heap_t	*heap = NULL;
    uint8_t	*p = NULL;
    size_t	nalloc, need;
    size_t      max_idx = 0;            /* The maximum index seen */
    H5HG_heap_t	*ret_value = NULL;      /* Return value */

    FUNC_ENTER_NOAPI(H5HG_load, NULL)

    /* check arguments */
    HDassert(f);
    HDassert(H5F_addr_defined(addr));
    HDassert(!udata1);
    HDassert(!udata2);

    /* Read the initial 4k page */
    if(NULL == (heap = H5FL_CALLOC(H5HG_heap_t)))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")
    heap->addr = addr;
    if(NULL == (heap->chunk = H5FL_BLK_MALLOC(gheap_chunk, (size_t)H5HG_MINSIZE)))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")
    if(H5F_block_read(f, H5FD_MEM_GHEAP, addr, (size_t)H5HG_MINSIZE, dxpl_id, heap->chunk) < 0)
	HGOTO_ERROR(H5E_HEAP, H5E_READERROR, NULL, "unable to read global heap collection")

    /* Magic number */
    if(HDmemcmp(heap->chunk, H5HG_MAGIC, (size_t)H5_SIZEOF_MAGIC))
	HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, NULL, "bad global heap collection signature")
    p = heap->chunk + H5_SIZEOF_MAGIC;

    /* Version */
    if(H5HG_VERSION != *p++)
	HGOTO_ERROR(H5E_HEAP, H5E_CANTLOAD, NULL, "wrong version number in global heap")

    /* Reserved */
    p += 3;

    /* Size */
    H5F_DECODE_LENGTH(f, p, heap->size);
    HDassert(heap->size >= H5HG_MINSIZE);

    /*
     * If we didn't read enough in the first try, then read the rest of the
     * collection now.
     */
    if(heap->size > H5HG_MINSIZE) {
	haddr_t next_addr = addr + (hsize_t)H5HG_MINSIZE;

	if(NULL == (heap->chunk = H5FL_BLK_REALLOC(gheap_chunk, heap->chunk, heap->size)))
	    HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")
	if(H5F_block_read(f, H5FD_MEM_GHEAP, next_addr, (heap->size - H5HG_MINSIZE), dxpl_id, heap->chunk + H5HG_MINSIZE) < 0)
	    HGOTO_ERROR(H5E_HEAP, H5E_READERROR, NULL, "unable to read global heap collection")
    } /* end if */

    /* Decode each object */
    p = heap->chunk + H5HG_SIZEOF_HDR(f);
    nalloc = H5HG_NOBJS(f, heap->size);
    if(NULL == (heap->obj = H5FL_SEQ_MALLOC(H5HG_obj_t, nalloc)))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")
    heap->obj[0].size = heap->obj[0].nrefs = 0;
    heap->obj[0].begin = NULL;

    heap->nalloc = nalloc;
    while(p < (heap->chunk + heap->size)) {
	if((p + H5HG_SIZEOF_OBJHDR(f)) > (heap->chunk + heap->size)) {
	    /*
	     * The last bit of space is too tiny for an object header, so we
	     * assume that it's free space.
	     */
	    HDassert(NULL == heap->obj[0].begin);
	    heap->obj[0].size = (heap->chunk + heap->size) - p;
	    heap->obj[0].begin = p;
	    p += heap->obj[0].size;
	} /* end if */
        else {
	    unsigned idx;
	    uint8_t *begin = p;

	    UINT16DECODE(p, idx);

            /* Check if we need more room to store heap objects */
            if(idx >= heap->nalloc) {
                size_t new_alloc;       /* New allocation number */
                H5HG_obj_t *new_obj;	/* New array of object descriptions */

                /* Determine the new number of objects to index */
                new_alloc = MAX(heap->nalloc * 2, (idx + 1));

                /* Reallocate array of objects */
                if(NULL == (new_obj = H5FL_SEQ_REALLOC(H5HG_obj_t, heap->obj, new_alloc)))
                    HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")

                /* Update heap information */
                heap->nalloc = new_alloc;
                heap->obj = new_obj;
            } /* end if */

	    UINT16DECODE(p, heap->obj[idx].nrefs);
	    p += 4; /*reserved*/
	    H5F_DECODE_LENGTH(f, p, heap->obj[idx].size);
	    heap->obj[idx].begin = begin;

	    /*
	     * The total storage size includes the size of the object header
	     * and is zero padded so the next object header is properly
	     * aligned. The last bit of space is the free space object whose
	     * size is never padded and already includes the object header.
	     */
	    if(idx > 0) {
		need = H5HG_SIZEOF_OBJHDR(f) + H5HG_ALIGN(heap->obj[idx].size);

                /* Check for "gap" in index numbers (caused by deletions) and fill in heap object values */
                if(idx > (max_idx + 1))
                    HDmemset(&heap->obj[max_idx + 1], 0, sizeof(H5HG_obj_t) * (idx - (max_idx + 1)));
                max_idx = idx;
	    } /* end if */
            else
		need = heap->obj[idx].size;
	    p = begin + need;
	} /* end else */
    } /* end while */
    HDassert(p == heap->chunk + heap->size);
    HDassert(H5HG_ISALIGNED(heap->obj[0].size));

    /* Set the next index value to use */
    if(max_idx > 0)
        heap->nused = max_idx + 1;
    else
        heap->nused = 1;

    /*
     * Add the new heap to the CWFS list, removing some other entry if
     * necessary to make room. We remove the right-most entry that has less
     * free space than this heap.
     */
    if(!f->shared->cwfs) {
        f->shared->cwfs = (H5HG_heap_t **)H5MM_malloc(H5HG_NCWFS * sizeof(H5HG_heap_t *));
        if(NULL == f->shared->cwfs)
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")
        f->shared->ncwfs = 1;
        f->shared->cwfs[0] = heap;
    } else if(H5HG_NCWFS == f->shared->ncwfs) {
        int i;          /* Local index variable */

        for(i = H5HG_NCWFS - 1; i >= 0; --i)
            if(f->shared->cwfs[i]->obj[0].size < heap->obj[0].size) {
                HDmemmove(f->shared->cwfs + 1, f->shared->cwfs, i * sizeof(H5HG_heap_t *));
                f->shared->cwfs[0] = heap;
                break;
            } /* end if */
    } else {
        HDmemmove(f->shared->cwfs + 1, f->shared->cwfs, f->shared->ncwfs * sizeof(H5HG_heap_t *));
        f->shared->ncwfs += 1;
        f->shared->cwfs[0] = heap;
    } /* end else */

    ret_value = heap;

done:
    if(!ret_value && heap)
        if(H5HG_dest(f, heap) < 0)
	    HDONE_ERROR(H5E_HEAP, H5E_CANTFREE, NULL, "unable to destroy global heap collection")

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5HG_load() */


/*-------------------------------------------------------------------------
 * Function:	H5HG_flush
 *
 * Purpose:	Flushes a global heap collection from memory to disk if it's
 *		dirty.  Optionally deletes teh heap from memory.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Friday, March 27, 1998
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5HG_flush(H5F_t *f, hid_t dxpl_id, hbool_t destroy, haddr_t addr, H5HG_heap_t *heap, unsigned UNUSED * flags_ptr)
{
    herr_t ret_value = SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5HG_flush, FAIL)

    /* Check arguments */
    HDassert(f);
    HDassert(H5F_addr_defined(addr));
    HDassert(H5F_addr_eq(addr, heap->addr));
    HDassert(heap);

    if(heap->cache_info.is_dirty) {
	if(H5F_block_write(f, H5FD_MEM_GHEAP, addr, heap->size, dxpl_id, heap->chunk) < 0)
	    HGOTO_ERROR(H5E_HEAP, H5E_WRITEERROR, FAIL, "unable to write global heap collection to file")
	heap->cache_info.is_dirty = FALSE;
    } /* end if */

    if(destroy)
        if(H5HG_dest(f, heap) < 0)
	    HGOTO_ERROR(H5E_HEAP, H5E_CANTFREE, FAIL, "unable to destroy global heap collection")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5HG_flush() */


/*-------------------------------------------------------------------------
 * Function:	H5HG_dest
 *
 * Purpose:	Destroys a global heap collection in memory
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Wednesday, January 15, 2003
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5HG_dest(H5F_t *f, H5HG_heap_t *heap)
{
    int	i;                              /* Local index variable */
    herr_t ret_value = SUCCEED;         /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5HG_dest)

    /* Check arguments */
    HDassert(heap);

    /* Verify that node is clean */
    HDassert(heap->cache_info.is_dirty == FALSE);

    /* If we're going to free the space on disk, the address must be valid */
    HDassert(!heap->cache_info.free_file_space_on_destroy || H5F_addr_defined(heap->cache_info.addr));

    /* Check for freeing file space for globalheap */
    if(heap->cache_info.free_file_space_on_destroy) {
        /* Release the space on disk */
        /* (XXX: Nasty usage of internal DXPL value! -QAK) */
        if(H5MF_xfree(f, H5FD_MEM_GHEAP, H5AC_dxpl_id, heap->cache_info.addr, (hsize_t)heap->size) < 0)
            HGOTO_ERROR(H5E_BTREE, H5E_CANTFREE, FAIL, "unable to free global heap")
    } /* end if */

    /* Remove heap from NCWFS array, if it's present */
    for(i = 0; i < f->shared->ncwfs; i++)
        if(f->shared->cwfs[i] == heap) {
            f->shared->ncwfs -= 1;
            HDmemmove(f->shared->cwfs + i, f->shared->cwfs + i + 1, (f->shared->ncwfs - i) * sizeof(H5HG_heap_t *));
            break;
        } /* end if */

    /* Release resources */
    if(heap->chunk)
        heap->chunk = H5FL_BLK_FREE(gheap_chunk, heap->chunk);
    if(heap->obj)
        heap->obj = H5FL_SEQ_FREE(H5HG_obj_t, heap->obj);
    (void)H5FL_FREE(H5HG_heap_t, heap);

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* H5HG_dest() */


/*-------------------------------------------------------------------------
 * Function:	H5HG_clear
 *
 * Purpose:	Mark a global heap in memory as non-dirty.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, March 20, 2003
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5HG_clear(H5F_t *f, H5HG_heap_t *heap, hbool_t destroy)
{
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_NOAPI_NOINIT(H5HG_clear)

    /* Check arguments */
    HDassert(heap);

    /* Mark heap as clean */
    heap->cache_info.is_dirty = FALSE;

    if(destroy)
        if(H5HG_dest(f, heap) < 0)
	    HGOTO_ERROR(H5E_HEAP, H5E_CANTFREE, FAIL, "unable to destroy global heap collection")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* H5HG_clear() */


/*-------------------------------------------------------------------------
 * Function:	H5HG_size
 *
 * Purpose:	Compute the size in bytes of the specified instance of
 *              H5HG_heap_t on disk, and return it in *len_ptr.  On failure,
 *              the value of *len_ptr is undefined.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	John Mainzer
 *              5/13/04
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5HG_size(const H5F_t UNUSED *f, const H5HG_heap_t *heap, size_t *size_ptr)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5HG_size)

    /* Check arguments */
    HDassert(heap);
    HDassert(size_ptr);

    *size_ptr = heap->size;

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* H5HG_size() */

