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

/* Programmer: 	Robb Matzke <matzke@llnl.gov>
 *	       	Wednesday, October  8, 1997
 *
 * Purpose:	Indexed (chunked) I/O functions.  The logical
 *		multi-dimensional data space is regularly partitioned into
 *		same-sized "chunks", the first of which is aligned with the
 *		logical origin.  The chunks are given a multi-dimensional
 *		index which is used as a lookup key in a B-tree that maps
 *		chunk index to disk address.  Each chunk can be compressed
 *		independently and the chunks may move around in the file as
 *		their storage requirements change.
 *
 * Cache:	Disk I/O is performed in units of chunks and H5MF_alloc()
 *		contains code to optionally align chunks on disk block
 *		boundaries for performance.
 *
 *		The chunk cache is an extendible hash indexed by a function
 *		of storage B-tree address and chunk N-dimensional offset
 *		within the dataset.  Collisions are not resolved -- one of
 *		the two chunks competing for the hash slot must be preempted
 *		from the cache.  All entries in the hash also participate in
 *		a doubly-linked list and entries are penalized by moving them
 *		toward the front of the list.  When a new chunk is about to
 *		be added to the cache the heap is pruned by preempting
 *		entries near the front of the list to make room for the new
 *		entry which is added to the end of the list.
 */

/****************/
/* Module Setup */
/****************/

#define H5B_PACKAGE		/*suppress error about including H5Bpkg	  */
#define H5D_PACKAGE		/*suppress error about including H5Dpkg	  */


/***********/
/* Headers */
/***********/
#include "H5private.h"		/* Generic Functions			*/
#include "H5Bpkg.h"		/* B-link trees				*/
#include "H5Dpkg.h"		/* Datasets				*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5Fprivate.h"		/* Files				*/
#include "H5FDprivate.h"	/* File drivers				*/
#include "H5FLprivate.h"	/* Free Lists                           */
#include "H5Iprivate.h"		/* IDs			  		*/
#include "H5MFprivate.h"	/* File space management		*/
#include "H5MMprivate.h"	/* Memory management			*/
#include "H5Oprivate.h"		/* Object headers		  	*/
#include "H5Pprivate.h"         /* Property lists                       */
#include "H5Sprivate.h"         /* Dataspaces                           */
#include "H5SLprivate.h"	/* Skip lists				*/
#include "H5Vprivate.h"		/* Vector and array functions		*/

/****************/
/* Local Macros */
/****************/

/*
 * Feature: If this constant is defined then every cache preemption and load
 *	    causes a character to be printed on the standard error stream:
 *
 *     `.': Entry was preempted because it has been completely read or
 *	    completely written but not partially read and not partially
 *	    written. This is often a good reason for preemption because such
 *	    a chunk will be unlikely to be referenced in the near future.
 *
 *     `:': Entry was preempted because it hasn't been used recently.
 *
 *     `#': Entry was preempted because another chunk collided with it. This
 *	    is usually a relatively bad thing.  If there are too many of
 *	    these then the number of entries in the cache can be increased.
 *
 *       c: Entry was preempted because the file is closing.
 *
 *	 w: A chunk read operation was eliminated because the library is
 *	    about to write new values to the entire chunk.  This is a good
 *	    thing, especially on files where the chunk size is the same as
 *	    the disk block size, chunks are aligned on disk block boundaries,
 *	    and the operating system can also eliminate a read operation.
 */

/*#define H5D_ISTORE_DEBUG */

/*
 * Given a B-tree node return the dimensionality of the chunks pointed to by
 * that node.
 */
#define H5D_ISTORE_NDIMS(X)	(((X)->sizeof_rkey-8)/8)

#define H5D_HASH(D,ADDR) H5F_addr_hash(ADDR,(D)->cache.chunk.nslots)

#define H5D_ISTORE_DEFAULT_SKIPLIST_HEIGHT     8

/******************/
/* Local Typedefs */
/******************/

/* Raw data chunks are cached.  Each entry in the cache is: */
typedef struct H5D_rdcc_ent_t {
    hbool_t	locked;		/*entry is locked in cache		*/
    hbool_t	dirty;		/*needs to be written to disk?		*/
    hsize_t	offset[H5O_LAYOUT_NDIMS]; /*chunk name			*/
    size_t	rd_count;	/*bytes remaining to be read		*/
    size_t	wr_count;	/*bytes remaining to be written		*/
    size_t	chunk_size;	/*size of a chunk			*/
    size_t	alloc_size;	/*amount allocated for the chunk	*/
    uint8_t	*chunk;		/*the unfiltered chunk data		*/
    unsigned	idx;		/*index in hash table			*/
    struct H5D_rdcc_ent_t *next;/*next item in doubly-linked list	*/
    struct H5D_rdcc_ent_t *prev;/*previous item in doubly-linked list	*/
} H5D_rdcc_ent_t;
typedef H5D_rdcc_ent_t *H5D_rdcc_ent_ptr_t; /* For free lists */

/*
 * Data exchange structure for indexed storage nodes.  This structure is
 * passed through the B-link tree layer to the methods for the objects
 * to which the B-link tree points for operations which require no
 * additional information.
 *
 * (Just an alias for the "common" info).
 */
typedef H5D_istore_bt_ud_common_t H5D_istore_ud0_t;

/* B-tree callback info for iteration to total allocated space */
typedef struct H5D_istore_it_ud1_t {
    H5D_istore_bt_ud_common_t common;           /* Common info for B-tree user data (must be first) */
    hsize_t		total_storage;	        /*output from iterator	*/
} H5D_istore_it_ud1_t;

/* B-tree callback info for iteration to dump node's info */
typedef struct H5D_istore_it_ud2_t {
    H5D_istore_bt_ud_common_t common;           /* Common info for B-tree user data (must be first) */
    FILE		*stream;		/*debug output stream	*/
    hbool_t             header_displayed;       /* Node's header is displayed? */
} H5D_istore_it_ud2_t;

/* B-tree callback info for iteration to prune chunks */
typedef struct H5D_istore_it_ud3_t {
    H5D_istore_bt_ud_common_t common;           /* Common info for B-tree user data (must be first) */
    const hsize_t	*dims;		        /* New dataset dimensions	*/
    const hsize_t       *down_chunks;           /* "down" size of number of chunks in each dimension */
    H5SL_t              *outside;               /* Skip list to hold chunks outside the new dimensions */
} H5D_istore_it_ud3_t;

/* B-tree callback info for iteration to copy data */
typedef struct H5D_istore_it_ud4_t {
    H5D_istore_bt_ud_common_t common;           /* Common info for B-tree user data (must be first) */
    H5F_t               *file_src;              /* Source file for copy */
    haddr_t             addr_dst;               /* Address of dest. B-tree */
    void                *buf;                   /* Buffer to hold chunk data for read/write */
    void                *bkg;                   /* Buffer for background information during type conversion */
    size_t              buf_size;               /* Buffer size */
    hbool_t             do_convert;             /* Whether to perform type conversions */

    /* needed for converting variable-length data */
    hid_t               tid_src;                /* Datatype ID for source datatype */
    hid_t               tid_dst;                /* Datatype ID for destination datatype */
    hid_t               tid_mem;                /* Datatype ID for memory datatype */
    H5T_t               *dt_src;                /* Source datatype */
    H5T_path_t          *tpath_src_mem;         /* Datatype conversion path from source file to memory */
    H5T_path_t          *tpath_mem_dst;         /* Datatype conversion path from memory to dest. file */
    void                *reclaim_buf;           /* Buffer for reclaiming data */
    size_t              reclaim_buf_size;       /* Reclaim buffer size */
    size_t              nelmts;                 /* Number of elements in buffer */
    H5S_t               *buf_space;             /* Dataspace describing buffer */

    /* needed for compressed variable-length data */
    H5O_pline_t         *pline;                 /* Filter pipeline */

    /* needed for copy object pointed by refs */
    H5F_t               *file_dst;              /* Destination file for copy */
    H5O_copy_t          *cpy_info;              /* Copy options */
} H5D_istore_it_ud4_t;

/* B-tree callback info for iteration to obtain chunk address and the index of the chunk for all chunks in the B-tree. */
typedef struct H5D_istore_it_ud5_t {
    H5D_istore_bt_ud_common_t common;           /* Common info for B-tree user data (must be first) */
    hsize_t             *down_chunks;
    haddr_t             *chunk_addr;
} H5D_istore_it_ud5_t;

/* Skip list node for storing chunks to remove during an iteration */
typedef struct H5D_istore_sl_ck_t {
    hsize_t             index;                  /* Index of chunk to remove (must be first) */
    H5D_istore_key_t	key;	                /* Chunk key */
} H5D_istore_sl_ck_t;

/* Skip list callback info when destroying list & removing chunks */
typedef struct H5D_istore_sl_rm_t {
    H5F_t               *f;                     /* Pointer to file for B-tree */
    hid_t               dxpl_id;                /* DXPL to use */
    const H5O_layout_t	*mesg;		        /* Layout message	*/
} H5D_istore_sl_rm_t;

/********************/
/* Local Prototypes */
/********************/

static void *H5D_istore_chunk_alloc(size_t size, const H5O_pline_t *pline);
static void *H5D_istore_chunk_xfree(void *chk, const H5O_pline_t *pline);
static herr_t H5D_istore_shared_create (const H5F_t *f, H5O_layout_t *layout);
static herr_t H5D_istore_shared_free (void *page);

/* B-tree iterator callbacks */
static int H5D_istore_iter_chunkmap(H5F_t *f, hid_t dxpl_id, const void *left_key, haddr_t addr,
                                 const void *right_key, void *_udata);
static int H5D_istore_iter_allocated(H5F_t *f, hid_t dxpl_id, const void *left_key, haddr_t addr,
				 const void *right_key, void *_udata);
static int H5D_istore_iter_dump(H5F_t *f, hid_t dxpl_id, const void *left_key, haddr_t addr,
				 const void *right_key, void *_udata);
static int H5D_istore_prune_check(H5F_t *f, hid_t dxpl_id, const void *_lt_key, haddr_t addr,
        const void *_rt_key, void *_udata);
static int H5D_istore_iter_copy(H5F_t *f, hid_t dxpl_id, const void *_lt_key, haddr_t addr,
                                 const void *_rt_key, void *_udata);

/* B-tree callbacks */
static H5RC_t *H5D_istore_get_shared(const H5F_t *f, const void *_udata);
static herr_t H5D_istore_new_node(H5F_t *f, hid_t dxpl_id, H5B_ins_t, void *_lt_key,
				  void *_udata, void *_rt_key,
				  haddr_t *addr_p /*out*/);
static int H5D_istore_cmp2(H5F_t *f, hid_t dxpl_id, void *_lt_key, void *_udata,
			    void *_rt_key);
static int H5D_istore_cmp3(H5F_t *f, hid_t dxpl_id, void *_lt_key, void *_udata,
			    void *_rt_key);
static herr_t H5D_istore_found(H5F_t *f, hid_t dxpl_id, haddr_t addr, const void *_lt_key,
			       void *_udata);
static H5B_ins_t H5D_istore_insert(H5F_t *f, hid_t dxpl_id, haddr_t addr, void *_lt_key,
				   hbool_t *lt_key_changed, void *_md_key,
				   void *_udata, void *_rt_key,
				   hbool_t *rt_key_changed,
				   haddr_t *new_node/*out*/);
static H5B_ins_t H5D_istore_remove( H5F_t *f, hid_t dxpl_id, haddr_t addr, void *_lt_key,
                  hbool_t *lt_key_changed, void *_udata, void *_rt_key,
                  hbool_t *rt_key_changed);
static herr_t H5D_istore_decode_key(const H5F_t *f, const H5B_t *bt, const uint8_t *raw,
				    void *_key);
static herr_t H5D_istore_encode_key(const H5F_t *f, const H5B_t *bt, uint8_t *raw,
				    void *_key);
static herr_t H5D_istore_debug_key(FILE *stream, H5F_t *f, hid_t dxpl_id,
                                int indent, int fwidth, const void *key,
                                    const void *udata);

/* inherits B-tree like properties from H5B */
H5B_class_t H5B_ISTORE[1] = {{
    H5B_ISTORE_ID,		/*id			*/
    sizeof(H5D_istore_key_t),	/*sizeof_nkey		*/
    H5D_istore_get_shared,	/*get_shared		*/
    H5D_istore_new_node,	/*new			*/
    H5D_istore_cmp2,		/*cmp2			*/
    H5D_istore_cmp3,		/*cmp3			*/
    H5D_istore_found,		/*found			*/
    H5D_istore_insert,		/*insert		*/
    FALSE,			/*follow min branch?	*/
    FALSE,			/*follow max branch?	*/
    H5D_istore_remove,          /*remove		*/
    H5D_istore_decode_key,	/*decode		*/
    H5D_istore_encode_key,	/*encode		*/
    H5D_istore_debug_key,	/*debug			*/
}};

/*********************/
/* Package Variables */
/*********************/

/*****************************/
/* Library Private Variables */
/*****************************/

/*******************/
/* Local Variables */
/*******************/

/* Declare a free list to manage H5F_rdcc_ent_t objects */
H5FL_DEFINE_STATIC(H5D_rdcc_ent_t);

/* Declare a free list to manage the H5F_rdcc_ent_ptr_t sequence information */
H5FL_SEQ_DEFINE_STATIC(H5D_rdcc_ent_ptr_t);

/* Declare a free list to manage the chunk sequence information */
H5FL_BLK_DEFINE_STATIC(chunk);

/* Declare a free list to manage the native key offset sequence information */
H5FL_SEQ_DEFINE_STATIC(size_t);

/* Declare a free list to manage the raw page information */
H5FL_BLK_DEFINE_STATIC(chunk_page);

/* Declare a free list to manage H5D_istore_sl_ck_t objects */
H5FL_DEFINE_STATIC(H5D_istore_sl_ck_t);


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_get_shared
 *
 * Purpose:	Returns the shared B-tree info for the specified UDATA.
 *
 * Return:	Success:	Pointer to the raw B-tree page for this dataset
 *
 *		Failure:	Can't fail
 *
 * Programmer:	Quincey Koziol
 *		Monday, July  5, 2004
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static H5RC_t *
H5D_istore_get_shared(const H5F_t UNUSED *f, const void *_udata)
{
    const H5D_istore_ud0_t *udata = (const H5D_istore_ud0_t *) _udata;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_get_shared)

    HDassert(udata);
    HDassert(udata->mesg);
    HDassert(udata->mesg->u.chunk.btree_shared);

    /* Increment reference count on B-tree info */
    H5RC_INC(udata->mesg->u.chunk.btree_shared);

    /* Return the pointer to the ref-count object */
    FUNC_LEAVE_NOAPI(udata->mesg->u.chunk.btree_shared)
} /* end H5D_istore_get_shared() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_decode_key
 *
 * Purpose:	Decodes a raw key into a native key for the B-tree
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Friday, October 10, 1997
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_decode_key(const H5F_t UNUSED *f, const H5B_t *bt, const uint8_t *raw, void *_key)
{
    H5D_istore_key_t	*key = (H5D_istore_key_t *) _key;
    H5B_shared_t        *shared;        /* Pointer to shared B-tree info */
    size_t		ndims;
    unsigned		u;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_decode_key)

    /* check args */
    HDassert(f);
    HDassert(bt);
    shared = (H5B_shared_t *)H5RC_GET_OBJ(bt->rc_shared);
    HDassert(shared);
    HDassert(raw);
    HDassert(key);
    ndims = H5D_ISTORE_NDIMS(shared);
    HDassert(ndims <= H5O_LAYOUT_NDIMS);

    /* decode */
    UINT32DECODE(raw, key->nbytes);
    UINT32DECODE(raw, key->filter_mask);
    for(u = 0; u < ndims; u++)
	UINT64DECODE(raw, key->offset[u]);

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* end H5D_istore_decode_key() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_encode_key
 *
 * Purpose:	Encode a key from native format to raw format.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Friday, October 10, 1997
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_encode_key(const H5F_t UNUSED *f, const H5B_t *bt, uint8_t *raw, void *_key)
{
    H5D_istore_key_t	*key = (H5D_istore_key_t *) _key;
    H5B_shared_t        *shared;        /* Pointer to shared B-tree info */
    size_t		ndims;
    unsigned		u;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_encode_key)

    /* check args */
    HDassert(f);
    HDassert(bt);
    shared = (H5B_shared_t *)H5RC_GET_OBJ(bt->rc_shared);
    HDassert(shared);
    HDassert(raw);
    HDassert(key);
    ndims = H5D_ISTORE_NDIMS(shared);
    HDassert(ndims <= H5O_LAYOUT_NDIMS);

    /* encode */
    UINT32ENCODE(raw, key->nbytes);
    UINT32ENCODE(raw, key->filter_mask);
    for(u = 0; u < ndims; u++)
	UINT64ENCODE(raw, key->offset[u]);

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* end H5D_istore_encode_key() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_debug_key
 *
 * Purpose:	Prints a key.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, April 16, 1998
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static herr_t
H5D_istore_debug_key(FILE *stream, H5F_t UNUSED *f, hid_t UNUSED dxpl_id, int indent, int fwidth,
		      const void *_key, const void *_udata)
{
    const H5D_istore_key_t	*key = (const H5D_istore_key_t *)_key;
    const H5D_istore_ud0_t	*udata = (const H5D_istore_ud0_t *)_udata;
    unsigned		u;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_debug_key)

    HDassert(key);

    HDfprintf(stream, "%*s%-*s %Zd bytes\n", indent, "", fwidth, "Chunk size:", key->nbytes);
    HDfprintf(stream, "%*s%-*s 0x%08x\n", indent, "", fwidth, "Filter mask:", key->filter_mask);
    HDfprintf(stream, "%*s%-*s {", indent, "", fwidth, "Logical offset:");
    for(u = 0; u < udata->mesg->u.chunk.ndims; u++)
        HDfprintf(stream, "%s%Hd", u?", ":"", key->offset[u]);
    HDfputs("}\n", stream);

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* end H5D_istore_debug_key() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_cmp2
 *
 * Purpose:	Compares two keys sort of like strcmp().  The UDATA pointer
 *		is only to supply extra information not carried in the keys
 *		(in this case, the dimensionality) and is not compared
 *		against the keys.
 *
 * Return:	Success:	-1 if LT_KEY is less than RT_KEY;
 *				1 if LT_KEY is greater than RT_KEY;
 *				0 if LT_KEY and RT_KEY are equal.
 *
 *		Failure:	FAIL (same as LT_KEY<RT_KEY)
 *
 * Programmer:	Robb Matzke
 *		Thursday, November  6, 1997
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static int
H5D_istore_cmp2(H5F_t UNUSED *f, hid_t UNUSED dxpl_id, void *_lt_key, void *_udata,
		void *_rt_key)
{
    H5D_istore_key_t	*lt_key = (H5D_istore_key_t *) _lt_key;
    H5D_istore_key_t	*rt_key = (H5D_istore_key_t *) _rt_key;
    H5D_istore_ud0_t	*udata = (H5D_istore_ud0_t *) _udata;
    int		ret_value;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_cmp2)

    HDassert(lt_key);
    HDassert(rt_key);
    HDassert(udata);
    HDassert(udata->mesg->u.chunk.ndims > 0 && udata->mesg->u.chunk.ndims <= H5O_LAYOUT_NDIMS);

    /* Compare the offsets but ignore the other fields */
    ret_value = H5V_vector_cmp_u(udata->mesg->u.chunk.ndims, lt_key->offset, rt_key->offset);

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_cmp2() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_cmp3
 *
 * Purpose:	Compare the requested datum UDATA with the left and right
 *		keys of the B-tree.
 *
 * Return:	Success:	negative if the min_corner of UDATA is less
 *				than the min_corner of LT_KEY.
 *
 *				positive if the min_corner of UDATA is
 *				greater than or equal the min_corner of
 *				RT_KEY.
 *
 *				zero otherwise.	 The min_corner of UDATA is
 *				not necessarily contained within the address
 *				space represented by LT_KEY, but a key that
 *				would describe the UDATA min_corner address
 *				would fall lexicographically between LT_KEY
 *				and RT_KEY.
 *
 *		Failure:	FAIL (same as UDATA < LT_KEY)
 *
 * Programmer:	Robb Matzke
 *		Wednesday, October  8, 1997
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static int
H5D_istore_cmp3(H5F_t UNUSED *f, hid_t UNUSED dxpl_id, void *_lt_key, void *_udata,
		void *_rt_key)
{
    H5D_istore_key_t	*lt_key = (H5D_istore_key_t *) _lt_key;
    H5D_istore_key_t	*rt_key = (H5D_istore_key_t *) _rt_key;
    H5D_istore_ud0_t	*udata = (H5D_istore_ud0_t *) _udata;
    int		ret_value = 0;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_cmp3)

    HDassert(lt_key);
    HDassert(rt_key);
    HDassert(udata);
    HDassert(udata->mesg->u.chunk.ndims > 0 && udata->mesg->u.chunk.ndims <= H5O_LAYOUT_NDIMS);

    /* Special case for faster checks on 1-D chunks */
    /* (Checking for ndims==2 because last dimension is the datatype size) */
    /* The additional checking for the right key is necessary due to the */
    /* slightly odd way the library initializes the right-most node in the */
    /* indexed storage B-tree... */
    /* (Dump the B-tree with h5debug to look at it) -QAK */
    if(udata->mesg->u.chunk.ndims == 2) {
        if(udata->offset[0] > rt_key->offset[0])
            ret_value = 1;
        else if(udata->offset[0] == rt_key->offset[0] &&
                udata->offset[1] >= rt_key->offset[1])
            ret_value = 1;
        else if(udata->offset[0] < lt_key->offset[0])
            ret_value = (-1);
    } /* end if */
    else {
        if(H5V_vector_ge_u(udata->mesg->u.chunk.ndims, udata->offset, rt_key->offset))
            ret_value = 1;
        else if(H5V_vector_lt_u(udata->mesg->u.chunk.ndims, udata->offset, lt_key->offset))
            ret_value = (-1);
    } /* end else */

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_cmp3() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_new_node
 *
 * Purpose:	Adds a new entry to an i-storage B-tree.  We can assume that
 *		the domain represented by UDATA doesn't intersect the domain
 *		already represented by the B-tree.
 *
 * Return:	Success:	Non-negative. The address of leaf is returned
 *				through the ADDR argument.  It is also added
 *				to the UDATA.
 *
 * 		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *		Tuesday, October 14, 1997
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_new_node(H5F_t *f, hid_t dxpl_id, H5B_ins_t op,
		    void *_lt_key, void *_udata, void *_rt_key,
		    haddr_t *addr_p/*out*/)
{
    H5D_istore_key_t	*lt_key = (H5D_istore_key_t *) _lt_key;
    H5D_istore_key_t	*rt_key = (H5D_istore_key_t *) _rt_key;
    H5D_istore_ud1_t	*udata = (H5D_istore_ud1_t *) _udata;
    unsigned		u;
    herr_t      ret_value = SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_new_node)

    /* check args */
    HDassert(f);
    HDassert(lt_key);
    HDassert(rt_key);
    HDassert(udata);
    HDassert(udata->common.mesg->u.chunk.ndims > 0 && udata->common.mesg->u.chunk.ndims < H5O_LAYOUT_NDIMS);
    HDassert(addr_p);

    /* Allocate new storage */
    HDassert(udata->nbytes > 0);
    H5_CHECK_OVERFLOW(udata->nbytes, size_t, hsize_t);
    if(HADDR_UNDEF == (*addr_p = H5MF_alloc(f, H5FD_MEM_DRAW, dxpl_id, (hsize_t)udata->nbytes)))
        HGOTO_ERROR(H5E_IO, H5E_CANTINIT, FAIL, "couldn't allocate new file storage")
    udata->addr = *addr_p;

    /*
     * The left key describes the storage of the UDATA chunk being
     * inserted into the tree.
     */
    lt_key->nbytes = udata->nbytes;
    lt_key->filter_mask = udata->filter_mask;
    for(u = 0; u < udata->common.mesg->u.chunk.ndims; u++)
        lt_key->offset[u] = udata->common.offset[u];

    /*
     * The right key might already be present.  If not, then add a zero-width
     * chunk.
     */
    if(H5B_INS_LEFT != op) {
        rt_key->nbytes = 0;
        rt_key->filter_mask = 0;
        for(u = 0; u < udata->common.mesg->u.chunk.ndims; u++) {
            HDassert(udata->common.offset[u] + udata->common.mesg->u.chunk.dim[u] >
                udata->common.offset[u]);
            rt_key->offset[u] = udata->common.offset[u] + udata->common.mesg->u.chunk.dim[u];
        } /* end if */
    } /* end if */

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_new_node() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_found
 *
 * Purpose:	This function is called when the B-tree search engine has
 *		found the leaf entry that points to a chunk of storage that
 *		contains the beginning of the logical address space
 *		represented by UDATA.  The LT_KEY is the left key (the one
 *		that describes the chunk) and RT_KEY is the right key (the
 *		one that describes the next or last chunk).
 *
 * Note:	It's possible that the chunk isn't really found.  For
 *		instance, in a sparse dataset the requested chunk might fall
 *		between two stored chunks in which case this function is
 *		called with the maximum stored chunk indices less than the
 *		requested chunk indices.
 *
 * Return:	Non-negative on success with information about the chunk
 *		returned through the UDATA argument. Negative on failure.
 *
 * Programmer:	Robb Matzke
 *		Thursday, October  9, 1997
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static herr_t
H5D_istore_found(H5F_t UNUSED *f, hid_t UNUSED dxpl_id, haddr_t addr, const void *_lt_key,
		 void *_udata)
{
    H5D_istore_ud1_t	   *udata = (H5D_istore_ud1_t *) _udata;
    const H5D_istore_key_t *lt_key = (const H5D_istore_key_t *) _lt_key;
    unsigned		u;
    herr_t      ret_value = SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_found)

    /* Check arguments */
    HDassert(f);
    HDassert(H5F_addr_defined(addr));
    HDassert(udata);
    HDassert(lt_key);

    /* Is this *really* the requested chunk? */
    for(u = 0; u < udata->common.mesg->u.chunk.ndims; u++)
        if(udata->common.offset[u] >= lt_key->offset[u] + udata->common.mesg->u.chunk.dim[u])
            HGOTO_DONE(FAIL)

    /* Initialize return values */
    HDassert(lt_key->nbytes > 0);
    udata->addr = addr;
    udata->nbytes = lt_key->nbytes;
    udata->filter_mask = lt_key->filter_mask;

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_found() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_insert
 *
 * Purpose:	This function is called when the B-tree insert engine finds
 *		the node to use to insert new data.  The UDATA argument
 *		points to a struct that describes the logical addresses being
 *		added to the file.  This function allocates space for the
 *		data and returns information through UDATA describing a
 *		file chunk to receive (part of) the data.
 *
 *		The LT_KEY is always the key describing the chunk of file
 *		memory at address ADDR. On entry, UDATA describes the logical
 *		addresses for which storage is being requested (through the
 *		`offset' and `size' fields). On return, UDATA describes the
 *		logical addresses contained in a chunk on disk.
 *
 * Return:	Success:	An insertion command for the caller, one of
 *				the H5B_INS_* constants.  The address of the
 *				new chunk is returned through the NEW_NODE
 *				argument.
 *
 *		Failure:	H5B_INS_ERROR
 *
 * Programmer:	Robb Matzke
 *		Thursday, October  9, 1997
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static H5B_ins_t
H5D_istore_insert(H5F_t *f, hid_t dxpl_id, haddr_t addr, void *_lt_key,
		  hbool_t *lt_key_changed,
		  void *_md_key, void *_udata, void *_rt_key,
		  hbool_t UNUSED *rt_key_changed,
		  haddr_t *new_node_p/*out*/)
{
    H5D_istore_key_t	*lt_key = (H5D_istore_key_t *) _lt_key;
    H5D_istore_key_t	*md_key = (H5D_istore_key_t *) _md_key;
    H5D_istore_key_t	*rt_key = (H5D_istore_key_t *) _rt_key;
    H5D_istore_ud1_t	*udata = (H5D_istore_ud1_t *) _udata;
    int		cmp;
    unsigned		u;
    H5B_ins_t		ret_value;

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_insert)

    /* check args */
    HDassert(f);
    HDassert(H5F_addr_defined(addr));
    HDassert(lt_key);
    HDassert(lt_key_changed);
    HDassert(md_key);
    HDassert(udata);
    HDassert(rt_key);
    HDassert(new_node_p);

    cmp = H5D_istore_cmp3(f, dxpl_id, lt_key, udata, rt_key);
    HDassert(cmp <= 0);

    if(cmp < 0) {
        /* Negative indices not supported yet */
        HGOTO_ERROR(H5E_STORAGE, H5E_UNSUPPORTED, H5B_INS_ERROR, "internal error")

    } else if(H5V_vector_eq_u(udata->common.mesg->u.chunk.ndims,
				udata->common.offset, lt_key->offset) &&
	       lt_key->nbytes > 0) {
        /*
         * Already exists.  If the new size is not the same as the old size
         * then we should reallocate storage.
         */
        if(lt_key->nbytes != udata->nbytes) {
/* Currently, the old chunk data is "thrown away" after the space is reallocated,
 * so avoid data copy in H5MF_realloc() call by just free'ing the space and
 * allocating new space.
 *
 * This should keep the file smaller also, by freeing the space and then
 * allocating new space, instead of vice versa (in H5MF_realloc).
 *
 * QAK - 11/19/2002
 */
#ifdef OLD_WAY
            if(HADDR_UNDEF == (*new_node_p = H5MF_realloc(f, H5FD_MEM_DRAW, addr,
                      (hsize_t)lt_key->nbytes, (hsize_t)udata->nbytes)))
                HGOTO_ERROR(H5E_STORAGE, H5E_NOSPACE, H5B_INS_ERROR, "unable to reallocate chunk storage")
#else /* OLD_WAY */
            H5_CHECK_OVERFLOW( lt_key->nbytes ,size_t, hsize_t);
            if(H5MF_xfree(f, H5FD_MEM_DRAW, dxpl_id, addr, (hsize_t)lt_key->nbytes)<0)
                HGOTO_ERROR(H5E_STORAGE, H5E_CANTFREE, H5B_INS_ERROR, "unable to free chunk")
            H5_CHECK_OVERFLOW(udata->nbytes ,size_t, hsize_t);
            if(HADDR_UNDEF == (*new_node_p = H5MF_alloc(f, H5FD_MEM_DRAW, dxpl_id, (hsize_t)udata->nbytes)))
                HGOTO_ERROR(H5E_STORAGE, H5E_NOSPACE, H5B_INS_ERROR, "unable to reallocate chunk")
#endif /* OLD_WAY */
            lt_key->nbytes = udata->nbytes;
            lt_key->filter_mask = udata->filter_mask;
            *lt_key_changed = TRUE;
            udata->addr = *new_node_p;
            ret_value = H5B_INS_CHANGE;
        } else {
            udata->addr = addr;
            ret_value = H5B_INS_NOOP;
        }

    } else if (H5V_hyper_disjointp(udata->common.mesg->u.chunk.ndims,
				   lt_key->offset, udata->common.mesg->u.chunk.dim,
				   udata->common.offset, udata->common.mesg->u.chunk.dim)) {
        HDassert(H5V_hyper_disjointp(udata->common.mesg->u.chunk.ndims,
				   rt_key->offset, udata->common.mesg->u.chunk.dim,
				   udata->common.offset, udata->common.mesg->u.chunk.dim));
        /*
         * Split this node, inserting the new new node to the right of the
         * current node.  The MD_KEY is where the split occurs.
         */
        md_key->nbytes = udata->nbytes;
        md_key->filter_mask = udata->filter_mask;
        for(u = 0; u < udata->common.mesg->u.chunk.ndims; u++) {
            HDassert(0 == udata->common.offset[u] % udata->common.mesg->u.chunk.dim[u]);
            md_key->offset[u] = udata->common.offset[u];
        } /* end for */

        /*
         * Allocate storage for the new chunk
         */
        H5_CHECK_OVERFLOW(udata->nbytes, size_t, hsize_t);
        if(HADDR_UNDEF == (*new_node_p = H5MF_alloc(f, H5FD_MEM_DRAW, dxpl_id, (hsize_t)udata->nbytes)))
            HGOTO_ERROR(H5E_STORAGE, H5E_NOSPACE, H5B_INS_ERROR, "file allocation failed")
        udata->addr = *new_node_p;
        ret_value = H5B_INS_RIGHT;

    } else {
        HGOTO_ERROR(H5E_IO, H5E_UNSUPPORTED, H5B_INS_ERROR, "internal error")
    }

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_insert() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_iter_allocated
 *
 * Purpose:	Simply counts the number of chunks for a dataset.
 *
 * Return:	Success:	Non-negative
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *              Wednesday, April 21, 1999
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static int
H5D_istore_iter_allocated (H5F_t UNUSED *f, hid_t UNUSED dxpl_id, const void *_lt_key, haddr_t UNUSED addr,
		    const void UNUSED *_rt_key, void *_udata)
{
    H5D_istore_it_ud1_t	*udata = (H5D_istore_it_ud1_t *)_udata;
    const H5D_istore_key_t	*lt_key = (const H5D_istore_key_t *)_lt_key;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_iter_allocated)

    udata->total_storage += lt_key->nbytes;

    FUNC_LEAVE_NOAPI(H5_ITER_CONT)
} /* H5D_istore_iter_allocated() */

/*-------------------------------------------------------------------------
 * Function:	H5D_istore_iter_chunkmap
 *
 * Purpose:	obtain chunk address and the corresponding index
 *
 * Return:	Success:	Non-negative
 *
 *		Failure:	Negative
 *
 * Programmer:	Kent Yang
 *              Tuesday, November 15, 2005
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static int
H5D_istore_iter_chunkmap (H5F_t UNUSED *f, hid_t UNUSED dxpl_id, const void *_lt_key, haddr_t addr,
		    const void UNUSED *_rt_key, void *_udata)
{
    H5D_istore_it_ud5_t	*udata = (H5D_istore_it_ud5_t *)_udata;
    const H5D_istore_key_t	*lt_key = (const H5D_istore_key_t *)_lt_key;
    unsigned       rank;
    hsize_t        chunk_index;
    int            ret_value = H5_ITER_CONT;     /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_iter_chunkmap)

    rank = udata->common.mesg->u.chunk.ndims - 1;

    if(H5V_chunk_index(rank, lt_key->offset, udata->common.mesg->u.chunk.dim, udata->down_chunks, &chunk_index) < 0)
       HGOTO_ERROR(H5E_DATASPACE, H5E_BADRANGE, FAIL, "can't get chunk index")

    udata->chunk_addr[chunk_index] = addr;

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* H5D_istore_iter_chunkmap() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_iter_dump
 *
 * Purpose:	If the UDATA.STREAM member is non-null then debugging
 *              information is written to that stream.
 *
 * Return:	Success:	Non-negative
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *              Wednesday, April 21, 1999
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static int
H5D_istore_iter_dump (H5F_t UNUSED *f, hid_t UNUSED dxpl_id, const void *_lt_key, haddr_t UNUSED addr,
		    const void UNUSED *_rt_key, void *_udata)
{
    H5D_istore_it_ud2_t	*udata = (H5D_istore_it_ud2_t *)_udata;
    const H5D_istore_key_t	*lt_key = (const H5D_istore_key_t *)_lt_key;
    unsigned		u;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_iter_dump)

    if(udata->stream) {
        if(!udata->header_displayed) {
            HDfprintf(udata->stream, "           Flags    Bytes     Address          Logical Offset\n");
            HDfprintf(udata->stream, "        ========== ======== ========== ==============================\n");

            /* Set flag that the headers has been printed */
            udata->header_displayed = TRUE;
        } /* end if */
        HDfprintf(udata->stream,     "        0x%08x %8Zu %10a [", lt_key->filter_mask, lt_key->nbytes, addr);
        for(u = 0; u < udata->common.mesg->u.chunk.ndims; u++)
            HDfprintf(udata->stream, "%s%Hd", (u ? ", " : ""), lt_key->offset[u]);
        HDfputs("]\n", udata->stream);
    } /* end if */

    FUNC_LEAVE_NOAPI(H5_ITER_CONT)
} /* H5D_istore_iter_dump() */


/*-------------------------------------------------------------------------
 * Function:    H5D_istore_iter_copy
 *
 * Purpose:     copy chunked raw data from source file and insert to the
 *              B-tree node in the destination file
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Peter Cao
 *              August 20, 2005
 *
 *-------------------------------------------------------------------------
 */
static int
H5D_istore_iter_copy(H5F_t *f_src, hid_t dxpl_id, const void *_lt_key,
    haddr_t addr_src, const void UNUSED *_rt_key, void *_udata)
{
    H5D_istore_it_ud4_t     *udata = (H5D_istore_it_ud4_t *)_udata;
    const H5D_istore_key_t  *lt_key = (const H5D_istore_key_t *)_lt_key;
    H5D_istore_ud1_t        udata_dst;                  /* User data about new destination chunk */
    hbool_t                 is_vlen = FALSE;
    hbool_t                 fix_ref = FALSE;

    /* General information about chunk copy */
    void                    *bkg = udata->bkg;
    void                    *buf = udata->buf;
    size_t                  buf_size = udata->buf_size;
    H5O_pline_t             *pline = udata->pline;

    /* needed for commpressed variable length data */
    hbool_t                 is_compressed = FALSE;
    H5Z_EDC_t               edc_read = H5Z_NO_EDC;
    size_t                  nbytes = lt_key->nbytes;
    H5Z_cb_t                cb_struct;

    int                     ret_value = H5_ITER_CONT; /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_iter_copy)

    /* Check parameter for type conversion */
    if(udata->do_convert) {
        if(H5T_detect_class(udata->dt_src, H5T_VLEN) > 0)
            is_vlen = TRUE;
        else if((H5T_get_class(udata->dt_src, FALSE) == H5T_REFERENCE) && (udata->file_src != udata->file_dst))
            fix_ref = TRUE;
        else
            HGOTO_ERROR(H5E_DATASET, H5E_CANTCOPY, FAIL, "unable to copy dataset elements")
    } /* end if */

    /* Check for filtered chunks */
    if(pline && pline->nused) {
        is_compressed = TRUE;
        cb_struct.func = NULL; /* no callback function when failed */
    } /* end if */

    /* Resize the buf if it is too small to hold the data */
    if(nbytes > buf_size) {
        void *new_buf;          /* New buffer for data */

        /* Re-allocate memory for copying the chunk */
        if(NULL == (new_buf = H5MM_realloc(udata->buf, nbytes)))
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, H5_ITER_ERROR, "memory allocation failed for raw data chunk")
        udata->buf = new_buf;
        if(udata->bkg) {
            if(NULL == (new_buf = H5MM_realloc(udata->bkg, nbytes)))
                HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, H5_ITER_ERROR, "memory allocation failed for raw data chunk")
            udata->bkg = new_buf;
            if(!udata->cpy_info->expand_ref)
                HDmemset((uint8_t *)udata->bkg + buf_size, 0, (size_t)(nbytes - buf_size));

            bkg = udata->bkg;
        } /* end if */

        buf = udata->buf;
        udata->buf_size = buf_size = nbytes;
    } /* end if */

    /* read chunk data from the source file */
    if(H5F_block_read(f_src, H5FD_MEM_DRAW, addr_src, nbytes, dxpl_id, buf) < 0)
        HGOTO_ERROR(H5E_IO, H5E_READERROR, H5_ITER_ERROR, "unable to read raw data chunk")

    /* Need to uncompress variable-length & reference data elements */
    if(is_compressed && (is_vlen || fix_ref)) {
        unsigned filter_mask = lt_key->filter_mask;

        if(H5Z_pipeline(pline, H5Z_FLAG_REVERSE, &filter_mask, edc_read, cb_struct, &nbytes, &buf_size, &buf) < 0)
            HGOTO_ERROR(H5E_PLINE, H5E_CANTFILTER, H5_ITER_ERROR, "data pipeline read failed")
    } /* end if */

    /* Perform datatype conversion, if necessary */
    if(is_vlen) {
        H5T_path_t              *tpath_src_mem = udata->tpath_src_mem;
        H5T_path_t              *tpath_mem_dst = udata->tpath_mem_dst;
        H5S_t                   *buf_space = udata->buf_space;
        hid_t                   tid_src = udata->tid_src;
        hid_t                   tid_dst = udata->tid_dst;
        hid_t                   tid_mem = udata->tid_mem;
        size_t                  nelmts = udata->nelmts;
        void                    *reclaim_buf = udata->reclaim_buf;
        size_t                  reclaim_buf_size = udata->reclaim_buf_size;

        /* Convert from source file to memory */
        if(H5T_convert(tpath_src_mem, tid_src, tid_mem, nelmts, (size_t)0, (size_t)0, buf, NULL, dxpl_id) < 0)
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, H5_ITER_ERROR, "datatype conversion failed")

        /* Copy into another buffer, to reclaim memory later */
        HDmemcpy(reclaim_buf, buf, reclaim_buf_size);

        /* Set background buffer to all zeros */
        HDmemset(bkg, 0, buf_size);

        /* Convert from memory to destination file */
        if(H5T_convert(tpath_mem_dst, tid_mem, tid_dst, nelmts, (size_t)0, (size_t)0, buf, bkg, dxpl_id) < 0)
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, H5_ITER_ERROR, "datatype conversion failed")

        /* Reclaim space from variable length data */
        if(H5D_vlen_reclaim(tid_mem, buf_space, H5P_DATASET_XFER_DEFAULT, reclaim_buf) < 0)
            HGOTO_ERROR(H5E_DATASET, H5E_BADITER, H5_ITER_ERROR, "unable to reclaim variable-length data")
    } /* end if */
    else if(fix_ref) {
        /* Check for expanding references */
        /* (background buffer has already been zeroed out, if not expanding) */
        if(udata->cpy_info->expand_ref) {
            size_t ref_count;

            /* Determine # of reference elements to copy */
            ref_count = nbytes / H5T_get_size(udata->dt_src);

            /* Copy the reference elements */
            if(H5O_copy_expand_ref(f_src, buf, dxpl_id, udata->file_dst, bkg, ref_count, H5T_get_ref_type(udata->dt_src), udata->cpy_info) < 0)
                HGOTO_ERROR(H5E_DATASET, H5E_CANTCOPY, FAIL, "unable to copy reference attribute")
        } /* end if */

        /* After fix ref, copy the new reference elements to the buffer to write out */
        HDmemcpy(buf, bkg, buf_size);
    } /* end if */

    /* Set up destination chunk callback information for insertion */
    udata_dst.common.mesg = udata->common.mesg;     /* Share this pointer for a short while */
    udata_dst.common.offset = lt_key->offset;
    udata_dst.nbytes = lt_key->nbytes;
    udata_dst.filter_mask = lt_key->filter_mask;
    udata_dst.addr = HADDR_UNDEF;

    /* Need to compress variable-length & reference data elements before writing to file */
    if(is_compressed && (is_vlen || fix_ref) ) {
        if(H5Z_pipeline(pline, 0, &(udata_dst.filter_mask), edc_read,
                cb_struct, &nbytes, &buf_size, &buf) < 0)
            HGOTO_ERROR(H5E_PLINE, H5E_CANTFILTER, H5_ITER_ERROR, "output pipeline failed")
        udata_dst.nbytes = nbytes;
	udata->buf = buf;
	udata->buf_size = buf_size;
    } /* end if */

    /* Insert chunk into the destination Btree */
    if(H5B_insert(udata->file_dst, dxpl_id, H5B_ISTORE, udata->addr_dst, &udata_dst) < 0)
        HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, H5_ITER_ERROR, "unable to allocate chunk")

    /* Write chunk data to destination file */
    HDassert(H5F_addr_defined(udata_dst.addr));
    if(H5F_block_write(udata->file_dst, H5FD_MEM_DRAW, udata_dst.addr, nbytes, dxpl_id, buf) < 0)
        HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, H5_ITER_ERROR, "unable to write raw data to file")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_iter_copy() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_cinfo_cache_reset
 *
 * Purpose:	Reset the cached chunk info
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              November 27, 2007
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_cinfo_cache_reset(H5D_chunk_cached_t *last)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_cinfo_cache_reset)

    /* Sanity check */
    HDassert(last);

    /* Indicate that the cached info is not valid */
    last->valid = FALSE;

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* H5D_istore_cinfo_cache_reset() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_cinfo_cache_update
 *
 * Purpose:	Update the cached chunk info
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              November 27, 2007
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_cinfo_cache_update(H5D_chunk_cached_t *last, const H5D_istore_ud1_t *udata)
{
    unsigned    u;                              /* Local index variable */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_cinfo_cache_update)

    /* Sanity check */
    HDassert(last);
    HDassert(udata);
    HDassert(udata->common.mesg);
    HDassert(udata->common.offset);

    /* Stored the information to cache */
    for(u = 0; u < udata->common.mesg->u.chunk.ndims; u++)
        last->offset[u] = udata->common.offset[u];
    last->nbytes = udata->nbytes;
    last->filter_mask = udata->filter_mask;
    last->addr = udata->addr;

    /* Indicate that the cached info is valid */
    last->valid = TRUE;

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* H5D_istore_cinfo_cache_update() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_cinfo_cache_found
 *
 * Purpose:	Look for chunk info in cache
 *
 * Return:	TRUE/FALSE/FAIL
 *
 * Programmer:	Quincey Koziol
 *              November 27, 2007
 *
 *-------------------------------------------------------------------------
 */
static hbool_t
H5D_istore_cinfo_cache_found(const H5D_chunk_cached_t *last, H5D_istore_ud1_t *udata)
{
    hbool_t ret_value = FALSE;          /* Return value */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_cinfo_cache_found)

    /* Sanity check */
    HDassert(last);
    HDassert(udata);
    HDassert(udata->common.mesg);
    HDassert(udata->common.offset);

    /* Check if the cached information is what is desired */
    if(last->valid) {
        unsigned    u;                      /* Local index variable */

        /* Check that the offset is the same */
        for(u = 0; u < udata->common.mesg->u.chunk.ndims; u++)
            if(last->offset[u] != udata->common.offset[u])
                HGOTO_DONE(FALSE)

        /* Retrieve the information from the cache */
        udata->nbytes = last->nbytes;
        udata->filter_mask = last->filter_mask;
        udata->addr = last->addr;

        /* Indicate that the data was found */
        HGOTO_DONE(TRUE)
    } /* end if */

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* H5D_istore_cinfo_cache_found() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_init
 *
 * Purpose:	Initialize the raw data chunk cache for a dataset.  This is
 *		called when the dataset is initialized.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Monday, May 18, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_init(const H5F_t *f, const H5D_t *dset)
{
    H5D_rdcc_t	*rdcc = &(dset->shared->cache.chunk);
    herr_t      ret_value = SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_init, FAIL)

    if(H5F_RDCC_NBYTES(f) > 0 && H5F_RDCC_NELMTS(f) > 0) {
        rdcc->nbytes = H5F_RDCC_NBYTES(f);
	rdcc->nslots = H5F_RDCC_NELMTS(f);
	rdcc->slot = H5FL_SEQ_CALLOC (H5D_rdcc_ent_ptr_t,rdcc->nslots);
	if(NULL==rdcc->slot)
	    HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed")

        /* Reset any cached chunk info for this dataset */
        H5D_istore_cinfo_cache_reset(&(rdcc->last));
    } /* end if */

    /* Allocate the shared structure */
    if(H5D_istore_shared_create(f, &dset->shared->layout) < 0)
	HGOTO_ERROR(H5E_RESOURCE, H5E_CANTINIT, FAIL, "can't create wrapper for shared B-tree info")
done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_init() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_flush_entry
 *
 * Purpose:	Writes a chunk to disk.  If RESET is non-zero then the
 *		entry is cleared -- it's slightly faster to flush a chunk if
 *		the RESET flag is turned on because it results in one fewer
 *		memory copy.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, May 21, 1998
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_flush_entry(const H5D_io_info_t *io_info, H5D_rdcc_ent_t *ent, hbool_t reset)
{
    void		*buf = NULL;	/*temporary buffer		*/
    size_t		alloc;		/*bytes allocated for BUF	*/
    hbool_t		point_of_no_return = FALSE;
    herr_t		ret_value = SUCCEED;	/*return value			*/

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_flush_entry)

    assert(io_info);
    assert(io_info->dset);
    assert(ent);
    assert(!ent->locked);

    buf = ent->chunk;
    if(ent->dirty) {
        H5D_istore_ud1_t 	udata;		/*pass through B-tree		*/

        /* Initial buffer size */
        alloc = ent->alloc_size;

        /* Set up user data for B-tree callbacks */
        udata.common.mesg = &io_info->dset->shared->layout;
        udata.common.offset = ent->offset;
        udata.filter_mask = 0;
        udata.nbytes = ent->chunk_size;
        udata.addr = HADDR_UNDEF;

        /* Should the chunk be filtered before writing it to disk? */
        if(io_info->dset->shared->dcpl_cache.pline.nused) {
            if(!reset) {
                /*
                 * Copy the chunk to a new buffer before running it through
                 * the pipeline because we'll want to save the original buffer
                 * for later.
                 */
                alloc = ent->chunk_size;
                if(NULL == (buf = H5MM_malloc(alloc)))
                    HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed for pipeline")
                HDmemcpy(buf, ent->chunk, ent->chunk_size);
            } /* end if */
            else {
                /*
                 * If we are reseting and something goes wrong after this
                 * point then it's too late to recover because we may have
                 * destroyed the original data by calling H5Z_pipeline().
                 * The only safe option is to continue with the reset
                 * even if we can't write the data to disk.
                 */
                point_of_no_return = TRUE;
                ent->chunk = NULL;
            } /* end else */
            if(H5Z_pipeline(&(io_info->dset->shared->dcpl_cache.pline), 0, &(udata.filter_mask), io_info->dxpl_cache->err_detect,
                     io_info->dxpl_cache->filter_cb, &(udata.nbytes), &alloc, &buf) < 0)
                HGOTO_ERROR(H5E_PLINE, H5E_CANTFILTER, FAIL, "output pipeline failed")
        } /* end if */

        /*
         * Create the chunk it if it doesn't exist, or reallocate the chunk if
         * its size changed.  Then write the data into the file.
         */
        if(H5B_insert(io_info->dset->oloc.file, io_info->dxpl_id, H5B_ISTORE, io_info->dset->shared->layout.u.chunk.addr, &udata)<0)
            HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to allocate chunk")
        if(H5F_block_write(io_info->dset->oloc.file, H5FD_MEM_DRAW, udata.addr, udata.nbytes, io_info->dxpl_id, buf) < 0)
            HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to write raw data to file")

        /* Cache the chunk's info, in case it's accessed again shortly */
        H5D_istore_cinfo_cache_update(&io_info->dset->shared->cache.chunk.last, &udata);

        /* Mark cache entry as clean */
        ent->dirty = FALSE;
#ifdef H5D_ISTORE_DEBUG
        io_info->dset->shared->cache.chunk.nflushes++;
#endif /* H5D_ISTORE_DEBUG */
    } /* end if */

    /* Reset, but do not free or removed from list */
    if(reset) {
        point_of_no_return = FALSE;
        if(buf == ent->chunk)
            buf = NULL;
        if(ent->chunk != NULL)
            ent->chunk = (uint8_t *)H5D_istore_chunk_xfree(ent->chunk, &(io_info->dset->shared->dcpl_cache.pline));
    } /* end if */

done:
    /* Free the temp buffer only if it's different than the entry chunk */
    if(buf != ent->chunk)
        H5MM_xfree(buf);

    /*
     * If we reached the point of no return then we have no choice but to
     * reset the entry.  This can only happen if RESET is true but the
     * output pipeline failed.  Do not free the entry or remove it from the
     * list.
     */
    if(ret_value < 0 && point_of_no_return) {
        if(ent->chunk)
            ent->chunk = (uint8_t *)H5D_istore_chunk_xfree(ent->chunk, &(io_info->dset->shared->dcpl_cache.pline));
    } /* end if */

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_flush_entry() */


/*-------------------------------------------------------------------------
 * Function:    H5D_istore_preempt
 *
 * Purpose:     Preempts the specified entry from the cache, flushing it to
 *              disk if necessary.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *              Thursday, May 21, 1998
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_preempt(const H5D_io_info_t *io_info, H5D_rdcc_ent_t * ent, hbool_t flush)
{
    H5D_rdcc_t *rdcc = &(io_info->dset->shared->cache.chunk);
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_preempt)

    assert(io_info);
    assert(ent);
    assert(!ent->locked);
    assert(ent->idx < rdcc->nslots);

    if(flush) {
	/* Flush */
	if(H5D_istore_flush_entry(io_info, ent, TRUE) < 0)
	    HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "cannot flush indexed storage buffer")
    } /* end if */
    else {
        /* Don't flush, just free chunk */
	if(ent->chunk != NULL)
	    ent->chunk = (uint8_t *)H5D_istore_chunk_xfree(ent->chunk, &(io_info->dset->shared->dcpl_cache.pline));
    } /* end else */

    /* Unlink from list */
    if(ent->prev)
	ent->prev->next = ent->next;
    else
	rdcc->head = ent->next;
    if(ent->next)
	ent->next->prev = ent->prev;
    else
	rdcc->tail = ent->prev;
    ent->prev = ent->next = NULL;

    /* Remove from cache */
    rdcc->slot[ent->idx] = NULL;
    ent->idx = UINT_MAX;
    rdcc->nbytes -= ent->chunk_size;
    --rdcc->nused;

    /* Free */
    H5FL_FREE(H5D_rdcc_ent_t, ent);

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_preempt() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_flush
 *
 * Purpose:	Writes all dirty chunks to disk and optionally preempts them
 *		from the cache.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, May 21, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_flush(H5D_t *dset, hid_t dxpl_id, unsigned flags)
{
    H5D_io_info_t io_info;              /* Temporary I/O info object */
    H5D_dxpl_cache_t _dxpl_cache;       /* Data transfer property cache buffer */
    H5D_dxpl_cache_t *dxpl_cache = &_dxpl_cache;   /* Data transfer property cache */
    H5D_rdcc_t *rdcc = &(dset->shared->cache.chunk);
    unsigned		nerrors = 0;
    H5D_rdcc_ent_t	*ent, *next;
    herr_t ret_value = SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_flush, FAIL)

    /* Fill the DXPL cache values for later use */
    if(H5D_get_dxpl_cache(dxpl_id, &dxpl_cache) < 0)
        HGOTO_ERROR(H5E_DATASET, H5E_CANTGET, FAIL, "can't fill dxpl cache")

    /* Construct dataset I/O info */
    H5D_BUILD_IO_INFO(&io_info, dset, dxpl_cache, dxpl_id, NULL);

    /* Loop over all entries in the chunk cache */
    for(ent = rdcc->head; ent; ent = next) {
	next = ent->next;
	if((flags & H5F_FLUSH_INVALIDATE)) {
	    if(H5D_istore_preempt(&io_info, ent, TRUE) < 0)
		nerrors++;
	} else {
	    if(H5D_istore_flush_entry(&io_info, ent, FALSE) < 0)
		nerrors++;
	}
    } /* end for */
    if(nerrors)
	HGOTO_ERROR(H5E_IO, H5E_CANTFLUSH, FAIL, "unable to flush one or more raw data chunks")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_flush() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_dest
 *
 * Purpose:	Destroy the entire chunk cache by flushing dirty entries,
 *		preempting all entries, and freeing the cache itself.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, May 21, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_dest (H5D_t *dset, hid_t dxpl_id)
{
    H5D_io_info_t io_info;              /* Temporary I/O info object */
    H5D_dxpl_cache_t _dxpl_cache;       /* Data transfer property cache buffer */
    H5D_dxpl_cache_t *dxpl_cache=&_dxpl_cache;   /* Data transfer property cache */
    H5D_rdcc_t		*rdcc = &(dset->shared->cache.chunk);
    int		nerrors=0;
    H5D_rdcc_ent_t	*ent=NULL, *next=NULL;
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_dest, FAIL)

    assert(dset);

    /* Fill the DXPL cache values for later use */
    if (H5D_get_dxpl_cache(dxpl_id,&dxpl_cache)<0)
        HGOTO_ERROR(H5E_DATASET, H5E_CANTGET, FAIL, "can't fill dxpl cache")

    /* Construct dataset I/O info */
    H5D_BUILD_IO_INFO(&io_info,dset,dxpl_cache,dxpl_id,NULL);

    /* Flush all the cached chunks */
    for (ent=rdcc->head; ent; ent=next) {
#ifdef H5D_ISTORE_DEBUG
	HDfputc('c', stderr);
	HDfflush(stderr);
#endif
	next = ent->next;
	if (H5D_istore_preempt(&io_info, ent, TRUE )<0)
	    nerrors++;
    }
    if (nerrors)
	HGOTO_ERROR(H5E_IO, H5E_CANTFLUSH, FAIL, "unable to flush one or more raw data chunks")

    if(rdcc->slot)
        H5FL_SEQ_FREE (H5D_rdcc_ent_ptr_t,rdcc->slot);
    HDmemset (rdcc, 0, sizeof(H5D_rdcc_t));

    /* Free the raw B-tree node buffer */
    if(dset->shared->layout.u.chunk.btree_shared==NULL)
        HGOTO_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "ref-counted page nil")
    if(H5RC_DEC(dset->shared->layout.u.chunk.btree_shared)<0)
	HGOTO_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "unable to decrement ref-counted page")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_dest() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_shared_create
 *
 * Purpose:	Create & initialize B-tree shared info
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Monday, September 27, 2004
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_shared_create (const H5F_t *f, H5O_layout_t *layout)
{
    H5B_shared_t *shared;               /* Shared B-tree node info */
    size_t	u;                      /* Local index variable */
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_shared_create)

    /* Allocate space for the shared structure */
    if(NULL==(shared=H5FL_MALLOC(H5B_shared_t)))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed for shared B-tree info")

    /* Set up the "global" information for this file's groups */
    shared->type= H5B_ISTORE;
    shared->two_k=2*H5F_KVALUE(f,H5B_ISTORE);
    shared->sizeof_rkey = 4 +			/*storage size		*/
                         4 +			/*filter mask		*/
                         layout->u.chunk.ndims*8;	/*dimension indices	*/
    assert(shared->sizeof_rkey);
    shared->sizeof_rnode = H5B_nodesize(f, shared, &shared->sizeof_keys);
    assert(shared->sizeof_rnode);
    if(NULL==(shared->page=H5FL_BLK_MALLOC(chunk_page,shared->sizeof_rnode)))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed for B-tree page")
#ifdef H5_CLEAR_MEMORY
HDmemset(shared->page, 0, shared->sizeof_rnode);
#endif /* H5_CLEAR_MEMORY */
    if(NULL==(shared->nkey=H5FL_SEQ_MALLOC(size_t,(size_t)(2*H5F_KVALUE(f,H5B_ISTORE)+1))))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed for B-tree page")

    /* Initialize the offsets into the native key buffer */
    for(u=0; u<(2*H5F_KVALUE(f,H5B_ISTORE)+1); u++)
        shared->nkey[u]=u*H5B_ISTORE[0].sizeof_nkey;

    /* Make shared B-tree info reference counted */
    if(NULL==(layout->u.chunk.btree_shared=H5RC_create(shared,H5D_istore_shared_free)))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "can't create ref-count wrapper for shared B-tree info")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_shared_create() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_shared_free
 *
 * Purpose:	Free B-tree shared info
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July  8, 2004
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_shared_free (void *_shared)
{
    H5B_shared_t *shared = (H5B_shared_t *)_shared;

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_shared_free)

    /* Free the raw B-tree node buffer */
    (void)H5FL_BLK_FREE(chunk_page, shared->page);

    /* Free the B-tree native key offsets buffer */
    H5FL_SEQ_FREE(size_t, shared->nkey);

    /* Free the shared B-tree info */
    H5FL_FREE(H5B_shared_t, shared);

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* end H5D_istore_shared_free() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_prune
 *
 * Purpose:	Prune the cache by preempting some things until the cache has
 *		room for something which is SIZE bytes.  Only unlocked
 *		entries are considered for preemption.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, May 21, 1998
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_prune (const H5D_io_info_t *io_info, size_t size)
{
    int		i, j, nerrors=0;
    const H5D_rdcc_t	*rdcc = &(io_info->dset->shared->cache.chunk);
    size_t		total = rdcc->nbytes;
    const int		nmeth=2;	/*number of methods		*/
    int		        w[1];		/*weighting as an interval	*/
    H5D_rdcc_ent_t	*p[2], *cur;	/*list pointers			*/
    H5D_rdcc_ent_t	*n[2];		/*list next pointers		*/
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_prune)

    /*
     * Preemption is accomplished by having multiple pointers (currently two)
     * slide down the list beginning at the head. Pointer p(N+1) will start
     * traversing the list when pointer pN reaches wN percent of the original
     * list.  In other words, preemption method N gets to consider entries in
     * approximate least recently used order w0 percent before method N+1
     * where 100% means tha method N will run to completion before method N+1
     * begins.  The pointers participating in the list traversal are each
     * given a chance at preemption before any of the pointers are advanced.
     */
    w[0] = (int)(rdcc->nused * H5F_RDCC_W0(io_info->dset->oloc.file));
    p[0] = rdcc->head;
    p[1] = NULL;

    while ((p[0] || p[1]) && rdcc->nbytes+size>total) {

	/* Introduce new pointers */
	for (i=0; i<nmeth-1; i++)
            if (0==w[i])
                p[i+1] = rdcc->head;

	/* Compute next value for each pointer */
	for (i=0; i<nmeth; i++)
            n[i] = p[i] ? p[i]->next : NULL;

	/* Give each method a chance */
	for (i=0; i<nmeth && rdcc->nbytes+size>total; i++) {
	    if (0==i && p[0] && !p[0]->locked &&
                    ((0==p[0]->rd_count && 0==p[0]->wr_count) ||
                     (0==p[0]->rd_count && p[0]->chunk_size==p[0]->wr_count) ||
                     (p[0]->chunk_size==p[0]->rd_count && 0==p[0]->wr_count))) {
		/*
		 * Method 0: Preempt entries that have been completely written
		 * and/or completely read but not entries that are partially
		 * written or partially read.
		 */
		cur = p[0];
#ifdef H5D_ISTORE_DEBUG
		HDputc('.', stderr);
		HDfflush(stderr);
#endif

	    } else if (1==i && p[1] && !p[1]->locked) {
		/*
		 * Method 1: Preempt the entry without regard to
		 * considerations other than being locked.  This is the last
		 * resort preemption.
		 */
		cur = p[1];
#ifdef H5D_ISTORE_DEBUG
		HDputc(':', stderr);
		HDfflush(stderr);
#endif

	    } else {
		/* Nothing to preempt at this point */
		cur= NULL;
	    }

	    if (cur) {
		for (j=0; j<nmeth; j++) {
		    if (p[j]==cur)
                        p[j] = NULL;
		    if (n[j]==cur)
                        n[j] = cur->next;
		}
		if (H5D_istore_preempt(io_info, cur, TRUE)<0)
                    nerrors++;
	    }
	}

	/* Advance pointers */
	for (i=0; i<nmeth; i++)
            p[i] = n[i];
	for (i=0; i<nmeth-1; i++)
            w[i] -= 1;
    }

    if (nerrors)
	HGOTO_ERROR(H5E_IO, H5E_CANTFLUSH, FAIL, "unable to preempt one or more raw data cache entry")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_prune() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_lock
 *
 * Purpose:	Return a pointer to a dataset chunk.  The pointer points
 *		directly into the chunk cache and should not be freed
 *		by the caller but will be valid until it is unlocked.  The
 *		input value IDX_HINT is used to speed up cache lookups and
 *		it's output value should be given to H5F_istore_unlock().
 *		IDX_HINT is ignored if it is out of range, and if it points
 *		to the wrong entry then we fall back to the normal search
 *		method.
 *
 *		If RELAX is non-zero and the chunk isn't in the cache then
 *		don't try to read it from the file, but just allocate an
 *		uninitialized buffer to hold the result.  This is intended
 *		for output functions that are about to overwrite the entire
 *		chunk.
 *
 * Return:	Success:	Ptr to a file chunk.
 *
 *		Failure:	NULL
 *
 * Programmer:	Robb Matzke
 *              Thursday, May 21, 1998
 *
 *-------------------------------------------------------------------------
 */
void *
H5D_istore_lock(const H5D_io_info_t *io_info, H5D_istore_ud1_t *udata,
    hbool_t relax, unsigned *idx_hint/*in,out*/)
{
    H5D_t *dset = io_info->dset;                /* Local pointer to the dataset info */
    const H5O_pline_t  *pline = &(dset->shared->dcpl_cache.pline);    /* I/O pipeline info */
    const H5O_layout_t *layout = &(dset->shared->layout);       /* Dataset layout */
    const H5O_fill_t    *fill = &(dset->shared->dcpl_cache.fill);    /* Fill value info */
    H5D_fill_buf_info_t fb_info;                /* Dataset's fill buffer info */
    hbool_t             fb_info_init = FALSE;   /* Whether the fill value buffer has been initialized */
    H5D_rdcc_t		*rdcc = &(dset->shared->cache.chunk);   /*raw data chunk cache*/
    H5D_rdcc_ent_t	*ent = NULL;		/*cache entry		*/
    unsigned		idx = 0;		/*hash index number	*/
    hbool_t		found = FALSE;		/*already in cache?	*/
    size_t		chunk_size;		/*size of a chunk	*/
    void		*chunk = NULL;		/*the file chunk	*/
    unsigned		u;			/*counters		*/
    void		*ret_value;	        /*return value		*/

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_lock)

    HDassert(io_info);
    HDassert(dset);
    HDassert(io_info->dxpl_cache);
    HDassert(io_info->store);
    HDassert(TRUE == H5P_isa_class(io_info->dxpl_id, H5P_DATASET_XFER));

    /* Get the chunk's size */
    HDassert(layout->u.chunk.size > 0);
    H5_ASSIGN_OVERFLOW(chunk_size, layout->u.chunk.size, hsize_t, size_t);

    /* Search for the chunk in the cache */
    if(rdcc->nslots > 0) {
        idx = H5D_HASH(dset->shared,io_info->store->chunk.index);
        ent = rdcc->slot[idx];

        if(ent)
            for(u = 0, found = TRUE; u < layout->u.chunk.ndims; u++)
                if(io_info->store->chunk.offset[u] != ent->offset[u]) {
                    found = FALSE;
                    break;
                } /* end if */
    } /* end if */

    if(found) {
        /*
         * Already in the cache.  Count a hit.
         */
#ifdef H5D_ISTORE_DEBUG
        rdcc->nhits++;
#endif /* H5D_ISTORE_DEBUG */
    } /* end if */
    else if(relax) {
        /*
         * Not in the cache, but we're about to overwrite the whole thing
         * anyway, so just allocate a buffer for it but don't initialize that
         * buffer with the file contents. Count this as a hit instead of a
         * miss because we saved ourselves lots of work.
         */
#ifdef H5D_ISTORE_DEBUG
        HDputc('w', stderr);
        HDfflush(stderr);
        rdcc->nhits++;
#endif
        if(NULL == (chunk = H5D_istore_chunk_alloc(chunk_size, pline)))
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed for raw data chunk")
	
        /* In the case that some dataset functions look through this data,
         * clear it to all 0s. */
        HDmemset(chunk, 0, chunk_size);

    } /* end if */
    else {
        H5D_istore_ud1_t tmp_udata;		/*B-tree pass-through	*/
        haddr_t chunk_addr;             /* Address of chunk on disk */

        if(udata!=NULL)
            chunk_addr = udata->addr;
        else {
            /* Point at temporary storage for B-tree pass through */
            udata = &tmp_udata;

            /*
             * Not in the cache.  Read it from the file and count this as a miss
             * if it's in the file or an init if it isn't.
             */
            chunk_addr = H5D_istore_get_addr(io_info, udata);
        } /* end else */

        if (H5F_addr_defined(chunk_addr)) {
            size_t		chunk_alloc = 0;		/*allocated chunk size	*/

            /*
             * The chunk exists on disk.
             */
            /* Chunk size on disk isn't [likely] the same size as the final chunk
             * size in memory, so allocate memory big enough. */
            chunk_alloc = udata->nbytes;
            if(NULL == (chunk = H5D_istore_chunk_alloc (chunk_alloc, pline)))
                HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed for raw data chunk")
            if(H5F_block_read(dset->oloc.file, H5FD_MEM_DRAW, chunk_addr, udata->nbytes, io_info->dxpl_id, chunk) < 0)
                HGOTO_ERROR(H5E_IO, H5E_READERROR, NULL, "unable to read raw data chunk")

            if(pline->nused)
                if(H5Z_pipeline(pline, H5Z_FLAG_REVERSE, &(udata->filter_mask), io_info->dxpl_cache->err_detect,
                        io_info->dxpl_cache->filter_cb, &(udata->nbytes), &chunk_alloc, &chunk) < 0)
                    HGOTO_ERROR(H5E_PLINE, H5E_CANTFILTER, NULL, "data pipeline read failed")
#ifdef H5D_ISTORE_DEBUG
            rdcc->nmisses++;
#endif /* H5D_ISTORE_DEBUG */
        } else {
            H5D_fill_value_t	fill_status;

#ifdef OLD_WAY
            /* Clear the error stack from not finding the chunk on disk */
            H5E_clear_stack(NULL);
#endif /* OLD_WAY */

            /* Chunk size on disk isn't [likely] the same size as the final chunk
             * size in memory, so allocate memory big enough. */
            if(NULL == (chunk = H5D_istore_chunk_alloc (chunk_size, pline)))
                HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed for raw data chunk")

            if(H5P_is_fill_value_defined(fill, &fill_status) < 0)
                HGOTO_ERROR(H5E_PLIST, H5E_CANTGET, NULL, "can't tell if fill value defined")

            if(fill->fill_time == H5D_FILL_TIME_ALLOC ||
                    (fill->fill_time == H5D_FILL_TIME_IFSET && fill_status == H5D_FILL_VALUE_USER_DEFINED)) {
                /*
                 * The chunk doesn't exist in the file.  Replicate the fill
                 * value throughout the chunk, if the fill value is defined.
                 */

                /* Initialize the fill value buffer */
                /* (use the compact dataset storage buffer as the fill value buffer) */
                if(H5D_fill_init(&fb_info, chunk, FALSE,
                        NULL, NULL, NULL, NULL,
                        &dset->shared->dcpl_cache.fill, dset->shared->type,
                        dset->shared->type_id, (size_t)0, chunk_size, io_info->dxpl_id) < 0)
                    HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, NULL, "can't initialize fill buffer info")
                fb_info_init = TRUE;

                /* Check for VL datatype & non-default fill value */
                if(fb_info.has_vlen_fill_type)
                    /* Fill the buffer with VL datatype fill values */
                    if(H5D_fill_refill_vl(&fb_info, fb_info.elmts_per_buf, io_info->dxpl_id) < 0)
                        HGOTO_ERROR(H5E_DATASET, H5E_CANTCONVERT, NULL, "can't refill fill value buffer")
            } /* end if */
#ifdef H5_CLEAR_MEMORY
            else
                HDmemset(chunk, 0, chunk_size);
#endif /* H5_CLEAR_MEMORY */
#ifdef H5D_ISTORE_DEBUG
            rdcc->ninits++;
#endif /* H5D_ISTORE_DEBUG */
        } /* end else */
    } /* end else */
    HDassert(found || chunk_size > 0);

    if(!found && rdcc->nslots > 0 && chunk_size <= rdcc->nbytes &&
            (!ent || !ent->locked)) {
        /*
         * Add the chunk to the cache only if the slot is not already locked.
         * Preempt enough things from the cache to make room.
         */
        if (ent) {
#ifdef H5D_ISTORE_DEBUG
            HDputc('#', stderr);
            HDfflush(stderr);
#endif
            if (H5D_istore_preempt(io_info, ent, TRUE)<0)
                HGOTO_ERROR(H5E_IO, H5E_CANTINIT, NULL, "unable to preempt chunk from cache")
        }
        if (H5D_istore_prune(io_info, chunk_size)<0)
            HGOTO_ERROR(H5E_IO, H5E_CANTINIT, NULL, "unable to preempt chunk(s) from cache")

        /* Create a new entry */
        ent = H5FL_MALLOC(H5D_rdcc_ent_t);
        ent->locked = 0;
        ent->dirty = FALSE;
        ent->chunk_size = chunk_size;
        ent->alloc_size = chunk_size;
        for (u=0; u<layout->u.chunk.ndims; u++)
            ent->offset[u] = io_info->store->chunk.offset[u];
        ent->rd_count = chunk_size;
        ent->wr_count = chunk_size;
        ent->chunk = (uint8_t*)chunk;

        /* Add it to the cache */
        assert(NULL==rdcc->slot[idx]);
        rdcc->slot[idx] = ent;
        ent->idx = idx;
        rdcc->nbytes += chunk_size;
        rdcc->nused++;

        /* Add it to the linked list */
        ent->next = NULL;
        if (rdcc->tail) {
            rdcc->tail->next = ent;
            ent->prev = rdcc->tail;
            rdcc->tail = ent;
        } else {
            rdcc->head = rdcc->tail = ent;
            ent->prev = NULL;
        }
        found = TRUE;
    } else if (!found) {
        /*
         * The chunk is larger than the entire cache so we don't cache it.
         * This is the reason all those arguments have to be repeated for the
         * unlock function.
         */
        ent = NULL;
        idx = UINT_MAX;

    } else {
        /*
         * The chunk is not at the beginning of the cache; move it backward
         * by one slot.  This is how we implement the LRU preemption
         * algorithm.
         */
        assert(ent);
        if (ent->next) {
            if (ent->next->next)
                ent->next->next->prev = ent;
            else
                rdcc->tail = ent;
            ent->next->prev = ent->prev;
            if (ent->prev)
                ent->prev->next = ent->next;
            else
                rdcc->head = ent->next;
            ent->prev = ent->next;
            ent->next = ent->next->next;
            ent->prev->next = ent;
        }
    }

    /* Lock the chunk into the cache */
    if (ent) {
        assert (!ent->locked);
        ent->locked = TRUE;
        chunk = ent->chunk;
    }

    if (idx_hint)
        *idx_hint = idx;

    /* Set return value */
    ret_value = chunk;

done:
    /* Release the fill buffer info, if it's been initialized */
    if(fb_info_init && H5D_fill_term(&fb_info) < 0)
        HDONE_ERROR(H5E_DATASET, H5E_CANTFREE, NULL, "Can't release fill buffer info")

    /* Release the chunk allocated, on error */
    if(!ret_value)
        if(chunk)
            chunk = H5D_istore_chunk_xfree(chunk, pline);

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_lock() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_unlock
 *
 * Purpose:	Unlocks a previously locked chunk. The LAYOUT, COMP, and
 *		OFFSET arguments should be the same as for H5F_rdcc_lock().
 *		The DIRTY argument should be set to non-zero if the chunk has
 *		been modified since it was locked. The IDX_HINT argument is
 *		the returned index hint from the lock operation and BUF is
 *		the return value from the lock.
 *
 *		The NACCESSED argument should be the number of bytes accessed
 *		for reading or writing (depending on the value of DIRTY).
 *		It's only purpose is to provide additional information to the
 *		preemption policy.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, May 21, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_unlock(const H5D_io_info_t *io_info,
    hbool_t dirty, unsigned idx_hint, void *chunk, size_t naccessed)
{
    const H5O_layout_t *layout=&(io_info->dset->shared->layout); /* Dataset layout */
    const H5D_rdcc_t	*rdcc = &(io_info->dset->shared->cache.chunk);
    H5D_rdcc_ent_t	*ent = NULL;
    unsigned		u;
    herr_t              ret_value=SUCCEED;      /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_unlock)

    assert(io_info);

    if (UINT_MAX==idx_hint) {
        /*
         * It's not in the cache, probably because it's too big.  If it's
         * dirty then flush it to disk.  In any case, free the chunk.
         * Note: we have to copy the layout and filter messages so we
         *	 don't discard the `const' qualifier.
         */
        if (dirty) {
            H5D_rdcc_ent_t x;

            HDmemset (&x, 0, sizeof x);
            x.dirty = TRUE;
            for (u=0; u<layout->u.chunk.ndims; u++)
                x.offset[u] = io_info->store->chunk.offset[u];
            assert(layout->u.chunk.size>0);
            H5_ASSIGN_OVERFLOW(x.chunk_size,layout->u.chunk.size,hsize_t,size_t);
            x.alloc_size = x.chunk_size;
            x.chunk = (uint8_t*)chunk;

            if (H5D_istore_flush_entry(io_info, &x, TRUE)<0)
                HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "cannot flush indexed storage buffer")
        } else {
            if(chunk)
                chunk=H5D_istore_chunk_xfree (chunk,&(io_info->dset->shared->dcpl_cache.pline));
        }
    } else {
        /* Sanity check */
	assert(idx_hint<rdcc->nslots);
	assert(rdcc->slot[idx_hint]);
	assert(rdcc->slot[idx_hint]->chunk==chunk);

        /*
         * It's in the cache so unlock it.
         */
        ent = rdcc->slot[idx_hint];
        assert (ent->locked);
        if (dirty) {
            ent->dirty = TRUE;
            ent->wr_count -= MIN (ent->wr_count, naccessed);
        } else {
            ent->rd_count -= MIN (ent->rd_count, naccessed);
        }
        ent->locked = FALSE;
    }

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_unlock() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_if_load
 *
 * Purpose:	A small internal function to if it's necessary to load the
 *              chunk into cache.
 *
 * Return:	TRUE or FALSE
 *
 * Programmer:	Raymond Lu
 *		17 July 2007
 *
 *-------------------------------------------------------------------------
 */
hbool_t 
H5D_istore_if_load(const H5D_io_info_t *io_info, haddr_t caddr)
{
    const H5D_t *dataset = io_info->dset;
    hbool_t ret_value;
 
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_if_load)

    HDassert(io_info);
    HDassert(dataset);

    /*
     * If the chunk is too large to load into the cache and it has no
     * filters in the pipeline (i.e. not compressed) and if the address
     * for the chunk has been defined, then don't load the chunk into the
     * cache, just write the data to it directly.
     *
     * If MPI based VFD is used, must bypass the
     * chunk-cache scheme because other MPI processes could be
     * writing to other elements in the same chunk.  Do a direct
     * write-through of only the elements requested.
     */
    if(dataset->shared->dcpl_cache.pline.nused==0 &&
            ((dataset->shared->layout.u.chunk.size > dataset->shared->cache.chunk.nbytes && caddr != HADDR_UNDEF)
#ifdef H5_HAVE_PARALLEL
            || (io_info->using_mpi_vfd && (H5F_ACC_RDWR & H5F_get_intent(dataset->oloc.file)))
#endif /* H5_HAVE_PARALLEL */
            )) {
        ret_value = FALSE;
    } else
        ret_value = TRUE;

    FUNC_LEAVE_NOAPI(ret_value)
}


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_readvv
 *
 * Purpose:	Reads a multi-dimensional buffer from (part of) an indexed raw
 *		storage array.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *		Wednesday, May  7, 2003
 *
 * Modification: 
 *              Raymond Lu
 *              20 July 2007
 *              Moved H5D_istore_lock and H5D_istore_unlock to H5D_chunk_read
 *              from this function to avoid frequent lock and unlock.
 *
 *-------------------------------------------------------------------------
 */
ssize_t
H5D_istore_readvv(const H5D_io_info_t *io_info,
    size_t chunk_max_nseq, size_t *chunk_curr_seq, size_t chunk_len_arr[], hsize_t chunk_offset_arr[],
    size_t mem_max_nseq, size_t *mem_curr_seq, size_t mem_len_arr[], hsize_t mem_offset_arr[],
    haddr_t chunk_addr, void *chunk, void *buf)
{
    H5D_t *dset=io_info->dset;          /* Local pointer to the dataset info */
    size_t		u;              /* Local index variables */
    ssize_t             ret_value;      /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_readvv, FAIL)

    /* Check args */
    HDassert(io_info);
    HDassert(dset && H5D_CHUNKED==dset->shared->layout.type);
    HDassert(dset->shared->layout.u.chunk.ndims>0 && dset->shared->layout.u.chunk.ndims<=H5O_LAYOUT_NDIMS);
    HDassert(io_info->dxpl_cache);
    HDassert(io_info->store);
    HDassert(chunk_len_arr);
    HDassert(chunk_offset_arr);
    HDassert(mem_len_arr);
    HDassert(mem_offset_arr);
    HDassert(buf);

    /* Get the address of this chunk on disk */
#ifdef QAK
HDfprintf(stderr,"%s: io_info->store->chunk.offset={",FUNC);
for(u=0; u<dset->shared->layout.u.chunk.ndims; u++)
    HDfprintf(stderr,"%Hd%s",io_info->store->chunk.offset[u],(u<(dset->shared->layout.u.chunk.ndims-1) ? ", " : "}\n"));

HDfprintf(stderr,"%s: chunk_addr=%a, chunk_size=%Zu\n",FUNC,chunk_addr,dset->shared->layout.u.chunk.size);
HDfprintf(stderr,"%s: chunk_len_arr[%Zu]=%Zu\n",FUNC,*chunk_curr_seq,chunk_len_arr[*chunk_curr_seq]);
HDfprintf(stderr,"%s: chunk_offset_arr[%Zu]=%Hu\n",FUNC,*chunk_curr_seq,chunk_offset_arr[*chunk_curr_seq]);
HDfprintf(stderr,"%s: mem_len_arr[%Zu]=%Zu\n",FUNC,*mem_curr_seq,mem_len_arr[*mem_curr_seq]);
HDfprintf(stderr,"%s: mem_offset_arr[%Zu]=%Hu\n",FUNC,*mem_curr_seq,mem_offset_arr[*mem_curr_seq]);
HDfprintf(stderr,"%s: buf=%p\n",FUNC,buf);
#endif /* QAK */

    /*
     * If the chunk is too large to load into the cache and it has no
     * filters in the pipeline (i.e. not compressed) and if the address
     * for the chunk has been defined, then don't load the chunk into the
     * cache, just read the data from it directly.
     *
     * If MPI based VFD is used, must bypass the
     * chunk-cache scheme because other MPI processes could be
     * writing to other elements in the same chunk.  Do a direct
     * read-through of only the elements requested.
     */
    if(!H5D_istore_if_load(io_info, chunk_addr)) {
        H5D_io_info_t chk_io_info;      /* Temporary I/O info object */
        H5D_storage_t chk_store;        /* Chunk storage information */

        /* Set up the storage information for the chunk */
        chk_store.contig.dset_addr=chunk_addr;
        chk_store.contig.dset_size=(hsize_t)dset->shared->layout.u.chunk.size;

        /* Set up new dataset I/O info */
        H5D_BUILD_IO_INFO(&chk_io_info,dset,io_info->dxpl_cache,io_info->dxpl_id,&chk_store);

        /* Do I/O directly on chunk without reading it into the cache */
        if ((ret_value=H5D_contig_readvv(&chk_io_info, chunk_max_nseq, chunk_curr_seq, chunk_len_arr, 
                chunk_offset_arr, mem_max_nseq, mem_curr_seq, mem_len_arr, mem_offset_arr, (haddr_t)0, NULL, buf))<0)
            HGOTO_ERROR(H5E_IO, H5E_READERROR, FAIL, "unable to read raw data to file")
    } /* end if */
    else {
        ssize_t         naccessed;      /* Number of bytes accessed in chunk */

        /* If the chunk address is not defined, check if the fill value is
         * undefined also.  If both situations hold, don't bother copying
         * values to the destination buffer, since they will just be
         * garbage.
         *
         * Ideally, this will eventually be checked at a higher level and
         * the entire I/O operation on the chunk will be skipped.  -QAK
         */
        if(!H5F_addr_defined(chunk_addr)) {
            H5D_rdcc_t		*rdcc = &(dset->shared->cache.chunk);/*raw data chunk cache*/
            hbool_t		found = FALSE;		/*already in cache?	*/

            /* Check if the chunk is in the cache (but hasn't been written to disk yet) */
            if(rdcc->nslots>0) {
                unsigned idx = H5D_HASH(dset->shared, io_info->store->chunk.index); /* Cache entry index */
                H5D_rdcc_ent_t	*ent = rdcc->slot[idx]; /* Cache entry */

                /* Potential match... */
                if(ent) {
                    for(u = 0, found = TRUE; u < dset->shared->layout.u.chunk.ndims; u++) {
                        if(io_info->store->chunk.offset[u] != ent->offset[u]) {
                            found = FALSE;
                            break;
                        } /* end if */
                    } /* end for */
                } /* end if */
            } /* end if */

            /* If the chunk is in the cache, then it must have valid data */
            if(!found) {
                const H5O_fill_t *fill = &(dset->shared->dcpl_cache.fill);    /* Fill value info */
                H5D_fill_value_t fill_status;

                /* Check if the fill value is defined */
                if(H5P_is_fill_value_defined(fill, &fill_status) < 0)
                    HGOTO_ERROR(H5E_PLIST, H5E_CANTGET, FAIL, "can't tell if fill value defined")

                /* If we are never to return fill values, or if we would return them
                 * but they aren't set, process the entire set of I/O vectors and
                 * get out now.
                 */
                if(fill->fill_time == H5D_FILL_TIME_NEVER ||
                        (fill->fill_time == H5D_FILL_TIME_IFSET && fill_status!=H5D_FILL_VALUE_USER_DEFINED)) {
                    size_t size;                /* Size of sequence in bytes */
                    size_t v;                   /* Local index variable */
                    ssize_t bytes_processed = 0;  /* Eventual return value */

                    /* Work through all the sequences */
                    for(u = *mem_curr_seq, v = *chunk_curr_seq; u < mem_max_nseq && v < chunk_max_nseq; ) {
                        /* Choose smallest buffer to write */
                        if(chunk_len_arr[v] < mem_len_arr[u])
                            size = chunk_len_arr[v];
                        else
                            size = mem_len_arr[u];

                        /* Update source information */
                        chunk_len_arr[v] -= size;
                        chunk_offset_arr[v] += size;
                        if(chunk_len_arr[v] == 0)
                            v++;

                        /* Update destination information */
                        mem_len_arr[u] -= size;
                        mem_offset_arr[u] += size;
                        if(mem_len_arr[u] == 0)
                            u++;

                        /* Increment number of bytes copied */
                        bytes_processed += (ssize_t)size;
                    } /* end for */

                    /* Update current sequence vectors */
                    *mem_curr_seq = u;
                    *chunk_curr_seq = v;

                    HGOTO_DONE(bytes_processed)
                } /* end if */
            } /* end if */
        } /* end if */

        /* Use the vectorized memory copy routine to do actual work */
        if((naccessed = H5V_memcpyvv(buf, mem_max_nseq, mem_curr_seq, mem_len_arr, mem_offset_arr, chunk, chunk_max_nseq, chunk_curr_seq, chunk_len_arr, chunk_offset_arr)) < 0)
            HGOTO_ERROR(H5E_IO, H5E_READERROR, FAIL, "vectorized memcpy failed")

        H5_CHECK_OVERFLOW(naccessed, ssize_t, size_t);

        /* Set return value */
        ret_value = naccessed;
    } /* end else */

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* H5D_istore_readvv() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_writevv
 *
 * Purpose:	Writes a multi-dimensional buffer to (part of) an indexed raw
 *		storage array.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *		Friday, May  2, 2003
 *
 * Modification: 
 *              Raymond Lu
 *              20 July 2007
 *              Moved H5D_istore_lock and H5D_istore_unlock to H5D_chunk_write
 *              from this function to avoid frequent lock and unlock.
 *
 *-------------------------------------------------------------------------
 */
ssize_t
H5D_istore_writevv(const H5D_io_info_t *io_info,
    size_t chunk_max_nseq, size_t *chunk_curr_seq, size_t chunk_len_arr[], hsize_t chunk_offset_arr[],
    size_t mem_max_nseq, size_t *mem_curr_seq, size_t mem_len_arr[], hsize_t mem_offset_arr[],
    haddr_t chunk_addr, void *chunk, const void *buf)
{
    H5D_t *dset = io_info->dset;          /* Local pointer to the dataset info */
    ssize_t             ret_value;      /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_writevv, FAIL)

    /* Check args */
    HDassert(io_info);
    HDassert(dset && H5D_CHUNKED==dset->shared->layout.type);
    HDassert(dset->shared->layout.u.chunk.ndims>0 && dset->shared->layout.u.chunk.ndims<=H5O_LAYOUT_NDIMS);
    HDassert(io_info->dxpl_cache);
    HDassert(io_info->store);
    HDassert(chunk_len_arr);
    HDassert(chunk_offset_arr);
    HDassert(mem_len_arr);
    HDassert(mem_offset_arr);
    HDassert(buf);

#ifdef QAK
{
size_t		u;              /* Local index variables */

HDfprintf(stderr,"%s: io_info->store->chunk.offset={",FUNC);
for(u=0; u<dset->shared->layout.u.chunk.ndims; u++)
    HDfprintf(stderr,"%Hd%s",io_info->store->chunk.offset[u],(u<(dset->shared->layout.u.chunk.ndims-1) ? ", " : "}\n"));

HDfprintf(stderr,"%s: chunk_addr=%a, chunk_size=%Zu\n",FUNC,chunk_addr,dset->shared->layout.u.chunk.size);
HDfprintf(stderr,"%s: chunk_len_arr[%Zu]=%Zu\n",FUNC,*chunk_curr_seq,chunk_len_arr[*chunk_curr_seq]);
HDfprintf(stderr,"%s: chunk_offset_arr[%Zu]=%Hu\n",FUNC,*chunk_curr_seq,chunk_offset_arr[*chunk_curr_seq]);
HDfprintf(stderr,"%s: mem_len_arr[%Zu]=%Zu\n",FUNC,*mem_curr_seq,mem_len_arr[*mem_curr_seq]);
HDfprintf(stderr,"%s: mem_offset_arr[%Zu]=%Hu\n",FUNC,*mem_curr_seq,mem_offset_arr[*mem_curr_seq]);
}
#endif /* QAK */

    /*
     * If the chunk is too large to load into the cache and it has no
     * filters in the pipeline (i.e. not compressed) and if the address
     * for the chunk has been defined, then don't load the chunk into the
     * cache, just write the data to it directly.
     *
     * If MPI based VFD is used, must bypass the
     * chunk-cache scheme because other MPI processes could be
     * writing to other elements in the same chunk.  Do a direct
     * write-through of only the elements requested.
     */
    if(!H5D_istore_if_load(io_info, chunk_addr)) {
        H5D_io_info_t chk_io_info;      /* Temporary I/O info object */
        H5D_storage_t chk_store;        /* Chunk storage information */

        /* Set up the storage information for the chunk */
        chk_store.contig.dset_addr=chunk_addr;
        chk_store.contig.dset_size=(hsize_t)dset->shared->layout.u.chunk.size;

        /* Set up new dataset I/O info */
        H5D_BUILD_IO_INFO(&chk_io_info,dset,io_info->dxpl_cache,io_info->dxpl_id,&chk_store);

        /* Do I/O directly on chunk without reading it into the cache */
        if((ret_value = H5D_contig_writevv(&chk_io_info, chunk_max_nseq, chunk_curr_seq, chunk_len_arr, chunk_offset_arr, mem_max_nseq, mem_curr_seq, mem_len_arr, mem_offset_arr, (haddr_t)0, NULL, buf)) < 0)
            HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to write raw data to file")
    } /* end if */
    else {
        ssize_t         naccessed;      /* Number of bytes accessed in chunk */

        /* Use the vectorized memory copy routine to do actual work */
        if((naccessed=H5V_memcpyvv(chunk,chunk_max_nseq,chunk_curr_seq,chunk_len_arr,chunk_offset_arr,buf,mem_max_nseq,mem_curr_seq,mem_len_arr,mem_offset_arr))<0)
            HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "vectorized memcpy failed")

        H5_CHECK_OVERFLOW(naccessed,ssize_t,size_t);

        /* Set return value */
        ret_value=naccessed;
    } /* end else */

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* H5D_istore_writevv() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_create
 *
 * Purpose:	Creates a new indexed-storage B-tree and initializes the
 *		istore struct with information about the storage.  The
 *		struct should be immediately written to the object header.
 *
 *		This function must be called before passing ISTORE to any of
 *		the other indexed storage functions!
 *
 * Return:	Non-negative on success (with the ISTORE argument initialized
 *		and ready to write to an object header). Negative on failure.
 *
 * Programmer:	Robb Matzke
 *		Tuesday, October 21, 1997
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_create(H5F_t *f, hid_t dxpl_id, H5O_layout_t *layout /*out */)
{
    H5D_istore_ud0_t	udata;
#ifndef NDEBUG
    unsigned			u;
#endif
    herr_t      ret_value = SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_create, FAIL)

    /* Check args */
    HDassert(f);
    HDassert(layout && H5D_CHUNKED == layout->type);
    HDassert(layout->u.chunk.ndims > 0 && layout->u.chunk.ndims <= H5O_LAYOUT_NDIMS);
#ifndef NDEBUG
    for(u = 0; u < layout->u.chunk.ndims; u++)
	HDassert(layout->u.chunk.dim[u] > 0);
#endif

    /* Initialize "user" data for B-tree callbacks, etc. */
    udata.mesg = layout;

    if(H5B_create(f, dxpl_id, H5B_ISTORE, &udata, &(layout->u.chunk.addr)/*out*/) < 0)
	HGOTO_ERROR(H5E_IO, H5E_CANTINIT, FAIL, "can't create B-tree")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_create() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_allocated
 *
 * Purpose:	Return the number of bytes allocated in the file for storage
 *		of raw data under the specified B-tree (ADDR is the address
 *		of the B-tree).
 *
 * Return:	Success:	Number of bytes stored in all chunks.
 *
 *		Failure:	0
 *
 * Programmer:	Robb Matzke
 *              Wednesday, April 21, 1999
 *
 *-------------------------------------------------------------------------
 */
hsize_t
H5D_istore_allocated(H5D_t *dset, hid_t dxpl_id)
{
    H5D_io_info_t io_info;              /* Temporary I/O info object */
    const H5D_rdcc_t   *rdcc = &(dset->shared->cache.chunk);	/*raw data chunk cache */
    H5D_rdcc_ent_t     *ent;    /*cache entry  */
    H5D_dxpl_cache_t _dxpl_cache;       /* Data transfer property cache buffer */
    H5D_dxpl_cache_t *dxpl_cache=&_dxpl_cache;   /* Data transfer property cache */
    H5D_istore_it_ud1_t	udata;
    hsize_t      ret_value;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_allocated, 0)

    HDassert(dset);

    /* Fill the DXPL cache values for later use */
    if (H5D_get_dxpl_cache(dxpl_id,&dxpl_cache)<0)
        HGOTO_ERROR(H5E_DATASET, H5E_CANTGET, 0, "can't fill dxpl cache")

    /* Construct dataset I/O info */
    H5D_BUILD_IO_INFO(&io_info,dset,dxpl_cache,dxpl_id,NULL);

    /* Search for cached chunks that haven't been written out */
    for(ent = rdcc->head; ent; ent = ent->next) {
        /* Flush the chunk out to disk, to make certain the size is correct later */
        if (H5D_istore_flush_entry(&io_info, ent, FALSE)<0)
            HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, 0, "cannot flush indexed storage buffer")
    } /* end for */

    HDmemset(&udata, 0, sizeof udata);
    udata.common.mesg = &dset->shared->layout;
    if(H5B_iterate(dset->oloc.file, dxpl_id, H5B_ISTORE, H5D_istore_iter_allocated, dset->shared->layout.u.chunk.addr, &udata) < 0)
        HGOTO_ERROR(H5E_IO, H5E_CANTINIT, 0, "unable to iterate over chunk B-tree")

    /* Set return value */
    ret_value = udata.total_storage;

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_allocated() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_chunkmap
 *
 * Purpose:     obtain the chunk address and corresponding chunk index
 *
 * Return:	Success:	Non-negative on succeed.
 *
 *		Failure:	negative value
 *
 * Programmer:  Kent Yang
 *              November 15, 2005
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_chunkmap(const H5D_io_info_t *io_info, haddr_t chunk_addr[],
    hsize_t down_chunks[])
{
    H5D_t *dset = io_info->dset;       /* Local pointer to dataset info */
    H5D_istore_it_ud5_t	udata;
    herr_t      ret_value = SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_chunkmap, FAIL)

    HDassert(dset);

    /* Set up user data for B-tree callback */
    HDmemset(&udata, 0, sizeof(udata));
    udata.common.mesg = &dset->shared->layout;
    udata.down_chunks = down_chunks;
    udata.chunk_addr  = chunk_addr;

    /* Build mapping of chunk addresses and indices */
    if(H5B_iterate(dset->oloc.file, io_info->dxpl_id, H5B_ISTORE, H5D_istore_iter_chunkmap, dset->shared->layout.u.chunk.addr, &udata) < 0)
        HGOTO_ERROR(H5E_IO, H5E_CANTINIT, FAIL, "unable to iterate over chunk B-tree")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_chunkmap() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_get_addr
 *
 * Purpose:	Get the file address of a chunk if file space has been
 *		assigned.  Save the retrieved information in the udata
 *		supplied.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Albert Cheng
 *              June 27, 1998
 *
 *-------------------------------------------------------------------------
 */
haddr_t
H5D_istore_get_addr(const H5D_io_info_t *io_info, H5D_istore_ud1_t *_udata)
{
    H5D_istore_ud1_t	tmp_udata;      /* Information about a chunk */
    H5D_istore_ud1_t	*udata;         /* Pointer to information about a chunk */
    haddr_t	ret_value;		/* Return value */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_get_addr)

    HDassert(io_info);
    HDassert(io_info->dset);
    HDassert(io_info->dset->shared->layout.u.chunk.ndims > 0);
    HDassert(io_info->store->chunk.offset);

    /* Check for udata struct to return */
    udata = (_udata != NULL ? _udata : &tmp_udata);

    /* Initialize the information about the chunk we are looking for */
    udata->common.mesg = &(io_info->dset->shared->layout);
    udata->common.offset = io_info->store->chunk.offset;
    udata->nbytes = 0;
    udata->filter_mask = 0;
    udata->addr = HADDR_UNDEF;

    /* Check for cached information */
    if(!H5D_istore_cinfo_cache_found(&io_info->dset->shared->cache.chunk.last, udata)) {
        /* Go get the chunk information */
        if(H5B_find(io_info->dset->oloc.file, io_info->dxpl_id, H5B_ISTORE, io_info->dset->shared->layout.u.chunk.addr, udata) < 0) {
            /* Note: don't push error on stack, leave that to next higher level,
             *      since many times the B-tree is searched in order to determine
             *      if a chunk exists in the B-tree or not. -QAK
             */
#ifdef OLD_WAY
            H5E_clear_stack(NULL);

            HGOTO_ERROR(H5E_BTREE, H5E_NOTFOUND, HADDR_UNDEF, "Can't locate chunk info")
#else /* OLD_WAY */
            /* Cache the fact that the chunk is not in the B-tree */
            H5D_istore_cinfo_cache_update(&io_info->dset->shared->cache.chunk.last, udata);

            HGOTO_DONE(HADDR_UNDEF)
#endif /* OLD_WAY */
        } /* end if */

        /* Cache the information retrieved */
        HDassert(H5F_addr_defined(udata->addr));
        H5D_istore_cinfo_cache_update(&io_info->dset->shared->cache.chunk.last, udata);
    } /* end else */

    /* Success!  Set the return value */
    ret_value = udata->addr;

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* H5D_istore_get_addr() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_chunk_alloc
 *
 * Purpose:	Allocate space for a chunk in memory.  This routine allocates
 *              memory space for non-filtered chunks from a block free list
 *              and uses malloc()/free() for filtered chunks.
 *
 * Return:	Pointer to memory for chunk on success/NULL on failure
 *
 * Programmer:	Quincey Koziol
 *              April 22, 2004
 *
 *-------------------------------------------------------------------------
 */
static void *
H5D_istore_chunk_alloc(size_t size, const H5O_pline_t *pline)
{
    void *ret_value = NULL;		/* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_chunk_alloc)

    HDassert(size);
    HDassert(pline);

    if(pline->nused > 0)
        ret_value = H5MM_malloc(size);
    else
        ret_value = H5FL_BLK_MALLOC(chunk, size);

    FUNC_LEAVE_NOAPI(ret_value)
} /* H5D_istore_chunk_alloc() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_chunk_xfree
 *
 * Purpose:	Free space for a chunk in memory.  This routine allocates
 *              memory space for non-filtered chunks from a block free list
 *              and uses malloc()/free() for filtered chunks.
 *
 * Return:	NULL (never fails)
 *
 * Programmer:	Quincey Koziol
 *              April 22, 2004
 *
 *-------------------------------------------------------------------------
 */
static void *
H5D_istore_chunk_xfree(void *chk, const H5O_pline_t *pline)
{
    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5D_istore_chunk_xfree)

    HDassert(pline);

    if(chk) {
        if(pline->nused > 0)
            H5MM_xfree(chk);
        else
            (void)H5FL_BLK_FREE(chunk, chk);
    } /* end if */

    FUNC_LEAVE_NOAPI(NULL)
} /* H5D_istore_chunk_xfree() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_allocate
 *
 * Purpose:	Allocate file space for all chunks that are not allocated yet.
 *		Return SUCCEED if all needed allocation succeed, otherwise
 *		FAIL.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Note:	Current implementation relies on cache_size being 0,
 *		thus no chunk is cached and written to disk immediately
 *		when a chunk is unlocked (via H5F_istore_unlock)
 *		This should be changed to do a direct flush independent
 *		of the cache value.
 *
 * Programmer:	Albert Cheng
 *		June 26, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_allocate(H5D_t *dset, hid_t dxpl_id, hbool_t full_overwrite)
{
    H5D_io_info_t io_info;      /* Dataset I/O info */
    H5D_storage_t store;        /* Dataset storage information */
    hsize_t	chunk_offset[H5O_LAYOUT_NDIMS]; /* Offset of current chunk */
    size_t	orig_chunk_size; /* Original size of chunk in bytes */
    unsigned    filter_mask = 0; /* Filter mask for chunks that have them */
    const H5O_layout_t *layout = &(dset->shared->layout);       /* Dataset layout */
    const H5O_pline_t *pline = &(dset->shared->dcpl_cache.pline);    /* I/O pipeline info */
    const H5O_fill_t *fill = &(dset->shared->dcpl_cache.fill);    /* Fill value info */
    H5D_fill_value_t fill_status; /* The fill value status */
    hbool_t     should_fill = FALSE; /* Whether fill values should be written */
    H5D_dxpl_cache_t _dxpl_cache;       /* Data transfer property cache buffer */
    H5D_dxpl_cache_t *dxpl_cache = &_dxpl_cache;   /* Data transfer property cache */
#ifdef H5_HAVE_PARALLEL
    MPI_Comm	mpi_comm = MPI_COMM_NULL;	/* MPI communicator for file */
    int         mpi_rank = (-1); /* This process's rank  */
    int         mpi_code;       /* MPI return code */
    hbool_t     blocks_written = FALSE; /* Flag to indicate that chunk was actually written */
    hbool_t     using_mpi = FALSE;    /* Flag to indicate that the file is being accessed with an MPI-capable file driver */
#endif /* H5_HAVE_PARALLEL */
    hbool_t	carry;          /* Flag to indicate that chunk increment carrys to higher dimension (sorta) */
    int         space_ndims;    /* Dataset's space rank */
    hsize_t     space_dim[H5O_LAYOUT_NDIMS];    /* Dataset's dataspace dimensions */
    H5D_fill_buf_info_t fb_info;        /* Dataset's fill buffer info */
    hbool_t     fb_info_init = FALSE;   /* Whether the fill value buffer has been initialized */
    hid_t       data_dxpl_id;           /* DXPL ID to use for raw data I/O operations */
    herr_t	ret_value = SUCCEED;	/* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_allocate, FAIL)

    /* Check args */
    HDassert(dset && H5D_CHUNKED == layout->type);
    HDassert(layout->u.chunk.ndims > 0 && layout->u.chunk.ndims <= H5O_LAYOUT_NDIMS);
    HDassert(H5F_addr_defined(layout->u.chunk.addr));
    HDassert(TRUE == H5P_isa_class(dxpl_id, H5P_DATASET_XFER));

    /* Retrieve the dataset dimensions */
    if((space_ndims = H5S_get_simple_extent_dims(dset->shared->space, space_dim, NULL)) < 0)
         HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "unable to get simple dataspace info")
    space_dim[space_ndims] = layout->u.chunk.dim[space_ndims];

#ifdef H5_HAVE_PARALLEL
    /* Retrieve MPI parameters */
    if(IS_H5FD_MPI(dset->oloc.file)) {
        /* Get the MPI communicator */
        if(MPI_COMM_NULL == (mpi_comm = H5F_mpi_get_comm(dset->oloc.file)))
            HGOTO_ERROR(H5E_INTERNAL, H5E_MPI, FAIL, "Can't retrieve MPI communicator")

        /* Get the MPI rank */
        if((mpi_rank = H5F_mpi_get_rank(dset->oloc.file)) < 0)
            HGOTO_ERROR(H5E_INTERNAL, H5E_MPI, FAIL, "Can't retrieve MPI rank")

        /* Set the MPI-capable file driver flag */
        using_mpi = TRUE;

        /* Use the internal "independent" DXPL */
        data_dxpl_id = H5AC_ind_dxpl_id;
    } /* end if */
    else {
#endif  /* H5_HAVE_PARALLEL */
        /* Use the DXPL we were given */
        data_dxpl_id = dxpl_id;
#ifdef H5_HAVE_PARALLEL
    } /* end else */
#endif  /* H5_HAVE_PARALLEL */

    /* Fill the DXPL cache values for later use */
    if(H5D_get_dxpl_cache(data_dxpl_id, &dxpl_cache) < 0)
        HGOTO_ERROR(H5E_DATASET, H5E_CANTGET, FAIL, "can't fill dxpl cache")

    /* Get original chunk size */
    H5_CHECK_OVERFLOW(layout->u.chunk.size, hsize_t, size_t);
    orig_chunk_size = (size_t)layout->u.chunk.size;

    /* Check the dataset's fill-value status */
    if(H5P_is_fill_value_defined(fill, &fill_status) < 0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTGET, FAIL, "can't tell if fill value defined")

    /* If we are filling the dataset on allocation or "if set" and
     * the fill value _is_ set, _and_ we are not overwriting the new blocks,
     * or if there are any pipeline filters defined,
     * set the "should fill" flag
     */
    if((!full_overwrite && (fill->fill_time == H5D_FILL_TIME_ALLOC ||
            (fill->fill_time == H5D_FILL_TIME_IFSET && fill_status == H5D_FILL_VALUE_USER_DEFINED)))
            || pline->nused > 0)
        should_fill = TRUE;

    /* Check if fill values should be written to chunks */
    if(should_fill) {
        /* Initialize the fill value buffer */
        /* (delay allocating fill buffer for VL datatypes until refilling) */
        /* (casting away const OK - QAK) */
        if(H5D_fill_init(&fb_info, NULL, (hbool_t)(pline->nused > 0),
                (H5MM_allocate_t)H5D_istore_chunk_alloc, (void *)pline,
                (H5MM_free_t)H5D_istore_chunk_xfree, (void *)pline,
                &dset->shared->dcpl_cache.fill, dset->shared->type,
                dset->shared->type_id, (size_t)0, orig_chunk_size, data_dxpl_id) < 0)
            HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "can't initialize fill buffer info")
        fb_info_init = TRUE;

       /* Check if there are filters which need to be applied to the chunk */
       /* (only do this in advance when the chunk info can be re-used (i.e.
        *      it doesn't contain any non-default VL datatype fill values)
        */
       if(!fb_info.has_vlen_fill_type && pline->nused > 0) {
           size_t buf_size = orig_chunk_size;

           /* Push the chunk through the filters */
           if(H5Z_pipeline(pline, 0, &filter_mask, dxpl_cache->err_detect, dxpl_cache->filter_cb, &orig_chunk_size, &buf_size, &fb_info.fill_buf) < 0)
               HGOTO_ERROR(H5E_PLINE, H5E_WRITEERROR, FAIL, "output pipeline failed")
       } /* end if */
    } /* end if */

    /* Set up dataset I/O info */
    store.chunk.offset = chunk_offset;
    H5D_BUILD_IO_INFO(&io_info, dset, dxpl_cache, data_dxpl_id, &store);

    /* Reset the chunk offset indices */
    HDmemset(chunk_offset, 0, (layout->u.chunk.ndims * sizeof(chunk_offset[0])));

    /* Loop over all chunks */
    carry = FALSE;
    while(!carry) {
        int i;                  /* Local index variable */

        /* Check if the chunk exists yet on disk */
        if(!H5F_addr_defined(H5D_istore_get_addr(&io_info, NULL))) {
            const H5D_rdcc_t *rdcc = &(dset->shared->cache.chunk);  /* Raw data chunk cache */
            H5D_rdcc_ent_t *ent;    /* Cache entry  */
            hbool_t chunk_exists;   /* Flag to indicate whether a chunk exists already */
            unsigned u;             /* Local index variable */

            /* Didn't find the chunk on disk */
            chunk_exists = FALSE;

            /* Look for chunk in cache */
            for(ent = rdcc->head; ent && !chunk_exists; ent = ent->next) {
                /* Assume a match */
                chunk_exists = TRUE;
                for(u = 0; u < layout->u.chunk.ndims; u++)
                    if(ent->offset[u] != chunk_offset[u]) {
                        chunk_exists = FALSE;       /* Reset if no match */
                        break;
                    } /* end if */
            } /* end for */

            /* Chunk wasn't in cache either, create it now */
            if(!chunk_exists) {
                H5D_istore_ud1_t udata;	/* B-tree pass-through for creating chunk */
                size_t	chunk_size;     /* Size of chunk in bytes, possibly filtered */

                /* Check for VL datatype & non-default fill value */
                if(fb_info_init && fb_info.has_vlen_fill_type) {
                    /* Sanity check */
                    HDassert(should_fill);

                    /* Fill the buffer with VL datatype fill values */
                    if(H5D_fill_refill_vl(&fb_info, fb_info.elmts_per_buf, data_dxpl_id) < 0)
                        HGOTO_ERROR(H5E_DATASET, H5E_CANTCONVERT, FAIL, "can't refill fill value buffer")

                    /* Check if there are filters which need to be applied to the chunk */
                    if(pline->nused > 0) {
                        size_t buf_size = orig_chunk_size;
                        size_t nbytes = fb_info.fill_buf_size;

                        /* Push the chunk through the filters */
                        if(H5Z_pipeline(pline, 0, &filter_mask, dxpl_cache->err_detect, dxpl_cache->filter_cb, &nbytes, &buf_size, &fb_info.fill_buf) < 0)
                            HGOTO_ERROR(H5E_PLINE, H5E_WRITEERROR, FAIL, "output pipeline failed")

                        /* Keep the number of bytes the chunk turned in to */
                        chunk_size = nbytes;
                    } /* end if */
                    else
                        chunk_size = (size_t)layout->u.chunk.size;
                } /* end if */
                else
                    chunk_size = orig_chunk_size;

                /* Initialize the chunk information */
                udata.common.mesg = layout;
                udata.common.offset = chunk_offset;
                udata.nbytes = chunk_size;
                udata.filter_mask = filter_mask;
                udata.addr = HADDR_UNDEF;

                /* Allocate the chunk with all processes */
                if(H5B_insert(dset->oloc.file, dxpl_id, H5B_ISTORE, layout->u.chunk.addr, &udata) < 0)
                    HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to allocate chunk")

                /* Check if fill values should be written to chunks */
                if(should_fill) {
                    /* Sanity check */
                    HDassert(fb_info_init);

#ifdef H5_HAVE_PARALLEL
                    /* Check if this file is accessed with an MPI-capable file driver */
                    if(using_mpi) {
                        /* Write the chunks out from only one process */
                        /* !! Use the internal "independent" DXPL!! -QAK */
                        if(H5_PAR_META_WRITE == mpi_rank)
                            if(H5F_block_write(dset->oloc.file, H5FD_MEM_DRAW, udata.addr, udata.nbytes, data_dxpl_id, fb_info.fill_buf) < 0)
                                HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to write raw data to file")

                        /* Indicate that blocks are being written */
                        blocks_written = TRUE;
                    } /* end if */
                    else {
#endif /* H5_HAVE_PARALLEL */
                        if(H5F_block_write(dset->oloc.file, H5FD_MEM_DRAW, udata.addr, udata.nbytes, data_dxpl_id, fb_info.fill_buf) < 0)
                            HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to write raw data to file")
#ifdef H5_HAVE_PARALLEL
                    } /* end else */
#endif /* H5_HAVE_PARALLEL */
                } /* end if */

                /* Release the fill buffer if we need to re-allocate it each time */
                if(fb_info_init && fb_info.has_vlen_fill_type && pline->nused > 0)
                    H5D_fill_release(&fb_info);
            } /* end if */
        } /* end if */

        /* Increment indices */
        carry = TRUE;
        for(i = (int)(space_ndims - 1); i >= 0; --i) {
            chunk_offset[i] += layout->u.chunk.dim[i];
            if(chunk_offset[i] >= space_dim[i])
                chunk_offset[i] = 0;
            else {
                carry = FALSE;
                break;
            } /* end else */
        } /* end for */
    } /* end while */

#ifdef H5_HAVE_PARALLEL
    /* Only need to block at the barrier if we actually initialized a chunk */
    /* using an MPI-capable file driver */
    if(using_mpi && blocks_written) {
        /* Wait at barrier to avoid race conditions where some processes are
         * still writing out chunks and other processes race ahead to read
         * them in, getting bogus data.
         */
        if(MPI_SUCCESS != (mpi_code = MPI_Barrier(mpi_comm)))
            HMPI_GOTO_ERROR(FAIL, "MPI_Barrier failed", mpi_code)
    } /* end if */
#endif /* H5_HAVE_PARALLEL */

    /* Reset any cached chunk info for this dataset */
    H5D_istore_cinfo_cache_reset(&dset->shared->cache.chunk.last);

done:
    /* Release the fill buffer info, if it's been initialized */
    if(fb_info_init && H5D_fill_term(&fb_info) < 0)
        HDONE_ERROR(H5E_DATASET, H5E_CANTFREE, FAIL, "Can't release fill buffer info")

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_allocate() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_prune_check
 *
 * Purpose:	Search for chunks that are no longer necessary in the B-tree.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Pedro Vicente, pvn@ncsa.uiuc.edu
 * 		March 26, 2002
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static int
H5D_istore_prune_check(H5F_t UNUSED *f, hid_t UNUSED dxpl_id,
    const void *_lt_key, haddr_t UNUSED addr, const void UNUSED *_rt_key,
    void *_udata)
{
    H5D_istore_it_ud3_t       *udata = (H5D_istore_it_ud3_t *)_udata;
    const H5D_istore_key_t       *lt_key = (const H5D_istore_key_t *)_lt_key;
    unsigned                rank;	/*current # of dimensions */
    unsigned                u;
    int                     ret_value = H5_ITER_CONT;       /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_prune_check)

    /* Figure out what chunks are no longer in use for the specified extent and release them */
    rank = udata->common.mesg->u.chunk.ndims - 1;
    for(u = 0; u < rank; u++)
        /* The LT_KEY is the left key (the one that describes the chunk). It points to a chunk of
         * storage that contains the beginning of the logical address space represented by UDATA.
         */
	if((hsize_t)lt_key->offset[u] > udata->dims[u]) {
            H5D_istore_sl_ck_t *sl_node;        /* Skip list node for chunk to remove */

            /* Allocate space for the shared structure */
            if(NULL == (sl_node = H5FL_MALLOC(H5D_istore_sl_ck_t)))
                HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, H5_ITER_ERROR, "memory allocation failed for shared B-tree info")

            /* Calculate the index of this chunk */
            if(H5V_chunk_index(rank, lt_key->offset, udata->common.mesg->u.chunk.dim, udata->down_chunks, &sl_node->index) < 0) {
                H5FL_FREE(H5D_istore_sl_ck_t, sl_node);
                HGOTO_ERROR(H5E_IO, H5E_BADRANGE, H5_ITER_ERROR, "can't get chunk index")
            } /* end if */

            /* Store the key for the chunk */
            sl_node->key = *lt_key;

            /* Insert the chunk description in the skip list */
            if(H5SL_insert(udata->outside, sl_node, &sl_node->index) < 0) {
                H5FL_FREE(H5D_istore_sl_ck_t, sl_node);
                HGOTO_ERROR(H5E_IO, H5E_CANTINSERT, H5_ITER_ERROR, "can't insert chunk into skip list")
            } /* end if */

            /* Break out of loop, we know the chunk is outside the current dimensions */
	    break;
	} /* end if */

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_prune_check() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_prune_remove
 *
 * Purpose:	Destroy a skip list node for "pruning" chunks, also removes
 *              the chunk from the B-tree.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol, koziol@hdfgroup.org
 * 		May 3, 2007
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5D_istore_prune_remove(void *item, void UNUSED *key, void *op_data)
{
    H5D_istore_sl_ck_t *sl_node = (H5D_istore_sl_ck_t *)item;           /* Temporary pointer to chunk to remove */
    H5D_istore_sl_rm_t *rm_info = (H5D_istore_sl_rm_t *)op_data;       /* Information needed for removing chunk from B-tree */
    H5D_istore_ud0_t bt_udata;                  /* User data for B-tree removal routine */
    herr_t ret_value = H5_ITER_CONT;            /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_prune_remove)

    /* Sanity checks */
    HDassert(sl_node);
    HDassert(rm_info);

    /* Initialize the user data for the B-tree callback */
    bt_udata.mesg = rm_info->mesg;
    bt_udata.offset = sl_node->key.offset;

    /* Remove */
    if(H5B_remove(rm_info->f, rm_info->dxpl_id, H5B_ISTORE, rm_info->mesg->u.chunk.addr, &bt_udata) < 0)
        HGOTO_ERROR(H5E_IO, H5E_CANTINIT, H5_ITER_ERROR, "unable to remove entry")

    /* Free the chunk checking node */
    H5FL_FREE(H5D_istore_sl_ck_t, sl_node);

done:
    FUNC_LEAVE_NOAPI(ret_value)
}   /* H5D_istore_prune_remove() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_prune_by_extent
 *
 * Purpose:	This function searches for chunks that are no longer necessary both in the
 *		raw data cache and in the B-tree.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Pedro Vicente, pvn@ncsa.uiuc.edu
 * Algorithm:	Robb Matzke
 * 		March 27, 2002
 *
 * The algorithm is:
 *
 *  For chunks that are no longer necessary:
 *
 *  1. Search in the raw data cache for each chunk
 *  2. If found then preempt it from the cache
 *  3. Search in the B-tree for each chunk
 *  4. If found then remove it from the B-tree and deallocate file storage for the chunk
 *
 * This example shows a 2d dataset of 90x90 with a chunk size of 20x20.
 *
 *
 *     0         20        40        60        80    90   100
 *    0 +---------+---------+---------+---------+-----+...+
 *      |:::::X::::::::::::::         :         :     |   :
 *      |:::::::X::::::::::::         :         :     |   :   Key
 *      |::::::::::X:::::::::         :         :     |   :   --------
 *      |::::::::::::X:::::::         :         :     |   :  +-+ Dataset
 *    20+::::::::::::::::::::.........:.........:.....+...:  | | Extent
 *      |         :::::X:::::         :         :     |   :  +-+
 *      |         :::::::::::         :         :     |   :
 *      |         :::::::::::         :         :     |   :  ... Chunk
 *      |         :::::::X:::         :         :     |   :  : : Boundary
 *    40+.........:::::::::::.........:.........:.....+...:  :.:
 *      |         :         :         :         :     |   :
 *      |         :         :         :         :     |   :  ... Allocated
 *      |         :         :         :         :     |   :  ::: & Filled
 *      |         :         :         :         :     |   :  ::: Chunk
 *    60+.........:.........:.........:.........:.....+...:
 *      |         :         :::::::X:::         :     |   :   X  Element
 *      |         :         :::::::::::         :     |   :      Written
 *      |         :         :::::::::::         :     |   :
 *      |         :         :::::::::::         :     |   :
 *    80+.........:.........:::::::::::.........:.....+...:   O  Fill Val
 *      |         :         :         :::::::::::     |   :      Explicitly
 *      |         :         :         ::::::X::::     |   :      Written
 *    90+---------+---------+---------+---------+-----+   :
 *      :         :         :         :::::::::::         :
 *   100:.........:.........:.........:::::::::::.........:
 *
 *
 * We have 25 total chunks for this dataset, 5 of which have space
 * allocated in the file because they were written to one or more
 * elements. These five chunks (and only these five) also have entries in
 * the storage B-tree for this dataset.
 *
 * Now lets say we want to shrink the dataset down to 70x70:
 *
 *
 *      0         20        40        60   70   80    90   100
 *    0 +---------+---------+---------+----+----+-----+...+
 *      |:::::X::::::::::::::         :    |    :     |   :
 *      |:::::::X::::::::::::         :    |    :     |   :    Key
 *      |::::::::::X:::::::::         :    |    :     |   :    --------
 *      |::::::::::::X:::::::         :    |    :     |   :   +-+ Dataset
 *    20+::::::::::::::::::::.........:....+....:.....|...:   | | Extent
 *      |         :::::X:::::         :    |    :     |   :   +-+
 *      |         :::::::::::         :    |    :     |   :
 *      |         :::::::::::         :    |    :     |   :   ... Chunk
 *      |         :::::::X:::         :    |    :     |   :   : : Boundary
 *    40+.........:::::::::::.........:....+....:.....|...:   :.:
 *      |         :         :         :    |    :     |   :
 *      |         :         :         :    |    :     |   :   ... Allocated
 *      |         :         :         :    |    :     |   :   ::: & Filled
 *      |         :         :         :    |    :     |   :   ::: Chunk
 *    60+.........:.........:.........:....+....:.....|...:
 *      |         :         :::::::X:::    |    :     |   :    X  Element
 *      |         :         :::::::::::    |    :     |   :       Written
 *      +---------+---------+---------+----+    :     |   :
 *      |         :         :::::::::::         :     |   :
 *    80+.........:.........:::::::::X:.........:.....|...:    O  Fill Val
 *      |         :         :         :::::::::::     |   :       Explicitly
 *      |         :         :         ::::::X::::     |   :       Written
 *    90+---------+---------+---------+---------+-----+   :
 *      :         :         :         :::::::::::         :
 *   100:.........:.........:.........:::::::::::.........:
 *
 *
 * That means that the nine chunks along the bottom and right side should
 * no longer exist. Of those nine chunks, (0,80), (20,80), (40,80),
 * (60,80), (80,80), (80,60), (80,40), (80,20), and (80,0), one is actually allocated
 * that needs to be released.
 * To release the chunks, we traverse the B-tree to obtain a list of unused
 * allocated chunks, and then call H5B_remove() for each chunk.
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_prune_by_extent(const H5D_io_info_t *io_info, const hsize_t *old_dims)
{
    H5D_t *dset = io_info->dset;          /* Local pointer to the dataset info */
    const H5D_rdcc_t       *rdcc = &(dset->shared->cache.chunk);	/*raw data chunk cache */
    H5D_rdcc_ent_t         *ent = NULL, *next = NULL;	/*cache entry  */
    H5D_istore_it_ud3_t     udata;	/*B-tree pass-through */
    H5D_istore_sl_rm_t      rm_info;    /* User data for skip list destroy callback */
    hsize_t                 curr_dims[H5O_LAYOUT_NDIMS];	/*current dataspace dimensions */
    hsize_t                 chunks[H5O_LAYOUT_NDIMS];	        /*current number of chunks in each dimension */
    hsize_t                 down_chunks[H5O_LAYOUT_NDIMS];      /* "down" size of number of elements in each dimension */
    unsigned                rank;	/* Current # of dimensions */
    unsigned                u;	        /* Local index variable */
    herr_t                  ret_value = SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_prune_by_extent, FAIL)

    /* Check args */
    HDassert(io_info);
    HDassert(dset && H5D_CHUNKED == dset->shared->layout.type);
    HDassert(dset->shared->layout.u.chunk.ndims > 0 && dset->shared->layout.u.chunk.ndims <= H5O_LAYOUT_NDIMS);
    HDassert(H5F_addr_defined(dset->shared->layout.u.chunk.addr));

    /* Go get the rank & dimensions */
    rank = dset->shared->layout.u.chunk.ndims - 1;
    if(H5S_get_simple_extent_dims(dset->shared->space, curr_dims, NULL) < 0)
	HGOTO_ERROR(H5E_DATASET, H5E_CANTGET, FAIL, "can't get dataset dimensions")

    /*-------------------------------------------------------------------------
     * Figure out what chunks are no longer in use for the specified extent
     * and release them from the linked list raw data cache
     *-------------------------------------------------------------------------
     */
    for(ent = rdcc->head; ent; ent = next) {
        /* Get pointer to next extry in cache, in case this one is evicted */
	next = ent->next;

        /* Check for chunk offset outside of new dimensions */
        for(u = 0; u < rank; u++)
            if((hsize_t)ent->offset[u] > curr_dims[u]) {
#ifdef H5D_ISTORE_DEBUG
                HDfputs("cache:remove:[", stderr);
                for(u = 0; u < rank; u++)
                    HDfprintf(stderr, "%s%Hd", (u ? ", " : ""), ent->offset[u]);
                HDfputs("]\n", stderr);
#endif

                /* Preempt the entry from the cache, but do not flush it to disk */
                if(H5D_istore_preempt(io_info, ent, FALSE) < 0)
                    HGOTO_ERROR(H5E_IO, H5E_CANTINIT, FAIL, "unable to preempt chunk")

                /* Break out of loop, chunk is evicted */
                break;
            } /* end if */
    } /* end for */

    /* Round up to the next integer # of chunks, to accomodate partial chunks */
    for(u = 0; u < rank; u++)
        chunks[u] = ((old_dims[u] + dset->shared->layout.u.chunk.dim[u]) - 1) / dset->shared->layout.u.chunk.dim[u];

    /* Get the "down" sizes for each dimension */
    if(H5V_array_down(rank, chunks, down_chunks) < 0)
        HGOTO_ERROR(H5E_IO, H5E_BADVALUE, FAIL, "can't compute 'down' sizes")

    /* Initialize the user data for the iteration */
    HDmemset(&udata, 0, sizeof udata);
    udata.common.mesg = &dset->shared->layout;
    udata.dims = curr_dims;
    udata.down_chunks = down_chunks;

    /* Initialize the skip list that will hold the chunks outside the dimensions */
    if(NULL == (udata.outside = H5SL_create(H5SL_TYPE_HSIZE, 0.5, (size_t)H5D_ISTORE_DEFAULT_SKIPLIST_HEIGHT)))
        HGOTO_ERROR(H5E_IO, H5E_CANTCREATE, FAIL, "can't create skip list for chunks outside new dimensions")

    /* Iterate over chunks in dataset, creating a list of chunks which are
     *  now completely outside the dataset's dimensions.
     *
     * Note: It would be more efficient to create a new B-tree routine that
     *          performed a "remove if" operation on the B-tree and remove all
     *          the chunks that were outside the dataset's dimensions through
     *          that routine.  However, that's a fair amount of work and it's
     *          unlikely that shrinking a dataset is a performance critical
     *          operation. - QAK
     */
    if(H5B_iterate(dset->oloc.file, io_info->dxpl_id, H5B_ISTORE, H5D_istore_prune_check, dset->shared->layout.u.chunk.addr, &udata) < 0)
	HGOTO_ERROR(H5E_IO, H5E_CANTINIT, FAIL, "unable to iterate over B-tree")

    /* Set up user data for skip list callback */
    rm_info.f = dset->oloc.file;
    rm_info.dxpl_id = io_info->dxpl_id;
    rm_info.mesg = &dset->shared->layout;

    /* Destroy the skip list, deleting the chunks in the callback */
    H5SL_destroy(udata.outside, H5D_istore_prune_remove, &rm_info);

    /* Reset any cached chunk info for this dataset */
    H5D_istore_cinfo_cache_reset(&dset->shared->cache.chunk.last);

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_prune_by_extent() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_remove
 *
 * Purpose:	Removes chunks that are no longer necessary in the B-tree.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer: Robb Matzke
 *             Pedro Vicente, pvn@ncsa.uiuc.edu
 * 		March 28, 2002
 *
 *-------------------------------------------------------------------------
 */
/* ARGSUSED */
static H5B_ins_t
H5D_istore_remove(H5F_t *f, hid_t dxpl_id, haddr_t addr, void *_lt_key /*in,out */ ,
	hbool_t *lt_key_changed /*out */ ,
	void UNUSED * _udata /*in,out */ ,
	void UNUSED * _rt_key /*in,out */ ,
	hbool_t *rt_key_changed /*out */ )
{
    H5D_istore_key_t    *lt_key = (H5D_istore_key_t *)_lt_key;
    H5B_ins_t ret_value=H5B_INS_REMOVE; /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5D_istore_remove)

    /* Remove raw data chunk from file */
    if(H5MF_xfree(f, H5FD_MEM_DRAW, dxpl_id, addr, (hsize_t)lt_key->nbytes)<0)
        HGOTO_ERROR(H5E_STORAGE, H5E_CANTFREE, H5B_INS_ERROR, "unable to free chunk")

    /* Mark keys as unchanged */
    *lt_key_changed = FALSE;
    *rt_key_changed = FALSE;

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_remove() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_initialize_by_extent
 *
 * Purpose:	This function searches for chunks that have to be initialized with the fill
 *		value both in the raw data cache and in the B-tree.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
 *		April 4, 2002
 *
 * Comments:
 *
 * (See the example of H5D_istore_prune_by_extent)
 * Next, there are seven chunks where the database extent boundary is
 * within the chunk. We find those seven just like we did with the previous nine.
 * Fot the ones that are allocated we initialize the part that lies outside the boundary
 * with the fill value.
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_initialize_by_extent(H5D_io_info_t *io_info)
{
    H5S_t           *space_chunk = NULL;	/* Dataspace for a chunk */
    H5D_t           *dset = io_info->dset;      /* Local pointer to the dataset info */
    const H5O_layout_t *layout = &(dset->shared->layout); /* Dataset's layout */
    H5D_storage_t    store;                     /* Dataset storage information */
    H5D_fill_buf_info_t fb_info;        /* Dataset's fill buffer info */
    hbool_t     fb_info_init = FALSE;   /* Whether the fill value buffer has been initialized */
    hsize_t          dset_dims[H5O_LAYOUT_NDIMS];   /* Current dataspace dimensions */
    hsize_t          chunk_dims[H5O_LAYOUT_NDIMS];  /* Current chunk dimensions */
    hsize_t          chunk_offset[H5O_LAYOUT_NDIMS]; /* Logical location of the chunks */
    hsize_t          hyper_start[H5O_LAYOUT_NDIMS];	/* Starting location of hyperslab */
    hsize_t          nchunks[H5O_LAYOUT_NDIMS]; /* Current number of chunks in each dimension */
    hsize_t          down_chunks[H5O_LAYOUT_NDIMS]; /* "down" size of number of elements in each dimension */
    hsize_t          bytes_per_chunk;	/* Bytes in chunk */
    int              srank;	        /* # of chunk dimensions (signed) */
    unsigned         rank;	        /* # of chunk dimensions */
    hbool_t	     carry;             /* Flag to indicate that chunk increment carrys to higher dimension (sorta) */
    unsigned         u;                 /* Local index variable */
    herr_t	     ret_value = SUCCEED;	/* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_initialize_by_extent, FAIL)

    /* Check args */
    HDassert(io_info);
    HDassert(io_info->dset && H5D_CHUNKED == layout->type);
    HDassert(layout->u.chunk.ndims > 0 && layout->u.chunk.ndims <= H5O_LAYOUT_NDIMS);
    HDassert(H5F_addr_defined(layout->u.chunk.addr));

    /* Go get the rank & dimensions */
    if((srank = H5S_get_simple_extent_dims(dset->shared->space, dset_dims, NULL)) < 0)
	HGOTO_ERROR(H5E_DATASET, H5E_CANTGET, FAIL, "can't get dataset dimensions")
    H5_ASSIGN_OVERFLOW(rank, srank, int, unsigned);

    /* Set size of lowest chunk dimension (the dataset element size) */
    dset_dims[rank] = layout->u.chunk.dim[rank];

    /* Compute the number of chunks in dataset & the # of bytes in a chunk */
    /* (round up to the next integer # of chunks, to accomodate partial chunks) */
    for(u = 0, bytes_per_chunk = layout->u.chunk.dim[rank]; u < rank; u++) {
        nchunks[u] = ((dset_dims[u] - 1) / layout->u.chunk.dim[u]) + 1;
        bytes_per_chunk *= layout->u.chunk.dim[u];
    } /* end for */

    /* Get the "down" sizes for each dimension */
    if(H5V_array_down(rank, nchunks, down_chunks) < 0)
        HGOTO_ERROR(H5E_INTERNAL, H5E_BADVALUE, FAIL, "can't compute 'down' sizes")

    /* Create a data space for a chunk & set the extent */
    for(u = 0; u < rank; u++)
	chunk_dims[u] = layout->u.chunk.dim[u];
    if(NULL == (space_chunk = H5S_create_simple(rank, chunk_dims, NULL)))
	HGOTO_ERROR(H5E_DATASPACE, H5E_CANTCREATE, FAIL, "can't create simple dataspace")

    /* Point to local dataset storage info */
    HDassert(io_info->store == NULL);       /* Make certain we aren't blowing anything away */
    io_info->store = &store;

    /* Reset hyperslab start array */
    HDmemset(hyper_start, 0, sizeof(hyper_start));

    /* Initialize current chunk offset to the origin (0, 0, 0, ...) */
    HDmemset(chunk_offset, 0, sizeof(chunk_offset));

    /* Loop over all chunks */
    carry = FALSE;
    while(!carry) {
        hbool_t found;	        /* Initialize this entry */
        int i;	                /* Local index variable */

	/*
	 * Figure out what chunks have to be initialized. These are the chunks where the dataspace
	 * extent boundary is within the chunk
	 */
        found = FALSE;
	for(u = 0; u < rank; u++)
	    if((chunk_offset[u] + layout->u.chunk.dim[u]) > dset_dims[u]) {
		found = TRUE;
		break;
	    } /* end if */

	if(found) {
            H5S_sel_iter_t chunk_iter;  /* Memory selection iteration info */
            hssize_t nelmts;            /* Number of data elements */
            hsize_t count[H5O_LAYOUT_NDIMS];	/* Element count of hyperslab */
            uint8_t *chunk;	        /* The file chunk  */
            unsigned idx_hint;	        /* Which chunk we're dealing with */
            hsize_t bytes_accessed;	/* Bytes accessed in chunk */

            /* Initialize the fill value buffer, if necessary */
            if(!fb_info_init) {
                if(H5D_fill_init(&fb_info, NULL, FALSE, NULL, NULL, NULL, NULL,
                        &dset->shared->dcpl_cache.fill,
                        dset->shared->type, dset->shared->type_id, (size_t)(bytes_per_chunk / layout->u.chunk.dim[rank]),
                        io_info->dxpl_cache->max_temp_buf, io_info->dxpl_id) < 0)
                    HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "can't initialize fill buffer info")
                fb_info_init = TRUE;
            } /* end if */

            /* Compute the # of elements to leave with existing value, in each dimension */
            for(u = 0; u < rank; u++)
		count[u] = MIN(layout->u.chunk.dim[u], (dset_dims[u] - chunk_offset[u]));

#ifdef H5D_ISTORE_DEBUG
	    HDfputs("cache:initialize:offset:[", stdout);
	    for(u = 0; u < rank; u++)
		HDfprintf(stdout, "%s%Hd", u ? ", " : "", chunk_offset[u]);
	    HDfputs("]", stdout);
	    HDfputs(":count:[", stdout);
	    for(u = 0; u < rank; u++)
		HDfprintf(stdout, "%s%Hd", u ? ", " : "", count[u]);
	    HDfputs("]\n", stdout);
#endif

            /* Select all elements in chunk, to begin with */
	    if(H5S_select_all(space_chunk, TRUE) < 0)
		HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to select space")

            /* "Subtract out" the elements to keep */
	    if(H5S_select_hyperslab(space_chunk, H5S_SELECT_NOTB, hyper_start, NULL, count, NULL) < 0)
		HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to select hyperslab")

            /* Calculate the index of this chunk */
            if(H5V_chunk_index(rank, chunk_offset, layout->u.chunk.dim, down_chunks, &store.chunk.index) < 0)
                HGOTO_ERROR(H5E_DATASPACE, H5E_BADRANGE, FAIL, "can't get chunk index")

            /* Lock the chunk into the cache, to get a pointer to the chunk buffer */
            store.chunk.offset = chunk_offset;
	    if(NULL == (chunk = (uint8_t *)H5D_istore_lock(io_info, NULL, FALSE, &idx_hint)))
		HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to read raw data chunk")


	    /* Fill the selection in the memory buffer */
            /* Use the size of the elements in the chunk directly instead of */
            /* relying on the fill.size, which might be set to 0 if there is */
            /* no fill-value defined for the dataset -QAK */

            /* Get the number of elements in the selection */
            nelmts = H5S_GET_SELECT_NPOINTS(space_chunk);
            HDassert(nelmts >= 0);
            H5_CHECK_OVERFLOW(nelmts, hssize_t, size_t);

            /* Check for VL datatype & non-default fill value */
            if(fb_info.has_vlen_fill_type)
                /* Re-fill the buffer to use for this I/O operation */
                if(H5D_fill_refill_vl(&fb_info, (size_t)nelmts, io_info->dxpl_id) < 0)
                    HGOTO_ERROR(H5E_DATASET, H5E_CANTCONVERT, FAIL, "can't refill fill value buffer")

            /* Create a selection iterator for scattering the elements to memory buffer */
            if(H5S_select_iter_init(&chunk_iter, space_chunk, layout->u.chunk.dim[rank]) < 0)
                HGOTO_ERROR(H5E_DATASET, H5E_CANTINIT, FAIL, "unable to initialize chunk selection information")

            /* Scatter the data into memory */
            if(H5D_select_mscat(fb_info.fill_buf, space_chunk, &chunk_iter, (size_t)nelmts, io_info->dxpl_cache, chunk/*out*/) < 0) {
                H5S_SELECT_ITER_RELEASE(&chunk_iter);
                HGOTO_ERROR(H5E_DATASET, H5E_READERROR, FAIL, "scatter failed")
            } /* end if */

            /* Release the selection iterator */
            if(H5S_SELECT_ITER_RELEASE(&chunk_iter) < 0)
                HGOTO_ERROR(H5E_DATASET, H5E_CANTFREE, FAIL, "Can't release selection iterator")


            /* The number of bytes accessed in the chunk */
            /* (i.e. the bytes replaced with fill values) */
            bytes_accessed = nelmts * layout->u.chunk.dim[rank];

            /* Release lock on chunk */
            H5_CHECK_OVERFLOW(bytes_accessed, hsize_t, size_t);
	    if(H5D_istore_unlock(io_info, TRUE, idx_hint, chunk, (size_t)bytes_accessed) < 0)
		HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "unable to unlock raw data chunk")
	} /* end if */

	/* Increment indices */
        carry = TRUE;
	for(i = (int)(rank - 1); i >= 0; --i) {
            chunk_offset[i] += layout->u.chunk.dim[i];
	    if(chunk_offset[i] >= dset_dims[i])
		chunk_offset[i] = 0;
	    else {
                carry = FALSE;
                break;
            } /* end else */
	} /* end for */
    } /* end while */

done:
    /* Release resources */
    if(space_chunk && H5S_close(space_chunk) < 0)
        HDONE_ERROR(H5E_DATASET, H5E_CLOSEERROR, FAIL, "unable to release dataspace")
    if(fb_info_init && H5D_fill_term(&fb_info) < 0)
        HDONE_ERROR(H5E_DATASET, H5E_CANTFREE, FAIL, "Can't release fill buffer info")

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_initialize_by_extent() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_delete
 *
 * Purpose:	Delete raw data storage for entire dataset (i.e. all chunks)
 *
 * Return:	Success:	Non-negative
 *		Failure:	negative
 *
 * Programmer:	Quincey Koziol
 *              Thursday, March 20, 2003
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_delete(H5F_t *f, hid_t dxpl_id, const H5O_layout_t *layout)
{
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_delete, FAIL)

    /* Check if the B-tree has been created in the file */
    if(H5F_addr_defined(layout->u.chunk.addr)) {
        H5O_layout_t tmp_layout=*layout;/* Local copy of layout info */
        H5D_istore_ud0_t	udata;  /* User data for B-tree iterator call */

        /* Set up user data for B-tree deletion */
        HDmemset(&udata, 0, sizeof udata);
        udata.mesg = &tmp_layout;

        /* Allocate the shared structure */
        if(H5D_istore_shared_create(f, &tmp_layout)<0)
            HGOTO_ERROR(H5E_RESOURCE, H5E_CANTINIT, FAIL, "can't create wrapper for shared B-tree info")

        /* Delete entire B-tree */
        if(H5B_delete(f, dxpl_id, H5B_ISTORE, tmp_layout.u.chunk.addr, &udata)<0)
            HGOTO_ERROR(H5E_IO, H5E_CANTDELETE, 0, "unable to delete chunk B-tree")

        /* Free the raw B-tree node buffer */
        if(tmp_layout.u.chunk.btree_shared==NULL)
            HGOTO_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "ref-counted page nil")
        if(H5RC_DEC(tmp_layout.u.chunk.btree_shared)<0)
            HGOTO_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "unable to decrement ref-counted page")
    } /* end if */

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_delete() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_update_cache
 *
 * Purpose:	Update any cached chunks index values after the dataspace
 *              size has changed
 *
 * Return:	Success:	Non-negative
 *		Failure:	negative
 *
 * Programmer:	Quincey Koziol
 *              Saturday, May 29, 2004
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_update_cache(H5D_t *dset, hid_t dxpl_id)
{
    H5D_io_info_t io_info;              /* Temporary I/O info object */
    H5D_rdcc_t         *rdcc = &(dset->shared->cache.chunk);	/*raw data chunk cache */
    H5D_rdcc_ent_t     *ent, *next;	/*cache entry  */
    H5D_rdcc_ent_t     *old_ent;	/* Old cache entry  */
    H5D_dxpl_cache_t _dxpl_cache;       /* Data transfer property cache buffer */
    H5D_dxpl_cache_t *dxpl_cache = &_dxpl_cache;   /* Data transfer property cache */
    unsigned            rank;	        /*current # of dimensions */
    hsize_t             curr_dims[H5O_LAYOUT_NDIMS];	/*current dataspace dimensions */
    hsize_t             chunks[H5O_LAYOUT_NDIMS];	        /*current number of chunks in each dimension */
    hsize_t             down_chunks[H5O_LAYOUT_NDIMS];   /* "down" size of number of elements in each dimension */
    unsigned            u;	        /*counters  */
    herr_t              ret_value = SUCCEED;      /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_update_cache, FAIL)

    /* Check args */
    HDassert(dset && H5D_CHUNKED == dset->shared->layout.type);
    HDassert(dset->shared->layout.u.chunk.ndims > 0 && dset->shared->layout.u.chunk.ndims <= H5O_LAYOUT_NDIMS);

    /* Get the rank */
    rank = dset->shared->layout.u.chunk.ndims-1;
    HDassert(rank > 0);

    /* 1-D dataset's chunks can't have their index change */
    if(rank == 1)
        HGOTO_DONE(SUCCEED)

    /* Go get the dimensions */
    if(H5S_get_simple_extent_dims(dset->shared->space, curr_dims, NULL) < 0)
	HGOTO_ERROR(H5E_DATASET, H5E_CANTGET, FAIL, "can't get dataset dimensions")

    /* Round up to the next integer # of chunks, to accomodate partial chunks */
    for(u = 0; u < rank; u++)
        chunks[u] = ((curr_dims[u] + dset->shared->layout.u.chunk.dim[u]) - 1) / dset->shared->layout.u.chunk.dim[u];

    /* Get the "down" sizes for each dimension */
    if(H5V_array_down(rank, chunks, down_chunks) < 0)
        HGOTO_ERROR(H5E_INTERNAL, H5E_BADVALUE, FAIL, "can't compute 'down' sizes")

    /* Fill the DXPL cache values for later use */
    if(H5D_get_dxpl_cache(dxpl_id, &dxpl_cache) < 0)
        HGOTO_ERROR(H5E_DATASET, H5E_CANTGET, FAIL, "can't fill dxpl cache")

    /* Construct dataset I/O info */
    H5D_BUILD_IO_INFO(&io_info, dset, dxpl_cache, dxpl_id, NULL);

    /* Recompute the index for each cached chunk that is in a dataset */
    for(ent = rdcc->head; ent; ent = next) {
        hsize_t             idx;        /* Chunk index */
        unsigned	    old_idx;	/* Previous index number	*/

        /* Get the pointer to the next cache entry */
        next = ent->next;

        /* Calculate the index of this chunk */
        if(H5V_chunk_index(rank,ent->offset,dset->shared->layout.u.chunk.dim,down_chunks,&idx)<0)
            HGOTO_ERROR(H5E_DATASPACE, H5E_BADRANGE, FAIL, "can't get chunk index")

        /* Compute the index for the chunk entry */
        old_idx=ent->idx;   /* Save for later */
        ent->idx=H5D_HASH(dset->shared,idx);

        if(old_idx != ent->idx) {
            /* Check if there is already a chunk at this chunk's new location */
            old_ent = rdcc->slot[ent->idx];
            if(old_ent != NULL) {
                HDassert(old_ent->locked == 0);

                /* Check if we are removing the entry we would walk to next */
                if(old_ent == next)
                    next = old_ent->next;

                /* Remove the old entry from the cache */
                if(H5D_istore_preempt(&io_info, old_ent, TRUE) < 0)
                    HGOTO_ERROR(H5E_IO, H5E_CANTFLUSH, FAIL, "unable to flush one or more raw data chunks")
            } /* end if */

            /* Insert this chunk into correct location in hash table */
            rdcc->slot[ent->idx] = ent;

            /* Null out previous location */
            rdcc->slot[old_idx] = NULL;
        } /* end if */
    } /* end for */

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5F_istore_update_cache() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_copy
 *
 * Purpose:	copy an indexed storage B-tree from SRC file to DST file.
 *
 * Return:	Non-negative on success (with the ISTORE argument initialized
 *		and ready to write to an object header). Negative on failure.
 *
 * Programmer:  Peter Cao
 *	        August 20, 2005
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_copy(H5F_t *f_src, H5O_layout_t *layout_src, H5F_t *f_dst,
    H5O_layout_t *layout_dst,  H5T_t *dt_src, H5O_copy_t *cpy_info, H5O_pline_t *pline, hid_t dxpl_id)
{
    H5D_istore_it_ud4_t    udata;
    H5T_path_t  *tpath_src_mem = NULL, *tpath_mem_dst = NULL;   /* Datatype conversion paths */
    hid_t       tid_src = -1;           /* Datatype ID for source datatype */
    hid_t       tid_dst = -1;           /* Datatype ID for destination datatype */
    hid_t       tid_mem = -1;           /* Datatype ID for memory datatype */
    size_t      buf_size;               /* Size of copy buffer */
    size_t      reclaim_buf_size;       /* Size of reclaim buffer */
    void       *buf = NULL;             /* Buffer for copying data */
    void       *bkg = NULL;             /* Buffer for background during type conversion */
    void       *reclaim_buf = NULL;     /* Buffer for reclaiming data */
    H5S_t      *buf_space = NULL;       /* Dataspace describing buffer */
    hid_t       sid_buf = -1;           /* ID for buffer dataspace */
    size_t      nelmts = 0;             /* Number of elements in buffer */
    hbool_t     do_convert = FALSE;     /* Indicate that type conversions should be performed */
    herr_t      ret_value = SUCCEED;    /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_copy, FAIL)

    /* Check args */
    HDassert(f_src);
    HDassert(f_dst);
    HDassert(layout_src && H5D_CHUNKED == layout_src->type);
    HDassert(layout_dst && H5D_CHUNKED == layout_dst->type);
    HDassert(dt_src);

    /* Create datatype ID for src datatype */
    if((tid_src = H5I_register(H5I_DATATYPE, dt_src)) < 0)
        HGOTO_ERROR(H5E_DATATYPE, H5E_CANTREGISTER, FAIL, "unable to register source file datatype")

    /* Create shared B-tree info for each file */
    if(H5D_istore_shared_create(f_src, layout_src) < 0)
        HGOTO_ERROR(H5E_RESOURCE, H5E_CANTINIT, FAIL, "can't create wrapper for shared B-tree info")
    if(H5D_istore_shared_create(f_dst, layout_dst) < 0)
        HGOTO_ERROR(H5E_RESOURCE, H5E_CANTINIT, FAIL, "can't create wrapper for shared B-tree info")

    /* Check if we need to create the B-tree in the dest. file */
    if(layout_dst->u.chunk.addr == HADDR_UNDEF) {
        /* Create the root of the B-tree that describes chunked storage */
        if(H5D_istore_create(f_dst, dxpl_id, layout_dst) < 0)
            HGOTO_ERROR(H5E_IO, H5E_CANTINIT, FAIL, "unable to initialize chunked storage")
    } /* end if */

    /* If there's a VLEN source datatype, set up type conversion information */
    if(H5T_detect_class(dt_src, H5T_VLEN) > 0) {
        H5T_t *dt_dst;              /* Destination datatype */
        H5T_t *dt_mem;              /* Memory datatype */
        size_t mem_dt_size;         /* Memory datatype size */
        size_t tmp_dt_size;         /* Temp. datatype size */
        size_t max_dt_size;         /* Max atatype size */
        hsize_t buf_dim;            /* Dimension for buffer */
        unsigned u;

        /* create a memory copy of the variable-length datatype */
        if(NULL == (dt_mem = H5T_copy(dt_src, H5T_COPY_TRANSIENT)))
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL, "unable to copy")
        if((tid_mem = H5I_register(H5I_DATATYPE, dt_mem)) < 0)
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTREGISTER, FAIL, "unable to register memory datatype")

        /* create variable-length datatype at the destinaton file */
        if(NULL == (dt_dst = H5T_copy(dt_src, H5T_COPY_TRANSIENT)))
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL, "unable to copy")
        if(H5T_set_loc(dt_dst, f_dst, H5T_LOC_DISK) < 0)
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL, "cannot mark datatype on disk")
        if((tid_dst = H5I_register(H5I_DATATYPE, dt_dst)) < 0)
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTREGISTER, FAIL, "unable to register destination file datatype")

        /* Set up the conversion functions */
        if(NULL == (tpath_src_mem = H5T_path_find(dt_src, dt_mem, NULL, NULL, dxpl_id, FALSE)))
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL, "unable to convert between src and mem datatypes")
        if(NULL == (tpath_mem_dst = H5T_path_find(dt_mem, dt_dst, NULL, NULL, dxpl_id, FALSE)))
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL, "unable to convert between mem and dst datatypes")

        /* Determine largest datatype size */
        if(0 == (max_dt_size = H5T_get_size(dt_src)))
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL, "unable to determine datatype size")
        if(0 == (mem_dt_size = H5T_get_size(dt_mem)))
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL, "unable to determine datatype size")
        max_dt_size = MAX(max_dt_size, mem_dt_size);
        if(0 == (tmp_dt_size = H5T_get_size(dt_dst)))
            HGOTO_ERROR(H5E_DATATYPE, H5E_CANTINIT, FAIL, "unable to determine datatype size")
        max_dt_size = MAX(max_dt_size, tmp_dt_size);

        /* Compute the number of elements per chunk */
        nelmts = 1;
        for(u = 0;  u < (layout_src->u.chunk.ndims - 1); u++)
            nelmts *= layout_src->u.chunk.dim[u];

        /* Create the space and set the initial extent */
        buf_dim = nelmts;
        if(NULL == (buf_space = H5S_create_simple((unsigned)1, &buf_dim, NULL)))
            HGOTO_ERROR(H5E_DATASPACE, H5E_CANTCREATE, FAIL, "can't create simple dataspace")

        /* Atomize */
        if((sid_buf = H5I_register(H5I_DATASPACE, buf_space)) < 0) {
            H5S_close(buf_space);
            HGOTO_ERROR(H5E_ATOM, H5E_CANTREGISTER, FAIL, "unable to register dataspace ID")
        } /* end if */

        /* Set initial buffer sizes */
        buf_size = nelmts * max_dt_size;
        reclaim_buf_size = nelmts * mem_dt_size;

        /* Allocate memory for reclaim buf */
        if(NULL == (reclaim_buf = H5MM_malloc(reclaim_buf_size)))
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed for raw data chunk")

        /* Indicate that type conversion should be performed */
        do_convert = TRUE;
    } /* end if */
    else {
        if(H5T_get_class(dt_src, FALSE) == H5T_REFERENCE) {
            /* Indicate that type conversion should be performed */
            do_convert = TRUE;
        } /* end if */

        buf_size = layout_src->u.chunk.size;
        reclaim_buf_size = 0;
    } /* end else */

    /* Set up conversion buffer, if appropriate */
    if(do_convert) {
        /* Allocate background memory for converting the chunk */
        if(NULL == (bkg = H5MM_malloc(buf_size)))
            HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed for raw data chunk")

        /* Check for reference datatype and no expanding references & clear background buffer */
        if(!cpy_info->expand_ref && 
                ((H5T_get_class(dt_src, FALSE) == H5T_REFERENCE) && (f_src != f_dst)))
            /* Reset value to zero */
            HDmemset(bkg, 0, buf_size);
    } /* end if */

    /* Allocate memory for copying the chunk */
    if(NULL == (buf = H5MM_malloc(buf_size)))
        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed for raw data chunk")

    /* Initialize the callback structure for the source */
    HDmemset(&udata, 0, sizeof udata);
    udata.common.mesg = layout_src;
    udata.file_src = f_src;
    udata.addr_dst = layout_dst->u.chunk.addr;
    udata.buf = buf;
    udata.bkg = bkg;
    udata.buf_size = buf_size;
    udata.tid_src = tid_src;
    udata.tid_mem = tid_mem;
    udata.tid_dst = tid_dst;
    udata.dt_src = dt_src;
    udata.do_convert = do_convert;
    udata.tpath_src_mem = tpath_src_mem;
    udata.tpath_mem_dst = tpath_mem_dst;
    udata.reclaim_buf = reclaim_buf;
    udata.reclaim_buf_size = reclaim_buf_size;
    udata.buf_space = buf_space;
    udata.nelmts = nelmts;
    udata.pline = pline;
    udata.file_dst = f_dst;
    udata.cpy_info = cpy_info;

    /* copy the chunked data by iteration */
    if(H5B_iterate(f_src, dxpl_id, H5B_ISTORE, H5D_istore_iter_copy, layout_src->u.chunk.addr, &udata) < 0)
        HGOTO_ERROR(H5E_IO, H5E_CANTINIT, FAIL, "unable to iterate over chunk B-tree")

    /* I/O buffers may have been re-allocated */
    buf = udata.buf;
    bkg = udata.bkg;

done:
    if(sid_buf > 0)
        if(H5I_dec_ref(sid_buf) < 0)
            HDONE_ERROR(H5E_DATASET, H5E_CANTFREE, FAIL, "Can't decrement temporary dataspace ID")
    if(tid_src > 0)
        if(H5I_dec_ref(tid_src) < 0)
            HDONE_ERROR(H5E_DATASET, H5E_CANTFREE, FAIL, "Can't decrement temporary datatype ID")
    if(tid_dst > 0)
        if(H5I_dec_ref(tid_dst) < 0)
            HDONE_ERROR(H5E_DATASET, H5E_CANTFREE, FAIL, "Can't decrement temporary datatype ID")
    if(tid_mem > 0)
        if(H5I_dec_ref(tid_mem) < 0)
            HDONE_ERROR(H5E_DATASET, H5E_CANTFREE, FAIL, "Can't decrement temporary datatype ID")
    if(buf)
        H5MM_xfree(buf);
    if(bkg)
        H5MM_xfree(bkg);
    if(reclaim_buf)
        H5MM_xfree(reclaim_buf);

    if(H5RC_DEC(layout_src->u.chunk.btree_shared) < 0)
        HDONE_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "unable to decrement ref-counted page")

    if(H5RC_DEC(layout_dst->u.chunk.btree_shared) < 0)
        HDONE_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "unable to decrement ref-counted page")

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_copy() */


/*-------------------------------------------------------------------------
 * Function:    H5D_istore_bh_size
 *
 * Purpose:     Retrieve the amount of B-tree storage for chunked dataset
 *
 * Return:      Success:        Non-negative
 *              Failure:        negative
 *
 * Programmer:  Vailin Choi
 *              June 8, 2007
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_bh_info(H5F_t *f, hid_t dxpl_id, H5O_layout_t *layout, hsize_t *btree_size)
{
    H5D_istore_it_ud1_t udata;                  /* User-data for loading istore nodes */
    H5B_info_ud_t       bh_udata;               /* User-data for B-tree size iteration */
    herr_t              ret_value = SUCCEED;

    FUNC_ENTER_NOAPI(H5D_istore_bh_info, FAIL)

    /* Check args */
    HDassert(f);
    HDassert(layout);
    HDassert(btree_size);

    /* Initialize the shared info for the B-tree traversal */
    if(H5D_istore_shared_create(f, layout) < 0)
        HGOTO_ERROR(H5E_RESOURCE, H5E_CANTINIT, FAIL, "can't create wrapper for shared B-tree info")

    /* Initialize istore node user-data */
    HDmemset(&udata, 0, sizeof udata);
    udata.common.mesg = layout;

    /* Iterate over B-tree, accumulating metadata size */
    bh_udata.udata = &udata;
    bh_udata.btree_size = btree_size;
    if(H5B_iterate_size(f, dxpl_id, H5B_ISTORE, NULL, layout->u.chunk.addr, &bh_udata) < 0)
        HGOTO_ERROR(H5E_BTREE, H5E_CANTINIT, FAIL, "unable to iterate over chunk B-tree")

done:
    if(layout->u.chunk.btree_shared == NULL)
	HDONE_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "ref-counted page nil")
    if(H5RC_DEC(layout->u.chunk.btree_shared) < 0)
	HDONE_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "unable to decrement ref-counted page")

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_bh_info() */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_dump_btree
 *
 * Purpose:	Prints information about the storage B-tree to the specified
 *		stream.
 *
 * Return:	Success:	Non-negative
 *
 *		Failure:	negative
 *
 * Programmer:	Robb Matzke
 *              Wednesday, April 28, 1999
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_dump_btree(H5F_t *f, hid_t dxpl_id, FILE *stream, unsigned ndims, haddr_t addr)
{
    H5O_layout_t        layout;
    H5D_istore_it_ud2_t	udata;
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_dump_btree, FAIL)

    HDmemset(&udata, 0, sizeof udata);
    layout.u.chunk.ndims = ndims;
    udata.common.mesg = &layout;
    udata.stream = stream;
    if(stream)
        HDfprintf(stream, "    Address: %a\n", addr);
    if(H5B_iterate(f, dxpl_id, H5B_ISTORE, H5D_istore_iter_dump, addr, &udata)<0)
        HGOTO_ERROR(H5E_IO, H5E_CANTINIT, 0, "unable to iterate over chunk B-tree")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_dump_btree() */

#ifdef H5D_ISTORE_DEBUG

/*-------------------------------------------------------------------------
 * Function:	H5D_istore_stats
 *
 * Purpose:	Print raw data cache statistics to the debug stream.  If
 *		HEADERS is non-zero then print table column headers,
 *		otherwise assume that the H5AC layer has already printed them.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, May 21, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_stats (H5D_t *dset, hbool_t headers)
{
    H5D_rdcc_t	*rdcc = &(dset->shared->cache.chunk);
    double	miss_rate;
    char	ascii[32];
    herr_t      ret_value=SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_stats, FAIL)

    if (!H5DEBUG(AC))
        HGOTO_DONE(SUCCEED)

    if (headers) {
        fprintf(H5DEBUG(AC), "H5D: raw data cache statistics\n");
        fprintf(H5DEBUG(AC), "   %-18s %8s %8s %8s %8s+%-8s\n",
            "Layer", "Hits", "Misses", "MissRate", "Inits", "Flushes");
        fprintf(H5DEBUG(AC), "   %-18s %8s %8s %8s %8s-%-8s\n",
            "-----", "----", "------", "--------", "-----", "-------");
    }

#ifdef H5AC_DEBUG
    if (H5DEBUG(AC)) headers = TRUE;
#endif

    if (headers) {
        if (rdcc->nhits>0 || rdcc->nmisses>0) {
            miss_rate = 100.0 * rdcc->nmisses /
                    (rdcc->nhits + rdcc->nmisses);
        } else {
            miss_rate = 0.0;
        }
        if (miss_rate > 100) {
            sprintf(ascii, "%7d%%", (int) (miss_rate + 0.5));
        } else {
            sprintf(ascii, "%7.2f%%", miss_rate);
        }

        fprintf(H5DEBUG(AC), "   %-18s %8u %8u %7s %8d+%-9ld\n",
            "raw data chunks", rdcc->nhits, rdcc->nmisses, ascii,
            rdcc->ninits, (long)(rdcc->nflushes)-(long)(rdcc->ninits));
    }

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_stats() */
#endif /* H5D_ISTORE_DEBUG */


/*-------------------------------------------------------------------------
 * Function:	H5D_istore_debug
 *
 * Purpose:	Debugs a B-tree node for indexed raw data storage.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *              Thursday, April 16, 1998
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5D_istore_debug(H5F_t *f, hid_t dxpl_id, haddr_t addr, FILE * stream, int indent,
		 int fwidth, unsigned ndims)
{
    H5O_layout_t        layout;
    H5D_istore_ud0_t	udata;          /* B-tree user data */
    herr_t      ret_value = SUCCEED;    /* Return value */

    FUNC_ENTER_NOAPI(H5D_istore_debug,FAIL)

    layout.u.chunk.ndims = ndims;

    /* Allocate the shared structure */
    if(H5D_istore_shared_create(f, &layout)<0)
	HGOTO_ERROR(H5E_RESOURCE, H5E_CANTINIT, FAIL, "can't create wrapper for shared B-tree info")

    /* Set up B-tree user data */
    HDmemset(&udata, 0, sizeof udata);
    udata.mesg = &layout;

    (void)H5B_debug(f, dxpl_id, addr, stream, indent, fwidth, H5B_ISTORE, &udata);

    /* Free the raw B-tree node buffer */
    if(layout.u.chunk.btree_shared==NULL)
        HGOTO_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "ref-counted page nil")
    if(H5RC_DEC(layout.u.chunk.btree_shared)<0)
	HGOTO_ERROR(H5E_IO, H5E_CANTFREE, FAIL, "unable to decrement ref-counted page")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5D_istore_debug() */

