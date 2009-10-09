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
 * Created:		H5Ocache.c
 *			Sep 28 2005
 *			Quincey Koziol <koziol@ncsa.uiuc.edu>
 *
 * Purpose:		Object header metadata cache virtual functions.
 *
 *-------------------------------------------------------------------------
 */

/****************/
/* Module Setup */
/****************/

#define H5O_PACKAGE		/*suppress error about including H5Opkg	  */


/***********/
/* Headers */
/***********/
#include "H5private.h"		/* Generic Functions			*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5FLprivate.h"	/* Free lists                           */
#include "H5MFprivate.h"	/* File memory management		*/
#include "H5Opkg.h"             /* Object headers			*/


/****************/
/* Local Macros */
/****************/

/* Set the object header size to speculatively read in */
/* (needs to be more than the object header prefix size to work at all and
 *      should be larger than the largest object type's default object header
 *      size to save the extra I/O operations) */
#define H5O_SPEC_READ_SIZE 512


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
static H5O_t *H5O_load(H5F_t *f, hid_t dxpl_id, haddr_t addr, const void *_udata1,
		       void *_udata2);
static herr_t H5O_flush(H5F_t *f, hid_t dxpl_id, hbool_t destroy, haddr_t addr, H5O_t *oh, unsigned UNUSED * flags_ptr);
static herr_t H5O_clear(H5F_t *f, H5O_t *oh, hbool_t destroy);
static herr_t H5O_size(const H5F_t *f, const H5O_t *oh, size_t *size_ptr);


/*********************/
/* Package Variables */
/*********************/


/*****************************/
/* Library Private Variables */
/*****************************/

/* Declare external the free list for H5O_unknown_t's */
H5FL_EXTERN(H5O_unknown_t);


/*******************/
/* Local Variables */
/*******************/

/* H5O inherits cache-like properties from H5AC */
const H5AC_class_t H5AC_OHDR[1] = {{
    H5AC_OHDR_ID,
    (H5AC_load_func_t)H5O_load,
    (H5AC_flush_func_t)H5O_flush,
    (H5AC_dest_func_t)H5O_dest,
    (H5AC_clear_func_t)H5O_clear,
    (H5AC_size_func_t)H5O_size,
}};


/*-------------------------------------------------------------------------
 * Function:	H5O_load
 *
 * Purpose:	Loads an object header from disk.
 *
 * Return:	Success:	Pointer to the new object header.
 *
 *		Failure:	NULL
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Aug  5 1997
 *
 *-------------------------------------------------------------------------
 */
static H5O_t *
H5O_load(H5F_t *f, hid_t dxpl_id, haddr_t addr, const void UNUSED * _udata1,
	 void UNUSED * _udata2)
{
    H5O_t	*oh = NULL;     /* Object header read in */
    uint8_t     read_buf[H5O_SPEC_READ_SIZE];       /* Buffer for speculative read */
    const uint8_t *p;           /* Pointer into buffer to decode */
    size_t	spec_read_size; /* Size of buffer to speculatively read in */
    size_t	prefix_size;    /* Size of object header prefix */
    unsigned	nmesgs;         /* Total # of messages in this object header */
    unsigned	curmesg = 0;    /* Current message being decoded in object header */
    unsigned    merged_null_msgs = 0;   /* Number of null messages merged together */
    haddr_t	chunk_addr;     /* Address of first chunk */
    size_t	chunk_size;     /* Size of first chunk */
    haddr_t     eoa;		/* Relative end of file address	*/
    H5O_t	*ret_value;     /* Return value */

    FUNC_ENTER_NOAPI(H5O_load, NULL)

    /* check args */
    HDassert(f);
    HDassert(H5F_addr_defined(addr));
    HDassert(!_udata1);
    HDassert(!_udata2);

    /* Make certain we don't speculatively read off the end of the file */
    if(HADDR_UNDEF == (eoa = H5F_get_eoa(f, H5FD_MEM_OHDR)))
        HGOTO_ERROR(H5E_OHDR, H5E_CANTGET, NULL, "unable to determine file size")

    /* Compute the size of the speculative object header buffer */
    H5_ASSIGN_OVERFLOW(spec_read_size, MIN(eoa - addr, H5O_SPEC_READ_SIZE), /* From: */ hsize_t, /* To: */ size_t);

    /* Attempt to speculatively read both object header prefix and first chunk */
    if(H5F_block_read(f, H5FD_MEM_OHDR, addr, spec_read_size, dxpl_id, read_buf) < 0)
	HGOTO_ERROR(H5E_OHDR, H5E_READERROR, NULL, "unable to read object header")
    p = read_buf;

    /* allocate ohdr and init chunk list */
    if(NULL == (oh = H5FL_CALLOC(H5O_t)))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")

    /* File-specific, non-stored information */
    oh->sizeof_size = H5F_SIZEOF_SIZE(f);
    oh->sizeof_addr = H5F_SIZEOF_ADDR(f);

    /* Check for magic number */
    /* (indicates version 2 or later) */
    if(!HDmemcmp(p, H5O_HDR_MAGIC, (size_t)H5_SIZEOF_MAGIC)) {
        /* Magic number */
        p += H5_SIZEOF_MAGIC;

        /* Version */
        oh->version = *p++;
        if(H5O_VERSION_2 != oh->version)
            HGOTO_ERROR(H5E_OHDR, H5E_VERSION, NULL, "bad object header version number")

        /* Flags */
        oh->flags = *p++;
        if(oh->flags & ~H5O_HDR_ALL_FLAGS)
            HGOTO_ERROR(H5E_OHDR, H5E_BADVALUE, NULL, "unknown object header status flag(s)")

        /* Number of messages (to allocate initially) */
        nmesgs = 1;

        /* Number of links to object (unless overridden by refcount message) */
        oh->nlink = 1;

        /* Time fields */
        if(oh->flags & H5O_HDR_STORE_TIMES) {
            UINT32DECODE(p, oh->atime);
            UINT32DECODE(p, oh->mtime);
            UINT32DECODE(p, oh->ctime);
            UINT32DECODE(p, oh->btime);
        } /* end if */
        else
            oh->atime = oh->mtime = oh->ctime = oh->btime = 0;

        /* Attribute fields */
        if(oh->flags & H5O_HDR_ATTR_STORE_PHASE_CHANGE) {
            UINT16DECODE(p, oh->max_compact);
            UINT16DECODE(p, oh->min_dense);
            if(oh->max_compact < oh->min_dense)
                HGOTO_ERROR(H5E_OHDR, H5E_VERSION, NULL, "bad object header attribute phase change values")
        } /* end if */
        else {
            oh->max_compact = H5O_CRT_ATTR_MAX_COMPACT_DEF;
            oh->min_dense = H5O_CRT_ATTR_MIN_DENSE_DEF;
        } /* end else */

        /* First chunk size */
        switch(oh->flags & H5O_HDR_CHUNK0_SIZE) {
            case 0:     /* 1 byte size */
                chunk_size = *p++;
                break;

            case 1:     /* 2 byte size */
                UINT16DECODE(p, chunk_size);
                break;

            case 2:     /* 4 byte size */
                UINT32DECODE(p, chunk_size);
                break;

            case 3:     /* 8 byte size */
                UINT64DECODE(p, chunk_size);
                break;

            default:
                HGOTO_ERROR(H5E_OHDR, H5E_BADVALUE, NULL, "bad size for chunk 0")
        } /* end switch */
        if(chunk_size > 0 && chunk_size < H5O_SIZEOF_MSGHDR_OH(oh))
            HGOTO_ERROR(H5E_OHDR, H5E_VERSION, NULL, "bad object header chunk size")
    } /* end if */
    else {
        /* Version */
        oh->version = *p++;
        if(H5O_VERSION_1 != oh->version)
            HGOTO_ERROR(H5E_OHDR, H5E_VERSION, NULL, "bad object header version number")

        /* Flags */
        oh->flags = H5O_CRT_OHDR_FLAGS_DEF;

        /* Reserved */
        p++;

        /* Number of messages */
        UINT16DECODE(p, nmesgs);

        /* Link count */
        UINT32DECODE(p, oh->nlink);

        /* Reset unused time fields */
        oh->atime = oh->mtime = oh->ctime = oh->btime = 0;

        /* Reset unused attribute fields */
        oh->max_compact = 0;
        oh->min_dense = 0;

        /* First chunk size */
        UINT32DECODE(p, chunk_size);
        if((nmesgs > 0 && chunk_size < H5O_SIZEOF_MSGHDR_OH(oh)) ||
                (nmesgs == 0 && chunk_size > 0))
            HGOTO_ERROR(H5E_OHDR, H5E_VERSION, NULL, "bad object header chunk size")

        /* Reserved, in version 1 */
        p += 4;
    } /* end else */

    /* Determine object header prefix length */
    prefix_size = (size_t)(p - read_buf);
    HDassert((size_t)prefix_size == (size_t)(H5O_SIZEOF_HDR(oh) - H5O_SIZEOF_CHKSUM_OH(oh)));

    /* Compute first chunk address */
    chunk_addr = addr + (hsize_t)prefix_size;

    /* Allocate the message array */
    oh->alloc_nmesgs = (nmesgs > 0) ? nmesgs : 1;
    if(NULL == (oh->mesg = H5FL_SEQ_MALLOC(H5O_mesg_t, oh->alloc_nmesgs)))
	HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")

    /* Read each chunk from disk */
    while(H5F_addr_defined(chunk_addr)) {
        unsigned chunkno;       /* Current chunk's index */
#ifndef NDEBUG
        unsigned nullcnt;       /* Count of null messages (for sanity checking gaps in chunks) */
#endif /* NDEBUG */
        uint8_t *eom_ptr;       /* Pointer to end of messages for a chunk */

	/* Increase chunk array size, if necessary */
	if(oh->nchunks >= oh->alloc_nchunks) {
	    unsigned na = MAX(H5O_NCHUNKS, oh->alloc_nchunks * 2);        /* Double # of chunks allocated */
	    H5O_chunk_t *x = H5FL_SEQ_REALLOC(H5O_chunk_t, oh->chunk, (size_t)na);

	    if(!x)
                HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")
	    oh->alloc_nchunks = na;
	    oh->chunk = x;
	} /* end if */

	/* Init the chunk data info */
	chunkno = oh->nchunks++;
	oh->chunk[chunkno].dirty = FALSE;
	oh->chunk[chunkno].gap = 0;
        if(chunkno == 0) {
            /* First chunk's 'image' includes room for the object header prefix */
            oh->chunk[0].addr = addr;
            oh->chunk[0].size = chunk_size + H5O_SIZEOF_HDR(oh);
        } /* end if */
        else {
            oh->chunk[chunkno].addr = chunk_addr;
            oh->chunk[chunkno].size = chunk_size;
        } /* end else */
	if(NULL == (oh->chunk[chunkno].image = H5FL_BLK_MALLOC(chunk_image, oh->chunk[chunkno].size)))
	    HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")

        /* Handle chunk 0 as special case */
        if(chunkno == 0) {
            /* Check for speculative read of first chunk containing all the data needed */
            if(spec_read_size >= oh->chunk[0].size)
                HDmemcpy(oh->chunk[0].image, read_buf, oh->chunk[0].size);
            else {
                /* Copy the object header prefix into chunk 0's image */
                HDmemcpy(oh->chunk[0].image, read_buf, prefix_size);

                /* Read the chunk raw data */
                /* (probably re-reads some data we already retrieved, but since
                 *      we have to do the I/O operation anyway, we might as
                 *      well avoid memcpy()ing the data in our buffer already)
                 */
                if(H5F_block_read(f, H5FD_MEM_OHDR, chunk_addr, (oh->chunk[0].size - prefix_size),
                        dxpl_id, (oh->chunk[0].image + prefix_size)) < 0)
                    HGOTO_ERROR(H5E_OHDR, H5E_READERROR, NULL, "unable to read object header data")
            } /* end else */

            /* Point into chunk image to decode */
            p = oh->chunk[0].image + prefix_size;
        } /* end if */
        else {
            /* Read the chunk raw data */
            if(H5F_block_read(f, H5FD_MEM_OHDR, chunk_addr, chunk_size, dxpl_id, oh->chunk[chunkno].image) < 0)
                HGOTO_ERROR(H5E_OHDR, H5E_READERROR, NULL, "unable to read object header data")

            /* Point into chunk image to decode */
            p = oh->chunk[chunkno].image;
        } /* end else */

        /* Check for magic # on chunks > 0 in later versions of the format */
        if(chunkno > 0 && oh->version > H5O_VERSION_1) {
            /* Magic number */
            if(HDmemcmp(p, H5O_CHK_MAGIC, (size_t)H5_SIZEOF_MAGIC))
                HGOTO_ERROR(H5E_OHDR, H5E_CANTLOAD, NULL, "wrong object header chunk signature")
            p += H5_SIZEOF_MAGIC;
        } /* end if */

	/* Decode messages from this chunk */
        eom_ptr = oh->chunk[chunkno].image + (oh->chunk[chunkno].size - H5O_SIZEOF_CHKSUM_OH(oh));
#ifndef NDEBUG
        nullcnt = 0;
#endif /* NDEBUG */
	while(p < eom_ptr) {
            unsigned    mesgno;         /* Current message to operate on */
            size_t	mesg_size;      /* Size of message read in */
            unsigned	id;             /* ID (type) of current message */
            uint8_t	flags;          /* Flags for current message */
            H5O_msg_crt_idx_t crt_idx = 0;  /* Creation index for current message */

            /* Decode message prefix info */

            /* Version # */
            if(oh->version == H5O_VERSION_1)
                UINT16DECODE(p, id)
            else
                id = *p++;

            /* Check for unknown message ID getting encoded in file */
            if(id == H5O_UNKNOWN_ID)
                HGOTO_ERROR(H5E_OHDR, H5E_CANTLOAD, NULL, "'unknown' message ID encoded in file?!?")

            /* Message size */
	    UINT16DECODE(p, mesg_size);
	    HDassert(mesg_size == H5O_ALIGN_OH(oh, mesg_size));

            /* Message flags */
	    flags = *p++;
            if(flags & ~H5O_MSG_FLAG_BITS)
                HGOTO_ERROR(H5E_OHDR, H5E_CANTLOAD, NULL, "unknown flag for message")
            if((flags & H5O_MSG_FLAG_SHARED) && (flags & H5O_MSG_FLAG_DONTSHARE))
                HGOTO_ERROR(H5E_OHDR, H5E_CANTLOAD, NULL, "bad flag combination for message")
            if((flags & H5O_MSG_FLAG_WAS_UNKNOWN) && (flags & H5O_MSG_FLAG_FAIL_IF_UNKNOWN))
                HGOTO_ERROR(H5E_OHDR, H5E_CANTLOAD, NULL, "bad flag combination for message")
            if((flags & H5O_MSG_FLAG_WAS_UNKNOWN) && !(flags & H5O_MSG_FLAG_MARK_IF_UNKNOWN))
                HGOTO_ERROR(H5E_OHDR, H5E_CANTLOAD, NULL, "bad flag combination for message")

            /* Reserved bytes/creation index */
            if(oh->version == H5O_VERSION_1)
                p += 3; /*reserved*/
            else {
                /* Only encode creation index if they are being tracked */
                if(oh->flags & H5O_HDR_ATTR_CRT_ORDER_TRACKED)
                    UINT16DECODE(p, crt_idx);
            } /* end else */

            /* Try to detect invalidly formatted object header message that
             *  extends past end of chunk.
             */
	    if(p + mesg_size > eom_ptr)
		HGOTO_ERROR(H5E_OHDR, H5E_CANTINIT, NULL, "corrupt object header")

#ifndef NDEBUG
            /* Increment count of null messages */
            if(H5O_NULL_ID == id)
                nullcnt++;
#endif /* NDEBUG */

            /* Check for combining two adjacent 'null' messages */
            if((H5F_INTENT(f) & H5F_ACC_RDWR) &&
                    H5O_NULL_ID == id && oh->nmesgs > 0 &&
                    H5O_NULL_ID == oh->mesg[oh->nmesgs - 1].type->id &&
                    oh->mesg[oh->nmesgs - 1].chunkno == chunkno) {

                /* Combine adjacent null messages */
                mesgno = oh->nmesgs - 1;
                oh->mesg[mesgno].raw_size += H5O_SIZEOF_MSGHDR_OH(oh) + mesg_size;
                oh->mesg[mesgno].dirty = TRUE;
                merged_null_msgs++;
            } /* end if */
            else {
                /* Check if we need to extend message table to hold the new message */
                if(oh->nmesgs >= oh->alloc_nmesgs)
                    if(H5O_alloc_msgs(oh, (size_t)1) < 0)
                        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "can't allocate more space for messages")

                /* Get index for message */
                mesgno = oh->nmesgs++;

                /* Initialize information about message */
                oh->mesg[mesgno].dirty = FALSE;
                oh->mesg[mesgno].flags = flags;
                oh->mesg[mesgno].crt_idx = crt_idx;
                oh->mesg[mesgno].native = NULL;
                oh->mesg[mesgno].raw = (uint8_t *)p;        /* Casting away const OK - QAK */
                oh->mesg[mesgno].raw_size = mesg_size;
                oh->mesg[mesgno].chunkno = chunkno;

                /* Point unknown messages at 'unknown' message class */
                /* (Usually from future versions of the library) */
                if(id >= NELMTS(H5O_msg_class_g) || NULL == H5O_msg_class_g[id]) {
                    H5O_unknown_t *unknown;     /* Pointer to "unknown" message info */

                    /* Allocate "unknown" message info */
                    if(NULL == (unknown = H5FL_MALLOC(H5O_unknown_t)))
                        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed")

                    /* Save the original message type ID */
                    *unknown = id;

                    /* Save 'native' form of continuation message */
                    oh->mesg[mesgno].native = unknown;

                    /* Set message to "unknown" class */
                    oh->mesg[mesgno].type = H5O_msg_class_g[H5O_UNKNOWN_ID];

                    /* Check for "fail if unknown" message flag */
                    if(flags & H5O_MSG_FLAG_FAIL_IF_UNKNOWN)
                        HGOTO_ERROR(H5E_OHDR, H5E_BADMESG, NULL, "unknown message with 'fail if unknown' flag found")
                    /* Check for "mark if unknown" message flag, etc. */
                    else if((flags & H5O_MSG_FLAG_MARK_IF_UNKNOWN) &&
                            !(flags & H5O_MSG_FLAG_WAS_UNKNOWN) &&
                            (H5F_INTENT(f) & H5F_ACC_RDWR)) {

                        /* Mark the message as "unknown" */
                        /* This is a bit aggressive, since the application may
                         * never change anything about the object (metadata or
                         * raw data), but we can sort out the finer details
                         * when/if we start using the flag - QAK
                         */
                        /* Also, it's possible that this functionality may not
                         * get invoked if the object header is brought into
                         * the metadata cache in some other "weird" way, like
                         * using H5Ocopy() - QAK
                         */
                        oh->mesg[mesgno].flags |= H5O_MSG_FLAG_WAS_UNKNOWN;

                        /* Mark the message and object header as dirty */
                        oh->mesg[mesgno].dirty = TRUE;
                        oh->cache_info.is_dirty = TRUE;
                    } /* end if */
                } /* end if */
                else
                    /* Set message class for "known" messages */
                    oh->mesg[mesgno].type = H5O_msg_class_g[id];
            } /* end else */

            /* Advance decode pointer past message */
            p += mesg_size;

            /* Check for 'gap' at end of chunk */
            if((eom_ptr - p) > 0 && (eom_ptr - p) < H5O_SIZEOF_MSGHDR_OH(oh)) {
                /* Gaps can only occur in later versions of the format */
                HDassert(oh->version > H5O_VERSION_1);

                /* Gaps should only occur in chunks with no null messages */
                HDassert(nullcnt == 0);

                /* Set gap information for chunk */
                oh->chunk[chunkno].gap = (eom_ptr - p);

                /* Increment location in chunk */
                p += oh->chunk[chunkno].gap;
            } /* end if */
	} /* end while */

        /* Check for correct checksum on chunks, in later versions of the format */
        if(oh->version > H5O_VERSION_1) {
            uint32_t stored_chksum;     /* Checksum from file */
            uint32_t computed_chksum;   /* Checksum computed in memory */

            /* Metadata checksum */
            UINT32DECODE(p, stored_chksum);

            /* Compute checksum on chunk */
            computed_chksum = H5_checksum_metadata(oh->chunk[chunkno].image, (oh->chunk[chunkno].size - H5O_SIZEOF_CHKSUM), 0);

            /* Verify checksum */
            if(stored_chksum != computed_chksum)
                HGOTO_ERROR(H5E_OHDR, H5E_BADVALUE, NULL, "incorrect metadata checksum for object header chunk")
        } /* end if */

        /* Sanity check */
        HDassert(p == oh->chunk[chunkno].image + oh->chunk[chunkno].size);

        /* Check for another chunk to read in & parse */
        for(chunk_addr = HADDR_UNDEF; !H5F_addr_defined(chunk_addr) && curmesg < oh->nmesgs; ++curmesg) {
            /* Check if next message to examine is a continuation message */
            if(H5O_CONT_ID == oh->mesg[curmesg].type->id) {
                H5O_cont_t *cont;
                unsigned ioflags = 0;   /* Flags for decode routine */

                /* Decode continuation message */
                cont = (H5O_cont_t *)(H5O_MSG_CONT->decode)(f, dxpl_id, NULL, 0, &ioflags, oh->mesg[curmesg].raw);
                cont->chunkno = oh->nchunks;	/*the next chunk to allocate */

                /* Save 'native' form of continuation message */
                oh->mesg[curmesg].native = cont;

                /* Set up to read in next chunk */
                chunk_addr = cont->addr;
                chunk_size = cont->size;

                /* Mark the object header as dirty if the message was changed by decoding */
                if((ioflags & H5O_DECODEIO_DIRTY) && (H5F_get_intent(f) & H5F_ACC_RDWR)) {
                    oh->mesg[curmesg].dirty = TRUE;
                    oh->cache_info.is_dirty = TRUE;
	            }
            } /* end if */
            /* Check if next message to examine is a ref. count message */
            else if(H5O_REFCOUNT_ID == oh->mesg[curmesg].type->id) {
                H5O_refcount_t *refcount;
                unsigned ioflags = 0;   /* Flags for decode routine */

                /* Decode ref. count message */
                HDassert(oh->version > H5O_VERSION_1);
                refcount = (H5O_refcount_t *)(H5O_MSG_REFCOUNT->decode)(f, dxpl_id, NULL, 0, &ioflags, oh->mesg[curmesg].raw);

                /* Save 'native' form of ref. count message */
                oh->mesg[curmesg].native = refcount;

                /* Set object header values */
                oh->has_refcount_msg = TRUE;
                oh->nlink = *refcount;

                /* Mark the object header as dirty if the message was changed by decoding */
                if((ioflags & H5O_DECODEIO_DIRTY) && (H5F_get_intent(f) & H5F_ACC_RDWR)) {
                    oh->mesg[curmesg].dirty = TRUE;
                    oh->cache_info.is_dirty = TRUE;
	            }
            } /* end if */
            /* Check if next message to examine is a link message */
            else if(H5O_LINK_ID == oh->mesg[curmesg].type->id) {
                /* Increment the count of link messages */
                oh->link_msgs_seen++;
            } /* end if */
            /* Check if next message to examine is an attribute message */
            else if(H5O_ATTR_ID == oh->mesg[curmesg].type->id) {
                /* Increment the count of attribute messages */
                oh->attr_msgs_seen++;
            } /* end if */
        } /* end for */
    } /* end while */

    /* Mark the object header dirty if we've merged a message */
    if(merged_null_msgs)
	oh->cache_info.is_dirty = TRUE;

/* Don't check for the incorrect # of object header messages bug unless we've
 * enabled strict format checking.  This allows for older files, created with
 * a version of the library that had a bug in tracking the correct # of header
 * messages to be read in without the library fussing about things. -QAK
 */
#ifdef H5_STRICT_FORMAT_CHECKS
    /* Sanity check for the correct # of messages in object header */
    if(oh->version == H5O_VERSION_1)
        if((oh->nmesgs + merged_null_msgs) != nmesgs)
            HGOTO_ERROR(H5E_OHDR, H5E_CANTLOAD, NULL, "corrupt object header - too few messages")
#else /* H5_STRICT_FORMAT_CHECKS */
    /* Check for incorrect # of messages in object header and if we have write
     * access on the file, flag the object header as dirty, so it gets fixed.
     */
    if(oh->version == H5O_VERSION_1)
        if((oh->nmesgs + merged_null_msgs) != nmesgs &&
                (H5F_get_intent(f) & H5F_ACC_RDWR))
            oh->cache_info.is_dirty = TRUE;
#endif /* H5_STRICT_FORMAT_CHECKS */

#ifdef H5O_DEBUG
H5O_assert(oh);
#endif /* H5O_DEBUG */

    /* Set return value */
    ret_value = oh;

done:
    /* Release the [possibly partially initialized] object header on errors */
    if(!ret_value && oh)
        if(H5O_dest(f, oh) < 0)
	    HDONE_ERROR(H5E_OHDR, H5E_CANTFREE, NULL, "unable to destroy object header data")

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5O_load() */


/*-------------------------------------------------------------------------
 * Function:	H5O_flush
 *
 * Purpose:	Flushes (and destroys) an object header.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		matzke@llnl.gov
 *		Aug  5 1997
 *
 * Changes:     JRM -- 8/21/06
 *              Added the flags_ptr parameter.  This parameter exists to
 *              allow the flush routine to report to the cache if the
 *              entry is resized or renamed as a result of the flush.
 *              *flags_ptr is set to H5C_CALLBACK__NO_FLAGS_SET on entry.
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5O_flush(H5F_t *f, hid_t dxpl_id, hbool_t destroy, haddr_t UNUSED addr, H5O_t *oh, unsigned UNUSED * flags_ptr)
{
    herr_t      ret_value = SUCCEED;       /* Return value */

    FUNC_ENTER_NOAPI(H5O_flush, FAIL)

    /* check args */
    HDassert(f);
    HDassert(H5F_addr_defined(addr));
    HDassert(oh);

    /* flush */
    if(oh->cache_info.is_dirty) {
        uint8_t	*p;             /* Pointer to object header prefix buffer */
        unsigned u;             /* Local index variable */

#ifdef H5O_DEBUG
H5O_assert(oh);
#endif /* H5O_DEBUG */

        /* Point to raw data 'image' for first chunk, which has room for the prefix */
        p = oh->chunk[0].image;

        /* Later versions of object header prefix have different format and
         * also require that chunk 0 always be updated, since the checksum
         * on the entire block of memory needs to be updated if anything is
         * modified */
        if(oh->version > H5O_VERSION_1) {
            uint64_t chunk0_size = oh->chunk[0].size - H5O_SIZEOF_HDR(oh);  /* Size of chunk 0's data */

            /* Verify magic number */
            HDassert(!HDmemcmp(p, H5O_HDR_MAGIC, H5_SIZEOF_MAGIC));
            p += H5_SIZEOF_MAGIC;

            /* Version */
            *p++ = oh->version;

            /* Flags */
            *p++ = oh->flags;

            /* Time fields */
            if(oh->flags & H5O_HDR_STORE_TIMES) {
                UINT32ENCODE(p, oh->atime);
                UINT32ENCODE(p, oh->mtime);
                UINT32ENCODE(p, oh->ctime);
                UINT32ENCODE(p, oh->btime);
            } /* end if */

            /* Attribute fields */
            if(oh->flags & H5O_HDR_ATTR_STORE_PHASE_CHANGE) {
                UINT16ENCODE(p, oh->max_compact);
                UINT16ENCODE(p, oh->min_dense);
            } /* end if */

            /* First chunk size */
            switch(oh->flags & H5O_HDR_CHUNK0_SIZE) {
                case 0:     /* 1 byte size */
                    HDassert(chunk0_size < 256);
                    *p++ = (uint8_t)chunk0_size;
                    break;

                case 1:     /* 2 byte size */
                    HDassert(chunk0_size < 65536);
                    UINT16ENCODE(p, chunk0_size);
                    break;

                case 2:     /* 4 byte size */
		    /* use <= 2**32 -1 to stay within 4 bytes integer range */
                    HDassert(chunk0_size <= 4294967295UL);
                    UINT32ENCODE(p, chunk0_size);
                    break;

                case 3:     /* 8 byte size */
                    UINT64ENCODE(p, chunk0_size);
                    break;

                default:
                    HGOTO_ERROR(H5E_OHDR, H5E_BADVALUE, FAIL, "bad size for chunk 0")
            } /* end switch */
        } /* end if */
        else {
            /* Version */
            *p++ = oh->version;

            /* Reserved */
            *p++ = 0;

            /* Number of messages */
#ifdef H5O_ENABLE_BAD_MESG_COUNT
            if(oh->store_bad_mesg_count)
                UINT16ENCODE(p, (oh->nmesgs - 1))
            else
#endif /* H5O_ENABLE_BAD_MESG_COUNT */
                UINT16ENCODE(p, oh->nmesgs);

            /* Link count */
            UINT32ENCODE(p, oh->nlink);

            /* First chunk size */
            UINT32ENCODE(p, (oh->chunk[0].size - H5O_SIZEOF_HDR(oh)));

            /* Zero to alignment */
            HDmemset(p, 0, (size_t)(H5O_SIZEOF_HDR(oh) - 12));
            p += (size_t)(H5O_SIZEOF_HDR(oh) - 12);
        } /* end else */
        HDassert((size_t)(p - oh->chunk[0].image) == (size_t)(H5O_SIZEOF_HDR(oh) - H5O_SIZEOF_CHKSUM_OH(oh)));

        /* Mark chunk 0 as dirty, since the object header prefix has been updated */
        /* (this could be more sophisticated and track whether any prefix fields
         *      have been changed, which could save I/O accesses if the
         *      messages in chunk 0 haven't changed - QAK)
         */
        HDassert(H5F_addr_eq(addr, oh->chunk[0].addr));
        oh->chunk[0].dirty = TRUE;

	/* Encode any dirty messages */
        if(H5O_flush_msgs(f, oh) < 0)
            HGOTO_ERROR(H5E_OHDR, H5E_CANTFLUSH, FAIL, "unable to flush object header messages")

	/* Write each chunk to disk, if it's dirty */
	for(u = 0; u < oh->nchunks; u++) {
            /* Sanity checks */
            if(oh->version > H5O_VERSION_1)
                /* Make certain the magic # is present */
                HDassert(!HDmemcmp(oh->chunk[u].image, (u == 0 ? H5O_HDR_MAGIC : H5O_CHK_MAGIC), H5_SIZEOF_MAGIC));
            else
                /* Gaps should never occur in version 1 of the format */
                HDassert(oh->chunk[u].gap == 0);

            /* Write out chunk, if it's dirty */
	    if(oh->chunk[u].dirty) {
                /* Extra work, for later versions of the format */
                if(oh->version > H5O_VERSION_1) {
                    uint32_t metadata_chksum;   /* Computed metadata checksum value */

                    /* Check for gap in chunk & zero it out */
                    if(oh->chunk[u].gap)
                        HDmemset((oh->chunk[u].image + oh->chunk[u].size) -
                            (H5O_SIZEOF_CHKSUM + oh->chunk[u].gap), 0, oh->chunk[u].gap);

                    /* Compute metadata checksum */
                    metadata_chksum = H5_checksum_metadata(oh->chunk[u].image, (oh->chunk[u].size - H5O_SIZEOF_CHKSUM), 0);

                    /* Metadata checksum */
                    p = oh->chunk[u].image + (oh->chunk[u].size - H5O_SIZEOF_CHKSUM);
                    UINT32ENCODE(p, metadata_chksum);
                } /* end if */

                /* Write the chunk out */
                HDassert(H5F_addr_defined(oh->chunk[u].addr));
                if(H5F_block_write(f, H5FD_MEM_OHDR, oh->chunk[u].addr,
                            oh->chunk[u].size, dxpl_id, oh->chunk[u].image) < 0)
                    HGOTO_ERROR(H5E_OHDR, H5E_WRITEERROR, FAIL, "unable to write object header chunk to disk")

                /* Mark chunk as clean now */
                oh->chunk[u].dirty = FALSE;
	    } /* end if */
	} /* end for */

        /* Mark object header as clean now */
	oh->cache_info.is_dirty = FALSE;
    } /* end if */

    /* Destroy the object header, if requested */
    if(destroy)
        if(H5O_dest(f,oh) < 0)
	    HGOTO_ERROR(H5E_OHDR, H5E_CANTFREE, FAIL, "unable to destroy object header data")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5O_flush() */


/*-------------------------------------------------------------------------
 * Function:	H5O_dest
 *
 * Purpose:	Destroys an object header.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *		koziol@ncsa.uiuc.edu
 *		Jan 15 2003
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5O_dest(H5F_t *f, H5O_t *oh)
{
    unsigned	u;                      /* Local index variable */
    herr_t ret_value = SUCCEED;         /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5O_dest)
#ifdef QAK
HDfprintf(stderr, "%s: oh->cache_info.addr = %a\n", FUNC, oh->cache_info.addr);
HDfprintf(stderr, "%s: oh->cache_info.free_file_space_on_destroy = %t\n", FUNC, oh->cache_info.free_file_space_on_destroy);
#endif /* QAK */

    /* check args */
    HDassert(oh);

    /* Verify that node is clean */
    HDassert(!oh->cache_info.is_dirty);

    /* If we're going to free the space on disk, the address must be valid */
    HDassert(!oh->cache_info.free_file_space_on_destroy || H5F_addr_defined(oh->cache_info.addr));

    /* destroy chunks */
    if(oh->chunk) {
        /* Check for releasing file space for object header */
        if(oh->cache_info.free_file_space_on_destroy) {
            /* Free main (first) object header "chunk" */
            /* (XXX: Nasty usage of internal DXPL value! -QAK) */
            if(H5MF_xfree(f, H5FD_MEM_OHDR, H5AC_dxpl_id, oh->chunk[0].addr, (hsize_t)oh->chunk[0].size) < 0)
                HGOTO_ERROR(H5E_OHDR, H5E_CANTFREE, FAIL, "unable to free object header")
        } /* end if */

        /* Release buffer for each chunk */
        for(u = 0; u < oh->nchunks; u++) {
            /* Verify that chunk is clean */
            HDassert(oh->chunk[u].dirty == 0);

            oh->chunk[u].image = H5FL_BLK_FREE(chunk_image, oh->chunk[u].image);
        } /* end for */

        /* Release array of chunk info */
        oh->chunk = (H5O_chunk_t *)H5FL_SEQ_FREE(H5O_chunk_t, oh->chunk);
    } /* end if */

    /* destroy messages */
    if(oh->mesg) {
        for(u = 0; u < oh->nmesgs; u++) {
            /* Verify that message is clean, unless it could have been marked
             * dirty by decoding */
#ifndef NDEBUG
            if(oh->ndecode_dirtied && oh->mesg[u].dirty)
                oh->ndecode_dirtied--;
            else
                HDassert(oh->mesg[u].dirty == 0);
#endif /* NDEBUG */

            H5O_msg_free_mesg(&oh->mesg[u]);
        } /* end for */

        /* Make sure we accounted for all the messages dirtied by decoding */
        HDassert(!oh->ndecode_dirtied);

        oh->mesg = (H5O_mesg_t *)H5FL_SEQ_FREE(H5O_mesg_t, oh->mesg);
    } /* end if */

    /* destroy object header */
    (void)H5FL_FREE(H5O_t, oh);

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5O_dest() */


/*-------------------------------------------------------------------------
 * Function:	H5O_clear
 *
 * Purpose:	Mark a object header in memory as non-dirty.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *		koziol@ncsa.uiuc.edu
 *		Mar 20 2003
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5O_clear(H5F_t *f, H5O_t *oh, hbool_t destroy)
{
    unsigned	u;      /* Local index variable */
    herr_t ret_value = SUCCEED;

    FUNC_ENTER_NOAPI_NOINIT(H5O_clear)

    /* check args */
    HDassert(oh);

    /* Mark chunks as clean */
    for(u = 0; u < oh->nchunks; u++)
        oh->chunk[u].dirty = FALSE;

    /* Mark messages as clean */
    for(u = 0; u < oh->nmesgs; u++)
        oh->mesg[u].dirty = FALSE;

#ifndef NDEBUG
        /* Reset the number of messages dirtied by decoding */
        oh->ndecode_dirtied = 0;
#endif /* NDEBUG */

    /* Mark whole header as clean */
    oh->cache_info.is_dirty = FALSE;

    if(destroy)
        if(H5O_dest(f, oh) < 0)
            HGOTO_ERROR(H5E_OHDR, H5E_CANTFREE, FAIL, "unable to destroy object header data")

done:
    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5O_clear() */


/*-------------------------------------------------------------------------
 * Function:	H5O_size
 *
 * Purpose:	Compute the size in bytes of the specified instance of
 *              H5O_t on disk, and return it in *len_ptr.  On failure,
 *              the value of *len_ptr is undefined.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	John Mainzer
 *		5/13/04
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5O_size(const H5F_t UNUSED *f, const H5O_t *oh, size_t *size_ptr)
{
    size_t	size;           /* Running sum of the object header's size */
    unsigned	u;              /* Local index variable */

    FUNC_ENTER_NOAPI_NOINIT_NOFUNC(H5O_size)

    /* check args */
    HDassert(oh);
    HDassert(size_ptr);

    /* Add sizes of all the chunks */
    /* (includes size of prefix, in chunk 0) */
    for(u = 0, size = 0; u < oh->nchunks; u++)
        size += oh->chunk[u].size;

    *size_ptr = size;

    FUNC_LEAVE_NOAPI(SUCCEED)
} /* H5O_size() */

