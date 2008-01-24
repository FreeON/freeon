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

/*-------------------------------------------------------------------------
 *
 * Created:             hdf5cache.c
 *                      Jul  9 1997
 *                      Robb Matzke <matzke@llnl.gov>
 *
 * Purpose:             Functions in this file implement a cache for
 *                      things which exist on disk.  All "things" associated
 *                      with a particular HDF file share the same cache; each
 *                      HDF file has it's own cache.
 *
 * Modifications:
 *
 *      Robb Matzke, 4 Aug 1997
 *      Added calls to H5E.
 *
 *      Quincey Koziol, 22 Apr 2000
 *      Turned on "H5AC_SORT_BY_ADDR"
 *
 *-------------------------------------------------------------------------
 */

#define H5F_PACKAGE		/*suppress error about including H5Fpkg	  */

#include "H5private.h"
#include "H5ACprivate.h"
#include "H5Dprivate.h"         /* Datasets */
#include "H5Eprivate.h"
#include "H5Fpkg.h"
#include "H5FLprivate.h"	/*Free Lists                              */
#include "H5Iprivate.h"         /* IDs */
#include "H5MMprivate.h"
#include "H5Pprivate.h"         /* Property lists */

/*
 * The MPIO & MPIPOSIX drivers are needed because there are places where we
 * check for the parallel I/O transfer mode.
 */
#include "H5FDmpio.h"
#include "H5FDmpiposix.h"

/*
 * Private file-scope variables.
 */
#define PABLO_MASK      H5AC_mask

/* Interface initialization */
static int             interface_initialize_g = 0;
#define INTERFACE_INIT H5AC_init_interface
static herr_t H5AC_init_interface(void);

/* Default dataset transfer property list for metadata I/O calls */
/* (Collective set, "block before metadata write" set and "library internal" set) */
/* (Global variable definition, declaration is in H5ACprivate.h also) */
hid_t H5AC_dxpl_id=(-1);

/* Private dataset transfer property list for metadata I/O calls */
/* (Collective set and "library internal" set) */
/* (Static variable definition) */
static hid_t H5AC_noblock_dxpl_id=(-1);

/* Dataset transfer property list for independent metadata I/O calls */
/* (just "library internal" set - i.e. independent transfer mode) */
/* (Global variable definition, declaration is in H5ACprivate.h also) */
hid_t H5AC_ind_dxpl_id=(-1);

static H5AC_t          *current_cache_g = NULL;         /*for sorting */

/* Declare a free list to manage the H5AC_t struct */
H5FL_DEFINE_STATIC(H5AC_t);

/* Declare a PQ free list to manage the cache mapping array information */
H5FL_ARR_DEFINE_STATIC(int,-1);

/* Declare a PQ free list to manage the cache slot array information */
H5FL_ARR_DEFINE_STATIC(H5AC_info_ptr_t,-1);

#ifdef H5AC_DEBUG
/* Declare a PQ free list to manage the protected slot array information */
H5FL_ARR_DEFINE_STATIC(H5AC_prot_t,-1);
#endif /* H5AC_DEBUG */


/*-------------------------------------------------------------------------
 * Function:	H5AC_init
 *
 * Purpose:	Initialize the interface from some other layer.
 *
 * Return:	Success:	non-negative
 *
 *		Failure:	negative
 *
 * Programmer:	Quincey Koziol
 *              Saturday, January 18, 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5AC_init(void)
{
    FUNC_ENTER(H5AC_init, FAIL);
    /* FUNC_ENTER() does all the work */
    FUNC_LEAVE(SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5AC_init_interface
 *
 * Purpose:	Initialize interface-specific information
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 18, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5AC_init_interface(void)
{
    H5P_t *new_plist;           /* New dataset transfer property list created */
#ifdef H5_HAVE_PARALLEL
    H5D_xfer_t	*xfer_parms;    /* Pointer to copy of default dataset transfer property */
#endif /* H5_HAVE_PARALLEL */

    FUNC_ENTER(H5AC_init_interface, FAIL);

    /* Copy the default dataset transfer properties for use in flush calls */
    if (NULL==(new_plist=H5P_copy(H5P_DATASET_XFER, &H5D_xfer_dflt)))
        HRETURN_ERROR(H5E_CACHE, H5E_CANTCOPY, FAIL, "unable to copy default dataset transfer property list");

#ifdef H5_HAVE_PARALLEL
    /* Set the [private] flag to indicate parallel I/O needs to block before a metadata write */
    xfer_parms = &(new_plist->u.dxfer);
    xfer_parms->block_before_meta_write=1;      /* Should be enabled by default... */
    xfer_parms->xfer_mode = H5FD_MPIO_COLLECTIVE;
    xfer_parms->library_internal = 1;           /* This is a special, library internal property list */
#endif /* H5_HAVE_PARALLEL */

    /* Get an ID for the H5AC dxpl */
    if ((H5AC_dxpl_id=H5P_create(H5P_DATASET_XFER, new_plist)) < 0)
        HRETURN_ERROR(H5E_CACHE, H5E_CANTREGISTER, FAIL, "unable to register property list");

    /* Copy the default dataset transfer properties for use in flush calls */
    if (NULL==(new_plist=H5P_copy(H5P_DATASET_XFER, &H5D_xfer_dflt)))
        HRETURN_ERROR(H5E_CACHE, H5E_CANTCOPY, FAIL, "unable to copy default dataset transfer property list");

#ifdef H5_HAVE_PARALLEL
    /* Reset the [private] flag to indicate parallel I/O doesn't need to block before a metadata write */
    xfer_parms = &(new_plist->u.dxfer);
    xfer_parms->block_before_meta_write=0;
    xfer_parms->xfer_mode = H5FD_MPIO_COLLECTIVE;
    xfer_parms->library_internal = 1;           /* This is a special, library internal property list */
#endif /* H5_HAVE_PARALLEL */

    /* Get an ID for the non-blocking, collective H5AC dxpl */
    if ((H5AC_noblock_dxpl_id=H5P_create(H5P_DATASET_XFER, new_plist)) < 0)
        HRETURN_ERROR(H5E_CACHE, H5E_CANTREGISTER, FAIL, "unable to register property list");

    /* Copy the default dataset transfer properties for use in independent calls */
    if (NULL==(new_plist=H5P_copy(H5P_DATASET_XFER, &H5D_xfer_dflt)))
        HRETURN_ERROR(H5E_CACHE, H5E_CANTCOPY, FAIL, "unable to copy default dataset transfer property list");

#ifdef H5_HAVE_PARALLEL
    /* Reset the [private] flag to indicate parallel I/O doesn't need to block before a metadata write */
    xfer_parms = &(new_plist->u.dxfer);
    xfer_parms->block_before_meta_write=0;
    xfer_parms->xfer_mode = H5FD_MPIO_INDEPENDENT;
    xfer_parms->library_internal = 1;           /* This is a special, library internal property list */
#endif /* H5_HAVE_PARALLEL */

    /* Get an ID for the non-blocking, independent H5AC dxpl */
    if ((H5AC_ind_dxpl_id=H5P_create(H5P_DATASET_XFER, new_plist)) < 0)
        HRETURN_ERROR(H5E_CACHE, H5E_CANTREGISTER, FAIL, "unable to register property list");

    FUNC_LEAVE(SUCCEED);
} /* end H5AC_init_interface() */


/*-------------------------------------------------------------------------
 * Function:	H5AC_term_interface
 *
 * Purpose:	Terminate this interface.
 *
 * Return:	Success:	Positive if anything was done that might
 *				affect other interfaces; zero otherwise.
 *
 * 		Failure:	Negative.
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 18, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
H5AC_term_interface(void)
{
    int		n=0;

    if (interface_initialize_g) {
        if(H5AC_dxpl_id>0 || H5AC_noblock_dxpl_id>0 || H5AC_ind_dxpl_id>0) {
            /* Indicate more work to do */
            n = 1; /* H5I */

            /* Close H5AC dxpl */
            if (H5I_dec_ref(H5AC_dxpl_id ) < 0
                    || H5I_dec_ref(H5AC_noblock_dxpl_id)<0
                    || H5I_dec_ref(H5AC_ind_dxpl_id)<0)
                H5E_clear(); /*ignore the error*/
            else {
                /* Reset static IDs */
                H5AC_dxpl_id=(-1);
                H5AC_noblock_dxpl_id=(-1);
                H5AC_ind_dxpl_id=(-1);

                /* Reset interface initialization flag */
                interface_initialize_g = 0;
            } /* end else */
        } /* end if */
        else
            /* Reset interface initialization flag */
            interface_initialize_g = 0;
    } /* end if */

    return(n);
} /* end H5AC_term_interface() */


/*-------------------------------------------------------------------------
 * Function:    H5AC_create
 *
 * Purpose:     Initialize the cache just after a file is opened.  The
 *              SIZE_HINT is the number of cache slots desired.  If you
 *              pass an invalid value then H5AC_NSLOTS is used.  You can
 *              turn off caching by using 1 for the SIZE_HINT value.
 *
 * Return:      Success:        Number of slots actually used.
 *
 *              Failure:        Negative
 *
 * Programmer:  Robb Matzke
 *              matzke@llnl.gov
 *              Jul  9 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
H5AC_create(H5F_t *f, int size_hint)
{
    H5AC_t                 *cache = NULL;
    FUNC_ENTER(H5AC_create, FAIL);

    assert(f);
    assert(NULL == f->shared->cache);
    if (size_hint < 1) size_hint = H5AC_NSLOTS;

    if (NULL==(f->shared->cache = cache = H5FL_ALLOC(H5AC_t,1))) {
	HRETURN_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL,
		       "memory allocation failed");
    }
    cache->nslots = size_hint;
    cache->slot = H5FL_ARR_ALLOC(H5AC_info_ptr_t,(hsize_t)cache->nslots,1);
    if (NULL==cache->slot) {
        f->shared->cache = H5FL_FREE (H5AC_t,f->shared->cache);
        HRETURN_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
    }
    cache->dslot = H5FL_ARR_ALLOC(H5AC_info_ptr_t,(hsize_t)cache->nslots,1);
    if (NULL==cache->dslot) {
        cache->slot = H5FL_ARR_FREE (H5AC_info_ptr_t,cache->slot);
        f->shared->cache = H5FL_FREE (H5AC_t,f->shared->cache);
        HRETURN_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
    } /* end if */
#ifdef H5AC_DEBUG
    if ((cache->prot = H5FL_ARR_ALLOC(H5AC_prot_t,(hsize_t)cache->nslots,1))==NULL) {
        cache->dslot = H5FL_ARR_FREE (H5AC_info_ptr_t,cache->dslot);
        cache->slot = H5FL_ARR_FREE (H5AC_info_ptr_t,cache->slot);
        f->shared->cache = H5FL_FREE (H5AC_t,f->shared->cache);
        HRETURN_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");
    }
#endif /* H5AC_DEBUG */

    FUNC_LEAVE(size_hint);
}

/*-------------------------------------------------------------------------
 * Function:    H5AC_dest
 *
 * Purpose:     Flushes all data to disk and destroys the cache.
 *              This function fails if any object are protected since the
 *              resulting file might not be consistent.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *              matzke@llnl.gov
 *              Jul  9 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5AC_dest(H5F_t *f, hid_t dxpl_id)
{
    H5AC_t                 *cache = NULL;
    FUNC_ENTER(H5AC_dest, FAIL);

    assert(f);
    assert(f->shared->cache);
    cache = f->shared->cache;

    if (H5AC_flush(f, dxpl_id, NULL, HADDR_UNDEF, TRUE) < 0) {
        HRETURN_ERROR(H5E_CACHE, H5E_CANTFLUSH, FAIL,
                      "unable to flush cache");
    }
#ifdef H5AC_DEBUG
    {
        unsigned i;
        for (i=0; i<cache->nslots; i++) {
            cache->prot[i].slot = H5MM_xfree(cache->prot[i].slot);
            cache->prot[i].aprots = 0;
            cache->prot[i].nprots = 0;
        }
        cache->prot = H5FL_ARR_FREE(H5AC_prot_t,cache->prot);
    }
#endif

    cache->dslot = H5FL_ARR_FREE (H5AC_info_ptr_t,cache->dslot);
    cache->slot = H5FL_ARR_FREE(H5AC_info_ptr_t,cache->slot);
    cache->nslots = 0;
    f->shared->cache = cache = H5FL_FREE(H5AC_t,cache);
    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:    H5AC_find
 *
 * Purpose:     Given an object type and the address at which that object
 *              is located in the file, return a pointer to the object.
 *              The optional UDATA1 and UDATA2 structures are passed down to
 *              the function that is responsible for loading the object into
 *              memory.
 *
 *              The returned pointer is guaranteed to be valid until the next
 *              call to an H5AC function (if you want a pointer which is valid
 *              indefinately then see H5AC_protect()).
 *
 *              If H5AC_DEBUG is defined then this function also
 *              checks that the requested object is not currently
 *              protected since it is illegal to modify a protected object
 *              except through the pointer returned by H5AC_protect().
 *
 * Return:      Success:        Pointer to the object.  The pointer is
 *                              valid until some other cache function
 *                              is called.
 *
 *              Failure:        NULL
 *
 * Programmer:  Robb Matzke
 *              matzke@llnl.gov
 *              Jul  9 1997
 *
 * Modifications:
 *
 *      Robb Matzke, 4 Aug 1997
 *      Fails immediately if the cached object is at the correct address
 *      but is of the wrong type.  This happens if the caller doesn't know
 *      what type of object is at the address and calls this function with
 *      various type identifiers until one succeeds (cf., the debugger).
 *
 *      Robb Matzke, 30 Oct 1997
 *      Keeps track of hits, misses, and flushes per object type so we have
 *      some cache performance diagnostics.
 *
 * 	Robb Matzke, 1999-07-27
 *	The ADDR argument is passed by value.
 *
 *-------------------------------------------------------------------------
 */
void *
H5AC_find(H5F_t *f, hid_t dxpl_id, const H5AC_class_t *type, haddr_t addr,
	    const void *udata1, void *udata2)
{
    unsigned                idx;
    void                   *thing = NULL;
    H5AC_flush_func_t       flush;
    H5AC_info_t           **info = NULL;
    H5AC_t                 *cache = NULL;
#ifdef H5_HAVE_PARALLEL
    H5AC_info_t           **dinfo = NULL;
#endif /* H5_HAVE_PARALLEL */

    FUNC_ENTER(H5AC_find, NULL);

    assert(f);
    assert(f->shared->cache);
    assert(type);
    assert(type->load);
    assert(type->flush);
    assert(H5F_addr_defined(addr));

    /* Get local pointers to the file's cache information */
    idx = H5AC_HASH(f, addr);
    cache = f->shared->cache;
    info = cache->slot + idx;

#ifdef H5_HAVE_PARALLEL
    /* If MPIO or MPIPOSIX is used, do special parallel I/O actions */
    if(IS_H5FD_MPIO(f) || IS_H5FD_MPIPOSIX(f)) {
        H5AC_dest_func_t        dest;

        /* Get local pointer to file's dirty cache information */
        dinfo = cache->dslot + idx;

        /* Check if the cache has 'held' information for this cache slot */
        if (*dinfo) {
            /* Sanity check that the 'clean' item is really clean */
            assert(*info);
            assert((*info)->dirty==0);

            /* Destroy 'current' information */
            dest = (*info)->type->dest;
            if ((dest)(f, (*info))<0)
                HRETURN_ERROR(H5E_CACHE, H5E_CANTFREE, NULL, "unable to free cached object");

            /* Restore 'held' information back to 'current' information */
            (*info)=(*dinfo);

            /* Clear 'held' information */
            (*dinfo)=NULL;

            cache->diagnostics[type->id].nrestores++;
        } /* end if */
    } /* end if */
#endif /* H5_HAVE_PARALLEL */

    /*
     * Return right away if the item is in the cache.
     */
    if ((*info) && H5F_addr_eq((*info)->addr, addr) && (*info)->type==type) {
        cache->diagnostics[type->id].nhits++;
        HRETURN(*info);
    }

    cache->diagnostics[type->id].nmisses++;

#ifdef H5AC_DEBUG
    /*
     * Check that the requested thing isn't protected, for protected things
     * can only be modified through the pointer already handed out by the
     * H5AC_protect() function.
     */
    {
        H5AC_prot_t            *prot = NULL;
        int                    i;

        prot = cache->prot + idx;
        for (i = 0; i < prot->nprots; i++) {
            assert(H5F_addr_ne(addr, prot->slot[i]->addr));
        }
    }
#endif /* H5AC_DEBUG */

    /*
     * Load a new thing.  If it can't be loaded, then return an error
     * without preempting anything.
     */
    if (NULL == (thing = (type->load)(f, dxpl_id, addr, udata1, udata2)))
        HRETURN_ERROR(H5E_CACHE, H5E_CANTLOAD, NULL, "unable to load object");

#ifdef H5_HAVE_PARALLEL
    /* If MPIO or MPIPOSIX is used, do special parallel I/O actions */
    if(IS_H5FD_MPIO(f) || IS_H5FD_MPIPOSIX(f)) {
        const H5D_xfer_t	*xfer_parms;

        /* Get the dataset transfer property list */
        if (H5P_DEFAULT == dxpl_id) {
            xfer_parms = &H5D_xfer_dflt;
        } else if (H5P_DATASET_XFER != H5P_get_class(dxpl_id) ||
                   NULL == (xfer_parms = H5I_object(dxpl_id))) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, NULL, "not xfer parms");
        }

        /* Make certain there is no 'held' info for this slot */
        assert((*dinfo)==NULL);

        /* Must be using collective I/O to flush metadata in parallel */
        if(xfer_parms->xfer_mode==H5FD_MPIO_INDEPENDENT) {
            /* Check if there is dirty metadata in this slot */
            if((*info) && (*info)->dirty) {
                /* 'Hold' the current metadata for later */
                (*dinfo)=(*info);

                /* Reset the 'current' metadata, so it doesn't get flushed */
                (*info)=NULL;

                cache->diagnostics[(*dinfo)->type->id].nholds++;
            } /* end if */
        } /* end else */
    } /* end if */
#endif /* H5_HAVE_PARALLEL */

    /*
     * Flush & destroy the previous cache entry if there is one.
     */
    if (*info) {
        H5AC_subid_t type_id=(*info)->type->id;  /* Remember this for later */

        flush = (*info)->type->flush;
        if ( (flush)(f, dxpl_id, TRUE, (*info)->addr, (*info)) < 0) {
            /*
             * The old thing could not be removed from the stack.
             * Release the new thing and fail.
             */
            if ((type->dest)(f, thing) < 0)
                HRETURN_ERROR(H5E_CACHE, H5E_CANTFLUSH, NULL, "unable to flush just-loaded object");
            HRETURN_ERROR(H5E_CACHE, H5E_CANTFLUSH, NULL, "unable to flush existing cached object");
        } /* end if */
        cache->diagnostics[type_id].nflushes++;
    } /* end if */

    /*
     * Make the cache point to the new thing.
     */
    (*info)=thing;
    (*info)->type = type;
    (*info)->addr = addr;
    assert((*info)->dirty==0);  /* Should be clean after being loaded */

    FUNC_LEAVE(thing);
}

/*-------------------------------------------------------------------------
 * Function:    H5AC_compare
 *
 * Purpose:     Compare two hash entries by address.  Unused entries are
 *              all equal to one another and greater than all used entries.
 *
 * Return:      Success:        -1, 0, 1
 *
 *              Failure:        never fails
 *
 * Programmer:  Robb Matzke
 *              matzke@llnl.gov
 *              Aug 12 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static int
H5AC_compare(const void *_a, const void *_b)
{
    int                    a = *((const int *) _a);
    int                    b = *((const int *) _b);

    assert(current_cache_g);

    if(NULL==current_cache_g->slot[a] || NULL == current_cache_g->slot[b]) {
        if(NULL==current_cache_g->slot[a]) {
            if (NULL == current_cache_g->slot[b]) {
                return 0;
            } else {
                return -1;
            }
        }
        else {
            return 1;
        }
    }
    else if (NULL == current_cache_g->slot[a]->type) {
        if (NULL == current_cache_g->slot[b]->type) {
            return 0;
        } else {
            return -1;
        }
    } else if (NULL == current_cache_g->slot[b]->type) {
        return 1;
    } else if (current_cache_g->slot[a]->addr < current_cache_g->slot[b]->addr) {
        return -1;
    } else if (current_cache_g->slot[a]->addr > current_cache_g->slot[b]->addr) {
        return 1;
    }
    return 0;
}

/*-------------------------------------------------------------------------
 * Function:    H5AC_flush
 *
 * Purpose:     Flushes (and destroys if DESTROY is non-zero) the specified
 *              entry from the cache.  If the entry TYPE is CACHE_FREE and
 *              ADDR is HADDR_UNDEF then all types of entries are
 *              flushed. If TYPE is CACHE_FREE and ADDR is defined then
 *              whatever is cached at ADDR is flushed.  Otherwise the thing
 *              at ADDR is flushed if it is the correct type.
 *
 *              If there are protected objects they will not be flushed.
 *              However, an attempt will be made to flush all non-protected
 *              items before this function returns failure.
 *
 * Return:      Non-negative on success/Negative on failure if there was a
 *              request to flush all items and something was protected.
 *
 * Programmer:  Robb Matzke
 *              matzke@llnl.gov
 *              Jul  9 1997
 *
 * Modifications:
 * 		Robb Matzke, 1999-07-27
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
herr_t
H5AC_flush(H5F_t *f, hid_t dxpl_id, const H5AC_class_t *type, haddr_t addr, hbool_t destroy)
{
    unsigned                   i;
    herr_t                  status;
    H5AC_flush_func_t       flush=NULL;
    H5AC_info_t           **info;
    int                   *map = NULL;
    unsigned                   nslots;
    H5AC_t                 *cache = NULL;

    FUNC_ENTER(H5AC_flush, FAIL);

    assert(f);
    assert(f->shared->cache);

    /* Get local information */
    cache = f->shared->cache;

    if (!H5F_addr_defined(addr)) {
        unsigned first_flush=1;     /* Indicate if this is the first flush */

        /*
         * Sort the cache entries by address since flushing them in
         * ascending order by address is much more efficient.
         */

    /* Allocate space for the mapping array */
    if (NULL==(map=H5FL_ARR_ALLOC(int,(hsize_t)cache->nslots,0)))
        HRETURN_ERROR (H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed");

#ifdef H5_HAVE_PARALLEL
        /* If MPIO or MPIPOSIX is used, do special parallel I/O actions */
        if(IS_H5FD_MPIO(f) || IS_H5FD_MPIPOSIX(f)) {
            H5AC_info_t       **dinfo;
            H5AC_subid_t        type_id;
#ifndef NDEBUG
            const H5D_xfer_t   *xfer_parms;

            /* Get the dataset transfer property list */
            if (H5P_DEFAULT == dxpl_id) {
                xfer_parms = &H5D_xfer_dflt;
            } else if (H5P_DATASET_XFER != H5P_get_class(dxpl_id) ||
                       NULL == (xfer_parms = H5I_object(dxpl_id))) {
                HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, NULL, "not xfer parms");
            }

            /* Sanity check transfer mode */
            assert(xfer_parms->xfer_mode==H5FD_MPIO_COLLECTIVE);
#endif /* NDEBUG */

            /* Create the mapping */
            for (i = nslots = 0; i < cache->nslots; i++) {
                info = cache->slot + i;
                dinfo = cache->dslot + i;

                /* Move dirty metadata from 'held' slots into 'regular' slots */
                if((*dinfo)!=NULL) {
                    H5AC_dest_func_t        dest;

                    /* Various sanity checks */
                    assert((*dinfo)->dirty);
                    assert((*info)!=NULL);
                    assert((*info)->dirty==0);

                    type_id=(*info)->type->id;  /* Remember this for later */

                    /* Destroy 'current' information */
                    dest = (*info)->type->dest;
                    if ((dest)(f, (*info))<0)
                        HRETURN_ERROR(H5E_CACHE, H5E_CANTFREE, NULL, "unable to free cached object");

                    /* Restore 'held' information back to 'current' information */
                    (*info)=(*dinfo);

                    /* Clear 'held' information */
                    (*dinfo)=NULL;

                    cache->diagnostics[type_id].nrestores++;
                } /* end if */
                if ((*info))
                    map[nslots++] = i;
            } /* end for */
        } /* end if */
        else {
#endif /* H5_HAVE_PARALLEL */

            /* Create the mapping */
            for (i = nslots = 0; i < cache->nslots; i++) {
                if (cache->slot[i]!=NULL)
                    map[nslots++] = i;
            }
#ifdef H5_HAVE_PARALLEL
        } /* end else */
#endif /* H5_HAVE_PARALLEL */
        assert(NULL == current_cache_g);
        current_cache_g = cache;
        HDqsort(map, nslots, sizeof(int), H5AC_compare);
        current_cache_g = NULL;
#ifndef NDEBUG
        for (i = 1; i < nslots; i++)
            assert(H5F_addr_lt(cache->slot[map[i - 1]]->addr, cache->slot[map[i]]->addr));
#endif

        /*
         * Look at all cache entries.
         */
        for (i = 0; i < nslots; i++) {
            info = cache->slot + map[i];
            assert(*info);
            if (!type || type == (*info)->type) {
                H5AC_subid_t type_id=(*info)->type->id;  /* Remember this for later */

                flush = (*info)->type->flush;

                /* Only block for all the processes on the first piece of metadata */
                if(first_flush && (*info)->dirty) {
                    status = (flush)(f, dxpl_id, destroy, (*info)->addr, (*info));
                    first_flush=0;
                } /* end if */
                else
                    status = (flush)(f, H5AC_noblock_dxpl_id, destroy, (*info)->addr, (*info));
                if (status < 0) {
                    map = H5FL_ARR_FREE(int,map);
                    HRETURN_ERROR(H5E_CACHE, H5E_CANTFLUSH, FAIL, "unable to flush cache");
                }
                cache->diagnostics[type_id].nflushes++;
                if (destroy)
                    (*info)= NULL;
            }
        }
        map = H5FL_ARR_FREE(int,map);

        /*
         * If there are protected object then fail.  However, everything
         * else should have been flushed.
         */
        if (cache->nprots > 0)
            HRETURN_ERROR(H5E_CACHE, H5E_PROTECT, FAIL, "cache has protected items");
    } else {
        i = H5AC_HASH(f, addr);
        info = cache->slot + i;
#ifdef H5_HAVE_PARALLEL
        /* If MPIO or MPIPOSIX is used, do special parallel I/O actions */
        if(IS_H5FD_MPIO(f) || IS_H5FD_MPIPOSIX(f)) {
            H5AC_info_t       **dinfo;
            H5AC_subid_t        type_id;
#ifndef NDEBUG
            const H5D_xfer_t   *xfer_parms;

            /* Get the dataset transfer property list */
            if (H5P_DEFAULT == dxpl_id) {
                xfer_parms = &H5D_xfer_dflt;
            } else if (H5P_DATASET_XFER != H5P_get_class(dxpl_id) ||
                       NULL == (xfer_parms = H5I_object(dxpl_id))) {
                HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, NULL, "not xfer parms");
            }

            /* Sanity check transfer mode */
            assert(xfer_parms->xfer_mode==H5FD_MPIO_COLLECTIVE);
#endif /* NDEBUG */

            dinfo = cache->dslot + i;

            /* Restore dirty metadata from 'held' slot to 'current' slot */
            if((*dinfo)!=NULL) {
                H5AC_dest_func_t        dest;

                /* Various sanity checks */
                assert((*dinfo)->dirty);
                assert((*info)!=NULL);
                assert((*info)->dirty==0);

                type_id=(*info)->type->id;  /* Remember this for later */

                /* Destroy 'current' information */
                dest = (*info)->type->dest;
                if ((dest)(f, (*info))<0)
                    HRETURN_ERROR(H5E_CACHE, H5E_CANTFREE, NULL, "unable to free cached object");

                /* Restore 'held' information back to 'current' information */
                (*info)=(*dinfo);

                /* Clear 'held' information */
                (*dinfo)=NULL;

                cache->diagnostics[type_id].nrestores++;
            } /* end if */

        } /* end if */
#endif /* H5_HAVE_PARALLEL */
        if ((*info) && (!type || (*info)->type == type) &&
                H5F_addr_eq((*info)->addr, addr)) {

            H5AC_subid_t type_id=(*info)->type->id;  /* Remember this for later */

            /*
             * Flush just this entry.
             */
            flush = (*info)->type->flush;
            if ( (flush)(f, dxpl_id, destroy, (*info)->addr, (*info)) < 0)
                HRETURN_ERROR(H5E_CACHE, H5E_CANTFLUSH, FAIL, "unable to flush object");
            cache->diagnostics[type_id].nflushes++;
            if (destroy)
                (*info)= NULL;
        }
    }

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:    H5AC_set
 *
 * Purpose:     Adds the specified thing to the cache.  The thing need not
 *              exist on disk yet, but it must have an address and disk
 *              space reserved.
 *
 *              If H5AC_DEBUG is defined then this function checks
 *              that the object being inserted isn't a protected object.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *              matzke@llnl.gov
 *              Jul  9 1997
 *
 * Modifications:
 * 		Robb Matzke, 1999-07-27
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
herr_t
H5AC_set(H5F_t *f, hid_t dxpl_id, const H5AC_class_t *type, haddr_t addr, void *thing)
{
    unsigned                   idx;
    H5AC_flush_func_t       flush=NULL;
    H5AC_info_t           **info = NULL;
    H5AC_t                 *cache = NULL;

    FUNC_ENTER(H5AC_set, FAIL);

    assert(f);
    assert(f->shared->cache);
    assert(type);
    assert(type->flush);
    assert(H5F_addr_defined(addr));
    assert(thing);

    /* Get local information */
    idx = H5AC_HASH(f, addr);
    cache = f->shared->cache;
    info = cache->slot + idx;

#ifdef H5AC_DEBUG
    {
        H5AC_prot_t            *prot = NULL;
        int                    i;

        prot = cache->prot + idx;
        for (i = 0; i < prot->nprots; i++) {
            assert(H5F_addr_ne(addr, prot->slot[i]->addr));
        }
    }
#endif

#ifdef H5_HAVE_PARALLEL
    /* If MPIO or MPIPOSIX is used, do special parallel I/O actions */
    if(IS_H5FD_MPIO(f) || IS_H5FD_MPIPOSIX(f)) {
        H5AC_info_t       **dinfo;
        H5AC_subid_t        type_id;
        const H5D_xfer_t   *xfer_parms;

        /* Get the dataset transfer property list */
        if (H5P_DEFAULT == dxpl_id) {
            xfer_parms = &H5D_xfer_dflt;
        } else if (H5P_DATASET_XFER != H5P_get_class(dxpl_id) ||
                   NULL == (xfer_parms = H5I_object(dxpl_id))) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, NULL, "not xfer parms");
        }

        /* Get pointer to 'held' information */
        dinfo = cache->dslot + idx;

        /* Sanity check transfer mode */
        if(xfer_parms->xfer_mode==H5FD_MPIO_COLLECTIVE) {
            /* Check for dirty metadata */
            if(*dinfo) {
                H5AC_dest_func_t        dest;

                /* Various sanity checks */
                assert((*dinfo)->dirty);
                assert((*info)!=NULL);
                assert((*info)->dirty==0);

                type_id=(*info)->type->id;  /* Remember this for later */

                /* Destroy 'current' information */
                dest = (*info)->type->dest;
                if ((dest)(f, (*info))<0)
                    HRETURN_ERROR(H5E_CACHE, H5E_CANTFREE, NULL, "unable to free cached object");

                /* Restore 'held' information back to 'current' information */
                (*info)=(*dinfo);

                /* Clear 'held' information */
                (*dinfo)=NULL;

                cache->diagnostics[type_id].nrestores++;
            } /* end if */
        } /* end if */
        else {
            /* Sanity check */
            assert((*dinfo)==NULL);
            assert(xfer_parms->xfer_mode==H5FD_MPIO_INDEPENDENT);

            /* Make certain there will be no write of dirty metadata */
            if((*info) && (*info)->dirty) {
                /* Sanity check new item */
                assert(((H5AC_info_t*)thing)->dirty==0);

                /* 'Hold' the current metadata for later */
                (*dinfo)=(*info);

                /* Reset the 'current' metadata, so it doesn't get flushed */
                (*info)=NULL;

                cache->diagnostics[(*dinfo)->type->id].nholds++;
            } /* end if */
        } /* end else */
    } /* end if */
#endif /* H5_HAVE_PARALLEL */

    /* Flush any existing item out of the cache entry */
    if ((*info)) {
        H5AC_subid_t type_id=(*info)->type->id;  /* Remember this for later */

        flush = (*info)->type->flush;
        if ( (flush)(f, dxpl_id, TRUE, (*info)->addr, (*info)) < 0)
            HRETURN_ERROR(H5E_CACHE, H5E_CANTFLUSH, FAIL, "unable to flush object");
        cache->diagnostics[type_id].nflushes++;
    }

    /* Cache this item */
    (*info)=thing;
    (*info)->type = type;
    (*info)->addr = addr;
    cache->diagnostics[type->id].ninits++;

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:    H5AC_rename
 *
 * Purpose:     Use this function to notify the cache that an object's
 *              file address changed.
 *
 *              If H5AC_DEBUG is defined then this function checks
 *              that the old and new addresses don't correspond to the
 *              address of a protected object.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *              matzke@llnl.gov
 *              Jul  9 1997
 *
 * Modifications:
 * 		Robb Matzke, 1999-07-27
 *		The OLD_ADDR and NEW_ADDR arguments are passed by value.
 *-------------------------------------------------------------------------
 */
herr_t
H5AC_rename(H5F_t *f, hid_t dxpl_id, const H5AC_class_t *type, haddr_t old_addr,
	    haddr_t new_addr)
{
    unsigned                   old_idx, new_idx;
    H5AC_flush_func_t       flush=NULL;
    H5AC_t                 *cache = NULL;
    H5AC_info_t           **new_info = NULL;
    H5AC_info_t           **old_info = NULL;

    FUNC_ENTER(H5AC_rename, FAIL);

    assert(f);
    assert(f->shared->cache);
    assert(type);

    /* Get local information */
    old_idx = H5AC_HASH(f, old_addr);
    new_idx = H5AC_HASH(f, new_addr);
    cache = f->shared->cache;
    new_info = cache->slot + new_idx;
    old_info = cache->slot + old_idx;

#ifdef H5AC_DEBUG
    {
        H5AC_prot_t            *prot = NULL;
        int                     i;

        prot = cache->prot + old_idx;
        for (i = 0; i < prot->nprots; i++) {
            assert(H5F_addr_ne(old_addr, prot->slot[i]->addr));
        }
        prot = cache->prot + new_idx;
        for (i = 0; i < prot->nprots; i++) {
            assert(H5F_addr_ne(new_addr, prot->slot[i]->addr));
        }
    }
#endif

    /*
     * We don't need to do anything if the object isn't cached or if the
     * new hash value is the same as the old one.
     */
    if (H5F_addr_ne((*old_info)->addr, old_addr) || (*old_info)->type!=type) {
        HRETURN(SUCCEED);
    }
    if (old_idx == new_idx) {
        (*old_info)->addr = new_addr;
        HRETURN(SUCCEED);
    }

#ifdef H5_HAVE_PARALLEL
    /* If MPIO or MPIPOSIX is used, do special parallel I/O actions */
    if(IS_H5FD_MPIO(f) || IS_H5FD_MPIPOSIX(f)) {
        H5AC_info_t       **new_dinfo;
        H5AC_subid_t        type_id;
        const H5D_xfer_t   *xfer_parms;

        /* Get the dataset transfer property list */
        if (H5P_DEFAULT == dxpl_id) {
            xfer_parms = &H5D_xfer_dflt;
        } else if (H5P_DATASET_XFER != H5P_get_class(dxpl_id) ||
                   NULL == (xfer_parms = H5I_object(dxpl_id))) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, NULL, "not xfer parms");
        }

        /* Get pointer to new 'held' information */
        new_dinfo = cache->dslot + new_idx;

        /* Sanity check transfer mode */
        if(xfer_parms->xfer_mode==H5FD_MPIO_COLLECTIVE) {
            /* Check for dirty metadata */
            if(*new_dinfo) {
                H5AC_dest_func_t        dest;

                /* Various sanity checks */
                assert((*new_dinfo)->dirty);
                assert((*new_info)!=NULL);
                assert((*new_info)->dirty==0);

                type_id=(*new_info)->type->id;  /* Remember this for later */

                /* Destroy 'current' information */
                dest = (*new_info)->type->dest;
                if ((dest)(f, (*new_info))<0)
                    HRETURN_ERROR(H5E_CACHE, H5E_CANTFREE, NULL, "unable to free cached object");

                /* Restore 'held' information back to 'current' information */
                (*new_info)=(*new_dinfo);

                /* Clear 'held' information */
                (*new_dinfo)=NULL;

                cache->diagnostics[type_id].nrestores++;
            } /* end if */
        } /* end if */
        else {
            /* Sanity check that there will be no write of dirty metadata */
            assert((*new_dinfo)==NULL);
            assert(xfer_parms->xfer_mode==H5FD_MPIO_INDEPENDENT);

            /* Make certain there will be no write of dirty metadata */
            if((*new_info) && (*new_info)->dirty) {
                /* Sanity check that we won't put two pieces of dirty metadata in same cache location */
                assert((*old_info)->dirty==0);

                /* 'Hold' the current metadata for later */
                (*new_dinfo)=(*new_info);

                /* Reset the 'current' metadata, so it doesn't get flushed */
                (*new_info)=NULL;

                cache->diagnostics[(*new_dinfo)->type->id].nholds++;
            } /* end if */
        } /* end else */
    } /* end if */
#endif /* H5_HAVE_PARALLEL */

    /*
     * Free the item from the destination cache line.
     */
    if ((*new_info)) {
        H5AC_subid_t type_id=(*new_info)->type->id;  /* Remember this for later */

        flush = (*new_info)->type->flush;
        if ( (flush)(f, dxpl_id, TRUE, (*new_info)->addr, (*new_info)) < 0)
            HRETURN_ERROR(H5E_CACHE, H5E_CANTFLUSH, FAIL, "unable to flush object");
        cache->diagnostics[type_id].nflushes++;
    }

    /*
     * Move the source to the destination (it might not be cached)
     */
    (*new_info)= (*old_info);
    (*new_info)->addr = new_addr;

#ifdef H5_HAVE_PARALLEL
    /* If MPIO or MPIPOSIX is used, do special parallel I/O actions */
    if(IS_H5FD_MPIO(f) || IS_H5FD_MPIPOSIX(f)) {
        H5AC_info_t       **old_dinfo;
        H5AC_subid_t        type_id;

        /* Get pointer to new 'held' information */
        old_dinfo = cache->dslot + old_idx;

        /* Check for 'held' metadata in old location & restore it, if so */
        if(*old_dinfo) {
            /* Sanity check */
            assert((*old_dinfo)->dirty);

            type_id=(*old_info)->type->id;  /* Remember this for later */

            /* Restore 'held' information back to 'current' information */
            (*old_info)=(*old_dinfo);

            /* Clear 'held' information */
            (*old_dinfo)=NULL;

            cache->diagnostics[type_id].nrestores++;
        } /* end if */
        else
            (*old_info)= NULL;
    } /* end if */
    else {
#endif /* H5_HAVE_PARALLEL */

        (*old_info)= NULL;
#ifdef H5_HAVE_PARALLEL
    } /* end else */
#endif /* H5_HAVE_PARALLEL */

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:    H5AC_protect
 *
 * Purpose:     Similar to H5AC_find() except the object is removed from
 *              the cache and given to the caller, preventing other parts
 *              of the program from modifying the protected object or
 *              preempting it from the cache.
 *
 *              The caller must call H5AC_unprotect() when finished with
 *              the pointer.
 *
 *              If H5AC_DEBUG is defined then we check that the
 *              requested object isn't already protected.
 *
 * Return:      Success:        Ptr to the object.
 *
 *              Failure:        NULL
 *
 * Programmer:  Robb Matzke
 *              matzke@llnl.gov
 *              Sep  2 1997
 *
 * Modifications:
 *		Robb Matzke, 1999-07-27
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
void *
H5AC_protect(H5F_t *f, hid_t dxpl_id, const H5AC_class_t *type, haddr_t addr,
	     const void *udata1, void *udata2)
{
    int                     idx;
    void                   *thing = NULL;
    H5AC_t                 *cache = NULL;
    H5AC_info_t           **info = NULL;

#ifdef H5AC_DEBUG
    H5AC_prot_t            *prot = NULL;

    static int ncalls = 0;
    if (0 == ncalls++) {
	if (H5DEBUG(AC)) {
	    fprintf(H5DEBUG(AC), "H5AC: debugging cache (expensive)\n");
	}
    }
#endif

    FUNC_ENTER(H5AC_protect, NULL);

    /* check args */
    assert(f);
    assert(f->shared->cache);
    assert(type);
    assert(type->load);
    assert(type->flush);
    assert(H5F_addr_defined(addr));

    /* Get local information */
    idx = H5AC_HASH(f, addr);
    cache = f->shared->cache;
    info = cache->slot + idx;
#ifdef H5AC_DEBUG
    prot = cache->prot + idx;
#endif /* H5AC_DEBUG */

#ifdef H5_HAVE_PARALLEL
    /* If MPIO or MPIPOSIX is used, do special parallel I/O actions */
    if(IS_H5FD_MPIO(f) || IS_H5FD_MPIPOSIX(f)) {
        H5AC_info_t       **dinfo;

        /* Get pointer to new 'held' information */
        dinfo = cache->dslot + idx;

        /* Check for 'held' metadata in location & handle it */
        if(*dinfo) {
            /* Sanity checks */
            assert((*dinfo)->dirty);
            assert((*info));
            assert((*info)->dirty==0);
            assert((*dinfo)->addr!=(*info)->addr);

            /* Is 'held' metadata the metadata we are looking for? */
            if (H5F_addr_eq((*dinfo)->addr, addr) && (*dinfo)->type==type) {
                /*
                 * The object is already cached; simply remove it from the cache.
                 */
                cache->diagnostics[(*dinfo)->type->id].nhits++;
                thing = (*dinfo);
                (*dinfo)->type = NULL;
                (*dinfo)->addr = HADDR_UNDEF;
                (*dinfo)= NULL;
            } /* end if */
            /* 'held' metadata isn't what we are looking for, but check for 'current' metadata */
            else {
                if(H5F_addr_eq((*info)->addr, addr) && (*info)->type==type) {
                    /*
                     * The object is already cached; remove it from the cache.
                     * and bring the 'held' object into the 'regular' information
                     */
                    cache->diagnostics[(*info)->type->id].nhits++;
                    thing = (*info);
                    (*info)->type = NULL;
                    (*info)->addr = HADDR_UNDEF;
                    (*info)= (*dinfo);
                    (*dinfo)= NULL;
                } /* end if */
            } /* end else */
        } /* end if */
    } /* end if */

    /* Check if we've already found the object to protect */
    if(thing==NULL) {
#endif /* H5_HAVE_PARALLEL */

        if ((*info) && H5F_addr_eq((*info)->addr, addr) && (*info)->type==type) {
            /*
             * The object is already cached; simply remove it from the cache.
             */
            cache->diagnostics[(*info)->type->id].nhits++;
            thing = (*info);
            (*info)->type = NULL;
            (*info)->addr = HADDR_UNDEF;
            (*info)= NULL;
        } /* end if */
        else {
#ifdef H5AC_DEBUG
            /*
             * Check that the requested thing isn't protected, for protected things
             * can only be modified through the pointer already handed out by the
             * H5AC_protect() function.
             */
            int                    i;

            for (i = 0; i < prot->nprots; i++)
                assert(H5F_addr_ne(addr, prot->slot[i]->addr));
#endif /* H5AC_DEBUG */

            /*
             * Load a new thing.  If it can't be loaded, then return an error
             * without preempting anything.
             */
            cache->diagnostics[type->id].nmisses++;
            if (NULL == (thing = (type->load)(f, dxpl_id, addr, udata1, udata2)))
                HRETURN_ERROR(H5E_CACHE, H5E_CANTLOAD, NULL, "unable to load object");
        } /* end else */
#ifdef H5_HAVE_PARALLEL
    } /* end if */
#endif /* H5_HAVE_PARALLEL */

#ifdef H5AC_DEBUG
    /*
     * Add the protected object to the protect debugging fields of the
     * cache.
     */
    if (prot->nprots >= prot->aprots) {
        size_t na = prot->aprots + 10;
        H5AC_info_t **x = H5MM_realloc(prot->slot,
                          na * sizeof(H5AC_info_t *));
        if (NULL==x)
            HRETURN_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");
        prot->aprots = (int)na;
        prot->slot = x;
    }
    prot->slot[prot->nprots]= thing;
    prot->slot[prot->nprots]->type = type;
    prot->slot[prot->nprots]->addr = addr;
    prot->nprots += 1;
#endif /* H5AC_DEBUG */

    cache->nprots += 1;
    FUNC_LEAVE(thing);
}

/*-------------------------------------------------------------------------
 * Function:    H5AC_unprotect
 *
 * Purpose:     This function should be called to undo the effect of
 *              H5AC_protect().  The TYPE and ADDR arguments should be the
 *              same as the corresponding call to H5AC_protect() and the
 *              THING argument should be the value returned by H5AC_protect().
 *
 *              If H5AC_DEBUG is defined then this function fails
 *              if the TYPE and ADDR arguments are not what was used when the
 *              object was protected or if the object was never protected.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *              matzke@llnl.gov
 *              Sep  2 1997
 *
 * Modifications:
 * 		Robb Matzke, 1999-07-27
 *		The ADDR argument is passed by value.
 *-------------------------------------------------------------------------
 */
herr_t
H5AC_unprotect(H5F_t *f, hid_t dxpl_id, const H5AC_class_t *type, haddr_t addr, void *thing)
{
    unsigned                   idx;
    H5AC_flush_func_t       flush=NULL;
    H5AC_t                 *cache = NULL;
    H5AC_info_t           **info = NULL;

    FUNC_ENTER(H5AC_unprotect, FAIL);

    /* check args */
    assert(f);
    assert(f->shared->cache);
    assert(type);
    assert(type->flush);
    assert(H5F_addr_defined(addr));
    assert(thing);

    /* Get local information */
    idx = H5AC_HASH(f, addr);
    cache = f->shared->cache;
    info = cache->slot + idx;

#ifdef H5_HAVE_PARALLEL
    /* If MPIO or MPIPOSIX is used, do special parallel I/O actions */
    if(IS_H5FD_MPIO(f) || IS_H5FD_MPIPOSIX(f)) {
        H5AC_info_t       **dinfo;
        H5AC_subid_t        type_id;
        const H5D_xfer_t   *xfer_parms;

        /* Get the dataset transfer property list */
        if (H5P_DEFAULT == dxpl_id) {
            xfer_parms = &H5D_xfer_dflt;
        } else if (H5P_DATASET_XFER != H5P_get_class(dxpl_id) ||
                   NULL == (xfer_parms = H5I_object(dxpl_id))) {
            HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, NULL, "not xfer parms");
        }

        /* Get pointer to 'held' information */
        dinfo = cache->dslot + idx;

        /* Sanity check transfer mode */
        if(xfer_parms->xfer_mode==H5FD_MPIO_COLLECTIVE) {
            /* Check for dirty metadata */
            if(*dinfo) {
                H5AC_dest_func_t        dest;

                /* Various sanity checks */
                assert((*dinfo)->dirty);
                assert((*info)!=NULL);
                assert((*info)->dirty==0);

                type_id=(*info)->type->id;  /* Remember this for later */

                /* Destroy 'current' information */
                dest = (*info)->type->dest;
                if ((dest)(f, (*info))<0)
                    HRETURN_ERROR(H5E_CACHE, H5E_CANTFREE, NULL, "unable to free cached object");

                /* Restore 'held' information back to 'current' information */
                (*info)=(*dinfo);

                /* Clear 'held' information */
                (*dinfo)=NULL;

                cache->diagnostics[type_id].nrestores++;
            } /* end if */
        } /* end if */
        else {
            /* Sanity check */
            assert((*dinfo)==NULL);
            assert(xfer_parms->xfer_mode==H5FD_MPIO_INDEPENDENT);

            /* Make certain there will be no write of dirty metadata */
            if((*info) && (*info)->dirty) {
                /* Sanity check new item */
                assert(((H5AC_info_t*)thing)->dirty==0);

                /* 'Hold' the current metadata for later */
                (*dinfo)=(*info);

                /* Reset the 'current' metadata, so it doesn't get flushed */
                (*info)=NULL;

                cache->diagnostics[(*dinfo)->type->id].nholds++;
            } /* end if */
        } /* end else */
    } /* end if */
#endif /* H5_HAVE_PARALLEL */

    /*
     * Flush any object already in the cache at that location.  It had
     * better not be another copy of the protected object.
     */
    if (*info) {
        H5AC_subid_t type_id=(*info)->type->id;  /* Remember this for later */

        assert(H5F_addr_ne((*info)->addr, addr));
        flush = (*info)->type->flush;
        if ( (flush)(f, dxpl_id, TRUE, (*info)->addr, (*info)) < 0)
            HRETURN_ERROR(H5E_CACHE, H5E_CANTFLUSH, FAIL, "unable to flush object");
        cache->diagnostics[type_id].nflushes++;
    }
#ifdef H5AC_DEBUG
    /*
     * Remove the object's protect data to indicate that it is no longer
     * protected.
     */
    {
        H5AC_prot_t            *prot = NULL;
        int                     found, i;

        prot = cache->prot + idx;
        for (i = 0, found = FALSE; i < prot->nprots && !found; i++) {
            if (H5F_addr_eq(addr, prot->slot[i]->addr)) {
                assert(prot->slot[i]->type == type);
                HDmemmove(prot->slot + i, prot->slot + i + 1,
                          ((prot->nprots - i) - 1) * sizeof(H5AC_info_t *));
                prot->nprots -= 1;
                found = TRUE;
            }
        }
        assert(found);
    }
#endif /* H5AC_DEBUG */

    /*
     * Insert the object back into the cache; it is no longer protected.
     */
    (*info)=thing;
    (*info)->type = type;
    (*info)->addr = addr;
    cache->nprots -= 1;

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:    H5AC_debug
 *
 * Purpose:     Prints debugging info about the cache.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Robb Matzke
 *              Thursday, October 30, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5AC_debug(H5F_t UNUSED *f)
{
#ifdef H5AC_DEBUG
    H5AC_subid_t            i;
    char                    s[32], ascii[32];
    H5AC_t                 *cache = f->shared->cache;
    double                  miss_rate;
#endif

    FUNC_ENTER(H5AC_debug, FAIL);

#ifdef H5AC_DEBUG
    if (H5DEBUG(AC)) {
	fprintf(H5DEBUG(AC), "H5AC: meta data cache statistics for file %s\n",
		f->name);
	fprintf(H5DEBUG(AC), "   %-18s %8s %8s %8s %8s+%-8s\n",
		"Layer", "Hits", "Misses", "MissRate", "Inits", "Flushes");
	fprintf(H5DEBUG(AC), "   %-18s %8s %8s %8s %8s-%-8s\n",
		"-----", "----", "------", "--------", "-----", "-------");

	for (i = H5AC_BT_ID; i < H5AC_NTYPES; i++) {

	    switch (i) {
	    case H5AC_BT_ID:
		HDstrcpy(s, "B-tree nodes");
		break;
	    case H5AC_SNODE_ID:
		HDstrcpy(s, "symbol table nodes");
		break;
	    case H5AC_LHEAP_ID:
		HDstrcpy (s, "local heaps");
		break;
	    case H5AC_GHEAP_ID:
		HDstrcpy (s, "global heaps");
		break;
	    case H5AC_OHDR_ID:
		HDstrcpy(s, "object headers");
		break;
	    default:
		sprintf(s, "unknown id %d", i);
	    }

	    if (cache->diagnostics[i].nhits>0 ||
		cache->diagnostics[i].nmisses>0) {
		miss_rate = 100.0 * cache->diagnostics[i].nmisses /
			    (cache->diagnostics[i].nhits+
			     cache->diagnostics[i].nmisses);
	    } else {
		miss_rate = 0.0;
	    }

	    if (miss_rate > 100) {
		sprintf(ascii, "%7d%%", (int) (miss_rate + 0.5));
	    } else {
		sprintf(ascii, "%7.2f%%", miss_rate);
	    }
	    fprintf(H5DEBUG(AC), "   %-18s %8u %8u %7s %8u%+-9ld\n", s,
		    cache->diagnostics[i].nhits,
		    cache->diagnostics[i].nmisses,
		    ascii,
		    cache->diagnostics[i].ninits,
		    ((long)(cache->diagnostics[i].nflushes) -
		     (long)(cache->diagnostics[i].ninits)));
	}
    }
#endif

    FUNC_LEAVE(SUCCEED);
}
