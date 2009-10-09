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

/* Programmer:  Quincey Koziol <koziol@ncsa.uiuc.edu>
 *              Monday, October 17, 2005
 *
 * Purpose:	Group testing functions.
 */

/****************/
/* Module Setup */
/****************/

#define H5G_PACKAGE		/*suppress error about including H5Gpkg	  */
#define H5G_TESTING		/*suppress warning about H5G testing funcs*/


/***********/
/* Headers */
/***********/
#include "H5private.h"		/* Generic Functions			*/
#include "H5Dprivate.h"		/* Datasets				*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5Gpkg.h"		/* Groups		  		*/
#include "H5HLprivate.h"	/* Local Heaps				*/
#include "H5Iprivate.h"		/* IDs			  		*/

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


/*********************/
/* Package Variables */
/*********************/


/*****************************/
/* Library Private Variables */
/*****************************/


/*******************/
/* Local Variables */
/*******************/


/*--------------------------------------------------------------------------
 NAME
    H5G_is_empty_test
 PURPOSE
    Determine whether a group contains no objects
 USAGE
    htri_t H5G_is_empty_test(gid)
        hid_t gid;              IN: group to check
 RETURNS
    Non-negative TRUE/FALSE on success, negative on failure
 DESCRIPTION
    Checks to see if the group has no link messages and no symbol table message
    and no "dense" link storage
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
    DO NOT USE THIS FUNCTION FOR ANYTHING EXCEPT TESTING
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
htri_t
H5G_is_empty_test(hid_t gid)
{
    H5G_t *grp = NULL;          /* Pointer to group */
    htri_t msg_exists = FALSE;  /* Indicate that a header message is present */
    htri_t linfo_exists = FALSE;/* Indicate that the 'link info' message is present */
    htri_t ret_value = TRUE;    /* Return value */

    FUNC_ENTER_NOAPI(H5G_is_empty_test, FAIL)

    /* Get group structure */
    if(NULL == (grp = (H5G_t *)H5I_object_verify(gid, H5I_GROUP)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a group")

    /* "New format" checks */

    /* Check if the group has any link messages */
    if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_LINK_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(msg_exists > 0) {
        /* Sanity check that new group format shouldn't have old messages */
        if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_STAB_ID, H5AC_dxpl_id)) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
        if(msg_exists > 0)
            HGOTO_ERROR(H5E_SYM, H5E_BADVALUE, FAIL, "both symbol table and link messages found")

        HGOTO_DONE(FALSE)
    } /* end if */

    /* Check for a link info message */
    if((linfo_exists = H5O_msg_exists(&(grp->oloc), H5O_LINFO_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(linfo_exists > 0) {
        H5O_linfo_t linfo;		/* Link info message */

        /* Sanity check that new group format shouldn't have old messages */
        if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_STAB_ID, H5AC_dxpl_id)) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
        if(msg_exists > 0)
            HGOTO_ERROR(H5E_SYM, H5E_BADVALUE, FAIL, "both symbol table and link info messages found")

        /* Get the link info */
        if(H5G_obj_get_linfo(&(grp->oloc), &linfo, H5AC_dxpl_id) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_BADMESG, FAIL, "can't get link info")

        /* Check for 'dense' link storage file addresses being defined */
        if(H5F_addr_defined(linfo.fheap_addr))
            HGOTO_DONE(FALSE)
        if(H5F_addr_defined(linfo.name_bt2_addr))
            HGOTO_DONE(FALSE)
        if(H5F_addr_defined(linfo.corder_bt2_addr))
            HGOTO_DONE(FALSE)

        /* Check for link count */
        if(linfo.nlinks > 0)
            HGOTO_DONE(FALSE)
    } /* end if */

    /* "Old format" checks */

    /* Check if the group has a symbol table message */
    if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_STAB_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(msg_exists > 0) {
        H5O_stab_t stab;        /* Info about local heap & B-tree */
        hsize_t nlinks;         /* Number of links in the group */

        /* Sanity check that old group format shouldn't have new messages */
        if(linfo_exists > 0)
            HGOTO_ERROR(H5E_SYM, H5E_BADVALUE, FAIL, "both symbol table and link info messages found")
        if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_GINFO_ID, H5AC_dxpl_id)) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
        if(msg_exists > 0)
            HGOTO_ERROR(H5E_SYM, H5E_BADVALUE, FAIL, "both symbol table and group info messages found")

        /* Get the B-tree & local heap info */
        if(NULL == H5O_msg_read(&(grp->oloc), H5O_STAB_ID, &stab, H5AC_dxpl_id))
            HGOTO_ERROR(H5E_SYM, H5E_NOTFOUND, FAIL, "unable to read symbol table message")

        /* Get the count of links in the group */
        if(H5G_stab_count(&(grp->oloc), &nlinks, H5AC_dxpl_id) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_NOTFOUND, FAIL, "unable to count links")

        /* Check for link count */
        if(nlinks > 0)
            HGOTO_DONE(FALSE)
    } /* end if */

done:
    FUNC_LEAVE_NOAPI(ret_value)
}   /* H5G_is_empty_test() */


/*--------------------------------------------------------------------------
 NAME
    H5G_has_links_test
 PURPOSE
    Determine whether a group contains link messages
 USAGE
    htri_t H5G_has_links_test(gid)
        hid_t gid;              IN: group to check
        unsigned *nmsgs;        OUT: # of link messages in header
 RETURNS
    Non-negative TRUE/FALSE on success, negative on failure
 DESCRIPTION
    Checks to see if the group has link messages and how many.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
    DO NOT USE THIS FUNCTION FOR ANYTHING EXCEPT TESTING
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
htri_t
H5G_has_links_test(hid_t gid, unsigned *nmsgs)
{
    H5G_t *grp = NULL;          /* Pointer to group */
    htri_t msg_exists = 0;      /* Indicate that a header message is present */
    htri_t ret_value = TRUE;    /* Return value */

    FUNC_ENTER_NOAPI(H5G_has_links_test, FAIL)

    /* Get group structure */
    if(NULL == (grp = (H5G_t *)H5I_object_verify(gid, H5I_GROUP)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a group")

    /* Check if the group has any link messages */
    if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_LINK_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(msg_exists == 0)
        HGOTO_DONE(FALSE)

    /* Check if the group has a symbol table message */
    if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_STAB_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(msg_exists > 0)
	HGOTO_ERROR(H5E_SYM, H5E_BADVALUE, FAIL, "both symbol table and link messages found")

    /* Check if we should retrieve the number of link messages */
    if(nmsgs) {
        int msg_count;     /* Number of messages of a type */

        /* Check how many link messages there are */
        if((msg_count = H5O_msg_count(&(grp->oloc), H5O_LINK_ID, H5AC_dxpl_id)) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTCOUNT, FAIL, "unable to count link messages")
        *nmsgs = (unsigned)msg_count;
    } /* end if */

done:
    FUNC_LEAVE_NOAPI(ret_value)
}   /* H5G_has_links_test() */


/*--------------------------------------------------------------------------
 NAME
    H5G_has_stab_test
 PURPOSE
    Determine whether a group contains a symbol table message
 USAGE
    htri_t H5G_has_stab_test(gid)
        hid_t gid;              IN: group to check
 RETURNS
    Non-negative TRUE/FALSE on success, negative on failure
 DESCRIPTION
    Checks to see if the group has a symbol table message.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
    DO NOT USE THIS FUNCTION FOR ANYTHING EXCEPT TESTING
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
htri_t
H5G_has_stab_test(hid_t gid)
{
    H5G_t *grp = NULL;          /* Pointer to group */
    htri_t msg_exists = 0;      /* Indicate that a header message is present */
    htri_t ret_value = TRUE;    /* Return value */

    FUNC_ENTER_NOAPI(H5G_has_stab_test, FAIL)

    /* Get group structure */
    if(NULL == (grp = (H5G_t *)H5I_object_verify(gid, H5I_GROUP)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a group")

    /* Check if the group has a symbol table message */
    if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_STAB_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(msg_exists == 0)
        HGOTO_DONE(FALSE)

    /* Check if the group has any link messages */
    if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_LINK_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(msg_exists > 0)
	HGOTO_ERROR(H5E_SYM, H5E_BADVALUE, FAIL, "both symbol table and link messages found")

done:
    FUNC_LEAVE_NOAPI(ret_value)
}   /* H5G_has_stab_test() */


/*--------------------------------------------------------------------------
 NAME
    H5G_is_new_dense_test
 PURPOSE
    Determine whether a group is in the "new" format and dense
 USAGE
    htri_t H5G_is_new_dense_test(gid)
        hid_t gid;              IN: group to check
 RETURNS
    Non-negative TRUE/FALSE on success, negative on failure
 DESCRIPTION
    Checks to see if the group is in the "new" format for groups (link messages/
    fractal heap+v2 B-tree) and if it is in "dense" storage form (ie. it has
    a name B-tree index).
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
    DO NOT USE THIS FUNCTION FOR ANYTHING EXCEPT TESTING
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
htri_t
H5G_is_new_dense_test(hid_t gid)
{
    H5G_t *grp = NULL;          /* Pointer to group */
    htri_t msg_exists = 0;      /* Indicate that a header message is present */
    htri_t ret_value = TRUE;    /* Return value */

    FUNC_ENTER_NOAPI(H5G_is_new_dense_test, FAIL)

    /* Get group structure */
    if(NULL == (grp = (H5G_t *)H5I_object_verify(gid, H5I_GROUP)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a group")

    /* Check if the group has a symbol table message */
    if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_STAB_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(msg_exists > 0)
        HGOTO_DONE(FALSE)

    /* Check if the group has any link messages */
    if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_LINK_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(msg_exists > 0)
        HGOTO_DONE(FALSE)

    /* Check if the group has link info message */
    if((msg_exists = H5O_msg_exists(&(grp->oloc), H5O_LINFO_ID, H5AC_dxpl_id)) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read object header")
    if(msg_exists > 0) {
        H5O_linfo_t linfo;		/* Link info message */

        /* Get the link info */
        if(H5G_obj_get_linfo(&(grp->oloc), &linfo, H5AC_dxpl_id) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_BADMESG, FAIL, "can't get link info")

        /* Check for 'dense' link storage file addresses being defined */
        if(!H5F_addr_defined(linfo.fheap_addr))
            HGOTO_DONE(FALSE)
        if(!H5F_addr_defined(linfo.name_bt2_addr))
            HGOTO_DONE(FALSE)
    } /* end if */

done:
    FUNC_LEAVE_NOAPI(ret_value)
}   /* H5G_is_new_dense_test() */


/*--------------------------------------------------------------------------
 NAME
    H5G_new_dense_info_test
 PURPOSE
    Retrieve information about the state of the new "dense" storage for groups
 USAGE
    herr_t H5G_new_dense_info_test(gid, name_count, corder_count)
        hid_t gid;              IN: group to check
        hsize_t *name_count;    OUT: Number of links in name index
        hsize_t *corder_count;  OUT: Number of links in creation order index
 RETURNS
    Non-negative on success, negative on failure
 DESCRIPTION
    Currently, just retrieves the number of links in each index and returns
    them.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
    DO NOT USE THIS FUNCTION FOR ANYTHING EXCEPT TESTING
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5G_new_dense_info_test(hid_t gid, hsize_t *name_count, hsize_t *corder_count)
{
    H5O_linfo_t linfo;		/* Link info message */
    H5G_t *grp = NULL;          /* Pointer to group */
    herr_t ret_value = SUCCEED; /* Return value */

    FUNC_ENTER_NOAPI(H5G_new_dense_info_test, FAIL)

    /* Get group structure */
    if(NULL == (grp = (H5G_t *)H5I_object_verify(gid, H5I_GROUP)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a group")

    /* Get the link info */
    if(H5G_obj_get_linfo(&(grp->oloc), &linfo, H5AC_dxpl_id) < 0)
        HGOTO_ERROR(H5E_SYM, H5E_BADMESG, FAIL, "can't get link info")

    /* Check for 'dense' link storage file addresses being defined */
    if(!H5F_addr_defined(linfo.fheap_addr))
        HGOTO_DONE(FAIL)
    if(!H5F_addr_defined(linfo.name_bt2_addr))
        HGOTO_DONE(FAIL)

    /* Retrieve # of records in name index */
    if(H5B2_get_nrec(grp->oloc.file, H5AC_dxpl_id, H5G_BT2_NAME, linfo.name_bt2_addr, name_count) < 0)
        HGOTO_ERROR(H5E_SYM, H5E_CANTCOUNT, FAIL, "unable to retrieve # of records from name index")

    /* Check if there is a creation order index */
    if(H5F_addr_defined(linfo.corder_bt2_addr)) {
        /* Retrieve # of records in creation order index */
        if(H5B2_get_nrec(grp->oloc.file, H5AC_dxpl_id, H5G_BT2_CORDER, linfo.corder_bt2_addr, corder_count) < 0)
            HGOTO_ERROR(H5E_SYM, H5E_CANTCOUNT, FAIL, "unable to retrieve # of records from creation order index")
    } /* end if */
    else
        *corder_count = 0;

done:
    FUNC_LEAVE_NOAPI(ret_value)
}   /* H5G_new_dense_info_test() */


/*--------------------------------------------------------------------------
 NAME
    H5G_lheap_size_test
 PURPOSE
    Determine the size of a local heap for a group
 USAGE
    herr_t H5G_lheap_size_test(gid, lheap_size)
        hid_t gid;              IN: group to check
        size_t *lheap_size;     OUT: Size of local heap
 RETURNS
    Non-negative on success, negative on failure
 DESCRIPTION
    Checks the size of the local heap for a group
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
    DO NOT USE THIS FUNCTION FOR ANYTHING EXCEPT TESTING
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5G_lheap_size_test(hid_t gid, size_t *lheap_size)
{
    H5G_t *grp = NULL;          /* Pointer to group */
    H5O_stab_t stab;		/* Symbol table message		*/
    herr_t ret_value = SUCCEED; /* Return value */

    FUNC_ENTER_NOAPI(H5G_lheap_size_test, FAIL)

    /* Get group structure */
    if(NULL == (grp = (H5G_t *)H5I_object_verify(gid, H5I_GROUP)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a group")

    /* Make certain the group has a symbol table message */
    if(NULL == H5O_msg_read(&(grp->oloc), H5O_STAB_ID, &stab, H5AC_dxpl_id))
	HGOTO_ERROR(H5E_SYM, H5E_CANTINIT, FAIL, "unable to read symbol table message")

    /* Check the size of the local heap for the group */
    if(H5HL_get_size(grp->oloc.file, H5AC_dxpl_id, stab.heap_addr, lheap_size) < 0)
	HGOTO_ERROR(H5E_SYM, H5E_CANTGETSIZE, FAIL, "can't query local heap size")

done:
    FUNC_LEAVE_NOAPI(ret_value)
}   /* H5G_lheap_size_test() */


/*--------------------------------------------------------------------------
 NAME
    H5G_user_path_test
 PURPOSE
    Retrieve the user path for an ID
 USAGE
    herr_t H5G_user_path_test(obj_id, user_path, user_path_len)
        hid_t obj_id;           IN: ID to check
        char *user_path;        OUT: Pointer to buffer for User path
        size_t *user_path_len;  OUT: Size of user path
        unsigned *obj_hidden;   OUT: Whether object is hidden
 RETURNS
    Non-negative on success, negative on failure
 DESCRIPTION
    Retrieves the user path for an ID.  A zero for the length is returned in
    the case of no user path.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
    DO NOT USE THIS FUNCTION FOR ANYTHING EXCEPT TESTING
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5G_user_path_test(hid_t obj_id, char *user_path, size_t *user_path_len, unsigned *obj_hidden)
{
    void *obj_ptr;              /* Pointer to object for ID */
    H5G_name_t *obj_path;       /* Pointer to group hier. path for obj */
    herr_t ret_value = SUCCEED; /* Return value */

    FUNC_ENTER_NOAPI(H5G_user_path_test, FAIL)

    /* Sanity check */
    HDassert(user_path_len);
    HDassert(obj_hidden);

    /* Get pointer to object for ID */
    if(NULL == (obj_ptr = H5I_object(obj_id)))
         HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "can't get object for ID")

    /* Get the symbol table entry */
    switch(H5I_get_type(obj_id)) {
        case H5I_GROUP:
            obj_path = H5G_nameof((H5G_t *)obj_ptr);
            break;

        case H5I_DATASET:
            obj_path = H5D_nameof((H5D_t *)obj_ptr);
            break;

        case H5I_DATATYPE:
            /* Avoid non-named datatypes */
            if(!H5T_is_named((H5T_t *)obj_ptr))
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a named datatype")

            obj_path = H5T_nameof((H5T_t *)obj_ptr);
            break;

        default:
            HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "unknown data object type")
    } /* end switch */
    HDassert(obj_path);

    /* Retrieve a copy of the user path and put it into the buffer */
    if(obj_path->user_path_r) {
        ssize_t len = H5RS_len(obj_path->user_path_r);

        /* Set the user path, if given */
        if(user_path)
            HDstrcpy(user_path, H5RS_get_str(obj_path->user_path_r));

        /* Set the length of the path */
        *user_path_len = (size_t)len;

        /* Set the user path hidden flag */
        *obj_hidden = obj_path->obj_hidden;
    } /* end if */
    else {
        *user_path_len = 0;
        *obj_hidden = 0;
    } /* end else */

done:
    FUNC_LEAVE_NOAPI(ret_value)
}   /* H5G_user_path_test() */


/*-------------------------------------------------------------------------
 * Function:	H5G_verify_cached_stab_test
 *
 * Purpose:     Check that a that the provided group entry contains a
 *              cached symbol table entry, that the entry matches that in
 *              the provided group's object header, and check that the
 *              addresses are valid.
 *
 * Return:	Success:        Non-negative
 *		Failure:	Negative
 *
 * Programmer:	Neil Fortner
 *	        Mar  31, 2009
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5G_verify_cached_stab_test(H5O_loc_t *grp_oloc, H5G_entry_t *ent)
{
    H5O_stab_t  stab;                   /* Symbol table */
    H5HL_t      *heap = NULL;           /* Pointer to local heap */
    herr_t	ret_value = SUCCEED;    /* Return value */

    FUNC_ENTER_NOAPI_NOINIT(H5G_verify_cached_stab_test)

    /* Verify that stab info is cached in ent */
    if(ent->type != H5G_CACHED_STAB)
        HGOTO_ERROR(H5E_SYM, H5E_BADTYPE, FAIL, "symbol table information is not cached")

    /* Read the symbol table message from the group */
    if(NULL == H5O_msg_read(grp_oloc, H5O_STAB_ID, &stab, H5AC_ind_dxpl_id))
        HGOTO_ERROR(H5E_SYM, H5E_BADMESG, FAIL, "unable to read symbol table message")

    /* Verify that the cached symbol table info matches the symbol table message
     * in the object header */
    if((ent->cache.stab.btree_addr != stab.btree_addr)
            || (ent->cache.stab.heap_addr != stab.heap_addr))
        HGOTO_ERROR(H5E_SYM, H5E_BADVALUE, FAIL, "cached stab info does not match object header")

    /* Verify that the btree address is valid */
    if(H5B_valid(grp_oloc->file, H5AC_ind_dxpl_id, H5B_SNODE, stab.btree_addr) < 0)
        HGOTO_ERROR(H5E_BTREE, H5E_NOTFOUND, FAIL, "b-tree address is invalid")

    /* Verify that the heap address is valid */
    if(NULL == (heap = H5HL_protect(grp_oloc->file, H5AC_ind_dxpl_id, stab.heap_addr, H5AC_READ)))
        HGOTO_ERROR(H5E_HEAP, H5E_NOTFOUND, FAIL, "heap address is invalid")

done:
    /* Release resources */
    if(heap && H5HL_unprotect(grp_oloc->file, H5AC_ind_dxpl_id, heap, stab.heap_addr) < 0)
        HDONE_ERROR(H5E_SYM, H5E_PROTECT, FAIL, "unable to unprotect symbol table heap")

    FUNC_LEAVE_NOAPI(ret_value)
} /* end H5G_verify_cached_stab_test() */

