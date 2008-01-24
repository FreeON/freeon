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

/* $Id: H5P.c,v 1.140.2.18 2002/12/04 05:30:41 matzke Exp $ */

/* Private header files */
#include "H5private.h"		/* Generic Functions			*/
#include "H5Eprivate.h"		/* Error handling		  	*/
#include "H5FDprivate.h"	/* File drivers				*/
#include "H5FLprivate.h"	/* Free Lists	  */
#include "H5MMprivate.h"	/* Memory management			*/
#include "H5Iprivate.h"		/* IDs			  		*/
#include "H5Pprivate.h"		/* Property lists		  	*/

#define PABLO_MASK	H5P_mask

/* Is the interface initialized? */
static int		interface_initialize_g = 0;
#define INTERFACE_INIT H5P_init_interface
static herr_t		H5P_init_interface(void);

/*
 * Predefined data types. These are initialized at runtime by 
 * H5P_init_interface() in this source file.
 */
hid_t H5P_NO_CLASS_g            = FAIL;
hid_t H5P_FILE_CREATE_g         = FAIL;
hid_t H5P_FILE_ACCESS_g         = FAIL;
hid_t H5P_DATASET_CREATE_g      = FAIL;
hid_t H5P_DATASET_XFER_g        = FAIL;
hid_t H5P_MOUNT_g               = FAIL;

/* Local static functions */
static H5P_genclass_t *H5P_create_class(H5P_genclass_t *par_class,
     const char *name, unsigned hashsize, unsigned internal,
     H5P_cls_create_func_t cls_create, void *create_data,
     H5P_cls_close_func_t cls_close, void *close_data);
static herr_t H5P_close_list(void *_plist);
static herr_t H5P_close_class(void *_pclass);

/* Declare a free list to manage the H5P_t struct */
H5FL_DEFINE_STATIC(H5P_t);


/*-------------------------------------------------------------------------
 * Function:	H5P_init
 *
 * Purpose:	Initialize the interface from some other layer.
 *
 * Return:	Success:	non-negative
 *
 *		Failure:	negative
 *
 * Programmer:	Quincey Koziol
 *              Saturday, March 4, 2000
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5P_init(void)
{
    FUNC_ENTER(H5P_init, FAIL);
    /* FUNC_ENTER() does all the work */
    FUNC_LEAVE(SUCCEED);
}

/*--------------------------------------------------------------------------
NAME
    H5P_xor_name -- Generate an xor'ed value for a string
USAGE
    unsigned H5P_xor_name(s)
        const char *s;  IN: String to operate over
RETURNS
    Always returns valid value
DESCRIPTION
    Generates an xor'ed value for a string
--------------------------------------------------------------------------*/
static unsigned
H5P_xor_name(const char *s)
{
    unsigned ret=0;
    unsigned char temp;

    if(s!=NULL)
        while(*s!='\0') {
            temp=(ret>>24)&0xff;
            ret <<= 8;
            ret |= temp;
            ret ^= *s++;
        }

    return(ret);
}   /* end H5P_xor_name() */

/*--------------------------------------------------------------------------
NAME
    H5P_hash_name -- Generate a hash value for a string
USAGE
    unsigned H5P_hash_name(s, hashsize)
        const char *s;  IN: String to operate over
        unsigned;          IN: Size of hash table to clip against
RETURNS
    Always returns valid value
DESCRIPTION
    Generates a hash location based on an xor'ed value for a string
--------------------------------------------------------------------------*/
static unsigned
H5P_hash_name(const char *s, unsigned hashsize)
{
    return(H5P_xor_name(s)%hashsize);
}   /* end H5P_hash_name() */

/*--------------------------------------------------------------------------
NAME
   H5P_init_interface -- Initialize interface-specific information
USAGE
    herr_t H5P_init_interface()
   
RETURNS
    Non-negative on success/Negative on failure
DESCRIPTION
    Initializes any interface-specific data or routines.

--------------------------------------------------------------------------*/
static herr_t
H5P_init_interface(void)
{
    H5P_genclass_t  *root_class;    /* Pointer to root property list class created */
    H5P_genclass_t  *pclass;        /* Pointer to property list class to create */
    herr_t		    ret_value = SUCCEED;
    int		    i;
    herr_t		    status;

    FUNC_ENTER(H5P_init_interface, FAIL);

    assert(H5P_NCLASSES <= H5I_TEMPLATE_MAX - H5I_TEMPLATE_0);

    /*
     * Initialize the mappings between property list classes and atom
     * groups. We keep the two separate because property list classes are
     * publicly visible but atom groups aren't.
     */
    for (i = 0; i < H5P_NCLASSES; i++) {
        status = H5I_init_group((H5I_type_t)(H5I_TEMPLATE_0 +i),
                    H5I_TEMPID_HASHSIZE, 0, (H5I_free_t)H5P_close);
        if (status < 0)
            ret_value = FAIL;
    }

    if (ret_value < 0) {
        HRETURN_ERROR(H5E_ATOM, H5E_CANTINIT, FAIL,
		      "unable to initialize atom group");
    }
    
    /*
     * Initialize the Generic Property class & object groups.
     */
    if (H5I_init_group(H5I_GENPROP_CLS, H5I_GENPROPCLS_HASHSIZE, 0, (H5I_free_t)H5P_close_class) < 0)
        HRETURN_ERROR(H5E_ATOM, H5E_CANTINIT, FAIL,
		      "unable to initialize atom group");
    if (H5I_init_group(H5I_GENPROP_LST, H5I_GENPROPOBJ_HASHSIZE, 0, (H5I_free_t)H5P_close_list) < 0)
        HRETURN_ERROR(H5E_ATOM, H5E_CANTINIT, FAIL,
		      "unable to initialize atom group");

    /* Create root property list class */

    /* Allocate the root class */
    if (NULL==(root_class = H5P_create_class (NULL,"none",H5P_NO_CLASS_HASH_SIZE,1,NULL,NULL,NULL,NULL)))
        HRETURN_ERROR (H5E_PLIST, H5E_CANTINIT, FAIL, "class initialization failed");

    /* Register the root class */
    if ((H5P_NO_CLASS_g = H5I_register (H5I_GENPROP_CLS, root_class))<0)
        HRETURN_ERROR (H5E_PLIST, H5E_CANTREGISTER, FAIL, "can't register property list class");

    /* Register the file creation and file access property classes */

    /* Allocate the file creation class */
    if (NULL==(pclass = H5P_create_class (root_class,"file create",H5P_FILE_CREATE_HASH_SIZE,1,NULL,NULL,NULL,NULL)))
        HRETURN_ERROR (H5E_PLIST, H5E_CANTINIT, FAIL, "class initialization failed");

    /* Register the dataset creation class */
    if ((H5P_FILE_CREATE_g = H5I_register (H5I_GENPROP_CLS, pclass))<0)
        HRETURN_ERROR (H5E_PLIST, H5E_CANTREGISTER, FAIL, "can't register property list class");

    /* Allocate the file access class */
    if (NULL==(pclass = H5P_create_class (root_class,"file access",H5P_FILE_ACCESS_HASH_SIZE,1,NULL,NULL,NULL,NULL)))
        HRETURN_ERROR (H5E_PLIST, H5E_CANTINIT, FAIL, "class initialization failed");

    /* Register the file access class */
    if ((H5P_DATASET_XFER_g = H5I_register (H5I_GENPROP_CLS, pclass))<0)
        HRETURN_ERROR (H5E_PLIST, H5E_CANTREGISTER, FAIL, "can't register property list class");

    /* Register the dataset creation and data xfer property classes */

    /* Allocate the dataset creation class */
    if (NULL==(pclass = H5P_create_class (root_class,"dataset create",H5P_DATASET_CREATE_HASH_SIZE,1,NULL,NULL,NULL,NULL)))
        HRETURN_ERROR (H5E_PLIST, H5E_CANTINIT, FAIL, "class initialization failed");

    /* Register the dataset creation class */
    if ((H5P_DATASET_CREATE_g = H5I_register (H5I_GENPROP_CLS, pclass))<0)
        HRETURN_ERROR (H5E_PLIST, H5E_CANTREGISTER, FAIL, "can't register property list class");

    /* Allocate the data xfer class */
    if (NULL==(pclass = H5P_create_class (root_class,"data xfer",H5P_DATASET_XFER_HASH_SIZE,1,NULL,NULL,NULL,NULL)))
        HRETURN_ERROR (H5E_PLIST, H5E_CANTINIT, FAIL, "class initialization failed");

    /* Register the data xfer class */
    if ((H5P_DATASET_XFER_g = H5I_register (H5I_GENPROP_CLS, pclass))<0)
        HRETURN_ERROR (H5E_PLIST, H5E_CANTREGISTER, FAIL, "can't register property list class");

/* When do the "basic" properties for each of the library classes get added? */
/* Who adds them? */

    FUNC_LEAVE(ret_value);
}

/*--------------------------------------------------------------------------
 NAME
    H5P_term_interface
 PURPOSE
    Terminate various H5P objects
 USAGE
    void H5P_term_interface()
 RETURNS
    Non-negative on success/Negative on failure
 DESCRIPTION
    Release the atom group and any other resources allocated.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
     Can't report errors...
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
int
H5P_term_interface(void)
{
    int	i, n=0;

    if (interface_initialize_g) {
/* Destroy HDF5 library property classes */
        for (i=0; i<H5P_NCLASSES; i++)
            n += H5I_nmembers((H5I_type_t)(H5I_TEMPLATE_0+i));
        n += H5I_nmembers(H5I_GENPROP_CLS);
        n += H5I_nmembers(H5I_GENPROP_LST);
        if (n) {
            for (i=0; i<H5P_NCLASSES; i++)
                H5I_clear_group((H5I_type_t)(H5I_TEMPLATE_0+i), FALSE);
            H5I_clear_group(H5I_GENPROP_CLS, FALSE);
            H5I_clear_group(H5I_GENPROP_LST, FALSE);
        } else {
            for (i=0; i<H5P_NCLASSES; i++) {
                H5I_destroy_group((H5I_type_t)(H5I_TEMPLATE_0 + i));
                n++; /*H5I*/
            }

            H5I_destroy_group(H5I_GENPROP_CLS);
            n++; /*H5I*/
            H5I_destroy_group(H5I_GENPROP_LST);
            n++; /*H5I*/

            interface_initialize_g = 0;
        }
    }
    return n;
}


/*-------------------------------------------------------------------------
 * Function:	H5Pcreate
 *
 * Purpose:	Creates a new property list by copying a default property
 *		list.
 *
 * Return:	Success:	A new copy of a default property list.
 *
 *		Failure:	NULL
 *
 * Programmer:	Unknown
 *
 * Modifications:
 *		Robb Matzke, 1999-08-18
 *		Rewritten in terms of H5P_copy() to fix memory leaks.
 *-------------------------------------------------------------------------
 */
hid_t
H5Pcreate(H5P_class_t type)
{
    hid_t	ret_value = FAIL;
    const void	*src = NULL;
    H5P_t	*new_plist = NULL;

    FUNC_ENTER(H5Pcreate, FAIL);
    H5TRACE1("i","p",type);

    /* Allocate a new property list and initialize it with default values */
    switch (type) {
        case H5P_FILE_CREATE:
            src = &H5F_create_dflt;
            break;
        case H5P_FILE_ACCESS:
            src = &H5F_access_dflt;
            break;
        case H5P_DATASET_CREATE:
            src = &H5D_create_dflt;
            break;
        case H5P_DATASET_XFER:
            src = &H5D_xfer_dflt;
            break;
        case H5P_MOUNT:
            src = &H5F_mount_dflt;
            break;
        default:
            HRETURN_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL,
                  "unknown property list class");
    }

    /* Copy the property list */
    if (NULL==(new_plist=H5P_copy(type, src))) {
        HRETURN_ERROR(H5E_PLIST, H5E_CANTINIT, FAIL,
		      "unable to copy default property list");
    }
    
    /* Atomize the new property list */
    if ((ret_value = H5P_create(type, new_plist)) < 0) {
        HRETURN_ERROR(H5E_ATOM, H5E_CANTINIT, FAIL,
		      "unable to register property list");
    }
    
    FUNC_LEAVE(ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5P_create
 *
 * Purpose:	Given a pointer to some property list struct, atomize the
 *		property list and return its ID. The property list memory is
 *		not copied, so the caller should not free it; it will be
 *		freed by H5P_release().
 *
 * Return:	Success:	A new property list ID.
 *
 *		Failure:	Negative
 *
 * Programmer:	Robb Matzke
 *		Wednesday, December  3, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
hid_t
H5P_create(H5P_class_t type, H5P_t *plist)
{
    hid_t	ret_value = FAIL;

    FUNC_ENTER(H5P_create, FAIL);

    /* check args */
    assert(type >= 0 && type < H5P_NCLASSES);
    assert(plist);

    /* Atomize the new property list */
    if ((ret_value=H5I_register((H5I_type_t)(H5I_TEMPLATE_0+type), plist))<0) {
        HRETURN_ERROR(H5E_ATOM, H5E_CANTINIT, FAIL,
		      "unable to register property list");
    }
    
    FUNC_LEAVE(ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pclose
 *
 * Purpose:	Release access to a property list object, PLIST_ID.
 *
 * Return:	Success:	non-negative
 *
 *		Failure:	negative
 *
 * Programmer:	Unknown
 *
 * Modifications:
 * 		Robb Matzke, 1999-08-03
 *		Attempting to close H5P_DEFAULT is no longer an error, but
 *		rather a no-op.
 *-------------------------------------------------------------------------
 */
herr_t
H5Pclose(hid_t plist_id)
{
    FUNC_ENTER(H5Pclose, FAIL);
    H5TRACE1("e","i",plist_id);

    /* Check arguments */
    if (plist_id==H5P_DEFAULT)
        HRETURN(SUCCEED);
    if (H5P_get_class (plist_id)<0) {
        HRETURN_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");
    }
	
    /* When the reference count reaches zero the resources are freed */
    if (H5I_dec_ref(plist_id) < 0)
        HRETURN_ERROR(H5E_ATOM, H5E_BADATOM, FAIL, "problem freeing property list");

    FUNC_LEAVE (SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5P_close
 *
 * Purpose:	Closes a property list and frees the memory associated with
 *		the property list.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Wednesday, February 18, 1998
 *
 * Modifications:
 *		Robb Matzke, 1999-08-03
 *		Modified to work with the virtual file layer.
 *-------------------------------------------------------------------------
 */
herr_t
H5P_close(void *_plist)
{
    H5P_t           *plist=(H5P_t *)_plist;
    H5F_access_t	*fa_list = &(plist->u.faccess);
    H5D_xfer_t		*dx_list = &(plist->u.dxfer);
    H5D_create_t	*dc_list = &(plist->u.dcreate);
    
    FUNC_ENTER (H5P_close, FAIL);

    /* Check args */
    if (!plist)
        HRETURN (SUCCEED);

    /* Some property lists may need to do special things */
    switch (plist->cls) {
        case H5P_FILE_ACCESS:
            if (fa_list->driver_id>=0) {
                H5FD_fapl_free(fa_list->driver_id, fa_list->driver_info);
                H5I_dec_ref(fa_list->driver_id);
                fa_list->driver_info = NULL;
                fa_list->driver_id = -1;
            }
            break;
        
        case H5P_FILE_CREATE:
            break;
        
        case H5P_DATASET_CREATE:
            H5O_reset(H5O_FILL, &(dc_list->fill));
            H5O_reset(H5O_EFL, &(dc_list->efl));
            H5O_reset(H5O_PLINE, &(dc_list->pline));
            break;

        case H5P_DATASET_XFER:
            if (dx_list->driver_id>=0) {
                H5FD_dxpl_free(dx_list->driver_id, dx_list->driver_info);
                H5I_dec_ref(dx_list->driver_id);
                dx_list->driver_info = NULL;
                dx_list->driver_id = -1;
            }
            break;

        case H5P_MOUNT:
            break;

        default:
            HRETURN_ERROR (H5E_ARGS, H5E_BADTYPE, FAIL,
                   "unknown property list class");
    }

    /* Return the property list to the free list */
    H5FL_FREE(H5P_t,plist);

    FUNC_LEAVE(SUCCEED);
}


/*-------------------------------------------------------------------------
 * Function:	H5Pget_class
 *
 * Purpose:	Returns the class identifier for a property list.
 *
 * Return:	Success:	A property list class
 *
 *		Failure:	H5P_NO_CLASS (-1)
 *
 * Programmer:	Robb Matzke
 *		Wednesday, December  3, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
H5P_class_t
H5Pget_class(hid_t plist_id)
{
    H5I_type_t		    group;
    H5P_class_t		    ret_value = H5P_NO_CLASS;

    FUNC_ENTER(H5Pget_class, H5P_NO_CLASS);
    H5TRACE1("p","i",plist_id);

    if ((group = H5I_get_type(plist_id)) < 0 ||
            group >= H5I_TEMPLATE_MAX ||
            group < H5I_TEMPLATE_0) {
        HRETURN_ERROR(H5E_ATOM, H5E_BADATOM, H5P_NO_CLASS,
		      "not a property list");
    }

    ret_value = (H5P_class_t)(group - H5I_TEMPLATE_0);

    FUNC_LEAVE(ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5P_get_class
 *
 * Purpose:	Internal function for getting the property list class.
 *
 * Return:	Success:	A property list class
 *
 *		Failure:	H5P_NO_CLASS (-1)
 *
 * Programmer:	Robb Matzke
 *              Tuesday, July  7, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
H5P_class_t
H5P_get_class(hid_t plist_id)
{
    H5I_type_t		    group;
    H5P_class_t		    ret_value = H5P_NO_CLASS;

    FUNC_ENTER(H5P_get_class, H5P_NO_CLASS);

    if ((group = H5I_get_type(plist_id)) < 0 ||
            group >= H5I_TEMPLATE_MAX ||
            group < H5I_TEMPLATE_0) {
        HRETURN_ERROR(H5E_ATOM, H5E_BADATOM, H5P_NO_CLASS,
		      "not a property list");
    }

    ret_value = (H5P_class_t)(group - H5I_TEMPLATE_0);
    FUNC_LEAVE(ret_value);
}
    

/*-------------------------------------------------------------------------
 * Function:	H5Pcopy
 *
 * Purpose:	Deep-copies a property list PLIST_ID.
 *
 * Return:	Success:	The ID of the new copy of the property list.
 *				The ID will be different than the input ID
 *				since the new ID refers to a completely
 *				separate copy of the the structure that the
 *				original ID points to.
 *
 *		Failure:	Negative
 *
 * Programmer:	Unknown
 *
 * Modifications:
 *		Robb Matzke, 1999-08-03
 *		If PLIST_ID is H5P_DEFAULT then we return H5P_DEFAULT.
 *-------------------------------------------------------------------------
 */
hid_t
H5Pcopy(hid_t plist_id)
{
    const void		   *plist = NULL;
    void		   *new_plist = NULL;
    H5P_class_t		    type;
    hid_t		    ret_value = FAIL;
    H5I_type_t		    group;

    FUNC_ENTER(H5Pcopy, FAIL);
    H5TRACE1("i","i",plist_id);

    if (H5P_DEFAULT==plist_id)
        HRETURN(H5P_DEFAULT);

    /* Check args */
    if (NULL == (plist = H5I_object(plist_id)) ||
	(type = H5P_get_class(plist_id)) < 0 ||
	(group = H5I_get_type(plist_id)) < 0) {
	HRETURN_ERROR(H5E_ATOM, H5E_BADATOM, FAIL,
		      "unable to unatomize property list");
    }

    /* Copy it */
    if (NULL==(new_plist=H5P_copy (type, plist))) {
	HRETURN_ERROR (H5E_INTERNAL, H5E_CANTINIT, FAIL,
		       "unable to copy property list");
    }

    /* Register the atom for the new property list */
    if ((ret_value = H5I_register(group, new_plist)) < 0) {
	HRETURN_ERROR(H5E_ATOM, H5E_CANTREGISTER, FAIL,
		      "unable to atomize property list pointer");
    }
    FUNC_LEAVE(ret_value);
}


/*-------------------------------------------------------------------------
 * Function:	H5P_copy
 *
 * Purpose:	Creates a new property list and initializes it with some
 *		other property list.
 *
 * Return:	Success:	Ptr to new property list
 *
 *		Failure:	NULL
 *
 * Programmer:	Robb Matzke
 *		Tuesday, February  3, 1998
 *
 * Modifications:
 *		Robb Matzke, 1999-08-03
 *		Modified to use the virtual file layer.
 *-------------------------------------------------------------------------
 */
void *
H5P_copy (H5P_class_t type, const void *src)
{
    size_t		size;
    H5P_t		*dst = NULL;
    const H5D_create_t	*dc_src = NULL;
    H5D_create_t	*dc_dst = NULL;
    H5F_access_t	*fa_dst = NULL;
    H5D_xfer_t		*dx_dst = NULL;
    
    FUNC_ENTER (H5P_copy, NULL);
    
    /* How big is the property list */
    switch (type) {
        case H5P_FILE_CREATE:
            size = sizeof(H5F_create_t);
            break;

        case H5P_FILE_ACCESS:
            size = sizeof(H5F_access_t);
            break;

        case H5P_DATASET_CREATE:
            size = sizeof(H5D_create_t);
            break;

        case H5P_DATASET_XFER:
            size = sizeof(H5D_xfer_t);
            break;

        case H5P_MOUNT:
            size = sizeof(H5F_mprop_t);
            break;

        default:
            HRETURN_ERROR(H5E_ARGS, H5E_BADRANGE, NULL,
                  "unknown property list class");
    }

    /* Create the new property list */
    if (NULL==(dst = H5FL_ALLOC(H5P_t,0))) {
        HRETURN_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL,
               "memory allocation failed");
    }

    /* Copy into new object */
    HDmemcpy(dst, src, size);

    /* Set the type of the property list */
    dst->cls=type;

    /* Deep-copy pointers */
    switch (type) {
        case H5P_FILE_CREATE:
            break;
        
        case H5P_FILE_ACCESS:
            fa_dst = (H5F_access_t*)dst;

            if (fa_dst->driver_id>=0) {
                H5I_inc_ref(fa_dst->driver_id);
                fa_dst->driver_info = H5FD_fapl_copy(fa_dst->driver_id,
                                fa_dst->driver_info);
            }
            break;
        
        case H5P_DATASET_CREATE:
            dc_src = (const H5D_create_t*)src;
            dc_dst = (H5D_create_t*)dst;

            /* Copy the fill value */
            if (NULL==H5O_copy(H5O_FILL, &(dc_src->fill), &(dc_dst->fill))) {
                HRETURN_ERROR(H5E_PLIST, H5E_CANTINIT, NULL,
                      "unabe to copy fill value message");
            }
            
            /* Copy the external file list */
            HDmemset(&(dc_dst->efl), 0, sizeof(dc_dst->efl));
            if (NULL==H5O_copy(H5O_EFL, &(dc_src->efl), &(dc_dst->efl))) {
                HRETURN_ERROR(H5E_PLIST, H5E_CANTINIT, NULL,
                      "unable to copy external file list message");
            }

            /* Copy the filter pipeline */
            if (NULL==H5O_copy(H5O_PLINE, &(dc_src->pline), &(dc_dst->pline))) {
                HRETURN_ERROR(H5E_PLIST, H5E_CANTINIT, NULL,
                      "unable to copy filter pipeline message");
            }
            break;
        
        case H5P_DATASET_XFER:
            dx_dst = (H5D_xfer_t*)dst;

            if (dx_dst->driver_id>=0) {
                H5I_inc_ref(dx_dst->driver_id);
                dx_dst->driver_info = H5FD_dxpl_copy(dx_dst->driver_id,
                                dx_dst->driver_info);
            }
            break;

        case H5P_MOUNT:
            /* Nothing to do */
            break;

        default:
            HRETURN_ERROR(H5E_ARGS, H5E_BADRANGE, NULL,
                  "unknown property list class");
    }

    FUNC_LEAVE (dst);
}

/*=========================================*/
/* Start of generic property list routines */
/*=========================================*/


/*--------------------------------------------------------------------------
 NAME
    H5P_copy_prop
 PURPOSE
    Internal routine to copy a property
 USAGE
    H5P_genprop_t *H5P_copy_prop(oprop)
        H5P_genprop_t *oprop;   IN: Pointer to property to copy
 RETURNS
    Returns a pointer to the newly created property on success,
        NULL on failure.
 DESCRIPTION
    Allocates memory and copies property information into a new property object.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static H5P_genprop_t *
H5P_copy_prop(H5P_genprop_t *oprop)
{
    H5P_genprop_t *prop=NULL;        /* Pointer to new property copied */
    H5P_genprop_t *ret_value=NULL;   /* Return value */

    FUNC_ENTER (H5P_copy_prop, NULL);

    assert(oprop);

    /* Allocate the new property */
    if (NULL==(prop = H5MM_malloc (sizeof(H5P_genprop_t))))
        HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");

    /* Copy basic property information */
    HDmemcpy(prop,oprop,sizeof(H5P_genprop_t));

    /* Duplicate name */
    prop->name = HDstrdup(oprop->name);

    /* Duplicate current value, if it exists */
    if(oprop->value!=NULL) {
        if (NULL==(prop->value = H5MM_malloc (prop->size)))
            HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");
        HDmemcpy(prop->value,oprop->value,prop->size);
    } /* end if */

    /* Duplicate default value, if it exists */
    if(oprop->def_value!=NULL) {
        if (NULL==(prop->def_value = H5MM_malloc (prop->size)))
            HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");
        HDmemcpy(prop->def_value,oprop->def_value,prop->size);
    } /* end if */

    /* Reset the link to the next property */
    prop->next=NULL;

    /* Set return value */
    ret_value=prop;

done:
    /* Free any resources allocated */
    if(ret_value==NULL) {
        if(prop!=NULL) {
            if(prop->name!=NULL)
                H5MM_xfree(prop->name);
            if(prop->value!=NULL)
                H5MM_xfree(prop->value);
            if(prop->def_value!=NULL)
                H5MM_xfree(prop->def_value);
            H5MM_xfree(prop);
        } /* end if */
    } /* end if */

    FUNC_LEAVE (ret_value);
}   /* H5P_copy_prop() */


/*--------------------------------------------------------------------------
 NAME
    H5P_create_prop
 PURPOSE
    Internal routine to create a new property
 USAGE
    H5P_genprop_t *H5P_create_prop(name,size,def_value,prp_create,prp_set,
            prp_get,prp_close)
        const char *name;       IN: Name of property to register
        size_t size;            IN: Size of property in bytes
        void *def_value;        IN: Pointer to buffer containing default value
                                    for property in newly created property lists
        H5P_prp_create_func_t prp_create;   IN: Function pointer to property
                                    creation callback
        H5P_prp_set_func_t prp_set; IN: Function pointer to property set callback
        H5P_prp_get_func_t prp_get; IN: Function pointer to property get callback
        H5P_prp_close_func_t prp_close; IN: Function pointer to property close
                                    callback
 RETURNS
    Returns a pointer to the newly created property on success,
        NULL on failure.
 DESCRIPTION
    Allocates memory and copies property information into a new property object.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static H5P_genprop_t *
H5P_create_prop(const char *name, size_t size, void *def_value, void *value,
    H5P_prp_create_func_t prp_create, H5P_prp_set_func_t prp_set,
    H5P_prp_get_func_t prp_get, H5P_prp_close_func_t prp_close)
{
    H5P_genprop_t *prop=NULL;        /* Pointer to new property copied */
    H5P_genprop_t *ret_value=NULL;   /* Return value */

    FUNC_ENTER (H5P_create_prop, NULL);

    assert(name);
    assert((size>0 && def_value!=NULL) || (size==0));

    /* Allocate the new property */
    if (NULL==(prop = H5MM_malloc (sizeof(H5P_genprop_t))))
        HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");

    /* Set the property initial values */
    prop->xor_val = H5P_xor_name(name); /* Generate the XOR'd value for the name */
    prop->name = HDstrdup(name); /* Duplicate name */
    prop->size=size;
    /* Duplicate value, if it exists */
    if(value!=NULL) {
        if (NULL==(prop->value = H5MM_malloc (prop->size)))
            HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");
        HDmemcpy(prop->value,value,prop->size);
    } /* end if */
    else
        prop->value=NULL;


    /* Duplicate default value, if it exists */
    if(def_value!=NULL) {
        if (NULL==(prop->def_value = H5MM_malloc (prop->size)))
            HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");
        HDmemcpy(prop->def_value,def_value,prop->size);
    } /* end if */
    else
        prop->def_value=NULL;

    /* Set the function pointers */
    prop->create=prp_create;
    prop->set=prp_set;
    prop->get=prp_get;
    prop->close=prp_close;

    /* Reset the link to the next property */
    prop->next=NULL;

    /* Set return value */
    ret_value=prop;

done:
    /* Free any resources allocated */
    if(ret_value==NULL) {
        if(prop!=NULL) {
            if(prop->name!=NULL)
                H5MM_xfree(prop->name);
            if(prop->value!=NULL)
                H5MM_xfree(prop->value);
            if(prop->def_value!=NULL)
                H5MM_xfree(prop->def_value);
            H5MM_xfree(prop);
        } /* end if */
    } /* end if */

    FUNC_LEAVE (ret_value);
}   /* H5P_create_prop() */


/*--------------------------------------------------------------------------
 NAME
    H5P_add_prop
 PURPOSE
    Internal routine to insert a property into a property hash table
 USAGE
    herr_t H5P_add_prop(hash, hashsize, prop)
        H5P_gen_prop_t *hash[]; IN/OUT: Pointer to array of properties for hash table
        unsigned hashsize;         IN: Size of hash table
        H5P_genprop_t *prop;    IN: Pointer to property to insert
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
    Inserts a property into a hash table of properties, using the hashed
    property name.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5P_add_prop(H5P_genprop_t *hash[], unsigned hashsize, H5P_genprop_t *prop)
{
    unsigned loc;                  /* Hash table location */
    herr_t ret_value=SUCCEED;   /* Return value */

    FUNC_ENTER (H5P_add_prop, FAIL);

    assert(hash);
    assert(hashsize>0);
    assert(prop);

    /* Get correct hash table location */
    loc=H5P_hash_name(prop->name,hashsize);
    
    /* Insert property into hash table */
    prop->next=hash[loc];
    hash[loc]=prop;

#ifdef LATER
done:
#endif /* LATER */
    FUNC_LEAVE (ret_value);
}   /* H5P_add_prop() */


/*--------------------------------------------------------------------------
 NAME
    H5P_find_prop
 PURPOSE
    Internal routine to check for a property in a hash table
 USAGE
    H5P_genprop_t *H5P_find_prop(hash, hashsize, name)
        H5P_genprop_t *hash[];  IN: Pointer to array of properties for hash table
        unsigned hashsize;         IN: Size of hash table
        const char *name;       IN: Name of property to check for
 RETURNS
    Returns pointer to property on success, NULL on failure.
 DESCRIPTION
        Checks for a property in a hash table of properties.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static H5P_genprop_t *
H5P_find_prop(H5P_genprop_t *hash[], unsigned hashsize, const char *name)
{
    H5P_genprop_t *ret_value;   /* Property pointer return value */
    unsigned loc;                  /* Hash table location */
    unsigned xor_val;              /* XOR'ed value of the name to look for */

    FUNC_ENTER (H5P_add_prop, NULL);

    assert(hash);
    assert(hashsize>0);
    assert(name);

    /* Get correct hash table location */
    loc=H5P_hash_name(name,hashsize);
    
    /* Get the XOR'ed value for the name to search for, to speed up comparisons */
    xor_val=H5P_xor_name(name);
    
    /* Locate property in list */
    ret_value=hash[loc];
    while(ret_value!=NULL) {
        /* Check for name matching */
        if(ret_value->xor_val==xor_val && HDstrcmp(ret_value->name,name)==0)
            break;
    } /* end while */

    FUNC_LEAVE (ret_value);
}   /* H5P_find_prop() */


/*--------------------------------------------------------------------------
 NAME
    H5P_free_prop
 PURPOSE
    Internal routine to destroy a property node
 USAGE
    herr_t H5P_free_prop(prop)
        H5P_genprop_t *prop;    IN: Pointer to property to destroy
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
    Releases all the memory for a property list.  Does _not_ call the
    properties 'close' callback, that should already have been done.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5P_free_prop(H5P_genprop_t *prop)
{
    herr_t ret_value=SUCCEED;   /* Return value */

    FUNC_ENTER (H5P_free_prop, FAIL);

    assert(prop);

    /* Release the property value and default value if they exist */
    if(prop->size>0) {
        if(prop->value)
            H5MM_xfree(prop->value);
        if(prop->def_value)
            H5MM_xfree(prop->def_value);
    } /* end if */
    H5MM_xfree(prop->name);
    H5MM_xfree(prop);

#ifdef LATER
done:
#endif /* LATER */
    FUNC_LEAVE (ret_value);
}   /* H5P_free_prop() */


/*--------------------------------------------------------------------------
 NAME
    H5P_free_all_prop
 PURPOSE
    Internal routine to remove all properties from a property hash table
 USAGE
    herr_t H5P_free_all_prop(hash, hashsize, make_cb)
        H5P_gen_prop_t *hash[]; IN/OUT: Pointer to array of properties for hash table
        unsigned hashsize;         IN: Size of hash table
        unsigned make_cb;          IN: Whether to make property callbacks or not
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Remove all the properties from a property list.  Calls the property
    'close' callback for each property removed.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5P_free_all_prop(H5P_genprop_t *hash[], unsigned hashsize, unsigned make_cb)
{
    H5P_genprop_t *tprop, *next;/* Temporary pointer to properties */
    unsigned u;                    /* Local index variable */
    herr_t ret_value=SUCCEED;   /* Return value */

    FUNC_ENTER (H5P_free_all_prop, FAIL);

    assert(hash);
    assert(hashsize>0);

    /* Work through all the properties... */
    for(u=0; u<hashsize; u++) {
        tprop=hash[u];
        while(tprop!=NULL) {
            /* Find next node right away, to avoid accessing the current node after it's been free'd */
            next=tprop->next;

            /* Call the close callback and ignore the return value, there's nothing we can do about it */
            if(make_cb && tprop->close!=NULL)
                (tprop->close)(tprop->name,&(tprop->value));

            /* Free the property, ignoring return value, nothing we can do */
            H5P_free_prop(tprop);

            tprop=next;
        } /* end while */
    } /* end for */

#ifdef LATER
done:
#endif /* LATER */
    FUNC_LEAVE (ret_value);
}   /* H5P_free_all_prop() */


/*--------------------------------------------------------------------------
 NAME
    H5P_access_class
 PURPOSE
    Internal routine to increment or decrement list & class dependancies on a
        property list class
 USAGE
    herr_t H5P_access_class(pclass,mod)
        H5P_genclass_t *pclass;     IN: Pointer to class to modify
        H5P_class_mod_t mod;        IN: Type of modification to class
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Increment/Decrement the class or list dependancies for a given class.
    This routine is the final arbiter on decisions about actually releasing a
    class in memory, such action is only taken when the reference counts for
    both dependent classes & lists reach zero.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t
H5P_access_class(H5P_genclass_t *pclass, H5P_class_mod_t mod)
{
    herr_t ret_value=SUCCEED;   /* Return value */

    FUNC_ENTER (H5P_access_class, FAIL);

    assert(pclass);
    assert(mod>H5P_MOD_ERR && mod<H5P_MOD_MAX);

    switch(mod) {
        case H5P_MOD_INC_CLS:        /* Increment the dependant class count*/
            pclass->classes++;
            break;

        case H5P_MOD_DEC_CLS:        /* Decrement the dependant class count*/
            pclass->classes--;
            break;

        case H5P_MOD_INC_LST:        /* Increment the dependant list count*/
            pclass->plists++;
            break;

        case H5P_MOD_DEC_LST:        /* Decrement the dependant list count*/
            pclass->plists--;
            break;

        case H5P_MOD_INC_REF:        /* Increment the ID reference count*/
            pclass->ref_count++;
            break;

        case H5P_MOD_DEC_REF:        /* Decrement the ID reference count*/
            pclass->ref_count--;
            break;

        case H5P_MOD_CHECK:         /* NOOP, just check if we can delete the class */
            break;
        
        case H5P_MOD_ERR:
        case H5P_MOD_MAX:
            assert(0 && "Invalid H5P class modification");
    } /* end switch */

    /* Check if we can release the class information now */
    if(pclass->deleted && pclass->plists==0 && pclass->classes==0 ) {
        assert(pclass->name);
        H5MM_xfree(pclass->name);

        /* Free the class properties without making callbacks */
        H5P_free_all_prop(pclass->props,pclass->hashsize,0);

        H5MM_xfree(pclass);
    } /* end if */

#ifdef LATER
done:
#endif /* LATER */
    FUNC_LEAVE (ret_value);
}   /* H5P_access_class() */


/*--------------------------------------------------------------------------
 NAME
    H5P_create_class
 PURPOSE
    Internal routine to create a new property list class.
 USAGE
    H5P_genclass_t H5P_create_class(par_class, name, hashsize, internal,
                cls_create, create_data, cls_close, close_data)
        H5P_genclass_t *par_class;  IN: Pointer to parent class
        const char *name;       IN: Name of class we are creating
        unsigned hashsize; IN: Number of buckets in hash table
        unsigned internal; IN: Whether this is an internal class or not
        H5P_cls_create_func_t;  IN: The callback function to call when each
                                    property list in this class is created.
        void *create_data;      IN: Pointer to user data to pass along to class
                                    creation callback.
        H5P_cls_close_func_t;   IN: The callback function to call when each
                                    property list in this class is closed.
        void *close_data;       IN: Pointer to user data to pass along to class
                                    close callback.
 RETURNS
    Returns a pointer to the newly created property list class on success,
        NULL on failure.
 DESCRIPTION
    Allocates memory and attaches a class to the property list class hierarchy.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static H5P_genclass_t *
H5P_create_class(H5P_genclass_t *par_class, const char *name, unsigned hashsize, unsigned internal,
    H5P_cls_create_func_t cls_create, void *create_data,
    H5P_cls_close_func_t cls_close, void *close_data
    )
{
    H5P_genclass_t *pclass;             /* Property list class created */
    H5P_genclass_t *ret_value=NULL;     /* return value */

    FUNC_ENTER (H5P_create_class, NULL);

    assert(name);
    /* Allow internal classes to break some rules */
    /* (This allows the root of the tree to be created with this routine -QAK) */
    if(!internal) {
        assert(par_class);
        assert(hashsize>0);
    }

    /* Allocate room for the class & it's hash table of properties */
    if (NULL==(pclass = H5MM_calloc (sizeof(H5P_genclass_t)+((hashsize-1)*sizeof(H5P_genprop_t *)))))
        HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL,"memory allocation failed");

    /* Set class state */
    pclass->parent = par_class;
    pclass->name = HDstrdup(name);
    pclass->nprops = 0;     /* Classes are created without properties initially */
    pclass->hashsize = hashsize;
    pclass->plists = 0;     /* No properties lists of this class yet */
    pclass->classes = 0;    /* No classes derived from this class yet */
    pclass->ref_count = 1;  /* This is the first reference to the new class */
    pclass->internal = internal;
    pclass->deleted = 0;    /* Not deleted yet... :-) */

    /* Set callback functions and pass-along data */
    pclass->create_func = cls_create;
    pclass->create_data = create_data;
    pclass->close_func = cls_close;
    pclass->close_data = close_data;

    /* Increment parent class's derived class value */
    if(par_class!=NULL)
        if(H5P_access_class(par_class,H5P_MOD_INC_CLS)<0)
            HGOTO_ERROR (H5E_PLIST, H5E_CANTINIT, NULL,"Can't increment parent class ref count");

    /* Set return value */
    ret_value=pclass;

done:
    /* Free any resources allocated */
    if(ret_value==NULL) {
        if(pclass!=NULL)
            H5MM_xfree(pclass);
    }

    FUNC_LEAVE (ret_value);
}   /* H5P_create_class() */


/*--------------------------------------------------------------------------
 NAME
    H5Pcreate_class
 PURPOSE
    Create a new property list class.
 USAGE
    hid_t H5Pcreate_class(parent, name, hashsize, cls_create, create_data,
                cls_close, close_data)
        hid_t parent;       IN: Property list class ID of parent class
        const char *name;   IN: Name of class we are creating
        unsigned hashsize;     IN: Number of buckets in hash table
        H5P_cls_create_func_t cls_create;   IN: The callback function to call
                                    when each property list in this class is
                                    created.
        void *create_data;  IN: Pointer to user data to pass along to class
                                    creation callback.
        H5P_cls_close_func_t cls_close;     IN: The callback function to call
                                    when each property list in this class is
                                    closed.
        void *close_data;   IN: Pointer to user data to pass along to class
                                    close callback.
 RETURNS
    Returns a valid property list class ID on success, NULL on failure.
 DESCRIPTION
    Allocates memory and attaches a class to the property list class hierarchy.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
hid_t
H5Pcreate_class(hid_t parent, const char *name, unsigned hashsize,
    H5P_cls_create_func_t cls_create, void *create_data,
    H5P_cls_close_func_t cls_close, void *close_data
    )
{
    H5P_genclass_t	*par_class = NULL;  /* Pointer to the parent class */
    H5P_genclass_t	*pclass = NULL;     /* Property list class created */
    hid_t	ret_value = FAIL;           /* Return value			*/

    FUNC_ENTER(H5Pcreate_class, FAIL);
    H5TRACE7("i","isIuxxxx",parent,name,hashsize,cls_create,create_data,
             cls_close,close_data);

    /* Check arguments. */
    if (H5P_DEFAULT!=parent && (H5I_GENPROP_CLS!=H5I_get_type(parent)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list class");
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid class name");
    if (hashsize==0)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "hashsize too small");
    if ((create_data!=NULL && cls_create==NULL) || (close_data!=NULL && cls_close==NULL))
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "data specified, but no callback provided");
    
    /* Get the pointer to the parent class */
    if(parent==H5P_DEFAULT)
        par_class=NULL;
    else if (NULL == (par_class = H5I_object(parent)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "can't retrieve parent class");

    /* Create the new property list class */
    if (NULL==(pclass=H5P_create_class(par_class, name, hashsize, 0, cls_create, create_data, cls_close, close_data)))
        HGOTO_ERROR(H5E_PLIST, H5E_CANTCREATE, FAIL, "unable to create property list class");

    /* Get an atom for the class */
    if ((ret_value = H5I_register(H5I_GENPROP_CLS, pclass))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to atomize property list class");

done:
    if (ret_value<0 && pclass)
        H5P_close_class(pclass);

    FUNC_LEAVE(ret_value);
}   /* H5Pcreate_class() */


/*--------------------------------------------------------------------------
 NAME
    H5P_create_list
 PURPOSE
    Internal routine to create a new property list of a property list class.
 USAGE
    H5P_genplist_t *H5P_create_list(class)
        H5P_genclass_t *class;  IN: Property list class create list from
 RETURNS
    Returns a pointer to the newly created property list on success,
        NULL on failure.
 DESCRIPTION
        Creates a property list of a given class.  If a 'create' callback
    exists for the property list class, it is called before the
    property list is passed back to the user.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
        If this routine is called from a library routine other than
    H5Pcreate_list, the calling routine is responsible for getting an ID for
    the property list and calling the class 'create' callback (if one exists)
    and also setting the "class_init" flag.
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static H5P_genplist_t *H5P_create_list(H5P_genclass_t *pclass)
{
    H5P_genclass_t *tclass=NULL;        /* Temporary class pointer */
    H5P_genplist_t *plist=NULL;         /* New property list created */
    H5P_genplist_t *ret_value=NULL;     /* return value */
    H5P_genprop_t *tmp;                 /* Temporary pointer to parent class properties */
    H5P_genprop_t *pcopy;               /* Copy of property to insert into class */
    unsigned u;                            /* Local index variable */

    FUNC_ENTER (H5P_create_list, NULL);

    assert(pclass);

    /* 
     * Create new property list object
     */

    /* Allocate room for the property list & it's hash table of properties */
    if (NULL==(plist = H5MM_calloc (sizeof(H5P_genplist_t)+((pclass->hashsize-1)*sizeof(H5P_genprop_t *)))))
        HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL,"memory allocation failed");

    /* Set class state */
    plist->pclass = pclass;
    plist->nprops = 0;      /* Initially the plist has the same number of properties as the class */
    plist->class_init = 0;  /* Initially, wait until the class callback finishes to set */

    /*
     * Copy class properties (up through list of parent classes also),
     * initialize each with default value & make property 'create' callback.
     */
    tclass=pclass;
    while(tclass!=NULL) {
        if(tclass->nprops>0) {
            /* Walk through the hash table */
            for(u=0; u<tclass->hashsize; u++) {
                tmp=tclass->props[u];
                /* Walk through the list of properties at each hash location */
                while(tmp!=NULL) {
                    /* Check for property already existing in list */
                    if(H5P_find_prop(plist->props,tclass->hashsize,tmp->name)==NULL) {
                        /* Make a copy of the class's property */
                        if((pcopy=H5P_copy_prop(tmp))==NULL)
                            HGOTO_ERROR (H5E_PLIST, H5E_CANTCOPY, NULL,"Can't copy property");

                        /* Create initial value from default value for non-zero sized properties */
                        if(pcopy->size>0) {
                            /* Properties from the class should have any values yet, but should have a default */
                            assert(pcopy->value==NULL);
                            assert(pcopy->def_value);

                            /* Allocate space for the property value & copy default value */
                            if (NULL==(pcopy->value = H5MM_malloc (pcopy->size)))
                                HGOTO_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");

                            HDmemcpy(pcopy->value,pcopy->def_value,pcopy->size);
                        } /* end if */

                        /* Call property creation callback, if it exists */
                        if(pcopy->create) {
                            if((pcopy->create)(pcopy->name,&(pcopy->value))<0) {
                                H5P_free_prop(pcopy);
                                HGOTO_ERROR (H5E_PLIST, H5E_CANTINIT, NULL,"Can't initialize property");
                            } /* end if */
                        } /* end if */

                        /* Insert the initialized property into the property list */
                        if(H5P_add_prop(plist->props,tclass->hashsize,pcopy)<0)
                            HGOTO_ERROR (H5E_PLIST, H5E_CANTINSERT, NULL,"Can't insert property into class");
                    } /* end if */

                    /* Go to next registered property in class */
                    tmp=tmp->next;
                } /* end while */
            } /* end for */
        } /* end if */

        /* Go up to parent class */
        tclass=tclass->parent;
    } /* end while */

    /* Increment the number of property lists derived from class */
    if(H5P_access_class(plist->pclass,H5P_MOD_INC_LST)<0)
        HGOTO_ERROR (H5E_PLIST, H5E_CANTINIT, NULL,"Can't increment class ref count");

    /* Set return value */
    ret_value=plist;

done:
    /* Release resources allocated on failure */
    if(ret_value==NULL) {
        if(plist!=NULL) {
            /* Close & free all the properties */
            H5P_free_all_prop(plist->props,pclass->hashsize,1);

            /* Decrement the number of property lists derived from the class */
            pclass->plists--;
        } /* end if */
    } /* end if */

    FUNC_LEAVE (ret_value);
}   /* H5P_create_list() */


/*--------------------------------------------------------------------------
 NAME
    H5Pcreate_list
 PURPOSE
    Routine to create a new property list of a property list class.
 USAGE
    hid_t H5Pcreate_list(cls_id)
        hid_t cls_id;       IN: Property list class create list from
 RETURNS
    Returns a valid property list ID on success, FAIL on failure.
 DESCRIPTION
        Creates a property list of a given class.  If a 'create' callback
    exists for the property list class, it is called before the
    property list is passed back to the user.  If 'create' callbacks exist for
    any individual properties in the property list, they are called before the
    class 'create' callback.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
hid_t H5Pcreate_list(hid_t cls_id)
{
    H5P_genclass_t	*pclass;   /* Property list class to modify */
    H5P_genplist_t	*plist=NULL;    /* Property list created */
    hid_t plist_id=FAIL;       /* Property list ID */
    hid_t ret_value=FAIL;      /* return value */

    FUNC_ENTER (H5Pcreate_list, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_CLS != H5I_get_type(cls_id) || NULL == (pclass = H5I_object(cls_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list class");

    /* Create the new property list */
    if ((plist=H5P_create_list(pclass))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTCREATE, FAIL, "unable to create property list");

    /* Get an atom for the property list */
    if ((plist_id = H5I_register(H5I_GENPROP_LST, plist))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to atomize property list");

    /* Call the class callback (if it exists) now that we have the property list ID */
    if(plist->pclass->create_func!=NULL) {
        if((plist->pclass->create_func)(plist_id,plist->pclass->create_data)<0) {
            /* Delete ID, ignore return value */
            H5I_remove(plist_id);
            HGOTO_ERROR (H5E_PLIST, H5E_CANTINIT, FAIL,"Can't initialize property");
        } /* end if */
    } /* end if */

    /* Set the class initialization flag */
    plist->class_init=1;

    /* Set the return value */
    ret_value=plist_id;

done:
    if (ret_value<0 && plist)
        H5P_close_list(plist);

    FUNC_LEAVE (ret_value);
}   /* H5Pcreate_list() */


/*--------------------------------------------------------------------------
 NAME
    H5P_register
 PURPOSE
    Internal routine to register a new property in a property list class.
 USAGE
    herr_t H5P_register(class, name, size, default, prp_create, prp_set, prp_get, prp_close)
        H5P_genclass_t *class;  IN: Property list class to close
        const char *name;       IN: Name of property to register
        size_t size;            IN: Size of property in bytes
        void *def_value;        IN: Pointer to buffer containing default value
                                    for property in newly created property lists
        H5P_prp_create_func_t prp_create;   IN: Function pointer to property
                                    creation callback
        H5P_prp_set_func_t prp_set; IN: Function pointer to property set callback
        H5P_prp_get_func_t prp_get; IN: Function pointer to property get callback
        H5P_prp_close_func_t prp_close; IN: Function pointer to property close
                                    callback
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Registers a new property with a property list class.  The property will
    exist in all property list objects of that class after this routine is
    finished.  The name of the property must not already exist.  The default
    property value must be provided and all new property lists created with this
    property will have the property value set to the default provided.  Any of
    the callback routines may be set to NULL if they are not needed.

        Zero-sized properties are allowed and do not store any data in the
    property list.  These may be used as flags to indicate the presence or
    absence of a particular piece of information.  The 'default' pointer for a
    zero-sized property may be set to NULL.  The property 'create' & 'close'
    callbacks are called for zero-sized properties, but the 'set' and 'get'
    callbacks are never called.

        The 'create' callback is called when a new property list with this
    property is being created.  H5P_prp_create_func_t is defined as:
        typedef herr_t (*H5P_prp_create_func_t)(hid_t prop_id, const char *name,
                void **initial_value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being created.
        const char *name;   IN: The name of the property being modified.
        void **initial_value;   IN/OUT: The initial value for the property being created.
                                (The 'default' value passed to H5Pregister)
    The 'create' routine may modify the value to be set and those changes will
    be stored as the initial value of the property.  If the 'create' routine
    returns a negative value, the new property value is not copied into the
    property and the property list creation routine returns an error value.

        The 'set' callback is called before a new value is copied into the
    property.  H5P_prp_set_func_t is defined as:
        typedef herr_t (*H5P_prp_set_func_t)(hid_t prop_id, const char *name,
            void **value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being modified.
        const char *name;   IN: The name of the property being modified.
        void *new_value;    IN/OUT: The value being set for the property.
    The 'set' routine may modify the value to be set and those changes will be
    stored as the value of the property.  If the 'set' routine returns a
    negative value, the new property value is not copied into the property and
    the property list set routine returns an error value.

        The 'get' callback is called before a value is retrieved from the
    property.  H5P_prp_get_func_t is defined as:
        typedef herr_t (*H5P_prp_get_func_t)(hid_t prop_id, const char *name,
            void *value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being queried.
        const char *name;   IN: The name of the property being queried.
        void *value;        IN/OUT: The value being retrieved for the property.
    The 'get' routine may modify the value to be retrieved and those changes
    will be returned to the calling function.  If the 'get' routine returns a
    negative value, the property value is returned and the property list get
    routine returns an error value.

        The 'close' callback is called when a property list with this
    property is being destroyed.  H5P_prp_close_func_t is defined as:
        typedef herr_t (*H5P_prp_close_func_t)(const char *name, void *value);
    where the parameters to the callback function are:
        const char *name;   IN: The name of the property being closed.
        void *value;        IN: The value of the property being closed.
    The 'close' routine may modify the value passed in, but the value is not
    used by the library when the 'close' routine returns.  If the
    'close' routine returns a negative value, the property list close
    routine returns an error value but the property list is still closed.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
        The 'set' callback function may be useful to range check the value being
    set for the property or may perform some tranformation/translation of the
    value set.  The 'get' callback would then [probably] reverse the
    transformation, etc.  A single 'get' or 'set' callback could handle
    multiple properties by performing different actions based on the property
    name or other properties in the property list.

        I would like to say "the property list is not closed" when a 'close'
    routine fails, but I don't think that's possible due to other properties in
    the list being successfully closed & removed from the property list.  I
    suppose that it would be possible to just remove the properties which have
    successful 'close' callbacks, but I'm not happy with the ramifications
    of a mangled, un-closable property list hanging around...  Any comments? -QAK

 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_register(H5P_genclass_t *pclass, const char *name, size_t size,
    void *def_value, H5P_prp_create_func_t prp_create, H5P_prp_set_func_t prp_set,
    H5P_prp_get_func_t prp_get, H5P_prp_close_func_t prp_close)
{
    H5P_genclass_t *new_class; /* New class pointer */
    H5P_genprop_t *tmp_prop;   /* Temporary property pointer */
    H5P_genprop_t *new_prop=NULL;   /* Temporary property pointer */
    H5P_genprop_t *pcopy;      /* Property copy */
    unsigned u;                   /* Local index variable */
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5P_register, FAIL);

    assert(pclass);
    assert(name);
    assert((size>0 && def_value!=NULL) || (size==0));

    /* Check for duplicate named properties */
    if((tmp_prop=H5P_find_prop(pclass->props,pclass->hashsize,name))!=NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_EXISTS, FAIL, "property already exists");

    /* Check if class needs to be split because property lists or classes have
     *  been created since the last modification was made to the class.
     */
    if(pclass->plists>0 || pclass->classes>0) {
        if((new_class=H5P_create_class(pclass->parent,pclass->name,pclass->hashsize,
                pclass->internal,pclass->create_func,pclass->create_data,
                pclass->close_func,pclass->close_data))==NULL)
            HGOTO_ERROR(H5E_PLIST, H5E_CANTCOPY, FAIL, "can't copy class");

        /* Walk through the hash table of the old class and copy properties */
        for(u=0; u<pclass->hashsize; u++) {
            tmp_prop=pclass->props[u];
            while(tmp_prop!=NULL) {
                /* Make a copy of the class's property */
                if((pcopy=H5P_copy_prop(tmp_prop))==NULL)
                    HGOTO_ERROR (H5E_PLIST, H5E_CANTCOPY, FAIL,"Can't copy property");

                /* Insert the property into the new property class */
                if(H5P_add_prop(new_class->props,pclass->hashsize,pcopy)<0)
                    HGOTO_ERROR (H5E_PLIST, H5E_CANTINSERT, FAIL,"Can't insert property into class");

                /* Go to next registered property in class */
                tmp_prop=tmp_prop->next;
            } /* end while */
        } /* end for */

        /* Use the new class instead of the old one */
        pclass=new_class;
    } /* end if */

    /* Create property object from parameters */
    if((new_prop=H5P_create_prop(name,size,def_value,NULL,prp_create,prp_set,prp_get,prp_close))==NULL)
        HGOTO_ERROR (H5E_PLIST, H5E_CANTCREATE, FAIL,"Can't create property");

    /* Insert property into property list class */
    if(H5P_add_prop(pclass->props,pclass->hashsize,new_prop)<0)
        HGOTO_ERROR (H5E_PLIST, H5E_CANTINSERT, FAIL,"Can't insert property into class");

    /* Increment property count for class */
    pclass->nprops++;

    /* Set return value */
    ret_value=SUCCEED;

done:
    if(ret_value==FAIL) {
        if(new_prop!=NULL) {
            if(new_prop->name!=NULL)
                H5MM_xfree(new_prop->name);
            if(new_prop->value!=NULL)
                H5MM_xfree(new_prop->value);
            if(new_prop->def_value!=NULL)
                H5MM_xfree(new_prop->def_value);
            H5MM_xfree(new_prop);
        } /* end if */
    }  /* end if */
    FUNC_LEAVE (ret_value);
}   /* H5P_register() */


/*--------------------------------------------------------------------------
 NAME
    H5Pregister
 PURPOSE
    Routine to register a new property in a property list class.
 USAGE
    herr_t H5Pregister(class, name, size, default, prp_create, prp_set, prp_get, prp_close)
        hid_t class;            IN: Property list class to close
        const char *name;       IN: Name of property to register
        size_t size;            IN: Size of property in bytes
        void *def_value;        IN: Pointer to buffer containing default value
                                    for property in newly created property lists
        H5P_prp_create_func_t prp_create;   IN: Function pointer to property
                                    creation callback
        H5P_prp_set_func_t prp_set; IN: Function pointer to property set callback
        H5P_prp_get_func_t prp_get; IN: Function pointer to property get callback
        H5P_prp_close_func_t prp_close; IN: Function pointer to property close
                                    callback
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Registers a new property with a property list class.  The property will
    exist in all property list objects of that class after this routine is
    finished.  The name of the property must not already exist.  The default
    property value must be provided and all new property lists created with this
    property will have the property value set to the default provided.  Any of
    the callback routines may be set to NULL if they are not needed.

        Zero-sized properties are allowed and do not store any data in the
    property list.  These may be used as flags to indicate the presence or
    absence of a particular piece of information.  The 'default' pointer for a
    zero-sized property may be set to NULL.  The property 'create' & 'close'
    callbacks are called for zero-sized properties, but the 'set' and 'get'
    callbacks are never called.

        The 'create' callback is called when a new property list with this
    property is being created.  H5P_prp_create_func_t is defined as:
        typedef herr_t (*H5P_prp_create_func_t)(hid_t prop_id, const char *name,
                void **initial_value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being created.
        const char *name;   IN: The name of the property being modified.
        void **initial_value;   IN/OUT: The initial value for the property being created.
                                (The 'default' value passed to H5Pregister)
    The 'create' routine may modify the value to be set and those changes will
    be stored as the initial value of the property.  If the 'create' routine
    returns a negative value, the new property value is not copied into the
    property and the property list creation routine returns an error value.

        The 'set' callback is called before a new value is copied into the
    property.  H5P_prp_set_func_t is defined as:
        typedef herr_t (*H5P_prp_set_func_t)(hid_t prop_id, const char *name,
            void **value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being modified.
        const char *name;   IN: The name of the property being modified.
        void **new_value;   IN/OUT: The value being set for the property.
    The 'set' routine may modify the value to be set and those changes will be
    stored as the value of the property.  If the 'set' routine returns a
    negative value, the new property value is not copied into the property and
    the property list set routine returns an error value.

        The 'get' callback is called before a value is retrieved from the
    property.  H5P_prp_get_func_t is defined as:
        typedef herr_t (*H5P_prp_get_func_t)(hid_t prop_id, const char *name,
            void *value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being queried.
        const char *name;   IN: The name of the property being queried.
        void **value;       IN/OUT: The value being retrieved for the property.
    The 'get' routine may modify the value to be retrieved and those changes
    will be returned to the calling function.  If the 'get' routine returns a
    negative value, the property value is returned and the property list get
    routine returns an error value.

        The 'close' callback is called when a property list with this
    property is being destroyed.  H5P_prp_close_func_t is defined as:
        typedef herr_t (*H5P_prp_close_func_t)(const char *name, void *value);
    where the parameters to the callback function are:
        const char *name;   IN: The name of the property being closed.
        void *value;        IN: The value of the property being closed.
    The 'close' routine may modify the value passed in, but the value is not
    used by the library when the 'close' routine returns.  If the
    'close' routine returns a negative value, the property list close
    routine returns an error value but the property list is still closed.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
        The 'set' callback function may be useful to range check the value being
    set for the property or may perform some tranformation/translation of the
    value set.  The 'get' callback would then [probably] reverse the
    transformation, etc.  A single 'get' or 'set' callback could handle
    multiple properties by performing different actions based on the property
    name or other properties in the property list.

        I would like to say "the property list is not closed" when a 'close'
    routine fails, but I don't think that's possible due to other properties in
    the list being successfully closed & removed from the property list.  I
    suppose that it would be possible to just remove the properties which have
    successful 'close' callbacks, but I'm not happy with the ramifications
    of a mangled, un-closable property list hanging around...  Any comments? -QAK

 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Pregister(hid_t cls_id, const char *name, size_t size, void *def_value,
    H5P_prp_create_func_t prp_create, H5P_prp_set_func_t prp_set,
    H5P_prp_get_func_t prp_get, H5P_prp_close_func_t prp_close)
{
    H5P_genclass_t	*pclass;   /* Property list class to modify */
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Pregister, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_CLS != H5I_get_type(cls_id) || NULL == (pclass = H5I_object(cls_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list class");
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid class name");
    if (size>0 && def_value==NULL)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "properties >0 size must have default");

    /* Create the new property list class */
    if ((ret_value=H5P_register(pclass,name,size,def_value,prp_create,prp_set,prp_get,prp_close))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to register property in class");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pregister() */


/*--------------------------------------------------------------------------
 NAME
    H5P_insert
 PURPOSE
    Internal routine to insert a new property in a property list.
 USAGE
    herr_t H5P_insert(plist, name, size, value, prp_set, prp_get, prp_close)
        H5P_genplist_t *plist;  IN: Property list to add property to
        const char *name;       IN: Name of property to add
        size_t size;            IN: Size of property in bytes
        void *value;            IN: Pointer to the value for the property
        H5P_prp_set_func_t prp_set; IN: Function pointer to property set callback
        H5P_prp_get_func_t prp_get; IN: Function pointer to property get callback
        H5P_prp_close_func_t prp_close; IN: Function pointer to property close
                                    callback
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Inserts a temporary property into a property list.  The property will
    exist only in this property list object.  The name of the property must not
    already exist.  The value must be provided unless the property is zero-
    sized.  Any of the callback routines may be set to NULL if they are not
    needed.

        Zero-sized properties are allowed and do not store any data in the
    property list.  These may be used as flags to indicate the presence or
    absence of a particular piece of information.  The 'value' pointer for a
    zero-sized property may be set to NULL.  The property 'close' callback is
    called for zero-sized properties, but the 'set' and 'get' callbacks are
    never called.

        The 'set' callback is called before a new value is copied into the
    property.  H5P_prp_set_func_t is defined as:
        typedef herr_t (*H5P_prp_set_func_t)(hid_t prop_id, const char *name,
            void **value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being modified.
        const char *name;   IN: The name of the property being modified.
        void **new_value;   IN/OUT: The value being set for the property.
    The 'set' routine may modify the value to be set and those changes will be
    stored as the value of the property.  If the 'set' routine returns a
    negative value, the new property value is not copied into the property and
    the property list set routine returns an error value.

        The 'get' callback is called before a value is retrieved from the
    property.  H5P_prp_get_func_t is defined as:
        typedef herr_t (*H5P_prp_get_func_t)(hid_t prop_id, const char *name,
            void *value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being queried.
        const char *name;   IN: The name of the property being queried.
        void **value;       IN/OUT: The value being retrieved for the property.
    The 'get' routine may modify the value to be retrieved and those changes
    will be returned to the calling function.  If the 'get' routine returns a
    negative value, the property value is returned and the property list get
    routine returns an error value.

        The 'close' callback is called when a property list with this
    property is being destroyed.  H5P_prp_close_func_t is defined as:
        typedef herr_t (*H5P_prp_close_func_t)(const char *name,
            void *value);
    where the parameters to the callback function are:
        const char *name;   IN: The name of the property being closed.
        void *value;        IN: The value of the property being closed.
    The 'close' routine may modify the value passed in, but the value is not
    used by the library when the 'close' routine returns.  If the
    'close' routine returns a negative value, the property list close
    routine returns an error value but the property list is still closed.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
        The 'set' callback function may be useful to range check the value being
    set for the property or may perform some tranformation/translation of the
    value set.  The 'get' callback would then [probably] reverse the
    transformation, etc.  A single 'get' or 'set' callback could handle
    multiple properties by performing different actions based on the property
    name or other properties in the property list.
        
        There is no 'create' callback routine for temporary property list
    objects, the initial value is assumed to have any necessary setup already
    performed on it.

        I would like to say "the property list is not closed" when a 'close'
    routine fails, but I don't think that's possible due to other properties in
    the list being successfully closed & removed from the property list.  I
    suppose that it would be possible to just remove the properties which have
    successful 'close' callbacks, but I'm not happy with the ramifications
    of a mangled, un-closable property list hanging around...  Any comments? -QAK

 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_insert(H5P_genplist_t *plist, const char *name, size_t size,
    void *value, H5P_prp_set_func_t prp_set, H5P_prp_get_func_t prp_get,
    H5P_prp_close_func_t prp_close)
{
    H5P_genprop_t *new_prop=NULL;   /* Temporary property pointer */
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5P_insert, FAIL);

    assert(plist);
    assert(name);
    assert((size>0 && value!=NULL) || (size==0));

    /* Check for duplicate named properties */
    if(H5P_find_prop(plist->props,plist->pclass->hashsize,name)!=NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_EXISTS, FAIL, "property already exists");

    /* Create property object from parameters */
    if((new_prop=H5P_create_prop(name,size,NULL,value,NULL,prp_set,prp_get,prp_close))==NULL)
        HGOTO_ERROR (H5E_PLIST, H5E_CANTCREATE, FAIL,"Can't create property");

    /* Insert property into property list class */
    if(H5P_add_prop(plist->props,plist->pclass->hashsize,new_prop)<0)
        HGOTO_ERROR (H5E_PLIST, H5E_CANTINSERT, FAIL,"Can't insert property into class");

    /* Increment property count for class */
    plist->nprops++;

    /* Set return value */
    ret_value=SUCCEED;

done:
    if(ret_value==FAIL) {
        if(new_prop!=NULL) {
            if(new_prop->name!=NULL)
                H5MM_xfree(new_prop->name);
            if(new_prop->value!=NULL)
                H5MM_xfree(new_prop->value);
            if(new_prop->def_value!=NULL)
                H5MM_xfree(new_prop->def_value);
            H5MM_xfree(new_prop);
        } /* end if */
    }  /* end if */
    FUNC_LEAVE (ret_value);
}   /* H5P_insert() */


/*--------------------------------------------------------------------------
 NAME
    H5Pinsert
 PURPOSE
    Routine to insert a new property in a property list.
 USAGE
    herr_t H5Pinsert(plist, name, size, value, prp_set, prp_get, prp_close)
        hid_t plist;            IN: Property list to add property to
        const char *name;       IN: Name of property to add
        size_t size;            IN: Size of property in bytes
        void *value;            IN: Pointer to the value for the property
        H5P_prp_set_func_t prp_set; IN: Function pointer to property set callback
        H5P_prp_get_func_t prp_get; IN: Function pointer to property get callback
        H5P_prp_close_func_t prp_close; IN: Function pointer to property close
                                    callback
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Inserts a temporary property into a property list.  The property will
    exist only in this property list object.  The name of the property must not
    already exist.  The value must be provided unless the property is zero-
    sized.  Any of the callback routines may be set to NULL if they are not
    needed.

        Zero-sized properties are allowed and do not store any data in the
    property list.  These may be used as flags to indicate the presence or
    absence of a particular piece of information.  The 'value' pointer for a
    zero-sized property may be set to NULL.  The property 'close' callback is
    called for zero-sized properties, but the 'set' and 'get' callbacks are
    never called.

        The 'set' callback is called before a new value is copied into the
    property.  H5P_prp_set_func_t is defined as:
        typedef herr_t (*H5P_prp_set_func_t)(hid_t prop_id, const char *name,
            void **value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being modified.
        const char *name;   IN: The name of the property being modified.
        void **new_value;   IN/OUT: The value being set for the property.
    The 'set' routine may modify the value to be set and those changes will be
    stored as the value of the property.  If the 'set' routine returns a
    negative value, the new property value is not copied into the property and
    the property list set routine returns an error value.

        The 'get' callback is called before a value is retrieved from the
    property.  H5P_prp_get_func_t is defined as:
        typedef herr_t (*H5P_prp_get_func_t)(hid_t prop_id, const char *name,
            void *value);
    where the parameters to the callback function are:
        hid_t prop_id;      IN: The ID of the property list being queried.
        const char *name;   IN: The name of the property being queried.
        void **value;       IN/OUT: The value being retrieved for the property.
    The 'get' routine may modify the value to be retrieved and those changes
    will be returned to the calling function.  If the 'get' routine returns a
    negative value, the property value is returned and the property list get
    routine returns an error value.

        The 'close' callback is called when a property list with this
    property is being destroyed.  H5P_prp_close_func_t is defined as:
        typedef herr_t (*H5P_prp_close_func_t)(const char *name, void *value);
    where the parameters to the callback function are:
        const char *name;   IN: The name of the property being closed.
        void *value;        IN: The value of the property being closed.
    The 'close' routine may modify the value passed in, but the value is not
    used by the library when the 'close' routine returns.  If the
    'close' routine returns a negative value, the property list close
    routine returns an error value but the property list is still closed.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
        The 'set' callback function may be useful to range check the value being
    set for the property or may perform some tranformation/translation of the
    value set.  The 'get' callback would then [probably] reverse the
    transformation, etc.  A single 'get' or 'set' callback could handle
    multiple properties by performing different actions based on the property
    name or other properties in the property list.
        
        There is no 'create' callback routine for temporary property list
    objects, the initial value is assumed to have any necessary setup already
    performed on it.

        I would like to say "the property list is not closed" when a 'close'
    routine fails, but I don't think that's possible due to other properties in
    the list being successfully closed & removed from the property list.  I
    suppose that it would be possible to just remove the properties which have
    successful 'close' callbacks, but I'm not happy with the ramifications
    of a mangled, un-closable property list hanging around...  Any comments? -QAK

 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Pinsert(hid_t plist_id, const char *name, size_t size, void *value,
    H5P_prp_set_func_t prp_set, H5P_prp_get_func_t prp_get,
    H5P_prp_close_func_t prp_close)
{
    H5P_genplist_t	*plist;    /* Property list to modify */
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Pinsert, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(plist_id) || NULL == (plist = H5I_object(plist_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid property name");
    if (size>0 && value==NULL)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "properties >0 size must have default");

    /* Create the new property list class */
    if ((ret_value=H5P_insert(plist,name,size,value,prp_set,prp_get,prp_close))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to register property in plist");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pinsert() */


/*--------------------------------------------------------------------------
 NAME
    H5P_set
 PURPOSE
    Internal routine to set a property's value in a property list.
 USAGE
    herr_t H5P_set(plist, name, value)
        hid_t plist_id;         IN: Property list to find property in
        const char *name;       IN: Name of property to set
        void *value;            IN: Pointer to the value for the property
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Sets a new value for a property in a property list.  The property name
    must exist or this routine will fail.  If there is a 'set' callback routine
    registered for this property, the 'value' will be passed to that routine and
    any changes to the 'value' will be used when setting the property value.
    The information pointed at by the 'value' pointer (possibly modified by the
    'set' callback) is copied into the property list value and may be changed
    by the application making the H5Pset call without affecting the property
    value.

        If the 'set' callback routine returns an error, the property value will
    not be modified.  This routine may not be called for zero-sized properties
    and will return an error in that case.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_set(hid_t plist_id, const char *name, void *value)
{
    H5P_genplist_t *plist;      /* Property list to modify */
    H5P_genprop_t *prop;        /* Temporary property pointer */
    herr_t ret_value=FAIL;      /* return value */

    FUNC_ENTER (H5P_set, FAIL);

    assert(name);
    assert(value);

    /* Argument checking */
    if (H5I_GENPROP_LST != H5I_get_type(plist_id) || NULL == (plist = H5I_object(plist_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");

    /* Find property */
    if((prop=H5P_find_prop(plist->props,plist->pclass->hashsize,name))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "property doesn't exist");

    /* Check for property size >0 */
    if(prop->size==0)
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "property has zero size");

    /* Make a copy of the value and pass to 'set' callback */
    if(prop->set!=NULL) {
        void *tmp_value;            /* Temporary value for property */

        /* Make a copy of the current value, in case the callback fails */
        if (NULL==(tmp_value=H5MM_malloc(prop->size)))
            HRETURN_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed temporary property value");
        HDmemcpy(tmp_value,prop->value,prop->size);

        /* Call user's callback */
        if((*(prop->set))(plist_id,name,tmp_value)<0)
            HRETURN_ERROR(H5E_PLIST, H5E_CANTINIT, FAIL, "can't set property value");

        /* Copy new [possibly unchanged] value into property value */
        HDmemcpy(prop->value,tmp_value,prop->size);

        /* Free the temporary value buffer */
        H5MM_xfree(tmp_value);
    } /* end if */
    /* No 'set' callback, just copy value */
    else
        HDmemcpy(prop->value,value,prop->size);

    /* Set return value */
    ret_value=SUCCEED;

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_set() */


/*--------------------------------------------------------------------------
 NAME
    H5Pset
 PURPOSE
    Routine to set a property's value in a property list.
 USAGE
    herr_t H5P_set(plist_id, name, value)
        hid_t plist_id;         IN: Property list to find property in
        const char *name;       IN: Name of property to set
        void *value;            IN: Pointer to the value for the property
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Sets a new value for a property in a property list.  The property name
    must exist or this routine will fail.  If there is a 'set' callback routine
    registered for this property, the 'value' will be passed to that routine and
    any changes to the 'value' will be used when setting the property value.
    The information pointed at by the 'value' pointer (possibly modified by the
    'set' callback) is copied into the property list value and may be changed
    by the application making the H5Pset call without affecting the property
    value.

        If the 'set' callback routine returns an error, the property value will
    not be modified.  This routine may not be called for zero-sized properties
    and will return an error in that case.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Pset(hid_t plist_id, const char *name, void *value)
{
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Pset, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(plist_id))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid property name");
    if (value==NULL)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalied property value");

    /* Create the new property list class */
    if ((ret_value=H5P_set(plist_id,name,value))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to set value in plist");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pset() */


/*--------------------------------------------------------------------------
 NAME
    H5P_exist_plist
 PURPOSE
    Internal routine to query the existance of a property in a property list.
 USAGE
    herr_t H5P_exist_plist(plist, name)
        H5P_genplist_t *plist;  IN: Property list to check
        const char *name;       IN: Name of property to check for
 RETURNS
    Success: Positive if the property exists in the property list, zero
            if the property does not exist.
    Failure: negative value
 DESCRIPTION
        This routine checks if a property exists within a property list.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static htri_t H5P_exist_plist(H5P_genplist_t *plist, const char *name)
{
    htri_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5P_exist_plist, FAIL);

    assert(plist);
    assert(name);

    /* Check for property in property list */
    if(H5P_find_prop(plist->props,plist->pclass->hashsize,name)==NULL)
        ret_value=0;
    else
        ret_value=1;

    FUNC_LEAVE (ret_value);
}   /* H5P_exist_plist() */


/*--------------------------------------------------------------------------
 NAME
    H5P_exist_pclass
 PURPOSE
    Internal routine to query the existance of a property in a property class.
 USAGE
    herr_t H5P_exist_pclass(pclass, name)
        H5P_genclass_t *pclass;  IN: Property class to check
        const char *name;       IN: Name of property to check for
 RETURNS
    Success: Positive if the property exists in the property list, zero
            if the property does not exist.
    Failure: negative value
 DESCRIPTION
        This routine checks if a property exists within a property list.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static htri_t H5P_exist_pclass(H5P_genclass_t *pclass, const char *name)
{
    htri_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5P_exist_pclass, FAIL);

    assert(pclass);
    assert(name);

    /* Check for property in property list */
    if(H5P_find_prop(pclass->props,pclass->hashsize,name)==NULL)
        ret_value=0;
    else
        ret_value=1;

    FUNC_LEAVE (ret_value);
}   /* H5P_exist_pclass() */


/*--------------------------------------------------------------------------
 NAME
    H5Pexist
 PURPOSE
    Routine to query the existance of a property in a property object.
 USAGE
    htri_t H5P_exist(id, name)
        hid_t id;           IN: Property object ID to check
        const char *name;   IN: Name of property to check for
 RETURNS
    Success: Positive if the property exists in the property object, zero
            if the property does not exist.
    Failure: negative value
 DESCRIPTION
        This routine checks if a property exists within a property list or
    class.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
htri_t H5Pexist(hid_t id, const char *name)
{
    H5P_genplist_t	*plist;    /* Property list to query */
    H5P_genclass_t	*pclass;   /* Property class to query */
    htri_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Pexist, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(id) && H5I_GENPROP_CLS != H5I_get_type(id))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property object");
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid property name");

    /* Check for the existance of the property in the list or class */
    if(H5I_GENPROP_LST == H5I_get_type(id)) {
        if (NULL == (plist = H5I_object(id)))
            HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");
        if ((ret_value=H5P_exist_plist(plist,name))<0)
            HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to existance of property in plist");
    } /* end if */
    else 
        if(H5I_GENPROP_CLS == H5I_get_type(id)) {
            if (NULL == (pclass = H5I_object(id)))
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property class");
            if ((ret_value=H5P_exist_pclass(pclass,name))<0)
                HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to existance of property in pclass");
        } /* end if */
        else
            HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property object");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pexist() */


/*--------------------------------------------------------------------------
 NAME
    H5P_get_size_plist
 PURPOSE
    Internal routine to query the size of a property in a property list.
 USAGE
    herr_t H5P_get_size_plist(plist, name)
        H5P_genplist_t *plist;  IN: Property list to check
        const char *name;       IN: Name of property to query
        size_t *size;           OUT: Size of property
 RETURNS
    Success: non-negative value
    Failure: negative value
 DESCRIPTION
        This routine retrieves the size of a property's value in bytes.  Zero-
    sized properties are allowed and return a value of 0.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_get_size_plist(H5P_genplist_t *plist, const char *name, size_t *size)
{
    H5P_genprop_t *prop;        /* Temporary property pointer */
    herr_t ret_value=SUCCEED;      /* return value */

    FUNC_ENTER (H5P_get_size_plist, FAIL);

    assert(plist);
    assert(name);
    assert(size);

    /* Find property */
    if((prop=H5P_find_prop(plist->props,plist->pclass->hashsize,name))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "property doesn't exist");

    /* Get property size */
    *size=prop->size;

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_get_size_plist() */


/*--------------------------------------------------------------------------
 NAME
    H5P_get_size_pclass
 PURPOSE
    Internal routine to query the size of a property in a property class.
 USAGE
    herr_t H5P_get_size_pclass(pclass, name)
        H5P_genclass_t *pclass; IN: Property class to check
        const char *name;       IN: Name of property to query
        size_t *size;           OUT: Size of property
 RETURNS
    Success: non-negative value
    Failure: negative value
 DESCRIPTION
        This routine retrieves the size of a property's value in bytes.  Zero-
    sized properties are allowed and return a value of 0.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_get_size_pclass(H5P_genclass_t *pclass, const char *name, size_t *size)
{
    H5P_genprop_t *prop;        /* Temporary property pointer */
    herr_t ret_value=SUCCEED;      /* return value */

    FUNC_ENTER (H5P_get_size_pclass, FAIL);

    assert(pclass);
    assert(name);
    assert(size);

    /* Find property */
    if((prop=H5P_find_prop(pclass->props,pclass->hashsize,name))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "property doesn't exist");

    /* Get property size */
    *size=prop->size;

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_get_size_pclass() */


/*--------------------------------------------------------------------------
 NAME
    H5Pget_size
 PURPOSE
    Routine to query the size of a property in a property list or class.
 USAGE
    herr_t H5Pget_size(id, name)
        hid_t id;               IN: ID of property list or class to check
        const char *name;       IN: Name of property to query
        size_t *size;           OUT: Size of property
 RETURNS
    Success: non-negative value
    Failure: negative value
 DESCRIPTION
        This routine retrieves the size of a property's value in bytes.  Zero-
    sized properties are allowed and return a value of 0.  This function works
    for both property lists and classes.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Pget_size(hid_t id, const char *name, size_t *size)
{
    H5P_genclass_t	*pclass;   /* Property class to query */
    H5P_genplist_t	*plist;    /* Property list to query */
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Pget_size, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(id) && H5I_GENPROP_CLS != H5I_get_type(id))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property object");
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid property name");
    if (size==NULL)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid property size");

    if (H5I_GENPROP_LST == H5I_get_type(id)) {
        if (NULL == (plist = H5I_object(id)))
            HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");

        /* Check the property size */
        if ((ret_value=H5P_get_size_plist(plist,name,size))<0)
            HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to query size in plist");
    } /* end if */
    else 
        if (H5I_GENPROP_CLS == H5I_get_type(id)) {
            if (NULL == (pclass = H5I_object(id)))
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");

            /* Check the property size */
            if ((ret_value=H5P_get_size_pclass(pclass,name,size))<0)
                HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to query size in plist");
        } /* end if */
        else 
            HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property object");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pget_size() */


/*--------------------------------------------------------------------------
 NAME
    H5P_get_class
 PURPOSE
    Internal routine to query the class of a generic property list
 USAGE
    H5P_genclass_t *H5P_get_class(plist)
        H5P_genplist_t *plist;    IN: Property list to check
 RETURNS
    Success: Pointer to the class for a property list
    Failure: NULL
 DESCRIPTION
    This routine retrieves a pointer to the class for a property list.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static H5P_genclass_t *H5P_get_class_new(H5P_genplist_t *plist)
{
    H5P_genclass_t *ret_value=NULL;      /* return value */

    FUNC_ENTER (H5P_get_class, NULL);

    assert(plist);

    /* Get property size */
    ret_value=plist->pclass;

    FUNC_LEAVE (ret_value);
}   /* H5P_get_class() */


/*--------------------------------------------------------------------------
 NAME
    H5Pget_class
 PURPOSE
    Routine to query the class of a generic property list
 USAGE
    hid_t H5Pget_class(plist_id)
        hid_t plist_id;         IN: Property list to query
 RETURNS
    Success: ID of class object
    Failure: negative
 DESCRIPTION
    This routine retrieves the class of a property list.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
    Change the name of this function to H5Pget_class (and remove old H5Pget_class)
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
hid_t H5Pget_class_new(hid_t plist_id)
{
    H5P_genplist_t	*plist;         /* Property list to query */
    H5P_genclass_t	*pclass=NULL;   /* Property list class */
    hid_t ret_value=FAIL;           /* return value */

    FUNC_ENTER (H5Pget_class, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(plist_id) || NULL == (plist = H5I_object(plist_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");

    /* Retrieve the property list class */
    if ((pclass=H5P_get_class_new(plist))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "unable to query class of property list");

    /* Increment the outstanding references to the class object */
    if(H5P_access_class(pclass,H5P_MOD_INC_REF)<0)
        HGOTO_ERROR (H5E_PLIST, H5E_CANTINIT, FAIL,"Can't increment class ID ref count");

    /* Get an atom for the class */
    if ((ret_value = H5I_register(H5I_GENPROP_CLS, pclass))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to atomize property list class");

done:
    if (ret_value<0 && pclass)
        H5P_close_class(pclass);

    FUNC_LEAVE (ret_value);
}   /* H5Pget_class_name() */


/*--------------------------------------------------------------------------
 NAME
    H5P_get_nprops_plist
 PURPOSE
    Internal routine to query the number of properties in a property list
 USAGE
    herr_t H5P_get_nprops_plist(plist, nprops)
        H5P_genplist_t *plist;  IN: Property list to check
        size_t *nprops;         OUT: Number of properties in the property list
 RETURNS
    Success: non-negative value
    Failure: negative value
 DESCRIPTION
        This routine retrieves the number of a properties in a property list.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_get_nprops_plist(H5P_genplist_t *plist, size_t *nprops)
{
    herr_t ret_value=SUCCEED;      /* return value */

    FUNC_ENTER (H5P_get_nprops_plist, FAIL);

    assert(plist);
    assert(nprops);

    /* Get property size */
    *nprops=plist->nprops;

    FUNC_LEAVE (ret_value);
}   /* H5P_get_nprops_plist() */


/*--------------------------------------------------------------------------
 NAME
    H5P_get_nprops_pclass
 PURPOSE
    Internal routine to query the number of properties in a property class
 USAGE
    herr_t H5P_get_nprops_pclass(pclass, nprops)
        H5P_genclass_t *pclass;  IN: Property class to check
        size_t *nprops;         OUT: Number of properties in the property list
 RETURNS
    Success: non-negative value
    Failure: negative value
 DESCRIPTION
        This routine retrieves the number of a properties in a property class.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_get_nprops_pclass(H5P_genclass_t *pclass, size_t *nprops)
{
    herr_t ret_value=SUCCEED;      /* return value */

    FUNC_ENTER (H5P_get_nprops_pclass, FAIL);

    assert(pclass);
    assert(nprops);

    /* Get property size */
    *nprops=pclass->nprops;

    FUNC_LEAVE (ret_value);
}   /* H5P_get_nprops_pclass() */


/*--------------------------------------------------------------------------
 NAME
    H5Pget_nprops
 PURPOSE
    Routine to query the size of a property in a property list or class.
 USAGE
    herr_t H5Pget_nprops(id, nprops)
        hid_t id;               IN: ID of Property list or class to check
        size_t *nprops;         OUT: Number of properties in the property object
 RETURNS
    Success: non-negative value
    Failure: negative value
 DESCRIPTION
        This routine retrieves the number of properties in a property list or
    class.  If a property class ID is given, the number of registered properties
    in the class is returned in NPROPS.  If a property list ID is given, the
    current number of properties in the list is returned in NPROPS.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Pget_nprops(hid_t id, size_t *nprops)
{
    H5P_genplist_t	*plist;    /* Property list to query */
    H5P_genclass_t	*pclass;   /* Property class to query */
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Pget_nprops, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(id) && H5I_GENPROP_CLS != H5I_get_type(id))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property object");
    if (nprops==NULL)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid property nprops pointer");

    if(H5I_GENPROP_LST == H5I_get_type(id)) {
        if (NULL == (plist = H5I_object(id)))
            HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");
        if ((ret_value=H5P_get_nprops_plist(plist,nprops))<0)
            HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to query # of properties in plist");
    } /* end if */
    else 
        if(H5I_GENPROP_CLS == H5I_get_type(id)) {
            if (NULL == (pclass = H5I_object(id)))
                HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property class");
            if ((ret_value=H5P_get_nprops_pclass(pclass,nprops))<0)
                HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to query # of properties in pclass");
        } /* end if */
        else
            HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property object");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pget_nprops() */


/*--------------------------------------------------------------------------
 NAME
    H5P_cmp_prop
 PURPOSE
    Internal routine to compare two generic properties
 USAGE
    int H5P_cmp_prop(prop1, prop2)
        H5P_genprop_t *prop1;    IN: 1st property to compare
        H5P_genprop_t *prop1;    IN: 2nd property to compare
 RETURNS
    Success: negative if prop1 "less" than prop2, positive if prop1 "greater"
        than prop2, zero if prop1 is "equal" to prop2
    Failure: can't fail
 DESCRIPTION
        This function compares two generic properties together to see if
    they are the same property.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static int H5P_cmp_prop(H5P_genprop_t *prop1, H5P_genprop_t *prop2)
{
    int cmp_value;             /* Value from comparison */
    int ret_value=0;         /* return value */

    FUNC_ENTER (H5P_cmp_prop, FAIL);

    assert(prop1);
    assert(prop2);

    /* Check the name */
    if((cmp_value=HDstrcmp(prop1->name,prop2->name))!=0)
        HGOTO_DONE(cmp_value);

    /* Check the size of properties */
    if(prop1->size < prop2->size) HGOTO_DONE(-1);
    if(prop1->size > prop2->size) HGOTO_DONE(1);

    /* Check if they both have values allocated (or not allocated) */
    if(prop1->value==NULL && prop2->value!=NULL) HGOTO_DONE(-1);
    if(prop1->value!=NULL && prop2->value==NULL) HGOTO_DONE(1);
    if(prop1->value!=NULL)
        if((cmp_value=HDmemcmp(prop1->value,prop2->value,prop1->size))!=0)
            HGOTO_DONE(cmp_value);

    /* Check if they both have default values allocated (or not allocated) */
    if(prop1->def_value==NULL && prop2->def_value!=NULL) HGOTO_DONE(-1);
    if(prop1->def_value!=NULL && prop2->def_value==NULL) HGOTO_DONE(1);
    if(prop1->def_value!=NULL)
        if((cmp_value=HDmemcmp(prop1->def_value,prop2->def_value,prop1->size))!=0)
            HGOTO_DONE(cmp_value);

    /* Check if they both have the same 'create' callback */
    if(prop1->create==NULL && prop2->create!=NULL) HGOTO_DONE(-1);
    if(prop1->create!=NULL && prop2->create==NULL) HGOTO_DONE(1);
    if(prop1->create!=prop2->create) HGOTO_DONE(-1);

    /* Check if they both have the same 'set' callback */
    if(prop1->set==NULL && prop2->set!=NULL) HGOTO_DONE(-1);
    if(prop1->set!=NULL && prop2->set==NULL) HGOTO_DONE(1);
    if(prop1->set!=prop2->set) HGOTO_DONE(-1);

    /* Check if they both have the same 'get' callback */
    if(prop1->get==NULL && prop2->get!=NULL) HGOTO_DONE(-1);
    if(prop1->get!=NULL && prop2->get==NULL) HGOTO_DONE(1);
    if(prop1->get!=prop2->get) HGOTO_DONE(-1);

    /* Check if they both have the same 'close' callback */
    if(prop1->close==NULL && prop2->close!=NULL) HGOTO_DONE(-1);
    if(prop1->close!=NULL && prop2->close==NULL) HGOTO_DONE(1);
    if(prop1->close!=prop2->close) HGOTO_DONE(-1);

    /* Don't check the 'next' field, they must be equal by now */

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_cmp_prop() */


/*--------------------------------------------------------------------------
 NAME
    H5P_cmp_class
 PURPOSE
    Internal routine to compare two generic property classes 
 USAGE
    int H5P_cmp_class(pclass1, pclass2)
        H5P_genclass_t *pclass1;    IN: 1st property class to compare
        H5P_genclass_t *pclass2;    IN: 2nd property class to compare
 RETURNS
    Success: negative if class1 "less" than class2, positive if class1 "greater"
        than class2, zero if class1 is "equal" to class2
    Failure: can't fail
 DESCRIPTION
        This function compares two generic property classes together to see if
    they are the same class.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static int H5P_cmp_class(H5P_genclass_t *pclass1, H5P_genclass_t *pclass2)
{
    H5P_genprop_t *tprop1, *tprop2;/* Temporary pointer to properties */
    unsigned u;                    /* Local index variable */
    int cmp_value;             /* Value from comparison */
    int ret_value=0;         /* return value */

    FUNC_ENTER (H5P_cmp_class, FAIL);

    assert(pclass1);
    assert(pclass2);

    /* Check the name */
    if((cmp_value=HDstrcmp(pclass1->name,pclass2->name))!=0)
        HGOTO_DONE(cmp_value);

    /* Check the number of properties */
    if(pclass1->nprops < pclass2->nprops) HGOTO_DONE(-1);
    if(pclass1->nprops > pclass2->nprops) HGOTO_DONE(1);

    /* Check the hashsize */
    if(pclass1->hashsize < pclass2->hashsize) HGOTO_DONE(-1);
    if(pclass1->hashsize > pclass2->hashsize) HGOTO_DONE(1);

    /* Check the number of property lists created from the class */
    if(pclass1->plists < pclass2->plists) HGOTO_DONE(-1);
    if(pclass1->plists > pclass2->plists) HGOTO_DONE(1);

    /* Check the number of classes derived from the class */
    if(pclass1->classes < pclass2->classes) HGOTO_DONE(-1);
    if(pclass1->classes > pclass2->classes) HGOTO_DONE(1);

    /* Check the number of ID references open on the class */
    if(pclass1->ref_count < pclass2->ref_count) HGOTO_DONE(-1);
    if(pclass1->ref_count > pclass2->ref_count) HGOTO_DONE(1);

    /* Check whether they are internal or not */
    if(pclass1->internal < pclass2->internal) HGOTO_DONE(-1);
    if(pclass1->internal > pclass2->internal) HGOTO_DONE(1);

    /* Check whether they are deleted or not */
    if(pclass1->deleted < pclass2->deleted) HGOTO_DONE(-1);
    if(pclass1->deleted > pclass2->deleted) HGOTO_DONE(1);

    /* Check whether they have creation callback functions & data */
    if(pclass1->create_func==NULL && pclass2->create_func!=NULL) HGOTO_DONE(-1);
    if(pclass1->create_func!=NULL && pclass2->create_func==NULL) HGOTO_DONE(1);
    if(pclass1->create_func!=pclass2->create_func) HGOTO_DONE(-1);
    if(pclass1->create_data < pclass2->create_data) HGOTO_DONE(-1);
    if(pclass1->create_data > pclass2->create_data) HGOTO_DONE(1);

    /* Check whether they have close callback functions & data */
    if(pclass1->close_func==NULL && pclass2->close_func!=NULL) HGOTO_DONE(-1);
    if(pclass1->close_func!=NULL && pclass2->close_func==NULL) HGOTO_DONE(1);
    if(pclass1->close_func!=pclass2->close_func) HGOTO_DONE(-1);
    if(pclass1->close_data < pclass2->close_data) HGOTO_DONE(-1);
    if(pclass1->close_data > pclass2->close_data) HGOTO_DONE(1);

    /* Cycle through the properties and compare them also */
    for(u=0; u<pclass1->hashsize; u++) {
        tprop1=pclass1->props[u];
        tprop2=pclass2->props[u];

        /* Check if they both have properties in this hash location */
        if(tprop1==NULL && tprop2!=NULL) HGOTO_DONE(-1);
        if(tprop1!=NULL && tprop2==NULL) HGOTO_DONE(1);

        /* Check the actual properties */
        while(tprop1!=NULL && tprop2!=NULL) {
            /* Compare the two properties */
            if((cmp_value=H5P_cmp_prop(tprop1,tprop2))!=0)
                HGOTO_DONE(cmp_value);

            /* Advance the pointers */
            tprop1=tprop1->next;
            tprop2=tprop2->next;

            /* Check if they both have properties in this location */
            if(tprop1==NULL && tprop2!=NULL) HGOTO_DONE(-1);
            if(tprop1!=NULL && tprop2==NULL) HGOTO_DONE(1);
        } /* end while */
    } /* end for */

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_cmp_class() */


/*--------------------------------------------------------------------------
 NAME
    H5P_cmp_plist
 PURPOSE
    Internal routine to compare two generic property lists 
 USAGE
    int H5P_cmp_plist(plist1, plist2)
        H5P_genplist_t *plist1;    IN: 1st property list to compare
        H5P_genplist_t *plist2;    IN: 2nd property list to compare
 RETURNS
    Success: negative if list1 "less" than list2, positive if list1 "greater"
        than list2, zero if list1 is "equal" to list2
    Failure: can't fail
 DESCRIPTION
        This function compares two generic property lists together to see if
    they are the same list.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static int H5P_cmp_plist(H5P_genplist_t *plist1, H5P_genplist_t *plist2)
{
    H5P_genprop_t *tprop1, *tprop2;/* Temporary pointer to properties */
    unsigned u;                    /* Local index variable */
    int cmp_value;             /* Value from comparison */
    int ret_value=0;         /* return value */

    FUNC_ENTER (H5P_cmp_plist, FAIL);

    assert(plist1);
    assert(plist2);

    /* Check the parent class */
    if(plist1->pclass < plist2->pclass) HGOTO_DONE(-1);
    if(plist1->pclass > plist2->pclass) HGOTO_DONE(1);

    /* Check the number of properties */
    if(plist1->nprops < plist2->nprops) HGOTO_DONE(-1);
    if(plist1->nprops > plist2->nprops) HGOTO_DONE(1);

    /* Check whether they've been initialized */
    if(plist1->class_init < plist2->class_init) HGOTO_DONE(-1);
    if(plist1->class_init > plist2->class_init) HGOTO_DONE(1);

    /* Cycle through the properties and compare them also */
    for(u=0; u<plist1->pclass->hashsize; u++) {
        tprop1=plist1->props[u];
        tprop2=plist2->props[u];

        /* Check if they both have properties in this hash location */
        if(tprop1==NULL && tprop2!=NULL) HGOTO_DONE(-1);
        if(tprop1!=NULL && tprop2==NULL) HGOTO_DONE(1);

        /* Check the actual properties */
        while(tprop1!=NULL && tprop2!=NULL) {
            /* Compare the two properties */
            if((cmp_value=H5P_cmp_prop(tprop1,tprop2))!=0)
                HGOTO_DONE(cmp_value);

            /* Advance the pointers */
            tprop1=tprop1->next;
            tprop2=tprop2->next;

            /* Check if they both have properties in this location */
            if(tprop1==NULL && tprop2!=NULL) HGOTO_DONE(-1);
            if(tprop1!=NULL && tprop2==NULL) HGOTO_DONE(1);
        } /* end while */
    } /* end for */

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_cmp_plist() */


/*--------------------------------------------------------------------------
 NAME
    H5Pequal
 PURPOSE
    Routine to query whether two property lists or two property classes are equal
 USAGE
    htri_t H5Pequal(id1, id2)
        hid_t id1;         IN: Property list or class ID to compare
        hid_t id2;         IN: Property list or class ID to compare
 RETURNS
    Success: TRUE if equal, FALSE if unequal
    Failure: negative
 DESCRIPTION
    Determines whether two property lists or two property classes are equal.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
htri_t H5Pequal(hid_t id1, hid_t id2)
{
    void *obj1, *obj2;          /* Property objects to compare */
    htri_t ret_value=FALSE;     /* return value */

    FUNC_ENTER (H5Pequal, FAIL);

    /* Check arguments. */
    if ((H5I_GENPROP_LST != H5I_get_type(id1) && H5I_GENPROP_CLS != H5I_get_type(id1))
            || (H5I_GENPROP_CLS != H5I_get_type(id2) && H5I_GENPROP_CLS != H5I_get_type(id2)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not property objects");
    if (H5I_get_type(id1) != H5I_get_type(id2))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not the same kind of property objects");
    if(NULL == (obj1 = H5I_object(id1)) || NULL == (obj2 = H5I_object(id2)))
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "property object doesn't exist");

    /* Compare property lists */
    if(H5I_GENPROP_LST == H5I_get_type(id1)) {
        if(H5P_cmp_plist(obj1,obj2)==0)
            ret_value=TRUE;
    } /* end if */
    /* Must be property classes */
    else {
        if(H5P_cmp_class(obj1,obj2)==0)
            ret_value=TRUE;
    } /* end else */

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pequal() */


/*--------------------------------------------------------------------------
 NAME
    H5P_iterate_props
 PURPOSE
    Internal routine to iterate over a hashtable of properties 
 USAGE
    herr_t H5P_iterate_props(id, hash, hashsize, idx, iter_func, iter_data)
        hid_t id;                   IN: ID of property object iterating over
        H5P_gen_prop_t *hash[];     IN: Pointer to array of properties for hash table
        unsigned hashsize;             IN: Size of hash table
        int *idx;                   IN/OUT: Index of the property to begin with
        H5P_iterate_t iter_func;    IN: Function pointer to function to be
                                        called with each property iterated over.
        void *iter_data;            IN/OUT: Pointer to iteration data from user
 RETURNS
    Success: Returns the return value of the last call to ITER_FUNC if it was
                non-zero, or zero if all properties have been processed.
    Failure: negative value
 DESCRIPTION
    This routine iterates over the properties in the property object specified
with ID.  For each property in the object, the ITER_DATA and some
additional information, specified below, are passed to the ITER_FUNC function.
The iteration begins with the IDX property in the object and the next element
to be processed by the operator is returned in IDX.  If IDX is NULL, then the
iterator starts at the first property; since no stopping point is returned in
this case, the iterator cannot be restarted if one of the calls to its operator
returns non-zero. 

The prototype for H5P_iterate_t is: 
    typedef herr_t (*H5P_iterate_t)(hid_t id, const char *name, void *iter_data); 
The operation receives the property list or class identifier for the object
being iterated over, ID, the name of the current property within the object,
NAME, and the pointer to the operator data passed in to H5Piterate, ITER_DATA. 

The return values from an operator are: 
    Zero causes the iterator to continue, returning zero when all properties
        have been processed. 
    Positive causes the iterator to immediately return that positive value,
        indicating short-circuit success. The iterator can be restarted at the
        index of the next property. 
    Negative causes the iterator to immediately return that value, indicating
        failure. The iterator can be restarted at the index of the next
        property.

H5Piterate assumes that the properties in the object identified by ID remains
unchanged through the iteration.  If the membership changes during the
iteration, the function's behavior is undefined. 

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static int H5P_iterate_props(hid_t id, H5P_genprop_t *hash[], unsigned hashsize, int *idx, H5P_iterate_t iter_func, void *iter_data)
{
    H5P_genprop_t *prop;        /* Temporary property pointer */
    unsigned u;                    /* Local index variable */
    int curr_idx=0;             /* Current iteration index */
    int ret_value=0;            /* Return value */

    FUNC_ENTER (H5P_iterate_props, FAIL);

    assert(hash);
    assert(hashsize>0);
    assert(idx);
    assert(iter_func);

    /* Cycle through the properties and compare them also */
    for(u=0; u<hashsize && ret_value==0; u++) {
        prop=hash[u];

        /* Check the actual properties */
        while(prop!=NULL && ret_value==0) {
            /* Check if we are at the object to start iterating over */
            if(curr_idx>=*idx)
                ret_value=(*iter_func)(id,prop->name,iter_data);

            /* Increment the iteration index if iteration function succeeded */
            if(ret_value==0)
                curr_idx++;

            /* Advance the pointer */
            prop=prop->next;
        } /* end while */
    } /* end for */

    /* Set the index we stopped at */
    *idx=curr_idx;

    FUNC_LEAVE (ret_value);
}   /* H5P_iterate_props() */


/*--------------------------------------------------------------------------
 NAME
    H5P_iterate_plist
 PURPOSE
    Internal routine to iterate over the properties in a property list
 USAGE
    herr_t H5P_iterate_plist(plist_id, idx, iter_func, iter_data)
        hid_t plist_id;             IN: ID of property list to iterate over
        int *idx;                   IN/OUT: Index of the property to begin with
        H5P_iterate_t iter_func;    IN: Function pointer to function to be
                                        called with each property iterated over.
        void *iter_data;            IN/OUT: Pointer to iteration data from user
 RETURNS
    Success: Returns the return value of the last call to ITER_FUNC if it was
                non-zero, or zero if all properties have been processed.
    Failure: negative value
 DESCRIPTION
    This routine iterates over the properties in the property object specified
with PLIST_ID.  For each property in the object, the ITER_DATA and some
additional information, specified below, are passed to the ITER_FUNC function.
The iteration begins with the IDX property in the object and the next element
to be processed by the operator is returned in IDX.  If IDX is NULL, then the
iterator starts at the first property; since no stopping point is returned in
this case, the iterator cannot be restarted if one of the calls to its operator
returns non-zero. 

The prototype for H5P_iterate_t is: 
    typedef herr_t (*H5P_iterate_t)(hid_t id, const char *name, void *iter_data); 
The operation receives the property list or class identifier for the object
being iterated over, ID, the name of the current property within the object,
NAME, and the pointer to the operator data passed in to H5Piterate, ITER_DATA. 

The return values from an operator are: 
    Zero causes the iterator to continue, returning zero when all properties
        have been processed. 
    Positive causes the iterator to immediately return that positive value,
        indicating short-circuit success. The iterator can be restarted at the
        index of the next property. 
    Negative causes the iterator to immediately return that value, indicating
        failure. The iterator can be restarted at the index of the next
        property.

H5Piterate assumes that the properties in the object identified by ID remains
unchanged through the iteration.  If the membership changes during the
iteration, the function's behavior is undefined. 

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static int H5P_iterate_plist(hid_t plist_id, int *idx, H5P_iterate_t iter_func, void *iter_data)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    int ret_value=0;            /* Return value */

    FUNC_ENTER (H5P_iterate_plist, FAIL);

    assert(idx);
    assert(iter_func);

    /* Get the property list object */
    if (H5I_GENPROP_LST != H5I_get_type(plist_id) || NULL == (plist = H5I_object(plist_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");

    /* Iterate through the properties in the property list */
    ret_value=H5P_iterate_props(plist_id, plist->props, plist->pclass->hashsize, idx, iter_func, iter_data);

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_iterate_plist() */


/*--------------------------------------------------------------------------
 NAME
    H5P_iterate_pclass
 PURPOSE
    Internal routine to iterate over the properties in a property class
 USAGE
    herr_t H5P_iterate_pclass(pclass_id, idx, iter_func, iter_data)
        hid_t pclass_id;            IN: ID of property class to iterate over
        int *idx;                   IN/OUT: Index of the property to begin with
        H5P_iterate_t iter_func;    IN: Function pointer to function to be
                                        called with each property iterated over.
        void *iter_data;            IN/OUT: Pointer to iteration data from user
 RETURNS
    Success: Returns the return value of the last call to ITER_FUNC if it was
                non-zero, or zero if all properties have been processed.
    Failure: negative value
 DESCRIPTION
    This routine iterates over the properties in the property object specified
with PCLASS_ID.  For each property in the object, the ITER_DATA and some
additional information, specified below, are passed to the ITER_FUNC function.
The iteration begins with the IDX property in the object and the next element
to be processed by the operator is returned in IDX.  If IDX is NULL, then the
iterator starts at the first property; since no stopping point is returned in
this case, the iterator cannot be restarted if one of the calls to its operator
returns non-zero. 

The prototype for H5P_iterate_t is: 
    typedef herr_t (*H5P_iterate_t)(hid_t id, const char *name, void *iter_data); 
The operation receives the property list or class identifier for the object
being iterated over, ID, the name of the current property within the object,
NAME, and the pointer to the operator data passed in to H5Piterate, ITER_DATA. 

The return values from an operator are: 
    Zero causes the iterator to continue, returning zero when all properties
        have been processed. 
    Positive causes the iterator to immediately return that positive value,
        indicating short-circuit success. The iterator can be restarted at the
        index of the next property. 
    Negative causes the iterator to immediately return that value, indicating
        failure. The iterator can be restarted at the index of the next
        property.

H5Piterate assumes that the properties in the object identified by ID remains
unchanged through the iteration.  If the membership changes during the
iteration, the function's behavior is undefined. 

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static int H5P_iterate_pclass(hid_t pclass_id, int *idx, H5P_iterate_t iter_func, void *iter_data)
{
    H5P_genclass_t *pclass;     /* Property list pointer */
    int ret_value=0;            /* Return value */

    FUNC_ENTER (H5P_iterate_pclass, FAIL);

    assert(idx);
    assert(iter_func);

    /* Get the property list object */
    if (H5I_GENPROP_CLS != H5I_get_type(pclass_id) || NULL == (pclass = H5I_object(pclass_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property class");

    /* Iterate through the properties in the property list */
    ret_value=H5P_iterate_props(pclass_id, pclass->props, pclass->hashsize, idx, iter_func, iter_data);

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_iterate_pclass() */


/*--------------------------------------------------------------------------
 NAME
    H5Piterate
 PURPOSE
    Routine to iterate over the properties in a property list or class
 USAGE
    int H5Piterate(pclass_id, idx, iter_func, iter_data)
        hid_t id;                   IN: ID of property object to iterate over
        int *idx;                   IN/OUT: Index of the property to begin with
        H5P_iterate_t iter_func;    IN: Function pointer to function to be
                                        called with each property iterated over.
        void *iter_data;            IN/OUT: Pointer to iteration data from user
 RETURNS
    Success: Returns the return value of the last call to ITER_FUNC if it was
                non-zero, or zero if all properties have been processed.
    Failure: negative value
 DESCRIPTION
    This routine iterates over the properties in the property object specified
with ID.  The properties in both property lists and classes may be iterated
over with this function.  For each property in the object, the ITER_DATA and
some additional information, specified below, are passed to the ITER_FUNC
function.  The iteration begins with the IDX property in the object and the
next element to be processed by the operator is returned in IDX.  If IDX is
NULL, then the iterator starts at the first property; since no stopping point
is returned in this case, the iterator cannot be restarted if one of the calls
to its operator returns non-zero. 

The prototype for H5P_iterate_t is: 
    typedef herr_t (*H5P_iterate_t)(hid_t id, const char *name, void *iter_data); 
The operation receives the property list or class identifier for the object
being iterated over, ID, the name of the current property within the object,
NAME, and the pointer to the operator data passed in to H5Piterate, ITER_DATA. 

The return values from an operator are: 
    Zero causes the iterator to continue, returning zero when all properties
        have been processed. 
    Positive causes the iterator to immediately return that positive value,
        indicating short-circuit success. The iterator can be restarted at the
        index of the next property. 
    Negative causes the iterator to immediately return that value, indicating
        failure. The iterator can be restarted at the index of the next
        property.

H5Piterate assumes that the properties in the object identified by ID remains
unchanged through the iteration.  If the membership changes during the
iteration, the function's behavior is undefined. 

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
int H5Piterate(hid_t id, int *idx, H5P_iterate_t iter_func, void *iter_data)
{
    int fake_idx=0;         /* Index when user doesn't provide one */
    int ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Piterate, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(id) && H5I_GENPROP_CLS != H5I_get_type(id))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property object");
    if (iter_func==NULL)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid iteration callback");

    if (H5I_GENPROP_LST == H5I_get_type(id)) {
        /* Iterate over a property list */
        if ((ret_value=H5P_iterate_plist(id,(idx ? idx : &fake_idx),iter_func,iter_data))<0)
            HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to iterate over list");
    } /* end if */
    else 
        if (H5I_GENPROP_CLS == H5I_get_type(id)) {
            /* Iterate over a property class */
            if ((ret_value=H5P_iterate_pclass(id,(idx ? idx : &fake_idx),iter_func,iter_data))<0)
                HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to iterate over class");
        } /* end if */
        else 
            HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property object");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Piterate() */


/*--------------------------------------------------------------------------
 NAME
    H5P_get
 PURPOSE
    Internal routine to query the value of a property in a property list.
 USAGE
    herr_t H5P_get(plist, name, value)
        H5P_genplist_t *plist;  IN: Property list to check
        const char *name;       IN: Name of property to query
        void *value;            OUT: Pointer to the buffer for the property value
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Retrieves a copy of the value for a property in a property list.  The
    property name must exist or this routine will fail.  If there is a
    'get' callback routine registered for this property, the copy of the
    value of the property will first be passed to that routine and any changes
    to the copy of the value will be used when returning the property value
    from this routine.
        If the 'get' callback routine returns an error, 'value' will not be
    modified and this routine will return an error.  This routine may not be
    called for zero-sized properties.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_get(hid_t plist_id, const char *name, void *value)
{
    H5P_genplist_t *plist;      /* Property list pointer */
    H5P_genprop_t *prop;        /* Temporary property pointer */
    herr_t ret_value=FAIL;      /* return value */

    FUNC_ENTER (H5P_get, FAIL);

    assert(name);
    assert(value);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(plist_id) || NULL == (plist = H5I_object(plist_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");

    /* Find property */
    if((prop=H5P_find_prop(plist->props,plist->pclass->hashsize,name))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "property doesn't exist");

    /* Check for property size >0 */
    if(prop->size==0)
        HGOTO_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "property has zero size");

    /* Make a copy of the value and pass to 'get' callback */
    if(prop->get!=NULL) {
        void *tmp_value;            /* Temporary value for property */

        /* Make a copy of the current value, in case the callback fails */
        if (NULL==(tmp_value=H5MM_malloc(prop->size)))
            HRETURN_ERROR(H5E_RESOURCE, H5E_NOSPACE, FAIL, "memory allocation failed temporary property value");
        HDmemcpy(tmp_value,prop->value,prop->size);

        /* Call user's callback */
        if((*(prop->get))(plist_id,name,tmp_value)<0)
            HRETURN_ERROR(H5E_PLIST, H5E_CANTINIT, FAIL, "can't get property value");

        /* Copy new [possibly unchanged] value into return value */
        HDmemcpy(value,tmp_value,prop->size);

        /* Free the temporary value buffer */
        H5MM_xfree(tmp_value);
    } /* end if */
    /* No 'get' callback, just copy value */
    else
        HDmemcpy(value,prop->value,prop->size);

    /* Set return value */
    ret_value=SUCCEED;

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_get() */


/*--------------------------------------------------------------------------
 NAME
    H5Pget
 PURPOSE
    Routine to query the value of a property in a property list.
 USAGE
    herr_t H5Pget(plist_id, name, value)
        hid_t plist_id;         IN: Property list to check
        const char *name;       IN: Name of property to query
        void *value;            OUT: Pointer to the buffer for the property value
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Retrieves a copy of the value for a property in a property list.  The
    property name must exist or this routine will fail.  If there is a
    'get' callback routine registered for this property, the copy of the
    value of the property will first be passed to that routine and any changes
    to the copy of the value will be used when returning the property value
    from this routine.
        If the 'get' callback routine returns an error, 'value' will not be
    modified and this routine will return an error.  This routine may not be
    called for zero-sized properties.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Pget(hid_t plist_id, const char *name, void * value)
{
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Pget, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(plist_id))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid property name");
    if (value==NULL)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalied property value");

    /* Create the new property list class */
    if ((ret_value=H5P_get(plist_id,name,value))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to query property value");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pget() */


/*--------------------------------------------------------------------------
 NAME
    H5P_remove
 PURPOSE
    Internal routine to remove a property from a property list.
 USAGE
    herr_t H5P_remove(plist, name)
        H5P_genplist_t *plist;  IN: Property list to modify
        const char *name;       IN: Name of property to remove
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Removes a property from a property list.  Both properties which were
    in existance when the property list was created (i.e. properties registered
    with H5Pregister) and properties added to the list after it was created
    (i.e. added with H5Pinsert) may be removed from a property list.
    Properties do not need to be removed a property list before the list itself
    is closed, they will be released automatically when H5Pclose is called.
    The 'close' callback for this property is called before the property is
    release, if the callback exists.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_remove(H5P_genplist_t *plist, const char *name)
{
    H5P_genprop_t *prop;        /* Temporary property pointer */
    H5P_genprop_t *tprop, *prev;/* Temporary pointer to properties */
    unsigned loc;                  /* Hash table location */
    herr_t ret_value=FAIL;      /* return value */

    FUNC_ENTER (H5P_remove, FAIL);

    assert(plist);
    assert(name);

    /* Find the property in the property list */
    if((prop=H5P_find_prop(plist->props,plist->pclass->hashsize,name))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "property doesn't exist");

    /* Pass value to 'close' callback, if it exists */
    if(prop->close!=NULL) {
        /* Call user's callback */
        if((*(prop->close))(name,prop->value)<0)
            HRETURN_ERROR(H5E_PLIST, H5E_CANTINIT, FAIL, "can't close property value");
    } /* end if */

    /* Get correct hash table location */
    loc=H5P_hash_name(name,plist->pclass->hashsize);
    
    /* Remove property from property list */
    /* Check if the property being removed is at the head of the list for a hash location */
    if(prop==plist->props[loc])
        plist->props[loc]=prop->next;
    else {
        /* Set up initial pointers */
        prev=tprop=plist->props[loc];
        tprop=tprop->next;
        while(tprop!=NULL) {
            if(tprop==prop) {
                /* Jump over the property we are deleting */
                prev->next=prop->next;

                /* Free the property, ignoring return value, nothing we can do */
                H5P_free_prop(prop);

                /* Break out of while loop */
                break;
            } /* end if */

            /* Move to the next nodes */
            prev=tprop;
            tprop=tprop->next;
        } /* end while */
    } /* end else */

    /* Decrement the number of properties in list */
    plist->nprops--;

    /* Set return value */
    ret_value=SUCCEED;

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_remove() */


/*--------------------------------------------------------------------------
 NAME
    H5Premove
 PURPOSE
    Routine to remove a property from a property list.
 USAGE
    herr_t H5Premove(plist_id, name)
        hid_t plist_id;         IN: Property list to modify
        const char *name;       IN: Name of property to remove
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Removes a property from a property list.  Both properties which were
    in existance when the property list was created (i.e. properties registered
    with H5Pregister) and properties added to the list after it was created
    (i.e. added with H5Pinsert) may be removed from a property list.
    Properties do not need to be removed a property list before the list itself
    is closed, they will be released automatically when H5Pclose is called.
    The 'close' callback for this property is called before the property is
    release, if the callback exists.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Premove(hid_t plist_id, const char *name)
{
    H5P_genplist_t	*plist;    /* Property list to modify */
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Premove, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(plist_id) || NULL == (plist = H5I_object(plist_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid property name");

    /* Create the new property list class */
    if ((ret_value=H5P_remove(plist,name))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to remove property");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Premove() */


/*--------------------------------------------------------------------------
 NAME
    H5P_unregister
 PURPOSE
    Internal routine to remove a property from a property list class.
 USAGE
    herr_t H5P_unregister(pclass, name)
        H5P_genclass_t *pclass; IN: Property list class to modify
        const char *name;       IN: Name of property to remove
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Removes a property from a property list class.  Future property lists
    created of that class will not contain this property.  Existing property
    lists containing this property are not affected.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_unregister(H5P_genclass_t *pclass, const char *name)
{
    H5P_genprop_t *prop;        /* Temporary property pointer */
    H5P_genprop_t *tprop, *prev;/* Temporary pointer to properties */
    unsigned loc;                  /* Hash table location */
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5P_unregister, FAIL);

    assert(pclass);
    assert(name);

    /* Find the property in the property list */
    if((prop=H5P_find_prop(pclass->props,pclass->hashsize,name))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "property doesn't exist");

    /* Get correct hash table location */
    loc=H5P_hash_name(name,pclass->hashsize);
    
    /* Remove property from property list */
    /* Check if the property being removed is at the head of the list for a hash location */
    if(prop==pclass->props[loc]) {
        pclass->props[loc]=prop->next;

        /* Free the property, ignoring return value, nothing we can do */
        H5P_free_prop(prop);
    } /* end if */
    else {
        /* Set up initial pointers */
        prev=tprop=pclass->props[loc];
        tprop=tprop->next;
        while(tprop!=NULL) {
            if(tprop==prop) {
                /* Jump over the property we are deleting */
                prev->next=prop->next;

                /* Free the property, ignoring return value, nothing we can do */
                H5P_free_prop(prop);

                /* Break out of while loop */
                break;
            } /* end if */

            /* Move to the next nodes */
            prev=tprop;
            tprop=tprop->next;
        } /* end while */
    } /* end else */

    /* Decrement the number of registered properties in class */
    pclass->nprops--;

    /* Set return value */
    ret_value=SUCCEED;

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_unregister() */


/*--------------------------------------------------------------------------
 NAME
    H5Punregister
 PURPOSE
    Routine to remove a property from a property list class.
 USAGE
    herr_t H5Punregister(pclass_id, name)
        hid_t pclass_id;         IN: Property list class to modify
        const char *name;       IN: Name of property to remove
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Removes a property from a property list class.  Future property lists
    created of that class will not contain this property.  Existing property
    lists containing this property are not affected.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Punregister(hid_t pclass_id, const char *name)
{
    H5P_genclass_t	*pclass;   /* Property list class to modify */
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5Punregister, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_CLS != H5I_get_type(pclass_id) || NULL == (pclass = H5I_object(pclass_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list class");
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "invalid property name");

    /* Remove the property list from class */
    if ((ret_value=H5P_unregister(pclass,name))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to remove property from class");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Punregister() */


/*--------------------------------------------------------------------------
 NAME
    H5P_close_list
 PURPOSE
    Internal routine to close a property list.
 USAGE
    herr_t H5P_close_list(plist)
        H5P_genplist_t *plist;  IN: Property list to close
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Closes a property list.  If a 'close' callback exists for the property
    list class, it is called before the property list is destroyed.  If 'close'
    callbacks exist for any individual properties in the property list, they are
    called after the class 'close' callback.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
        The property list class 'close' callback routine is not called from
    here, it must have been check for and called properly prior to this routine
    being called
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_close_list(void *_plist)
{
    H5P_genplist_t *plist=(H5P_genplist_t *)_plist;
    herr_t ret_value=FAIL;     /* return value */

    FUNC_ENTER (H5P_close_list, FAIL);

    assert(plist);

    /* Make calls to any property close callbacks which exist */
    H5P_free_all_prop(plist->props,plist->pclass->hashsize,1);

    /* Decrement class's dependant property list value! */
    if(H5P_access_class(plist->pclass,H5P_MOD_DEC_LST)<0)
        HGOTO_ERROR (H5E_PLIST, H5E_CANTINIT, FAIL, "Can't decrement class ref count");

    /* Destroy property list object */
    H5MM_xfree(plist);

    ret_value=SUCCEED;     /* return value */

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_close_list() */


/*--------------------------------------------------------------------------
 NAME
    H5Pclose_list
 PURPOSE
    Routine to close a property list.
 USAGE
    herr_t H5Pclose_list(plist_id)
        hid_t plist_id;       IN: Property list to close
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
        Closes a property list.  If a 'close' callback exists for the property
    list class, it is called before the property list is destroyed.  If 'close'
    callbacks exist for any individual properties in the property list, they are
    called after the class 'close' callback.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t H5Pclose_list(hid_t plist_id)
{
    H5P_genplist_t	*plist;    /* Property list created */
    herr_t ret_value=FAIL;      /* return value */

    FUNC_ENTER (H5Pclose_list, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_LST != H5I_get_type(plist_id) || NULL == (plist = H5I_remove(plist_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list");

    /* Make call to property list class close callback, if needed */
    if(plist->class_init!=0 && plist->pclass->close_func!=NULL) {
        /* Call user's "close" callback function, ignoring return value */
        (plist->pclass->close_func)(plist_id,plist->pclass->close_data);
    } /* end if */

    /* Close the property list */
    if ((ret_value=H5P_close_list(plist)) < 0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTFREE, FAIL, "can't close");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pclose_list() */


/*--------------------------------------------------------------------------
 NAME
    H5P_get_class_name
 PURPOSE
    Internal routine to query the name of a generic property list class
 USAGE
    char *H5P_get_class_name(pclass)
        H5P_genpclass_t *pclass;    IN: Property list to check
 RETURNS
    Success: Pointer to a malloc'ed string containing the class name
    Failure: NULL
 DESCRIPTION
        This routine retrieves the name of a generic property list class.
    The pointer to the name must be free'd by the user for successful calls.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static char *H5P_get_class_name(H5P_genclass_t *pclass)
{
    char *ret_value=NULL;      /* return value */

    FUNC_ENTER (H5P_get_class_name, NULL);

    assert(pclass);

    /* Get property size */
    ret_value=HDstrdup(pclass->name);

    FUNC_LEAVE (ret_value);
}   /* H5P_get_class_name() */


/*--------------------------------------------------------------------------
 NAME
    H5Pget_class_name
 PURPOSE
    Internal routine to query the name of a generic property list class
 USAGE
    char *H5Pget_class_name(pclass_id)
        hid_t pclass_id;         IN: Property class to query
 RETURNS
    Success: Pointer to a malloc'ed string containing the class name
    Failure: NULL
 DESCRIPTION
        This routine retrieves the name of a generic property list class.
    The pointer to the name must be free'd by the user for successful calls.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
char *H5Pget_class_name(hid_t pclass_id)
{
    H5P_genclass_t	*pclass;    /* Property class to query */
    char *ret_value=NULL;       /* return value */

    FUNC_ENTER (H5Pget_class_name, NULL);

    /* Check arguments. */
    if (H5I_GENPROP_CLS != H5I_get_type(pclass_id) || NULL == (pclass = H5I_object(pclass_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, NULL, "not a property class");

    /* Create the new property list class */
    if ((ret_value=H5P_get_class_name(pclass))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, NULL, "unable to query name of class");

done:
    FUNC_LEAVE (ret_value);
}   /* H5Pget_class_name() */


/*--------------------------------------------------------------------------
 NAME
    H5P_get_class_parent
 PURPOSE
    Internal routine to query the parent class of a generic property class
 USAGE
    H5P_genclass_t *H5P_get_class_parent(pclass)
        H5P_genpclass_t *pclass;    IN: Property class to check
 RETURNS
    Success: Pointer to the parent class of a property class
    Failure: NULL
 DESCRIPTION
    This routine retrieves a pointer to the parent class for a property class.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static H5P_genclass_t *H5P_get_class_parent(H5P_genclass_t *pclass)
{
    H5P_genclass_t *ret_value=NULL;      /* return value */

    FUNC_ENTER (H5P_get_class_parent, NULL);

    assert(pclass);

    /* Get property size */
    ret_value=pclass->parent;

    FUNC_LEAVE (ret_value);
}   /* H5P_get_class_parent() */


/*--------------------------------------------------------------------------
 NAME
    H5Pget_class_parent
 PURPOSE
    routine to query the parent class of a generic property class
 USAGE
    hid_t H5Pget_class_parent(pclass_id)
        hid_t pclass_id;         IN: Property class to query
 RETURNS
    Success: ID of parent class object
    Failure: NULL
 DESCRIPTION
    This routine retrieves an ID for the parent class of a property class.

 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
hid_t H5Pget_class_parent(hid_t pclass_id)
{
    H5P_genclass_t	*pclass;    /* Property class to query */
    H5P_genclass_t	*parent=NULL;   /* Parent's property class */
    hid_t ret_value=FAIL;       /* return value */

    FUNC_ENTER (H5Pget_class_parent, FAIL);

    /* Check arguments. */
    if (H5I_GENPROP_CLS != H5I_get_type(pclass_id) || NULL == (pclass = H5I_object(pclass_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property class");

    /* Retrieve the property class's parent */
    if ((parent=H5P_get_class_parent(pclass))==NULL)
        HGOTO_ERROR(H5E_PLIST, H5E_NOTFOUND, FAIL, "unable to query class of property list");

    /* Increment the outstanding references to the class object */
    if(H5P_access_class(parent,H5P_MOD_INC_REF)<0)
        HGOTO_ERROR (H5E_PLIST, H5E_CANTINIT, FAIL,"Can't increment class ID ref count");

    /* Get an atom for the class */
    if ((ret_value = H5I_register(H5I_GENPROP_CLS, parent))<0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTREGISTER, FAIL, "unable to atomize property list class");

done:
    if (ret_value<0 && parent)
        H5P_close_class(parent);

    FUNC_LEAVE (ret_value);
}   /* H5Pget_class_parent() */


/*--------------------------------------------------------------------------
 NAME
    H5P_close_class
 PURPOSE
    Internal routine to close a property list class.
 USAGE
    herr_t H5P_close_class(class)
        H5P_genclass_t *class;  IN: Property list class to close
 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
    Releases memory and de-attach a class from the property list class hierarchy.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
static herr_t H5P_close_class(void *_pclass)
{
    H5P_genclass_t *pclass=(H5P_genclass_t *)_pclass;
    herr_t ret_value = FAIL;     /* return value */

    FUNC_ENTER (H5P_close_class, FAIL);

    assert(pclass);

    /* Decrement the reference count & check if the object should go away */
    if(H5P_access_class(pclass,H5P_MOD_DEC_REF)<0)
        HGOTO_ERROR (H5E_PLIST, H5E_NOTFOUND, FAIL, "Can't decrement ID ref count");

    /* Check if the object should go away */
    if(pclass->ref_count==0) {
        /* Decrement parent class's dependant property class value! */
        if(pclass->parent)
            if (H5P_access_class(pclass->parent, H5P_MOD_DEC_CLS) < 0)
                HGOTO_ERROR (H5E_PLIST, H5E_NOTFOUND, FAIL,"Can't decrement class ref count");
        
        /* Mark class as deleted */
        pclass->deleted = 1;

        /* Check dependancies on this class, deleting it if allowed */
        if (H5P_access_class(pclass, H5P_MOD_CHECK) < 0)
            HGOTO_ERROR (H5E_PLIST, H5E_NOTFOUND, FAIL,"Can't check class ref count");
    } /* end if */

    /* Set return value */
    ret_value = SUCCEED;

done:
    FUNC_LEAVE (ret_value);
}   /* H5P_close_class() */


/*--------------------------------------------------------------------------
 NAME
    H5Pclose_class
 PURPOSE
    Close a property list class.
 USAGE
    herr_t H5Pclose_class(cls_id)
        hid_t cls_id;       IN: Property list class ID to class

 RETURNS
    Returns non-negative on success, negative on failure.
 DESCRIPTION
    Releases memory and de-attach a class from the property list class hierarchy.
 GLOBAL VARIABLES
 COMMENTS, BUGS, ASSUMPTIONS
 EXAMPLES
 REVISION LOG
--------------------------------------------------------------------------*/
herr_t
H5Pclose_class(hid_t cls_id)
{
    H5P_genclass_t	*pclass;        /* Property list class created */
    hid_t	ret_value = SUCCEED;    /* Return value			*/

    FUNC_ENTER(H5Pclose_class, FAIL);
    H5TRACE1("e","i",cls_id);

    /* Check arguments */
    if (H5I_GENPROP_CLS != H5I_get_type(cls_id) || NULL == (pclass = H5I_remove(cls_id)))
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not a property list class");

    /* Delete the property list class */
    if (H5P_close_class(pclass) < 0)
        HGOTO_ERROR(H5E_PLIST, H5E_CANTFREE, FAIL, "can't close");

done:
    FUNC_LEAVE(ret_value);
}   /* H5Pclose_class() */

