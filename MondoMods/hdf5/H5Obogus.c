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
 * Created:             H5Obogus.c
 *                      Jan 21 2003
 *                      Quincey Koziol <koziol@ncsa.uiuc.edu>
 *
 * Purpose:             "bogus" message.  This message is guaranteed to never
 *                      be found in a valid HDF5 file and is only used to
 *                      generate a test file which verifies the library's
 *                      correct operation when parsing unknown object header
 *                      messages.
 *
 * Modifications:       
 *
 *-------------------------------------------------------------------------
 */
#include "H5private.h"
#include "H5Eprivate.h"
#include "H5MMprivate.h"
#include "H5Oprivate.h"

#ifdef H5O_ENABLE_BOGUS
#define PABLO_MASK      H5O_bogus_mask

/* PRIVATE PROTOTYPES */
static void *H5O_bogus_decode(H5F_t *f, const uint8_t *p, H5O_shared_t *sh);
static herr_t H5O_bogus_encode(H5F_t *f, uint8_t *p, const void *_mesg);
static size_t H5O_bogus_size(H5F_t *f, const void *_mesg);
static herr_t H5O_bogus_debug(H5F_t *f, const void *_mesg, FILE * stream,
			     int indent, int fwidth);

/* This message derives from H5O */
const H5O_class_t H5O_BOGUS[1] = {{
    H5O_BOGUS_ID,            	/*message id number             */
    "bogus",                 	/*message name for debugging    */
    0,     	                /*native message size           */
    H5O_bogus_decode,        	/*decode message                */
    H5O_bogus_encode,        	/*encode message                */
    NULL,          	        /*copy the native value         */
    H5O_bogus_size,          	/*raw message size              */
    NULL,         	        /*free internal memory          */
    NULL,		        /*free method			*/
    NULL,		    	/*get share method		*/
    NULL,			/*set share method		*/
    H5O_bogus_debug,         	/*debug the message             */
}};

/* Interface initialization */
static int interface_initialize_g = 0;
#define INTERFACE_INIT  NULL


/*-------------------------------------------------------------------------
 * Function:    H5O_bogus_decode
 *
 * Purpose:     Decode a "bogus" message and return a pointer to a new
 *              native message struct.
 *
 * Return:      Success:        Ptr to new message in native struct.
 *
 *              Failure:        NULL
 *
 * Programmer:  Quincey Koziol
 *              koziol@ncsa.uiuc.edu
 *              Jan 21 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void *
H5O_bogus_decode(H5F_t UNUSED *f, const uint8_t *p,
		H5O_shared_t UNUSED *sh)
{
    H5O_bogus_t             *mesg;

    FUNC_ENTER(H5O_bogus_decode, NULL);

    /* check args */
    assert(f);
    assert(p);
    assert(!sh);

    /* Allocate the bogus message */
    if (NULL==(mesg = H5MM_calloc(sizeof(H5O_bogus_t))))
	HRETURN_ERROR (H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");

    /* decode */
    UINT32DECODE(p, mesg->u);

    /* Validate the bogus info */
    if(mesg->u!=H5O_BOGUS_VALUE) {
        H5MM_xfree(mesg);
	HRETURN_ERROR (H5E_OHDR, H5E_BADVALUE, NULL, "invalid bogus value :-)");
    } /* end if */

    FUNC_LEAVE(mesg);
} /* end H5O_bogus_decode() */


/*-------------------------------------------------------------------------
 * Function:    H5O_bogus_encode
 *
 * Purpose:     Encodes a "bogus" message.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Quincey Koziol
 *              koziol@ncsa.uiuc.edu
 *              Jan 21 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5O_bogus_encode(H5F_t UNUSED *f, uint8_t *p, const void UNUSED *mesg)
{
    FUNC_ENTER(H5O_bogus_encode, FAIL);

    /* check args */
    assert(f);
    assert(p);
    assert(mesg);

    /* encode */
    UINT32ENCODE(p, H5O_BOGUS_VALUE);

    FUNC_LEAVE(SUCCEED);
} /* end H5O_bogus_encode() */


/*-------------------------------------------------------------------------
 * Function:    H5O_bogus_size
 *
 * Purpose:     Returns the size of the raw message in bytes not
 *              counting the message typ or size fields, but only the data
 *              fields.  This function doesn't take into account
 *              alignment.
 *
 * Return:      Success:        Message data size in bytes w/o alignment.
 *
 *              Failure:        Negative
 *
 * Programmer:  Quincey Koziol
 *              koziol@ncsa.uiuc.edu
 *              Jan 21 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static size_t
H5O_bogus_size(H5F_t UNUSED *f, const void UNUSED *mesg)
{
    FUNC_ENTER(H5O_bogus_size, 0);

    /* check args */
    assert(f);

    FUNC_LEAVE(4);
} /* end H5O_bogus_size() */


/*-------------------------------------------------------------------------
 * Function:    H5O_bogus_debug
 *
 * Purpose:     Prints debugging info for the message.
 *
 * Return:      Non-negative on success/Negative on failure
 *
 * Programmer:  Quincey Koziol
 *              koziol@ncsa.uiuc.edu
 *              Jan 21 2003
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5O_bogus_debug(H5F_t UNUSED *f, const void *_mesg, FILE *stream,
	       int indent, int fwidth)
{
    const H5O_bogus_t	*mesg = (const H5O_bogus_t *)_mesg;

    FUNC_ENTER(H5O_name_debug, FAIL);

    /* check args */
    assert(f);
    assert(mesg);
    assert(stream);
    assert(indent >= 0);
    assert(fwidth >= 0);

    fprintf(stream, "%*s%-*s `%u'\n", indent, "", fwidth,
            "Bogus Value:", mesg->u);

    FUNC_LEAVE(SUCCEED);
}
#endif /* H5O_ENABLE_BOGUS */

