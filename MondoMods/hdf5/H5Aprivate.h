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
 * This file contains private information about the H5D module
 */
#ifndef _H5Aprivate_H
#define _H5Aprivate_H

#include "H5Apublic.h"
#include "H5Gprivate.h"

#define H5A_RESERVED_ATOMS  0
typedef struct H5A_t H5A_t;

/* Private headers needed by this file */

/* Functions defined in H5A.c */
H5_DLL H5A_t * H5A_copy(const H5A_t *old_attr);
H5_DLL herr_t H5A_close(H5A_t *attr);
H5_DLL H5G_entry_t *H5A_entof(H5A_t *attr);

#endif
