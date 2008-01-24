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
 * Programmer:  Robb Matzke <matzke@llnl.gov>
 *              Monday, July 26, 1999
 */
#ifndef _H5FDprivate_H
#define _H5FDprivate_H

#include "H5FDpublic.h"

#define H5FD_has_cmp(id) (NULL!=H5FD_get_class(id)->cmp)

H5_DLL int H5FD_term_interface(void);
H5_DLL H5FD_class_t *H5FD_get_class(hid_t id);
H5_DLL hsize_t H5FD_sb_size(H5FD_t *file);
H5_DLL herr_t H5FD_sb_encode(H5FD_t *file, char *name/*out*/, uint8_t *buf);
H5_DLL herr_t H5FD_sb_decode(H5FD_t *file, const char *name, const uint8_t *buf);
H5_DLL void *H5FD_fapl_get(H5FD_t *file);
H5_DLL void *H5FD_fapl_copy(hid_t driver_id, const void *fapl);
H5_DLL herr_t H5FD_fapl_free(hid_t driver_id, void *fapl);
H5_DLL void *H5FD_dxpl_copy(hid_t driver_id, const void *dxpl);
H5_DLL herr_t H5FD_dxpl_free(hid_t driver_id, void *dxpl);
H5_DLL H5FD_t *H5FD_open(const char *name, unsigned flags, hid_t fapl_id,
		  haddr_t maxaddr);
H5_DLL herr_t H5FD_close(H5FD_t *file);
H5_DLL int H5FD_cmp(const H5FD_t *f1, const H5FD_t *f2);
H5_DLL int H5FD_query(const H5FD_t *f, unsigned long *flags/*out*/);
H5_DLL haddr_t H5FD_alloc(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id, hsize_t size);
H5_DLL herr_t H5FD_free(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id, haddr_t addr, hsize_t size);
H5_DLL haddr_t H5FD_realloc(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id, haddr_t old_addr,
		     hsize_t old_size, hsize_t new_size);
H5_DLL haddr_t H5FD_get_eoa(H5FD_t *file);
H5_DLL herr_t H5FD_set_eoa(H5FD_t *file, haddr_t addr);
H5_DLL haddr_t H5FD_get_eof(H5FD_t *file);
H5_DLL herr_t H5FD_read(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id, haddr_t addr, hsize_t size,
		 void *buf/*out*/);
H5_DLL herr_t H5FD_write(H5FD_t *file, H5FD_mem_t type, hid_t dxpl_id, haddr_t addr, hsize_t size,
		  const void *buf);
H5_DLL herr_t H5FD_flush(H5FD_t *file, hid_t dxpl_id);
H5_DLL herr_t H5FD_get_fileno(const H5FD_t *file, unsigned long *filenum);

#endif /* !_H5FDprivate_H */
