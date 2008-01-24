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
 *              Thursday, April 16, 1998
 */
#ifndef _H5Zpublic_H
#define _H5Zpublic_H

/*
 * Filter identifiers.  Values 0 through 255 are for filters defined by the
 * HDF5 library.  Values 256 through 511 are available for testing new
 * filters. Subsequent values should be obtained from the HDF5 development
 * team at hdf5dev@ncsa.uiuc.edu.  These values will never change because they
 * appear in the HDF5 files.
 */
typedef int H5Z_filter_t;
#define H5Z_FILTER_ERROR	(-1)	/*no filter			*/
#define H5Z_FILTER_NONE		0	/*reserved indefinitely		*/
#define H5Z_FILTER_DEFLATE	1 	/*deflation like gzip	     	*/
#define H5Z_FILTER_MAX		65535	/*maximum filter id		*/

/* Flags for filter definition */
#define H5Z_FLAG_DEFMASK	0x00ff	/*definition flag mask		*/
#define H5Z_FLAG_OPTIONAL	0x0001	/*filter is optional		*/

/* Additional flags for filter invocation */
#define H5Z_FLAG_INVMASK	0xff00	/*invocation flag mask		*/
#define H5Z_FLAG_REVERSE	0x0100	/*reverse direction; read	*/

#ifdef __cplusplus
extern "C" {
#endif

/*
 * A filter gets definition flags and invocation flags (defined above), the
 * client data array and size defined when the filter was added to the
 * pipeline, the size in bytes of the data on which to operate, and pointers
 * to a buffer and its allocated size.
 *
 * The filter should store the result in the supplied buffer if possible,
 * otherwise it can allocate a new buffer, freeing the original.  The
 * allocated size of the new buffer should be returned through the BUF_SIZE
 * pointer and the new buffer through the BUF pointer.
 *
 * The return value from the filter is the number of bytes in the output
 * buffer. If an error occurs then the function should return zero and leave
 * all pointer arguments unchanged.
 */
typedef size_t (*H5Z_func_t)(unsigned int flags, size_t cd_nelmts,
			     const unsigned int cd_values[], size_t nbytes,
			     size_t *buf_size, void **buf);

H5_DLL herr_t H5Zregister(H5Z_filter_t id, const char *comment,
			   H5Z_func_t filter);

size_t H5Z_filter_deflate(unsigned flags, size_t cd_nelmts,
			  const unsigned cd_values[], size_t nbytes,
			  size_t *buf_size, void **buf);

#ifdef __cplusplus
}
#endif
#endif

