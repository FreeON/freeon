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
 * Programmer: Robb Matzke <matzke@llnl.gov>
 *	       Friday, October 10, 1997
 */

#include "H5private.h"
#include "H5Eprivate.h"
#include "H5Oprivate.h"
#include "H5Vprivate.h"

#define H5V_HYPER_NDIMS H5O_LAYOUT_NDIMS
#define PABLO_MASK	H5V_mask
static int		interface_initialize_g = 0;
#define INTERFACE_INIT	NULL


/*-------------------------------------------------------------------------
 * Function:	H5V_stride_optimize1
 *
 * Purpose:	Given a stride vector which references elements of the
 *		specified size, optimize the dimensionality, the stride
 *		vector, and the element size to minimize the dimensionality
 *		and the number of memory accesses.
 *
 *		All arguments are passed by reference and their values may be
 *		modified by this function.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Saturday, October 11, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5V_stride_optimize1(unsigned *np/*in,out*/, hsize_t *elmt_size/*in,out*/,
		     hsize_t *size, hssize_t *stride1)
{
    FUNC_ENTER(H5V_stride_optimize1, FAIL);

    /*
     * This has to be true because if we optimize the dimensionality down to
     * zero we still must make one reference.
     */
    assert(1 == H5V_vector_reduce_product(0, NULL));

    /*
     * Combine adjacent memory accesses
     */
    while (*np && stride1[*np-1]>0 &&
           (hsize_t)(stride1[*np-1])==*elmt_size) {
        *elmt_size *= size[*np-1];
        if (--*np) {
            stride1[*np-1] += size[*np] * stride1[*np];
        }
    }

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:	H5V_stride_optimize2
 *
 * Purpose:	Given two stride vectors which reference elements of the
 *		specified size, optimize the dimensionality, the stride
 *		vectors, and the element size to minimize the dimensionality
 *		and the number of memory accesses.
 *
 *		All arguments are passed by reference and their values may be
 *		modified by this function.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Saturday, October 11, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5V_stride_optimize2(unsigned *np/*in,out*/, hsize_t *elmt_size/*in,out*/,
		     hsize_t *size, hssize_t *stride1, hssize_t *stride2)
{
    FUNC_ENTER(H5V_stride_optimize2, FAIL);

    /*
     * This has to be true because if we optimize the dimensionality down to
     * zero we still must make one reference.
     */
    assert(1 == H5V_vector_reduce_product(0, NULL));
    assert (*elmt_size>0);

    /*
     * Combine adjacent memory accesses
     */

    /* Unroll loop for common cases */
    switch(*np) {
        case 1: /* For 0-D datasets (dunno if this ever gets used...) */
            if((hsize_t)(stride1[0]) == *elmt_size &&
                   (hsize_t)(stride2[0]) == *elmt_size) {
                *elmt_size *= size[0];
                --*np;  /* *np decrements to a value of 0 now */
            } /* end if */
            break;

        case 2: /* For 1-D datasets */
            if((hsize_t)(stride1[1]) == *elmt_size &&
                   (hsize_t)(stride2[1]) == *elmt_size) {
                *elmt_size *= size[1];
                --*np;  /* *np decrements to a value of 1 now */
                stride1[0] += size[1] * stride1[1];
                stride2[0] += size[1] * stride2[1];

                if((hsize_t)(stride1[0]) == *elmt_size &&
                       (hsize_t)(stride2[0]) == *elmt_size) {
                    *elmt_size *= size[0];
                    --*np;  /* *np decrements to a value of 0 now */
                } /* end if */
            } /* end if */
            break;

        case 3: /* For 2-D datasets */
            if((hsize_t)(stride1[2]) == *elmt_size &&
                   (hsize_t)(stride2[2]) == *elmt_size) {
                *elmt_size *= size[2];
                --*np;  /* *np decrements to a value of 2 now */
                stride1[1] += size[2] * stride1[2];
                stride2[1] += size[2] * stride2[2];

                if((hsize_t)(stride1[1]) == *elmt_size &&
                       (hsize_t)(stride2[1]) == *elmt_size) {
                    *elmt_size *= size[1];
                    --*np;  /* *np decrements to a value of 1 now */
                    stride1[0] += size[1] * stride1[1];
                    stride2[0] += size[1] * stride2[1];

                    if((hsize_t)(stride1[0]) == *elmt_size &&
                           (hsize_t)(stride2[0]) == *elmt_size) {
                        *elmt_size *= size[0];
                        --*np;  /* *np decrements to a value of 0 now */
                    } /* end if */
                } /* end if */
            } /* end if */
            break;

        case 4: /* For 3-D datasets */
            if((hsize_t)(stride1[3]) == *elmt_size &&
                   (hsize_t)(stride2[3]) == *elmt_size) {
                *elmt_size *= size[3];
                --*np;  /* *np decrements to a value of 3 now */
                stride1[2] += size[3] * stride1[3];
                stride2[2] += size[3] * stride2[3];

                if((hsize_t)(stride1[2]) == *elmt_size &&
                       (hsize_t)(stride2[2]) == *elmt_size) {
                    *elmt_size *= size[2];
                    --*np;  /* *np decrements to a value of 2 now */
                    stride1[1] += size[2] * stride1[2];
                    stride2[1] += size[2] * stride2[2];

                    if((hsize_t)(stride1[1]) == *elmt_size &&
                           (hsize_t)(stride2[1]) == *elmt_size) {
                        *elmt_size *= size[1];
                        --*np;  /* *np decrements to a value of 1 now */
                        stride1[0] += size[1] * stride1[1];
                        stride2[0] += size[1] * stride2[1];

                        if((hsize_t)(stride1[0]) == *elmt_size &&
                               (hsize_t)(stride2[0]) == *elmt_size) {
                            *elmt_size *= size[0];
                            --*np;  /* *np decrements to a value of 0 now */
                        } /* end if */
                    } /* end if */
                } /* end if */
            } /* end if */
            break;

        default:
            while (*np &&
                   (hsize_t)(stride1[*np-1]) == *elmt_size &&
                   (hsize_t)(stride2[*np-1]) == *elmt_size) {
                *elmt_size *= size[*np-1];
                if (--*np) {
                    stride1[*np-1] += size[*np] * stride1[*np];
                    stride2[*np-1] += size[*np] * stride2[*np];
                }
            }
            break;
    } /* end switch */

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:	H5V_hyper_stride
 *
 * Purpose:	Given a description of a hyperslab, this function returns
 *		(through STRIDE[]) the byte strides appropriate for accessing
 *		all bytes of the hyperslab and the byte offset where the
 *		striding will begin.  The SIZE can be passed to the various
 *		stride functions.
 *
 *		The dimensionality of the whole array, the hyperslab, and the
 *		returned stride array is N.  The whole array dimensions are
 *		TOTAL_SIZE and the hyperslab is at offset OFFSET and has
 *		dimensions SIZE.
 *
 *		The stride and starting point returned will cause the
 *		hyperslab elements to be referenced in C order.
 *
 * Return:	Success:	Byte offset from beginning of array to start
 *				of striding.
 *
 *		Failure:	abort() -- should never fail
 *
 * Programmer:	Robb Matzke
 *		Saturday, October 11, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
hsize_t
H5V_hyper_stride(unsigned n, const hsize_t *size,
		 const hsize_t *total_size, const hssize_t *offset,
		 hssize_t *stride/*out*/)
{
    hsize_t	    skip;	/*starting point byte offset		*/
    hsize_t	    acc;	/*accumulator				*/
    hsize_t     tmp;
    int		i;		/*counter				*/

    FUNC_ENTER(H5V_hyper_stride, (HDabort(), 0));

    assert(n <= H5V_HYPER_NDIMS);
    assert(size);
    assert(total_size);
    assert(stride);

    /* init */
    stride[n-1] = 1;
    skip = offset ? offset[n-1] : 0;

    switch(n) {
        case 2: /* 1-D dataset */
            tmp = total_size[1] - size[1];
            assert (tmp<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
            stride[0] = (hssize_t)tmp; /*overflow checked*/
            acc = total_size[1];
            skip += acc * (offset ? offset[0] : 0);
            break;

        case 3: /* 2-D dataset */
            tmp = total_size[2] - size[2];
            assert (tmp<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
            stride[1] = (hssize_t)tmp; /*overflow checked*/
            acc = total_size[2];
            skip += acc * (offset ? offset[1] : 0);

            tmp = acc * (total_size[1] - size[1]);
            assert (tmp<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
            stride[0] = (hssize_t)tmp; /*overflow checked*/
            acc *= total_size[1];
            skip += acc * (offset ? offset[0] : 0);
            break;

        case 4: /* 3-D dataset */
            tmp = total_size[3] - size[3];
            assert (tmp<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
            stride[2] = (hssize_t)tmp; /*overflow checked*/
            acc = total_size[3];
            skip += acc * (offset ? offset[2] : 0);

            tmp = acc * (total_size[2] - size[2]);
            assert (tmp<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
            stride[1] = (hssize_t)tmp; /*overflow checked*/
            acc *= total_size[2];
            skip += acc * (offset ? offset[1] : 0);

            tmp = acc * (total_size[1] - size[1]);
            assert (tmp<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
            stride[0] = (hssize_t)tmp; /*overflow checked*/
            acc *= total_size[1];
            skip += acc * (offset ? offset[0] : 0);
            break;

        default:
            /* others */
            for (i=(int)(n-2), acc=1; i>=0; --i) {
                tmp = acc * (total_size[i+1] - size[i+1]);
                assert (tmp<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                stride[i] = (hssize_t)tmp; /*overflow checked*/
                acc *= total_size[i+1];
                skip += acc * (offset ? offset[i] : 0);
            }
            break;
    } /* end switch */

    FUNC_LEAVE(skip);
}


/*-------------------------------------------------------------------------
 * Function:	H5V_hyper_eq
 *
 * Purpose:	Determines whether two hyperslabs are equal.  This function
 *		assumes that both hyperslabs are relative to the same array,
 *		for if not, they could not possibly be equal.
 *
 * Return:	Success:	TRUE if the hyperslabs are equal (that is,
 *				both refer to exactly the same elements of an
 *				array)
 *
 *				FALSE otherwise.
 *
 *		Failure:	TRUE the rank is zero or if both hyperslabs
 *				are of zero size.
 *
 * Programmer:	Robb Matzke
 *		Friday, October 17, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
htri_t
H5V_hyper_eq(int n,
	     const hssize_t *offset1, const hsize_t *size1,
	     const hssize_t *offset2, const hsize_t *size2)
{
    hsize_t	nelmts1 = 1, nelmts2 = 1;
    int	i;

    if (n <= 0) return TRUE;

    for (i=0; i<n; i++) {
	if ((offset1 ? offset1[i] : 0) != (offset2 ? offset2[i] : 0)) {
	    return FALSE;
	}
	if ((size1 ? size1[i] : 0) != (size2 ? size2[i] : 0)) {
	    return FALSE;
	}
	if (0 == (nelmts1 *= (size1 ? size1[i] : 0))) return FALSE;
	if (0 == (nelmts2 *= (size2 ? size2[i] : 0))) return FALSE;
    }
    return TRUE;
}

/*-------------------------------------------------------------------------
 * Function:	H5V_hyper_disjointp
 *
 * Purpose:	Determines if two hyperslabs are disjoint.
 *
 * Return:	Success:	FALSE if they are not disjoint.
 *				TRUE if they are disjoint.
 *
 *		Failure:	A hyperslab of zero size is disjoint from all
 *				other hyperslabs.
 *
 * Programmer:	Robb Matzke
 *		Thursday, October 16, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
htri_t
H5V_hyper_disjointp(unsigned n,
		    const hssize_t *offset1, const hsize_t *size1,
		    const hssize_t *offset2, const hsize_t *size2)
{
    unsigned	u;

    if (!n || !size1 || !size2)	return TRUE;

    for (u=0; u<n; u++) {
        assert (size1[u]<HSSIZET_MAX);
        assert (size2[u]<HSSIZET_MAX);

        if (0==size1[u] || 0==size2[u])
            return TRUE;
        if (((offset1?offset1[u]:0) < (offset2?offset2[u]:0) &&
             ((offset1?offset1[u]:0) + (hssize_t)size1[u] <=
              (offset2?offset2[u]:0))) ||
            ((offset2?offset2[u]:0) < (offset1?offset1[u]:0) &&
             ((offset2?offset2[u]:0) + (hssize_t)size2[u] <=
              (offset1?offset1[u]:0)))) {
            return TRUE;
        }
    }
    return FALSE;
}

/*-------------------------------------------------------------------------
 * Function:	H5V_hyper_fill
 *
 * Purpose:	Similar to memset() except it operates on hyperslabs...
 *
 *		Fills a hyperslab of array BUF with some value VAL.  BUF
 *		is treated like a C-order array with N dimensions where the
 *		size of each dimension is TOTAL_SIZE[].	 The hyperslab which
 *		will be filled with VAL begins at byte offset OFFSET[] from
 *		the minimum corner of BUF and continues for SIZE[] bytes in
 *		each dimension.
 *		
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Friday, October 10, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5V_hyper_fill(unsigned n, const hsize_t *_size,
	       const hsize_t *total_size, const hssize_t *offset, void *_dst,
	       unsigned fill_value)
{
    uint8_t	*dst = (uint8_t*)_dst;	/*cast for ptr arithmetic	*/
    hsize_t	size[H5V_HYPER_NDIMS];	/*a modifiable copy of _size	*/
    hssize_t	dst_stride[H5V_HYPER_NDIMS]; /*destination stride info  */
    hsize_t	dst_start;		/*byte offset to start of stride*/
    hsize_t	elmt_size = 1;		/*bytes per element		*/
    herr_t	status;			/*function return status	*/
#ifndef NDEBUG
    unsigned	u;
#endif

    FUNC_ENTER(H5V_hyper_fill, FAIL);

    /* check args */
    assert(n > 0 && n <= H5V_HYPER_NDIMS);
    assert(_size);
    assert(total_size);
    assert(dst);
#ifndef NDEBUG
    for (u = 0; u < n; u++) {
        assert(_size[u] > 0);
        assert(total_size[u] > 0);
    }
#endif

    /* Copy the size vector so we can modify it */
    H5V_vector_cpy(n, size, _size);

    /* Compute an optimal destination stride vector */
    dst_start = H5V_hyper_stride(n, size, total_size, offset, dst_stride);
    H5V_stride_optimize1(&n, &elmt_size, size, dst_stride);

    /* Copy */
    status = H5V_stride_fill(n, elmt_size, size, dst_stride, dst+dst_start,
			     fill_value);

    FUNC_LEAVE(status);
}


/*-------------------------------------------------------------------------
 * Function:	H5V_hyper_copy
 *
 * Purpose:	Copies a hyperslab from the source to the destination.
 *
 *		A hyperslab is a logically contiguous region of
 *		multi-dimensional size SIZE of an array whose dimensionality
 *		is N and whose total size is DST_TOTAL_SIZE or SRC_TOTAL_SIZE.
 *		The minimum corner of the hyperslab begins at a
 *		multi-dimensional offset from the minimum corner of the DST
 *		(destination) or SRC (source) array.  The sizes and offsets
 *		are assumed to be in C order, that is, the first size/offset
 *		varies the slowest while the last varies the fastest in the
 *		mapping from N-dimensional space to linear space.  This
 *		function assumes that the array elements are single bytes (if
 *		your array has multi-byte elements then add an additional
 *		dimension whose size is that of your element).
 *
 *		The SRC and DST array may be the same array, but the results
 *		are undefined if the source hyperslab overlaps the
 *		destination hyperslab.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Friday, October 10, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5V_hyper_copy(unsigned n, const hsize_t *_size,

	       /*destination*/
	       const hsize_t *dst_size, const hssize_t *dst_offset,
	       void *_dst,

	       /*source*/
	       const hsize_t *src_size, const hssize_t *src_offset,
	       const void *_src)
{
    const uint8_t *src = (const uint8_t*)_src;	/*cast for ptr arithmtc */
    uint8_t	*dst = (uint8_t*) _dst;		/*cast for ptr arithmtc */
    hsize_t	size[H5V_HYPER_NDIMS];		/*a modifiable _size	*/
    hssize_t	src_stride[H5V_HYPER_NDIMS];	/*source stride info	*/
    hssize_t	dst_stride[H5V_HYPER_NDIMS];	/*dest stride info	*/
    hsize_t	dst_start, src_start;		/*offset to start at	*/
    hsize_t	elmt_size = 1;			/*element size in bytes */
    hsize_t tmp1;
    hsize_t tmp2;
    herr_t	status;				/*return status		*/
#ifndef NDEBUG		
    unsigned	u;
#endif

    FUNC_ENTER(H5V_hyper_copy, FAIL);

    /* check args */
    assert(n > 0 && n <= H5V_HYPER_NDIMS);
    assert(_size);
    assert(dst_size);
    assert(src_size);
    assert(dst);
    assert(src);
#ifndef NDEBUG
    for (u = 0; u < n; u++) {
        assert(_size[u] > 0);
        assert(dst_size[u] > 0);
        assert(src_size[u] > 0);
    }
#endif

    /* Copy the size vector so we can modify it */
    H5V_vector_cpy(n, size, _size);

    /* Compute stride vectors for source and destination */
#ifdef NO_INLINED_CODE
    dst_start = H5V_hyper_stride(n, size, dst_size, dst_offset, dst_stride);
    src_start = H5V_hyper_stride(n, size, src_size, src_offset, src_stride);
#else /* NO_INLINED_CODE */
    /* in-line version of two calls to H5V_hyper_stride() */
    {
        hsize_t	    dst_acc;	/*accumulator				*/
        hsize_t	    src_acc;	/*accumulator				*/
        int        ii;		    /*counter				*/

        /* init */
        dst_stride[n-1] = 1;
        src_stride[n-1] = 1;
        dst_start = dst_offset ? dst_offset[n-1] : 0;
        src_start = src_offset ? src_offset[n-1] : 0;

        /* Unroll loop for common cases */
        switch(n) {
            case 2:
                tmp1 = (dst_size[1] - size[1]);
                tmp2 = (src_size[1] - size[1]);
                assert (tmp1<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                assert (tmp2<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                dst_stride[0] = (hssize_t)tmp1; /*overflow checked*/
                src_stride[0] = (hssize_t)tmp2; /*overflow checked*/
                dst_acc = dst_size[1];
                src_acc = src_size[1];
                dst_start += dst_acc * (dst_offset ? dst_offset[0] : 0);
                src_start += src_acc * (src_offset ? src_offset[0] : 0);
                break;

            case 3:
                tmp1 = (dst_size[2] - size[2]);
                tmp2 = (src_size[2] - size[2]);
                assert (tmp1<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                assert (tmp2<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                dst_stride[1] = (hssize_t)tmp1; /*overflow checked*/
                src_stride[1] = (hssize_t)tmp2; /*overflow checked*/
                dst_acc = dst_size[2];
                src_acc = src_size[2];
                dst_start += dst_acc * (dst_offset ? dst_offset[1] : 0);
                src_start += src_acc * (src_offset ? src_offset[1] : 0);

                tmp1 = dst_acc * (dst_size[1] - size[1]);
                tmp2 = src_acc * (src_size[1] - size[1]);
                assert (tmp1<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                assert (tmp2<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                dst_stride[0] = (hssize_t)tmp1; /*overflow checked*/
                src_stride[0] = (hssize_t)tmp2; /*overflow checked*/
                dst_acc *= dst_size[1];
                src_acc *= src_size[1];
                dst_start += dst_acc * (dst_offset ? dst_offset[0] : 0);
                src_start += src_acc * (src_offset ? src_offset[0] : 0);
                break;

            case 4:
                tmp1 = (dst_size[3] - size[3]);
                tmp2 = (src_size[3] - size[3]);
                assert (tmp1<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                assert (tmp2<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                dst_stride[2] = (hssize_t)tmp1; /*overflow checked*/
                src_stride[2] = (hssize_t)tmp2; /*overflow checked*/
                dst_acc = dst_size[3];
                src_acc = src_size[3];
                dst_start += dst_acc * (dst_offset ? dst_offset[2] : 0);
                src_start += src_acc * (src_offset ? src_offset[2] : 0);

                tmp1 = dst_acc * (dst_size[2] - size[2]);
                tmp2 = src_acc * (src_size[2] - size[2]);
                assert (tmp1<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                assert (tmp2<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                dst_stride[1] = (hssize_t)tmp1; /*overflow checked*/
                src_stride[1] = (hssize_t)tmp2; /*overflow checked*/
                dst_acc *= dst_size[2];
                src_acc *= src_size[2];
                dst_start += dst_acc * (dst_offset ? dst_offset[1] : 0);
                src_start += src_acc * (src_offset ? src_offset[1] : 0);

                tmp1 = dst_acc * (dst_size[1] - size[1]);
                tmp2 = src_acc * (src_size[1] - size[1]);
                assert (tmp1<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                assert (tmp2<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                dst_stride[0] = (hssize_t)tmp1; /*overflow checked*/
                src_stride[0] = (hssize_t)tmp2; /*overflow checked*/
                dst_acc *= dst_size[1];
                src_acc *= src_size[1];
                dst_start += dst_acc * (dst_offset ? dst_offset[0] : 0);
                src_start += src_acc * (src_offset ? src_offset[0] : 0);
                break;

            default:
                /* others */
                for (ii=(int)(n-2), dst_acc=1, src_acc=1; ii>=0; --ii) {
                    tmp1 = dst_acc * (dst_size[ii+1] - size[ii+1]);
                    tmp2 = src_acc * (src_size[ii+1] - size[ii+1]);
                    assert (tmp1<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                    assert (tmp2<((hsize_t)1<<(8*sizeof(hssize_t)-1)));
                    dst_stride[ii] = (hssize_t)tmp1; /*overflow checked*/
                    src_stride[ii] = (hssize_t)tmp2; /*overflow checked*/
                    dst_acc *= dst_size[ii+1];
                    src_acc *= src_size[ii+1];
                    dst_start += dst_acc * (dst_offset ? dst_offset[ii] : 0);
                    src_start += src_acc * (src_offset ? src_offset[ii] : 0);
                }
                break;
        } /* end switch */
    }
#endif /* NO_INLINED_CODE */

    /* Optimize the strides as a pair */
    H5V_stride_optimize2(&n, &elmt_size, size, dst_stride, src_stride);

    /* Perform the copy in terms of stride */
    status = H5V_stride_copy(n, elmt_size, size,
             dst_stride, dst+dst_start, src_stride, src+src_start);

    FUNC_LEAVE(status);
}

/*-------------------------------------------------------------------------
 * Function:	H5V_stride_fill
 *
 * Purpose:	Fills all bytes of a hyperslab with the same value using
 *		memset().
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Saturday, October 11, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5V_stride_fill(unsigned n, hsize_t elmt_size, const hsize_t *size,
		const hssize_t *stride, void *_dst, unsigned fill_value)
{
    uint8_t	*dst = (uint8_t*)_dst; 	/*cast for ptr arithmetic	*/
    hsize_t	idx[H5V_HYPER_NDIMS]; 	/*1-origin indices		*/
    hsize_t	nelmts;			/*number of elements to fill	*/
    hsize_t	i;			/*counter			*/
    int	j;			/*counter			*/
    hbool_t	carry;			/*subtraction carray value	*/

    FUNC_ENTER(H5V_stride_fill, FAIL);
    assert (elmt_size < SIZET_MAX);

    H5V_vector_cpy(n, idx, size);
    nelmts = H5V_vector_reduce_product(n, size);
    for (i=0; i<nelmts; i++) {
        /* Copy an element */
        HDmemset(dst, (signed)fill_value, (size_t)elmt_size);

        /* Decrement indices and advance pointer */
        for (j=(int)(n-1), carry=TRUE; j>=0 && carry; --j) {
            dst += stride[j];

            if (--idx[j])
                carry = FALSE;
            else
                idx[j] = size[j];
        }
    }

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:	H5V_stride_copy
 *
 * Purpose:	Uses DST_STRIDE and SRC_STRIDE to advance through the arrays
 *		DST and SRC while copying bytes from SRC to DST.  This
 *		function minimizes the number of calls to memcpy() by
 *		combining various strides, but it will never touch memory
 *		outside the hyperslab defined by the strides.
 *
 * Note:	If the src_stride is all zero and elmt_size is one, then it's
 *		probably more efficient to use H5V_stride_fill() instead.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Saturday, October 11, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5V_stride_copy(unsigned n, hsize_t elmt_size, const hsize_t *size,
		const hssize_t *dst_stride, void *_dst,
		const hssize_t *src_stride, const void *_src)
{
    uint8_t	*dst = (uint8_t*)_dst;		/*cast for ptr arithmetic*/
    const uint8_t *src = (const uint8_t*) _src;	/*cast for ptr arithmetic*/
    hsize_t	idx[H5V_HYPER_NDIMS];		/*1-origin indices	*/
    hsize_t	nelmts;				/*num elements to copy	*/
    hsize_t	i;				/*counter		*/
    int	j;				/*counters		*/
    hbool_t	carry;				/*carray for subtraction*/

    FUNC_ENTER(H5V_stride_copy, FAIL);
    assert (elmt_size<SIZET_MAX);

    if (n) {
        H5V_vector_cpy(n, idx, size);
        nelmts = H5V_vector_reduce_product(n, size);
        for (i=0; i<nelmts; i++) {

            /* Copy an element */
	    H5_CHECK_OVERFLOW(elmt_size,hsize_t,size_t);
            HDmemcpy(dst, src, (size_t)elmt_size);

            /* Decrement indices and advance pointers */
            for (j=(int)(n-1), carry=TRUE; j>=0 && carry; --j) {
                src += src_stride[j];
                dst += dst_stride[j];

                if (--idx[j])
                    carry = FALSE;
                else
                    idx[j] = size[j];
            }
        }
    } else {
        H5_CHECK_OVERFLOW(elmt_size,hsize_t,size_t); /*check for overflow*/
        HDmemcpy (dst, src, (size_t)elmt_size);
        HRETURN (SUCCEED);
    }


    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:	H5V_stride_copy2
 *
 * Purpose:	Similar to H5V_stride_copy() except the source and
 *		destination each have their own dimensionality and size and
 *		we copy exactly NELMTS elements each of size ELMT_SIZE.	 The
 *		size counters wrap if NELMTS is more than a size counter.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Robb Matzke
 *		Saturday, October 11, 1997
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5V_stride_copy2(hsize_t nelmts, hsize_t elmt_size,

		 /* destination */
		 int dst_n, const hsize_t *dst_size,
		 const hssize_t *dst_stride,
		 void *_dst,

		 /* source */
		 int src_n, const hsize_t *src_size,
		 const hssize_t *src_stride,
		 const void *_src)
{
    uint8_t	*dst = (uint8_t *) _dst;
    const uint8_t *src = (const uint8_t *) _src;
    hsize_t	dst_idx[H5V_HYPER_NDIMS];
    hsize_t	src_idx[H5V_HYPER_NDIMS];
    hsize_t	i;
    int	j;
    hbool_t	carry;

    FUNC_ENTER(H5V_stride_copy2, FAIL);
    assert (elmt_size < SIZET_MAX);

    H5V_vector_cpy(dst_n, dst_idx, dst_size);
    H5V_vector_cpy(src_n, src_idx, src_size);

    for (i=0; i<nelmts; i++) {

	/* Copy an element */
	HDmemcpy(dst, src, (size_t)elmt_size);

	/* Decrement indices and advance pointers */
	for (j=dst_n-1, carry=TRUE; j>=0 && carry; --j) {
	    dst += dst_stride[j];
	    if (--dst_idx[j]) carry = FALSE;
	    else dst_idx[j] = dst_size[j];
	}
	for (j=src_n-1, carry=TRUE; j>=0 && carry; --j) {
	    src += src_stride[j];
	    if (--src_idx[j]) carry = FALSE;
	    else src_idx[j] = src_size[j];
	}
    }

    FUNC_LEAVE(SUCCEED);
}

/*-------------------------------------------------------------------------
 * Function:	H5V_array_fill
 *
 * Purpose:	Fills all bytes of an array with the same value using
 *		memset(). Increases amount copied by power of two until the
 *		halfway point is crossed, then copies the rest in one swoop.
 *
 * Return:	Non-negative on success/Negative on failure
 *
 * Programmer:	Quincey Koziol
 *		Thursday, June 18, 1998
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5V_array_fill(void *_dst, const void *src, size_t size, size_t count)
{
    size_t      copy_size;          /* size of the buffer to copy	*/
    size_t      copy_items;         /* number of items currently copying*/
    size_t      items_left;         /* number of items left to copy 	*/
    uint8_t     *dst=(uint8_t*)_dst;/* alias for pointer arithmetic	*/

    FUNC_ENTER(H5V_array_fill, FAIL);
    assert (dst);
    assert (src);
    assert (size < SIZET_MAX && size > 0);
    assert (count < SIZET_MAX && count > 0);

    HDmemcpy(dst, src, size);   /* copy first item */

    /* Initialize counters, etc. while compensating for first element copied */
    copy_size = size;
    copy_items = 1;
    items_left = count - 1;
    dst += size;

    /* copy until we've copied at least half of the items */
    while (items_left >= copy_items)
    {
        HDmemcpy(dst, _dst, copy_size);   /* copy the current chunk */
        dst += copy_size;     /* move the offset for the next chunk */
        items_left -= copy_items;   /* decrement the number of items left */

        copy_size *= 2;     /* increase the size of the chunk to copy */
        copy_items *= 2;    /* increase the count of items we are copying */
    }   /* end while */
    if (items_left > 0)   /* if there are any items left to copy */
        HDmemcpy(dst, _dst, items_left * size);

    FUNC_LEAVE(SUCCEED);
}   /* H5V_array_fill() */


/*-------------------------------------------------------------------------
 * Function:	H5V_array_offset
 *
 * Purpose:	Given a coordinate description of a location in an array, this
 *      function returns the byte offset of the coordinate.
 *
 *		The dimensionality of the whole array, the hyperslab, and the
 *		returned stride array is N.  The whole array dimensions are
 *		TOTAL_SIZE and the coordinate is at offset OFFSET.
 *
 * Return:	Success: Byte offset from beginning of array to start
 *				of striding.
 *
 *		Failure: abort() -- should never fail
 *
 * Programmer:	Quincey Koziol
 *		Tuesday, June 22, 1999
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
hsize_t
H5V_array_offset(unsigned n, const hsize_t *total_size, const hssize_t *offset)
{
    hsize_t	    skip;	/*starting point byte offset		*/
    hsize_t	    acc;	/*accumulator				*/
    int	    i;		/*counter				*/

    FUNC_ENTER(H5V_array_stride, (HDabort(), 0));

    assert(n <= H5V_HYPER_NDIMS);
    assert(total_size);
    assert(offset);

    /* others */
    for (i=(int)(n-1), acc=1, skip=0; i>=0; --i) {
        skip += acc * offset[i];
        acc *= total_size[i];
    }

    FUNC_LEAVE(skip);
}

