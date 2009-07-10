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

#ifndef H5DIFF_H__
#define H5DIFF_H__

#include "hdf5.h"
#include "h5trav.h"

/*-------------------------------------------------------------------------
 * command line options
 *-------------------------------------------------------------------------
 */

typedef struct {
    int    m_quiet;   /* quiet mode: no output at all */
    int    m_report;  /* report mode: print the data */
    int    m_verbose; /* verbose mode: print the data, list of objcets, warnings */
    int    d;         /* delta, absolute value to compare */
    double delta;     /* delta value */
    int    p;         /* relative error to compare*/
    double percent;   /* relative error value */
    int    n;         /* count, compare up to count */
    hsize_t count;    /* count value */
    int    err_stat;  /* an error ocurred (1, error, 0, no error) */
    int    cmn_objs;  /* do we have common objects */
    int    not_cmp;   /* are the objects comparable */
    int    contents;  /* equal contents */
    int    do_nans;   /* consider Nans while diffing floats */
    int    m_list_not_cmp;   /* list not comparable messages */
} diff_opt_t;


/*-------------------------------------------------------------------------
 * public functions
 *-------------------------------------------------------------------------
 */

#ifdef __cplusplus
extern "C" {
#endif

hsize_t  h5diff(const char *fname1,
                const char *fname2,
                const char *objname1,
                const char *objname2,
                diff_opt_t *options);

#ifdef H5_HAVE_PARALLEL
void phdiff_dismiss_workers(void);
void print_manager_output(void);
#endif

#ifdef __cplusplus
}
#endif



/*-------------------------------------------------------------------------
 * private functions
 *-------------------------------------------------------------------------
 */


hsize_t diff_dataset( hid_t file1_id,
                      hid_t file2_id,
                      const char *obj1_name,
                      const char *obj2_name,
                      diff_opt_t *options);

hsize_t diff_datasetid( hid_t dset1_id,
                        hid_t dset2_id,
                        const char *obj1_name,
                        const char *obj2_name,
                        diff_opt_t *options);

hsize_t diff( hid_t      file1_id,
              const char *path1,
              hid_t      file2_id,
              const char *path2,
              diff_opt_t *options,
              h5trav_type_t  type );

hsize_t diff_compare( hid_t file1_id,
                      const char *file1_name,
                      const char *obj1_name,
                      trav_info_t *info1,
                      hid_t file2_id,
                      const char *file2_name,
                      const char *obj2_name,
                      trav_info_t *info2,
                      diff_opt_t *options );

hsize_t diff_match( hid_t file1_id,
                    trav_info_t *info1,
                    hid_t file2_id,
                    trav_info_t *info2,
                    diff_opt_t *options );

hsize_t diff_array( void *_mem1,
                    void *_mem2,
                    hsize_t nelmts,
                    hsize_t hyper_start,
                    int rank,
                    hsize_t *dims,
                    diff_opt_t *options,
                    const char *name1,
                    const char *name2,
                    hid_t m_type,
                    hid_t container1_id,
                    hid_t container2_id); /* dataset where the reference came from*/


int diff_can_type( hid_t       f_type1, /* file data type */
                   hid_t       f_type2, /* file data type */
                   int         rank1,
                   int         rank2,
                   hsize_t     *dims1,
                   hsize_t     *dims2,
                   hsize_t     *maxdim1,
                   hsize_t     *maxdim2,
                   const char  *obj1_name,
                   const char  *obj2_name,
                   diff_opt_t  *options,
                   int         is_compound);


hsize_t diff_attr(hid_t loc1_id,
                  hid_t loc2_id,
                  const char *path1,
                  const char *path2,
                  diff_opt_t *options);


/*-------------------------------------------------------------------------
 * utility functions
 *-------------------------------------------------------------------------
 */

void        print_found(hsize_t nfound);
void        parallel_print(const char* format, ... );
void        print_type(hid_t type);
const char* diff_basename(const char *name);
const char* get_type(h5trav_type_t type);
const char* get_class(H5T_class_t tclass);
const char* get_sign(H5T_sign_t sign);
void        print_dimensions (int rank, hsize_t *dims);
int         print_objname(diff_opt_t *options, hsize_t nfound);
void        do_print_objname (const char *OBJ, const char *path1, const char *path2);



hsize_t diff_datum(void       *_mem1,
                   void       *_mem2,
                   hid_t      m_type,
                   hsize_t    i,
                   int        rank,
                   hsize_t    *dims,
                   hsize_t    *acc,
                   hsize_t    *pos,
                   diff_opt_t *options,
                   const char *obj1,
                   const char *obj2,
                   hid_t      container1_id,
                   hid_t      container2_id, /*where the reference came from*/
                   int        *ph);           /*print header */



hsize_t diff_float(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_double(unsigned char *mem1,
                    unsigned char *mem2,
                    hsize_t       nelmts,
                    hsize_t       hyper_start,
                    int           rank,
                    hsize_t       *dims,
                    hsize_t       *acc,
                    hsize_t       *pos,
                    diff_opt_t    *options,
                    const char    *obj1,
                    const char    *obj2,
                    int           *ph);

#if H5_SIZEOF_LONG_DOUBLE !=0

hsize_t diff_ldouble(unsigned char *mem1,
                     unsigned char *mem2,
                     hsize_t       nelmts,
                     hsize_t       hyper_start,
                     int           rank,
                     hsize_t       *dims,
                     hsize_t       *acc,
                     hsize_t       *pos,
                     diff_opt_t    *options,
                     const char    *obj1,
                     const char    *obj2,
                     int           *ph);

#endif

hsize_t diff_schar(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_uchar(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_short(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_ushort(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_int(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_uint(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_long(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_ulong(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_llong(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

hsize_t diff_ullong(unsigned char *mem1,
                   unsigned char *mem2,
                   hsize_t       nelmts,
                   hsize_t       hyper_start,
                   int           rank,
                   hsize_t       *dims,
                   hsize_t       *acc,
                   hsize_t       *pos,
                   diff_opt_t    *options,
                   const char    *obj1,
                   const char    *obj2,
                   int           *ph);

#endif  /* H5DIFF_H__ */

