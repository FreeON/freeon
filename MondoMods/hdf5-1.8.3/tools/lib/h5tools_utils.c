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

/*
 * Programmer:  Bill Wendling <wendling@ncsa.uiuc.edu>
 *              Tuesday, 6. March 2001
 */

/*
 * Portions of this work are derived from _Obfuscated C and Other Mysteries_,
 * by Don Libes, copyright (c) 1993 by John Wiley & Sons, Inc.
 */

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

#include "h5tools_utils.h"
#include "H5private.h"
#include "h5trav.h"

/* global variables */
int   nCols = 80;

/* ``get_option'' variables */
int         opt_err = 1;    /*get_option prints errors if this is on */
int         opt_ind = 1;    /*token pointer                          */
const char *opt_arg;        /*flag argument (or value)               */

/* local functions */
static void init_table(table_t **tbl);
#ifdef H5DUMP_DEBUG
static void dump_table(char* tablename, table_t *table);
#endif  /* H5DUMP_DEBUG */
static void add_obj(table_t *table, haddr_t objno, const char *objname, hbool_t recorded);


/*-------------------------------------------------------------------------
 * Function:    error_msg
 *
 * Purpose: Print a nicely formatted error message to stderr flushing the
 *              stdout stream first.
 *
 * Return:  Nothing
 *
 * Programmer:  Bill Wendling
 *              Tuesday, 20. February 2001
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
error_msg(const char *progname, const char *fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    HDfflush(stdout);
    HDfprintf(stderr, "%s error: ", progname);
    HDvfprintf(stderr, fmt, ap);

    va_end(ap);
}


/*-------------------------------------------------------------------------
 * Function:    warn_msg
 *
 * Purpose: Print a nicely formatted warning message to stderr flushing
 *              the stdout stream first.
 *
 * Return:  Nothing
 *
 * Programmer:  Bill Wendling
 *              Tuesday, 20. February 2001
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
warn_msg(const char *progname, const char *fmt, ...)
{
    va_list ap;

    va_start(ap, fmt);
    HDfflush(stdout);
    HDfprintf(stderr, "%s warning: ", progname);
    HDvfprintf(stderr, fmt, ap);
    va_end(ap);
}


/*-------------------------------------------------------------------------
 * Function:    get_option
 *
 * Purpose: Determine the command-line options a user specified. We can
 *      accept both short and long type command-lines.
 *
 * Return:  Success:    The short valued "name" of the command line
 *              parameter or EOF if there are no more
 *              parameters to process.
 *
 *      Failure:    A question mark.
 *
 * Programmer:  Bill Wendling
 *              Friday, 5. January 2001
 *
 * Modifications: Pedro Vicente
 *                October, 27 2008
 * Wilcard "*" argument type
 *
 *-------------------------------------------------------------------------
 */
int
get_option(int argc, const char **argv, const char *opts, const struct long_options *l_opts)
{
    static int sp = 1;    /* character index in current token */
    int opt_opt = '?';    /* option character passed back to user */

    if (sp == 1) {
        /* check for more flag-like tokens */
        if (opt_ind >= argc || argv[opt_ind][0] != '-' || argv[opt_ind][1] == '\0') {
            return EOF;
        } else if (HDstrcmp(argv[opt_ind], "--") == 0) {
            opt_ind++;
            return EOF;
        }
    }

    if (sp == 1 && argv[opt_ind][0] == '-' && argv[opt_ind][1] == '-') {
        /* long command line option */
        const char *arg = &argv[opt_ind][2];
        int i;

        for (i = 0; l_opts && l_opts[i].name; i++) {
            size_t len = HDstrlen(l_opts[i].name);

            if (HDstrncmp(arg, l_opts[i].name, len) == 0) {
                /* we've found a matching long command line flag */
                opt_opt = l_opts[i].shortval;

                if (l_opts[i].has_arg != no_arg) {
                    if (arg[len] == '=') {
                        opt_arg = &arg[len + 1];
                    } else if (opt_ind < (argc - 1) && argv[opt_ind + 1][0] != '-') {
                        opt_arg = argv[++opt_ind];
                    } else if (l_opts[i].has_arg == require_arg) {
                        if (opt_err)
                            HDfprintf(stderr,
                                    "%s: option required for \"--%s\" flag\n",
                                    argv[0], arg);

                        opt_opt = '?';
                    }
                } else {
                    if (arg[len] == '=') {
                        if (opt_err)
                            HDfprintf(stderr,
                                    "%s: no option required for \"%s\" flag\n",
                                    argv[0], arg);

                        opt_opt = '?';
                    }

                    opt_arg = NULL;
                }

                break;
            }
        }

        if (l_opts[i].name == NULL) {
            /* exhausted all of the l_opts we have and still didn't match */
            if (opt_err)
                HDfprintf(stderr, "%s: unknown option \"%s\"\n", argv[0], arg);

            opt_opt = '?';
        }

        opt_ind++;
        sp = 1;
    } else {
        register char *cp;    /* pointer into current token */

        /* short command line option */
        opt_opt = argv[opt_ind][sp];

        if (opt_opt == ':' || (cp = strchr(opts, opt_opt)) == 0) {
            if (opt_err)
                HDfprintf(stderr, "%s: unknown option \"%c\"\n",
                        argv[0], opt_opt);

            /* if no chars left in this token, move to next token */
            if (argv[opt_ind][++sp] == '\0') {
                opt_ind++;
                sp = 1;
            }

            return '?';
        }

        if (*++cp == ':') {
            /* if a value is expected, get it */
            if (argv[opt_ind][sp + 1] != '\0') {
                /* flag value is rest of current token */
                opt_arg = &argv[opt_ind++][sp + 1];
            } else if (++opt_ind >= argc) {
                if (opt_err)
                    HDfprintf(stderr,
                            "%s: value expected for option \"%c\"\n",
                            argv[0], opt_opt);

                opt_opt = '?';
            } else {
                /* flag value is next token */
                opt_arg = argv[opt_ind++];
            }

            sp = 1;
        } 
        
        /* wildcard argument */
        else if (*cp == '*')
        {
            /* check the next argument */
            opt_ind++;
            /* we do have an extra argument, check if not last */
            if ( argv[opt_ind][0] != '-' && (opt_ind+1) < argc )
            {
                opt_arg = argv[opt_ind++];
            }
            else
            {
                opt_arg = NULL;
            }
        }
        
        else 
        {
            /* set up to look at next char in token, next time */
            if (argv[opt_ind][++sp] == '\0') {
                /* no more in current token, so setup next token */
                opt_ind++;
                sp = 1;
                
                
            }

            opt_arg = NULL;
        }
    }

    /* return the current flag character found */
    return opt_opt;
}


/*-------------------------------------------------------------------------
 * Function:    indentation
 *
 * Purpose:     Print spaces for indentation
 *
 * Return:      void
 *
 * Programmer:  Ruey-Hsia Li
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
indentation(int x)
{
    if (x < nCols) {
        while (x-- > 0)
            printf(" ");
    } else {
        HDfprintf(stderr, "error: the indentation exceeds the number of cols.\n");
        exit(1);
    }
}


/*-------------------------------------------------------------------------
 * Function:    print_version
 *
 * Purpose:     Print the program name and the version information which is
 *      defined the same as the HDF5 library version.
 *
 * Return:      void
 *
 * Programmer:  unknown
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
print_version(const char *progname)
{
    printf("%s: Version %u.%u.%u%s%s\n",
           progname, H5_VERS_MAJOR, H5_VERS_MINOR, H5_VERS_RELEASE,
           H5_VERS_SUBRELEASE[0] ? "-" : "", H5_VERS_SUBRELEASE);
}


/*-------------------------------------------------------------------------
 * Function:    init_table
 *
 * Purpose:     allocate and initialize tables for shared groups, datasets,
 *              and committed types
 *
 * Return:      void
 *
 * Programmer:  Ruey-Hsia Li
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void
init_table(table_t **tbl)
{
    table_t *table = HDmalloc(sizeof(table_t));

    table->size = 20;
    table->nobjs = 0;
    table->objs = HDmalloc(table->size * sizeof(obj_t));

    *tbl = table;
}


/*-------------------------------------------------------------------------
 * Function:    free_table
 *
 * Purpose:     free tables for shared groups, datasets,
 *              and committed types
 *
 * Return:      void
 *
 * Programmer:  Paul Harten
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
free_table(table_t *table)
{
    unsigned u;         /* Local index value */

    /* Free the names for the objects in the table */
    for(u = 0; u < table->nobjs; u++)
        if(table->objs[u].objname)
            HDfree(table->objs[u].objname);

    HDfree(table->objs);
}

#ifdef H5DUMP_DEBUG

/*-------------------------------------------------------------------------
 * Function:    dump_table
 *
 * Purpose:     display the contents of tables for debugging purposes
 *
 * Return:      void
 *
 * Programmer:  Ruey-Hsia Li
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void
dump_table(char* tablename, table_t *table)
{
    unsigned u;

    printf("%s: # of entries = %d\n", tablename,table->nobjs);
    for (u = 0; u < table->nobjs; u++)
    HDfprintf(stdout,"%a %s %d %d\n", table->objs[u].objno,
           table->objs[u].objname,
           table->objs[u].displayed, table->objs[u].recorded);
}


/*-------------------------------------------------------------------------
 * Function:    dump_tables
 *
 * Purpose:     display the contents of tables for debugging purposes
 *
 * Return:      void
 *
 * Programmer:  Ruey-Hsia Li
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void
dump_tables(find_objs_t *info)
{
    dump_table("group_table", info->group_table);
    dump_table("dset_table", info->dset_table);
    dump_table("type_table", info->type_table);
}
#endif  /* H5DUMP_DEBUG */


/*-------------------------------------------------------------------------
 * Function:    search_obj
 *
 * Purpose:     search the object specified by objno in the table
 *
 * Return:      Success:    an integer, the location of the object
 *
 *              Failure:    FAIL   if object is not found
 *
 * Programmer:  Ruey-Hsia Li
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
obj_t *
search_obj(table_t *table, haddr_t objno)
{
    unsigned u;

    for(u = 0; u < table->nobjs; u++)
        if(table->objs[u].objno == objno)
            return &(table->objs[u]);

    return NULL;
}


/*-------------------------------------------------------------------------
 * Function:    find_objs_cb
 *
 * Purpose:     Callback to find objects, committed types and store them in tables
 *
 * Return:      Success:    SUCCEED
 *
 *              Failure:    FAIL
 *
 * Programmer:  Ruey-Hsia Li
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
find_objs_cb(const char *name, const H5O_info_t *oinfo, const char *already_seen,
    void *op_data)
{
    find_objs_t *info = (find_objs_t*)op_data;
    herr_t ret_value = 0;

    switch(oinfo->type) {
        case H5O_TYPE_GROUP:
            if(NULL == already_seen)
                add_obj(info->group_table, oinfo->addr, name, TRUE);
            break;

        case H5O_TYPE_DATASET:
            if(NULL == already_seen) {
                hid_t dset;

                /* Add the dataset to the list of objects */
                add_obj(info->dset_table, oinfo->addr, name, TRUE);

                /* Check for a dataset that uses a named datatype */
                if((dset = H5Dopen2(info->fid, name, H5P_DEFAULT)) >= 0) {
                    hid_t type = H5Dget_type(dset);

                    if(H5Tcommitted(type) > 0) {
                        H5O_info_t type_oinfo;

                        H5Oget_info(type, &type_oinfo);
                        if(search_obj(info->type_table, type_oinfo.addr) == NULL)
                            add_obj(info->type_table, type_oinfo.addr, name, FALSE);
                    } /* end if */

                    H5Tclose(type);
                    H5Dclose(dset);
                } /* end if */
                else
                    ret_value = FAIL;
            } /* end if */
            break;

        case H5O_TYPE_NAMED_DATATYPE:
            if(NULL == already_seen) {
                obj_t *found_obj;

                if((found_obj = search_obj(info->type_table, oinfo->addr)) == NULL)
                    add_obj(info->type_table, oinfo->addr, name, TRUE);
                else {
                    /* Use latest version of name */
                    HDfree(found_obj->objname);
                    found_obj->objname = HDstrdup(name);

                    /* Mark named datatype as having valid name */
                    found_obj->recorded = TRUE;
                } /* end else */
            } /* end if */
            break;

        default:
            break;
    } /* end switch */

    return ret_value;
}


/*-------------------------------------------------------------------------
 * Function:    init_objs
 *
 * Purpose:     Initialize tables for groups, datasets & named datatypes
 *
 * Return:      Success:    SUCCEED
 *
 *              Failure:    FAIL
 *
 * Programmer:  Ruey-Hsia Li
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
init_objs(hid_t fid, find_objs_t *info, table_t **group_table,
    table_t **dset_table, table_t **type_table)
{
    /* Initialize the tables */
    init_table(group_table);
    init_table(dset_table);
    init_table(type_table);

    /* Init the find_objs_t */
    info->fid = fid;
    info->group_table = *group_table;
    info->type_table = *type_table;
    info->dset_table = *dset_table;

    /* Find all shared objects */
    return(h5trav_visit(fid, "/", TRUE, TRUE, find_objs_cb, NULL, info));
}


/*-------------------------------------------------------------------------
 * Function:    add_obj
 *
 * Purpose:     add a shared object to the table
 *              realloc the table if necessary
 *
 * Return:      void
 *
 * Programmer:  Ruey-Hsia Li
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void
add_obj(table_t *table, haddr_t objno, const char *objname, hbool_t record)
{
    unsigned u;

    /* See if we need to make table larger */
    if(table->nobjs == table->size) {
        table->size *= 2;
        table->objs = HDrealloc(table->objs, table->size * sizeof(table->objs[0]));
    } /* end if */

    /* Increment number of objects in table */
    u = table->nobjs++;

    /* Set information about object */
    table->objs[u].objno = objno;
    table->objs[u].objname = HDstrdup(objname);
    table->objs[u].recorded = record;
    table->objs[u].displayed = 0;
}


#ifndef H5_HAVE_TMPFILE
/*-------------------------------------------------------------------------
 * Function:    tmpfile
 *
 * Purpose:     provide tmpfile() function when it is not supported by the
 *              system.  Always return NULL for now.
 *
 * Return:      a stream description when succeeds.
 *              NULL if fails.
 *
 * Programmer:  Albert Cheng, 2005/8/9
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
FILE *
tmpfile(void)
{
    return NULL;
}

#endif
