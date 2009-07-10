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

#include <stdlib.h>
#include <string.h>
#include "h5diff.h"
#include "h5diff_common.h"
#include "h5tools_utils.h"

static int check_n_input( const char* );
static int check_p_input( const char* );
static int check_d_input( const char* );


/* module-scoped variables */
const char  *progname = "h5diff";

/*
 * Command-line options: The user can specify short or long-named
 * parameters.
 */
static const char *s_opts = "hVrvqn:d:p:Nc";
static struct long_options l_opts[] = {
    { "help", no_arg, 'h' },
    { "version", no_arg, 'V' },
    { "report", no_arg, 'r' },
    { "verbose", no_arg, 'v' },
    { "quiet", no_arg, 'q' },
    { "count", require_arg, 'n' },
    { "delta", require_arg, 'd' },
    { "relative", require_arg, 'p' },
    { "nan", no_arg, 'N' },
    { "compare", no_arg, 'c' },
    { NULL, 0, '\0' }
};


/*-------------------------------------------------------------------------
 * Function: parse_command_line
 *
 * Purpose: parse command line input
 *
 *-------------------------------------------------------------------------
 */

void parse_command_line(int argc,
                        const char* argv[],
                        const char** fname1,
                        const char** fname2,
                        const char** objname1,
                        const char** objname2,
                        diff_opt_t* options)
{

    int opt;

    /* process the command-line */
    memset(options, 0, sizeof (diff_opt_t));

    /* assume equal contents initially */
    options->contents = 1;

    /* NaNs are handled by default */
    options->do_nans = 1;

    /* parse command line options */
    while ((opt = get_option(argc, argv, s_opts, l_opts)) != EOF)
    {
        switch ((char)opt)
        {
        default:
            usage();
            h5diff_exit(EXIT_FAILURE);
        case 'h':
            usage();
            h5diff_exit(EXIT_SUCCESS);
        case 'V':
            print_version(progname);
            h5diff_exit(EXIT_SUCCESS);
        case 'v':
            options->m_verbose = 1;
            break;
        case 'q':
            /* use quiet mode; supress the message "0 differences found" */
            options->m_quiet = 1;
            break;
        case 'r':
            options->m_report = 1;
            break;
        case 'd':
            options->d=1;

            if ( check_d_input( opt_arg )==-1)
            {
                printf("<-d %s> is not a valid option\n", opt_arg );
                usage();
                h5diff_exit(EXIT_FAILURE);
            }
            options->delta = atof( opt_arg );
            break;

        case 'p':

            options->p=1;
            if ( check_p_input( opt_arg )==-1)
            {
                printf("<-p %s> is not a valid option\n", opt_arg );
                usage();
                h5diff_exit(EXIT_FAILURE);
            }
            options->percent = atof( opt_arg );
            break;

        case 'n':

            options->n=1;
            if ( check_n_input( opt_arg )==-1)
            {
                printf("<-n %s> is not a valid option\n", opt_arg );
                usage();
                h5diff_exit(EXIT_FAILURE);
            }
            options->count = atol( opt_arg );

            break;

        case 'N':
            options->do_nans = 0;
            break;
        case 'c':
            options->m_list_not_cmp = 1;
            break;
        }
    }

    /* check for file names to be processed */
    if (argc <= opt_ind || argv[ opt_ind + 1 ] == NULL)
    {
        error_msg(progname, "missing file names\n");
        usage();
        h5diff_exit(EXIT_FAILURE);
    }

    *fname1 = argv[ opt_ind ];
    *fname2 = argv[ opt_ind + 1 ];
    *objname1 = argv[ opt_ind + 2 ];

    if ( *objname1 == NULL )
    {
        *objname2 = NULL;
        return;
    }

    if ( argv[ opt_ind + 3 ] != NULL)
    {
        *objname2 = argv[ opt_ind + 3 ];
    }
    else
    {
        *objname2 = *objname1;
    }


}


/*-------------------------------------------------------------------------
 * Function: print_info
 *
 * Purpose: print several information messages after h5diff call
 *
 *-------------------------------------------------------------------------
 */

 void  print_info(diff_opt_t* options)
 {
     if (options->m_quiet || options->err_stat )
         return;

     if (options->cmn_objs==0)
     {
         printf("No common objects found. Files are not comparable.\n");
         if (!options->m_verbose)
             printf("Use -v for a list of objects.\n");
     }

     if (options->not_cmp==1)
     {
         if ( options->m_list_not_cmp == 0 )
         {
             printf("--------------------------------\n");
             printf("Some objects are not comparable\n");
             printf("--------------------------------\n");
             printf("Use -c for a list of objects.\n");
         }
        
             
     }

 }

/*-------------------------------------------------------------------------
 * Function: check_n_input
 *
 * Purpose: check for valid input
 *
 * Return: 1 for ok, -1 for fail
 *
 * Programmer: Pedro Vicente, pvn@ncsa.uiuc.edu
 *
 * Date: May 9, 2003
 *
 * Comments:
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static int
check_n_input( const char *str )
{
    unsigned i;
    char c;

    for ( i = 0; i < strlen(str); i++)
    {
        c = str[i];
        if ( i==0 )
        {
            if ( c < 49 || c > 57  ) /* ascii values between 1 and 9 */
                return -1;
        }
        else
            if ( c < 48 || c > 57  ) /* 0 also */
                return -1;
    }
    return 1;
}

/*-------------------------------------------------------------------------
 * Function: check_p_input
 *
 * Purpose: check for a valid p option input
 *
 * Return: 1 for ok, -1 for fail
 *
 * Date: May 9, 2003
 *
 * Comments:
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static int
check_p_input( const char *str )
{
    double x;

    /*
    the atof return value on a hexadecimal input is different
    on some systems; we do a character check for this
    */
    if (strlen(str)>2 && str[0]=='0' && str[1]=='x')
        return -1;

    x=atof(str);
    if (x<=0)
        return -1;

    return 1;
}

/*-------------------------------------------------------------------------
 * Function: check_d_input
 *
 * Purpose: check for a valid d option input
 *
 * Return: 1 for ok, -1 for fail
 *
 * Date: November 11, 2007
 *
 * Comments:
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static int
check_d_input( const char *str )
{
    double x;

    /*
    the atof return value on a hexadecimal input is different
    on some systems; we do a character check for this
    */
    if (strlen(str)>2 && str[0]=='0' && str[1]=='x')
        return -1;

    x=atof(str);
    if (x <=0)
        return -1;

    return 1;
}

/*-------------------------------------------------------------------------
 * Function: usage
 *
 * Purpose: print a usage message
 *
 * Return: void
 *
 *-------------------------------------------------------------------------
 */

void usage(void)
{
 printf("usage: h5diff [OPTIONS] file1 file2 [obj1[obj2]] \n");
 printf("  file1                    File name of the first HDF5 file\n");
 printf("  file2                    File name of the second HDF5 file\n");
 printf("  [obj1]                   Name of an HDF5 object, in absolute path\n");
 printf("  [obj2]                   Name of an HDF5 object, in absolute path\n");

 printf("  OPTIONS\n");

 printf("   -h, --help              Print a usage message and exit\n");
 printf("   -V, --version           Print version number and exit\n");
 printf("   -r, --report            Report mode. Print differences\n");
 printf("   -v, --verbose           Verbose mode. Print differences, list of objects\n");
 printf("   -q, --quiet             Quiet mode. Do not do output\n");
 printf("   -c, --compare           List objects that are not comparable\n");
 printf("   -N, --nan               Avoid NaNs detection\n");

 printf("   -n C, --count=C         Print differences up to C number\n");
 printf("   -d D, --delta=D         Print difference when greater than limit D\n");
 printf("   -p R, --relative=R      Print difference when greater than relative limit R\n");


 printf("\n");

 printf("  C - is a positive integer\n");
 printf("  D - is a positive number. Compare criteria is |a - b| > D\n");
 printf("  R - is a positive number. Compare criteria is |(b-a)/a| > R\n");

 printf("\n");

 printf(" Modes of output:\n");
 printf("\n");
 printf("  Default mode: print the number of differences found and where they occured\n");
 printf("  -r Report mode: print the above plus the differences\n");
 printf("  -v Verbose mode: print the above plus a list of objects and warnings\n");
 printf("  -q Quiet mode: do not print output\n");

 printf("\n");

 printf(" Compare criteria\n");
 printf("\n");
 printf(" If no objects [obj1[obj2]] are specified, h5diff only compares objects\n");
 printf("   with the same absolute path in both files\n");
 printf("\n");

 printf(" The compare criteria is:\n");
 printf("   1) datasets: numerical array differences 2) groups: name string difference\n");
 printf("   3) datatypes: the return value of H5Tequal 4) links: name string difference\n");
 printf("   of the linked value\n");

 printf("\n");

 printf(" Return exit code:\n");
 printf("\n");
 printf("  1 if differences found, 0 if no differences, 2 if error\n");

 printf("\n");

 printf(" Examples of use:\n");
 printf("\n");
 printf(" 1) h5diff file1 file2 /g1/dset1 /g1/dset2\n");
 printf("\n");
 printf("    Compares object '/g1/dset1' in file1 with '/g1/dset2' in file2\n");
 printf("\n");
 printf(" 2) h5diff file1 file2 /g1/dset1\n");
 printf("\n");
 printf("    Compares object '/g1/dset1' in both files\n");
 printf("\n");
 printf(" 3) h5diff file1 file2\n");
 printf("\n");
 printf("    Compares all objects in both files\n");
 printf("\n");
 printf(" Note)  file1 and file2 can be the same file. Use\n");
 printf("\n");
 printf("    h5diff file1 file1 /g1/dset1 /g1/dset2\n");
 printf("\n");
 printf("    to compare '/g1/dset1' and '/g1/dset2' in the same file\n");
 printf("\n");


}


/*-------------------------------------------------------------------------
 * Function: h5diff_exit
 *
 * Purpose: dismiss phdiff worker processes and exit
 *
 * Return: none
 *
 * Programmer: Albert Cheng
 * Date: Feb 6, 2005
 *
 * Comments:
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
void h5diff_exit(int status)
{
#ifdef H5_HAVE_PARALLEL
    /* if in parallel mode, dismiss workers, close down MPI, then exit */
    if((g_nTasks > 1) && g_Parallel) {
        phdiff_dismiss_workers();
        MPI_Barrier(MPI_COMM_WORLD);
    }
    if(g_Parallel)
        MPI_Finalize();
#endif
    exit(status);
}

