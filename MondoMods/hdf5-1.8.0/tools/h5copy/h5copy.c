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


#include "h5tools.h"
#include "h5tools_utils.h"
#include <string.h>
#include <stdlib.h>

const char *progname="h5copy";
int   d_status;

/* command-line options: short and long-named parameters */
static const char *s_opts = "d:f:hi:o:ps:vV";
static struct long_options l_opts[] = {
    { "destination", require_arg, 'd' },
    { "flag", require_arg, 'f' },
    { "help", no_arg, 'h' },
    { "input", require_arg, 'i' },
    { "output", require_arg, 'o' },
    { "parents", no_arg, 'p' },
    { "source", require_arg, 's' },
    { "verbose", no_arg, 'v' },
    { "version", no_arg, 'V' },
    { NULL, 0, '\0' }
};

/*-------------------------------------------------------------------------
 * Function:    leave
 *
 * Purpose:     Shutdown MPI & HDF5 and call exit()
 *
 * Return:      Does not return
 *
 * Programmer:  Quincey Koziol
 *              Saturday, 31. January 2004
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void
leave(int ret)
{
 h5tools_close();
 
 exit(ret);
}


/*-------------------------------------------------------------------------
 * Function: usage
 *
 * Purpose: Prints a usage message on stderr and then returns.
 *
 * Return: void
 *
 * Programmer: Pedro Vicente Nunes, 7/8/2006
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void
usage (void)
{
 fprintf(stdout, "\
usage: h5copy [OPTIONS] [OBJECTS...]\n\
   OBJECTS\n\
      -i, --input        input file name\n\
      -o, --output       output file name\n\
      -s, --source       source object name\n\
      -d, --destination  destination object name\n\
   OPTIONS\n\
      -h, --help         Print a usage message and exit\n\
      -p, --parents      No error if existing, make parent groups as needed\n\
      -v, --verbose      Print information about OBJECTS and OPTIONS\n\
      -V, --version      Print version number and exit\n\
      -f, --flag         Flag type\n\n\
      Flag type is one of the following strings:\n\n\
      shallow     Copy only immediate members for groups\n\
      soft        Expand soft links into new objects\n\
      ext         Expand external links into new objects\n\
      ref         Copy objects that are pointed by references\n\
      noattr      Copy object without copying attributes\n\
      allflags    Switches all flags from the default to the non-default setting\n\n\
      These flag types correspond to the following API symbols\n\n\
      H5O_COPY_SHALLOW_HIERARCHY_FLAG\n\
      H5O_COPY_EXPAND_SOFT_LINK_FLAG\n\
      H5O_COPY_EXPAND_EXT_LINK_FLAG\n\
      H5O_COPY_EXPAND_REFERENCE_FLAG\n\
      H5O_COPY_WITHOUT_ATTR_FLAG\n\
      H5O_COPY_ALL\n");
}



/*-------------------------------------------------------------------------
 * Function: parse_flag
 *
 * Purpose: read the flag -f STRING
 *
 * STRING is one of the following (API symbol and description)
 *
 * shallow  H5O_COPY_SHALLOW_HIERARCHY_FLAG:  Copy only immediate members for groups
 * soft     H5O_COPY_EXPAND_SOFT_LINK_FLAG:  Expand soft links into new objects
 * ext      H5O_COPY_EXPAND_EXT_LINK_FLAG: Expand external links into new objects
 * ref      H5O_COPY_EXPAND_OBJ_REFERENCE_FLAG: Copy objects that are pointed by references
 * noattr   H5O_COPY_WITHOUT_ATTR_FLAG Copy object without copying attributes 
 * allflags Switches all flags from the default to the non-default setting 
 *
 * Return: Success:    SUCCEED
 *         Failure:    FAIL
 *
 * Programmer: Pedro Vicente Nunes, 7/8/2006
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */


static int parse_flag(const char* str_flag, unsigned *flag)
{
 unsigned fla=0;

 if (strcmp(str_flag,"shallow")==0)
 {
  fla = H5O_COPY_SHALLOW_HIERARCHY_FLAG;
 } 
 else  if (strcmp(str_flag,"soft")==0)
 {
  fla = H5O_COPY_EXPAND_SOFT_LINK_FLAG;
 }
 else  if (strcmp(str_flag,"ext")==0)
 {
  fla = H5O_COPY_EXPAND_EXT_LINK_FLAG;
 }
 else  if (strcmp(str_flag,"ref")==0)
 {
  fla = H5O_COPY_EXPAND_REFERENCE_FLAG;
 }
 else  if (strcmp(str_flag,"noattr")==0)
 {
  fla = H5O_COPY_WITHOUT_ATTR_FLAG;
 }
 else  if (strcmp(str_flag,"allflags")==0)
 {
  fla = H5O_COPY_ALL;
 }
 else  if (strcmp(str_flag,"nullmsg")==0)
 {
  fla = H5O_COPY_PRESERVE_NULL_FLAG;
 }
 else
 {
  error_msg(progname, "Error in input flag\n");
  return -1;
 }

 *flag = fla;

 return 0;
}

/*-------------------------------------------------------------------------
 * Function: main
 *
 * Purpose: main program
 *
 * Programmer: Pedro Vicente Nunes
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */

int
main (int argc, const char *argv[])
{
 hid_t        fid_src=-1;
 hid_t        fid_dst=-1;
 char         *fname_src=NULL;
 char         *fname_dst=NULL;
 char         *oname_src=NULL;
 char         *oname_dst=NULL;
 unsigned     flag=0;
 unsigned     verbose=0;
 unsigned     parents=0;
 hid_t        ocpl_id = (-1);          /* Object copy property list */
 hid_t        lcpl_id = (-1);          /* Link creation property list */
 char         str_flag[20];
 int          opt;

 /* initialize h5tools lib */
 h5tools_init();

 /* Check for no command line parameters */
 if(argc == 1) {
   usage();
   leave(EXIT_SUCCESS);
 } /* end if */

 /* parse command line options */
 while ((opt = get_option(argc, argv, s_opts, l_opts)) != EOF) 
 {
  switch ((char)opt) 
  {
  case 'd':
   oname_dst = strdup(opt_arg);
   break;
   
  case 'f':
   /* validate flag */
   if (parse_flag(opt_arg,&flag)<0)
   {
    usage();
    leave(EXIT_FAILURE);
   }
   strcpy(str_flag,opt_arg);
   break;

  case 'h':
   usage();
   leave(EXIT_SUCCESS);
   break;

  case 'i':
   fname_src = strdup(opt_arg);
   break;

  case 'o':
   fname_dst = strdup(opt_arg);
   break;

  case 'p':
   parents = 1;
   break;

  case 's':
   oname_src = strdup(opt_arg);
   break;

  case 'V':
   print_version(progname);
   leave(EXIT_SUCCESS);
   break;

  case 'v':
   verbose = 1;
   break;

  default:
   usage();
   leave(EXIT_FAILURE);
  }
 }

/*-------------------------------------------------------------------------
 * check for missing file/object names
 *-------------------------------------------------------------------------*/

 if (fname_src==NULL) 
 {
  error_msg(progname, "Input file name missing\n");
  usage();
  leave(EXIT_FAILURE);
 }

 if (fname_dst==NULL) 
 {
  error_msg(progname, "Output file name missing\n");
  usage();
  leave(EXIT_FAILURE);
 }

 if (oname_src==NULL) 
 {
  error_msg(progname, "Source object name missing\n");
  usage();
  leave(EXIT_FAILURE);
 }

 if (oname_dst==NULL) 
 {
  error_msg(progname, "Destination object name missing\n");
  usage();
  leave(EXIT_FAILURE);
 }


/*-------------------------------------------------------------------------
 * open input file
 *-------------------------------------------------------------------------*/
 
  fid_src = h5tools_fopen(fname_src, H5F_ACC_RDONLY, H5P_DEFAULT, NULL, NULL, 0);

/*-------------------------------------------------------------------------
 * test for error in opening input file
 *-------------------------------------------------------------------------*/
 if (fid_src==-1)
 {
  error_msg(progname, "Could not open input file <%s>...Exiting\n", fname_src);
  if (fname_src)
   free(fname_src);
  leave(EXIT_FAILURE);
 }

/*-------------------------------------------------------------------------
 * open output file
 *-------------------------------------------------------------------------*/

    /* Attempt to open an existing HDF5 file first */
    fid_dst = h5tools_fopen(fname_dst, H5F_ACC_RDWR, H5P_DEFAULT, NULL, NULL, 0);

    /* If we couldn't open an existing file, try creating file */
    /* (use "EXCL" instead of "TRUNC", so we don't blow away existing non-HDF5 file) */
    if(fid_dst < 0)
        fid_dst = H5Fcreate(fname_dst, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

/*-------------------------------------------------------------------------
 * test for error in opening output file
 *-------------------------------------------------------------------------*/
 if (fid_dst==-1)
 {
  error_msg(progname, "Could not open output file <%s>...Exiting\n", fname_dst);
  if (fname_src)
   free(fname_src);
  if (fname_dst)
   free(fname_dst);
  leave(EXIT_FAILURE);
 }
  
/*-------------------------------------------------------------------------
 * print some info
 *-------------------------------------------------------------------------*/
 
 if (verbose)
 {
  printf("Copying file <%s> and object <%s> to file <%s> and object <%s>\n",
   fname_src,
   oname_src,
   fname_dst,
   oname_dst);
  if (flag)
   printf("Using %s flag\n", str_flag);
 }

 
/*-------------------------------------------------------------------------
 * create property lists for copy
 *-------------------------------------------------------------------------*/
 
 /* create property to pass copy options */
 if ( (ocpl_id = H5Pcreate(H5P_OBJECT_COPY)) < 0) 
  goto error;

 /* set options for object copy */
 if (flag)
 {
  if ( H5Pset_copy_object(ocpl_id, flag) < 0) 
   goto error;
 }

    /* Create link creation property list */
    if((lcpl_id = H5Pcreate(H5P_LINK_CREATE)) < 0) {
        error_msg(progname, "Could not create link creation property list\n");
        goto error;
    } /* end if */

    /* Check for creating intermediate groups */
    if(parents) {
        /* Set the intermediate group creation property */
        if(H5Pset_create_intermediate_group(lcpl_id, 1) < 0) {
            error_msg(progname, "Could not set property for creating parent groups\n");
            goto error;
        } /* end if */

        /* Display some output if requested */
        if(verbose)
            printf("%s: Creating parent groups\n", progname);
    } /* end if */

/*-------------------------------------------------------------------------
 * do the copy
 *-------------------------------------------------------------------------*/

  if (H5Ocopy(fid_src,          /* Source file or group identifier */
              oname_src,        /* Name of the source object to be copied */
              fid_dst,          /* Destination file or group identifier  */
              oname_dst,        /* Name of the destination object  */
              ocpl_id,          /* Object copy property list */
              lcpl_id)<0)       /* Link creation property list */
              goto error;
 
 /* close propertis */
 if(H5Pclose(ocpl_id)<0)
  goto error;
 if(H5Pclose(lcpl_id)<0)
  goto error;

 /* close files */
 if (H5Fclose(fid_src)<0)
  goto error;
 if (H5Fclose(fid_dst)<0)
  goto error;

 if (fname_src)
  free(fname_src);
 if (fname_dst)
  free(fname_dst);
 if (oname_dst)
  free(oname_dst);
 if (oname_src)
  free(oname_src);

 h5tools_close();
 
 return 0;

error:
 printf("Error in copy...Exiting\n");
 H5E_BEGIN_TRY {
  H5Pclose(ocpl_id);
  H5Pclose(lcpl_id);
  H5Fclose(fid_src);
  H5Fclose(fid_dst);
  } H5E_END_TRY;
 if (fname_src)
  free(fname_src);
 if (fname_dst)
  free(fname_dst);
 if (oname_dst)
  free(oname_dst);
 if (oname_src)
  free(oname_src);

 h5tools_close();
 
 return 1;
}

