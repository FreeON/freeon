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
 * Programmer:  Quincey Koziol <koziol@ncsa.uiuc.ed>
 *              Thursday, July 11, 2002
 *
 * Purpose:	This is a "combination" MPI-2 and posix I/O driver.
 *              It uses MPI for coodinating the actions of several processes
 *              and posix I/O calls to do the actual I/O to the disk.
 *
 *              This driver was derived from the H5FDmpio.c driver and may
 *              share bugs/quirks/etc.
 *
 * Limitations:
 *              There is no "collective" I/O mode with this driver.
 *
 *              This will almost certainly _not_ work correctly for files
 *              accessed on distributed parallel systems with the file located
 *              on a non-parallel filesystem.
 *
 */
#include "H5private.h"		/*library functions			*/
#include "H5Eprivate.h"		/*error handling			*/
#include "H5Fprivate.h"		/*files					*/
#include "H5FDprivate.h"	/*file driver			        */
#include "H5FDmpiposix.h"       /* MPI/posix I/O file driver            */
#include "H5Iprivate.h"		/* IDs			  		*/
#include "H5MMprivate.h"        /*memory allocation                     */
#include "H5Pprivate.h"		/*property lists			*/

/* Features:
 *   USE_GPFS_HINTS -- issue gpfs_fcntl() calls to hopefully improve
 *                     performance when accessing files on a GPFS
 *                     file system.
 *
 *   REPORT_IO      -- if set then report all POSIX file calls to stderr.
 *
 */
#ifdef USE_GPFS_HINTS
#   include <gpfs_fcntl.h>
#endif

/*
 * The driver identification number, initialized at runtime if H5_HAVE_PARALLEL
 * is defined. This allows applications to still have the H5FD_MPIPOSIX
 * "constants" in their source code (it also makes this file strictly ANSI
 * compliant when H5_HAVE_PARALLEL isn't defined)
 */
static hid_t H5FD_MPIPOSIX_g = 0;

#ifdef H5_HAVE_PARALLEL

/* File operations */
#define OP_UNKNOWN	0
#define OP_READ		1
#define OP_WRITE	2

/*
 * The description of a file belonging to this driver.
 * The EOF value
 * is only used just after the file is opened in order for the library to
 * determine whether the file is empty, truncated, or okay. The MPIPOSIX driver
 * doesn't bother to keep it updated since it's an expensive operation.
 */
typedef struct H5FD_mpiposix_t {
    H5FD_t	pub;		/*public stuff, must be first		*/
    int		fd;		/*the unix file handle		        */
    MPI_Comm	comm;		/*communicator				*/
    int         mpi_rank;       /* This process's rank                  */
    int         mpi_size;       /* Total number of processes            */
    int         mpi_round;      /* Current round robin process (for metadata I/O) */
    haddr_t	eof;		/*end-of-file marker			*/
    haddr_t	eoa;		/*end-of-address marker			*/
    haddr_t	last_eoa;	/* Last known end-of-address marker	*/
    haddr_t	pos;		/* Current file I/O position	        */
    int		op;		/* Last file I/O operation		*/
    hsize_t	naccess;	/* Number of (write) accesses to file   */
    size_t      blksize;        /* Block size of file system            */
#ifndef WIN32
    /*
     * On most systems the combination of device and i-node number uniquely
     * identify a file.
     */
    dev_t	device;		/*file device number		*/
    ino_t	inode;		/*file i-node number		*/
#else
    /*
     * On WIN32 the low-order word of a unique identifier associated with the
     * file and the volume serial number uniquely identify a file. This number
     * (which, both? -rpm) may change when the system is restarted or when the
     * file is opened. After a process opens a file, the identifier is
     * constant until the file is closed. An application can use this
     * identifier and the volume serial number to determine whether two
     * handles refer to the same file.
     */
    int fileindexlo;
    int fileindexhi;
#endif
} H5FD_mpiposix_t;

/*
 * This driver supports systems that have the lseek64() function by defining
 * some macros here so we don't have to have conditional compilations later
 * throughout the code.
 *
 * file_offset_t:	The datatype for file offsets, the second argument of
 *			the lseek() or lseek64() call.
 *
 * file_seek:		The function which adjusts the current file position,
 *			either lseek() or lseek64().
 */
/* adding for windows NT file system support. */

#ifdef H5_HAVE_LSEEK64
#   define file_offset_t	off64_t
#   define file_seek		lseek64
#   define file_truncate	ftruncate64
#elif defined (WIN32) && !defined(__MWERKS__)
# /*MSVC*/
#   define file_offset_t __int64
#   define file_seek _lseeki64
#   define file_truncate	_ftruncatei64
#else
#   define file_offset_t	off_t
#   define file_seek		HDlseek
#   define file_truncate	HDftruncate
#endif

/*
 * These macros check for overflow of various quantities.  These macros
 * assume that file_offset_t is signed and haddr_t and size_t are unsigned.
 * 
 * ADDR_OVERFLOW:	Checks whether a file address of type `haddr_t'
 *			is too large to be represented by the second argument
 *			of the file seek function.
 *
 * SIZE_OVERFLOW:	Checks whether a buffer size of type `hsize_t' is too
 *			large to be represented by the `size_t' type.
 *
 * REGION_OVERFLOW:	Checks whether an address and size pair describe data
 *			which can be addressed entirely by the second
 *			argument of the file seek function.
 */
#define MAXADDR (((haddr_t)1<<(8*sizeof(file_offset_t)-1))-1)
#define ADDR_OVERFLOW(A)	(HADDR_UNDEF==(A) ||			      \
				 ((A) & ~(haddr_t)MAXADDR))
#define SIZE_OVERFLOW(Z)	((Z) & ~(hsize_t)MAXADDR)
#define REGION_OVERFLOW(A,Z)	(ADDR_OVERFLOW(A) || SIZE_OVERFLOW(Z) ||      \
				 sizeof(file_offset_t)<sizeof(size_t) ||      \
                                 HADDR_UNDEF==(A)+(Z) ||		      \
				 (file_offset_t)((A)+(Z))<(file_offset_t)(A))

/* Callbacks */
static void *H5FD_mpiposix_fapl_get(H5FD_t *_file);
static H5FD_t *H5FD_mpiposix_open(const char *name, unsigned flags, hid_t fapl_id,
			      haddr_t maxaddr);
static herr_t H5FD_mpiposix_close(H5FD_t *_file);
static int H5FD_mpiposix_cmp(const H5FD_t *_f1, const H5FD_t *_f2);
static herr_t H5FD_mpiposix_query(const H5FD_t *_f1, unsigned long *flags);
static haddr_t H5FD_mpiposix_get_eoa(H5FD_t *_file);
static herr_t H5FD_mpiposix_set_eoa(H5FD_t *_file, haddr_t addr);
static haddr_t H5FD_mpiposix_get_eof(H5FD_t *_file);
static herr_t H5FD_mpiposix_read(H5FD_t *_file, H5FD_mem_t type, hid_t fapl_id, haddr_t addr,
			     hsize_t size, void *buf);
static herr_t H5FD_mpiposix_write(H5FD_t *_file, H5FD_mem_t type, hid_t fapl_id, haddr_t addr,
			      hsize_t size, const void *buf);
static herr_t H5FD_mpiposix_flush(H5FD_t *_file, hid_t dxpl_id);

/* MPIPOSIX-specific file access properties */
typedef struct H5FD_mpiposix_fapl_t {
    MPI_Comm		comm;		/*communicator			*/
} H5FD_mpiposix_fapl_t;

/* The MPIPOSIX file driver information */
static const H5FD_class_t H5FD_mpiposix_g = {
    "mpiposix",					/*name			*/
    MAXADDR,					/*maxaddr		*/
    NULL,					/*sb_size		*/
    NULL,					/*sb_encode		*/
    NULL,					/*sb_decode		*/
    sizeof(H5FD_mpiposix_fapl_t),		/*fapl_size		*/
    H5FD_mpiposix_fapl_get,			/*fapl_get		*/
    NULL,					/*fapl_copy		*/
    NULL, 					/*fapl_free		*/
    0,						/*dxpl_size		*/
    NULL,					/*dxpl_copy		*/
    NULL,					/*dxpl_free		*/
    H5FD_mpiposix_open,				/*open			*/
    H5FD_mpiposix_close,			/*close			*/
    H5FD_mpiposix_cmp,			        /*cmp			*/
    H5FD_mpiposix_query,		        /*query			*/
    NULL,					/*alloc			*/
    NULL,					/*free			*/
    H5FD_mpiposix_get_eoa,			/*get_eoa		*/
    H5FD_mpiposix_set_eoa, 			/*set_eoa		*/
    H5FD_mpiposix_get_eof,			/*get_eof		*/
    H5FD_mpiposix_read,				/*read			*/
    H5FD_mpiposix_write,			/*write			*/
    H5FD_mpiposix_flush,			/*flush			*/
    H5FD_FLMAP_SINGLE,				/*fl_map		*/
};

/* Global var to allow elimination of redundant metadata writes
 * to be controlled by the value of an environment variable. */
/* Use the elimination by default unless this is the Intel Red machine */
#ifndef __PUMAGON__
hbool_t	H5_mpiposix_1_metawrite_g = TRUE;
#else
hbool_t	H5_mpiposix_1_metawrite_g = FALSE;
#endif

/* Interface initialization */
#define PABLO_MASK	H5FD_mpiposix_mask
#define INTERFACE_INIT	H5FD_mpiposix_init
static int interface_initialize_g = 0;


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_init
 *
 * Purpose:	Initialize this driver by registering the driver with the
 *		library.
 *
 * Return:	Success:	The driver ID for the mpiposix driver.
 *
 *		Failure:	Negative.
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
hid_t
H5FD_mpiposix_init(void)
{
    FUNC_ENTER(H5FD_mpiposix_init, FAIL);

    if (H5I_VFL!=H5Iget_type(H5FD_MPIPOSIX_g))
        H5FD_MPIPOSIX_g = H5FDregister(&H5FD_mpiposix_g);

    FUNC_LEAVE(H5FD_MPIPOSIX_g);
} /* end H5FD_mpiposix_init() */


/*-------------------------------------------------------------------------
 * Function:	H5Pset_fapl_mpiposix
 *
 * Purpose:	Store the user supplied MPI communicator COMM in
 *		the file access property list FAPL_ID which can then be used
 *		to create and/or open the file.  This function is available
 *		only in the parallel HDF5 library and is not collective.
 *
 *		COMM is the MPI communicator to be used for file open as
 *		defined in MPI_FILE_OPEN of MPI-2. This function does not
 *		make a duplicated communicator. Any modification to COMM
 *		after this function call returns may have undetermined effect
 *		on the access property list. Users should not modify the
 *		communicator while it is defined in a property list.
 *
 * Return:	Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:	Quincey Koziol
 *		Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pset_fapl_mpiposix(hid_t fapl_id, MPI_Comm comm)
{
    herr_t ret_value=FAIL;
    H5FD_mpiposix_fapl_t	fa;
    
    FUNC_ENTER(H5Pset_fapl_mpiposix, FAIL);
    H5TRACE2("e","iMc",fapl_id,comm);

    /* Check arguments */
    if (H5P_FILE_ACCESS!=H5Pget_class(fapl_id))
        HRETURN_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a fapl");
#ifdef LATER
#warning "We need to verify that COMM contains sensible information."
#endif

    /* Initialize driver specific properties */
    fa.comm = comm;

    ret_value= H5Pset_driver(fapl_id, H5FD_MPIPOSIX, &fa);

    FUNC_LEAVE(ret_value);
} /* end H5Pset_fapl_mpiposix() */


/*-------------------------------------------------------------------------
 * Function:	H5Pget_fapl_mpiposix
 *
 * Purpose:	If the file access property list is set to the H5FD_MPIPOSIX
 *		driver then this function returns the MPI communicator and
 *		information through the COMM pointer.
 *
 * Return:	Success:	Non-negative with the communicator and
 *				information returned through the COMM 
 *				argument if non-null. This piece of
 *				information is copied and is therefore
 *				valid only until the file access property
 *				list is modified or closed.
 *
 * 		Failure:	Negative
 *
 * Programmer:	Quincey Koziol
 *		Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
herr_t
H5Pget_fapl_mpiposix(hid_t fapl_id, MPI_Comm *comm/*out*/)
{
    H5FD_mpiposix_fapl_t	*fa;
    
    FUNC_ENTER(H5Pget_fapl_mpiposix, FAIL);
    H5TRACE2("e","ix",fapl_id,comm);

    if (H5P_FILE_ACCESS!=H5Pget_class(fapl_id))
        HRETURN_ERROR(H5E_PLIST, H5E_BADTYPE, FAIL, "not a fapl");
    if (H5FD_MPIPOSIX!=H5P_get_driver(fapl_id))
        HRETURN_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "incorrect VFL driver");
    if (NULL==(fa=H5Pget_driver_info(fapl_id)))
        HRETURN_ERROR(H5E_PLIST, H5E_BADVALUE, FAIL, "bad VFL driver info");

    if (comm) *comm = fa->comm;

    FUNC_LEAVE(SUCCEED);
} /* end H5Pget_fapl_mpiposix() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_communicator
 *
 * Purpose:	Returns the MPI communicator for the file.
 *
 * Return:	Success:	The communicator
 *
 *		Failure:	NULL
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
MPI_Comm
H5FD_mpiposix_communicator(H5FD_t *_file)
{
    H5FD_mpiposix_t *file = (H5FD_mpiposix_t*)_file;

    FUNC_ENTER(H5FD_mpiposix_communicator, MPI_COMM_NULL);
    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);

    FUNC_LEAVE(file->comm);
} /* end H5FD_mpi_posix_communicator() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_mpi_rank
 *
 * Purpose:	Returns the MPI rank for a process
 *
 * Return:	Success: non-negative
 *		Failure: negative
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
H5FD_mpiposix_mpi_rank(H5FD_t *_file)
{
    H5FD_mpiposix_t *file = (H5FD_mpiposix_t*)_file;

    FUNC_ENTER(H5FD_mpiposix_mpi_rank, FAIL);
    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);

    FUNC_LEAVE(file->mpi_rank);
} /* end H5FD_mpiposix_mpi_rank() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_mpi_size
 *
 * Purpose:	Returns the number of MPI processes
 *
 * Return:	Success: non-negative
 *		Failure: negative
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
int
H5FD_mpiposix_mpi_size(H5FD_t *_file)
{
    H5FD_mpiposix_t *file = (H5FD_mpiposix_t*)_file;

    FUNC_ENTER(H5FD_mpiposix_mpi_rank, FAIL);
    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);

    FUNC_LEAVE(file->mpi_size);
} /* end H5FD_mpiposix_mpi_size() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_fapl_get
 *
 * Purpose:	Returns a file access property list which could be used to
 *		create another file the same as this one.
 *
 * Return:	Success:	Ptr to new file access property list with all
 *				fields copied from the file pointer.
 *
 *		Failure:	NULL
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static void *
H5FD_mpiposix_fapl_get(H5FD_t *_file)
{
    H5FD_mpiposix_t	*file = (H5FD_mpiposix_t*)_file;
    H5FD_mpiposix_fapl_t *fa = NULL;

    FUNC_ENTER(H5FD_mpiposix_fapl_get, NULL);
    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);

    if (NULL==(fa=H5MM_calloc(sizeof(H5FD_mpiposix_fapl_t))))
        HRETURN_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");

    /* These should be copied. --QAK, 2002-07-11 */
    fa->comm = file->comm;

    FUNC_LEAVE(fa);
} /* end H5FD_mpiposix_fapl_get() */


/*-------------------------------------------------------------------------
 * Function:    H5FD_mpiposix_open
 *
 * Purpose:     Opens a file with name NAME.  The FLAGS are a bit field with
 *		purpose similar to the second argument of open(2) and which
 *		are defined in H5Fpublic.h. The file access property list
 *		FAPL_ID contains the properties driver properties and MAXADDR
 *		is the largest address which this file will be expected to
 *		access.  This is collective.
 *
 * Return:      Success:        A new file pointer.
 *              Failure:        NULL
 *
 * Programmer:  Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static H5FD_t *
H5FD_mpiposix_open(const char *name, unsigned flags, hid_t fapl_id,
	       haddr_t maxaddr)
{
    H5FD_mpiposix_t		*file=NULL;     /* New MPIPOSIX file struct */
    int                         o_flags;        /* Flags for file open call */
    int			        fd=(-1);        /* File handle for file opened */
    int				mpi_rank;       /* MPI rank of this process */
    int				mpi_size;       /* Total number of MPI processes */
    int				mpi_code;	/* mpi return code */
    const H5FD_mpiposix_fapl_t	*fa=NULL;       /* MPIPOSIX file access property list information */
    H5FD_mpiposix_fapl_t	_fa;            /* Private copy of default file access property list information */
    h5_stat_t                   sb;             /* Portable 'stat' struct */
#ifdef WIN32
    HFILE filehandle;
    struct _BY_HANDLE_FILE_INFORMATION fileinfo;
    int results;   
#endif
    H5FD_t                     *ret_value=NULL; /* Return value */

    FUNC_ENTER(H5FD_mpiposix_open, NULL);

    /* Check arguments */
    if (!name || !*name)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, NULL, "invalid file name");
    if (0==maxaddr || HADDR_UNDEF==maxaddr)
        HGOTO_ERROR(H5E_ARGS, H5E_BADRANGE, NULL, "bogus maxaddr");
    if (ADDR_OVERFLOW(maxaddr))
        HGOTO_ERROR(H5E_ARGS, H5E_OVERFLOW, NULL, "bogus maxaddr");

    /* Obtain a pointer to mpiposix-specific file access properties */
    if (H5P_DEFAULT==fapl_id || H5FD_MPIPOSIX!=H5P_get_driver(fapl_id)) {
	_fa.comm = MPI_COMM_SELF; /*default*/
	fa = &_fa;
    } /* end if */
    else {
	fa = H5Pget_driver_info(fapl_id);
	assert(fa);
    } /* end else */

    /* Get the MPI rank of this process and the total number of processes */
    if (MPI_SUCCESS != (mpi_code=MPI_Comm_rank (fa->comm, &mpi_rank)))
        HMPI_GOTO_ERROR(NULL, "MPI_Comm_rank failed", mpi_code);
    if (MPI_SUCCESS != (mpi_code=MPI_Comm_size (fa->comm, &mpi_size)))
        HMPI_GOTO_ERROR(NULL, "MPI_Comm_size failed", mpi_code);

    /* Build the open flags */
    o_flags = (H5F_ACC_RDWR & flags) ? O_RDWR : O_RDONLY;

    /* Only set the creation flag(s) for process 0 */
    if(mpi_rank==0) {
        if (H5F_ACC_TRUNC & flags)
            o_flags |= O_TRUNC;
        if (H5F_ACC_CREAT & flags)
            o_flags |= O_CREAT;
        if (H5F_ACC_EXCL & flags)
            o_flags |= O_EXCL;
    } /* end if */

    /* Process 0 opens (or creates) the file while the rest of the
     * processes wait.  Then process 0 signals the other processes and they
     * open (never create) the file and all processes proceed.
     */
    /* Process 0 opens (or creates) file and broadcasts result to other processes */
    if(mpi_rank==0) {
        /* Open the file */
        fd=HDopen(name, o_flags, 0666);
    } /* end if */

    /* Broadcast the results of the open() from process 0 */
    /* This is necessary because of the "tentative open" code in H5F_open()
     * where the file is attempted to be opened with different flags from the
     * user's, in order to check for the file's existence, etc.  Here, process 0
     * gets different flags from the other processes (since it is in charge of
     * creating the file, if necessary) and can fail in situations where the
     * other process's file opens would succeed, so allow the other processes
     * to check for that situation and bail out now also. - QAK
     */
    if (MPI_SUCCESS != (mpi_code= MPI_Bcast(&fd, sizeof(int), MPI_BYTE, 0, fa->comm)))
        HMPI_GOTO_ERROR(NULL, "MPI_Bcast failed", mpi_code);

    /* If the file open on process 0 failed, bail out on all processes now */
    if(fd<0)
        HGOTO_ERROR(H5E_FILE, H5E_CANTOPENFILE, NULL, "unable to open file");

    /* Other processes (non 0) wait for broadcast result from process 0 and then open file */
    if(mpi_rank!=0) {
        /* Open the file */
        if ((fd=HDopen(name, o_flags, 0666))<0)
            HGOTO_ERROR(H5E_FILE, H5E_CANTOPENFILE, NULL, "unable to open file");
    } /* end if */

    /* Process 0 fstat()s the file and broadcasts the results to the other processes */
    if(mpi_rank==0) {
        /* Get the stat information */
        if (HDfstat(fd, &sb)<0)
            HGOTO_ERROR(H5E_FILE, H5E_BADFILE, NULL, "unable to fstat file");
    } /* end if */

    /* Broadcast the results of the fstat() from process 0 */
    if (MPI_SUCCESS != (mpi_code= MPI_Bcast(&sb, sizeof(h5_stat_t), MPI_BYTE, 0, fa->comm)))
        HMPI_GOTO_ERROR(NULL, "MPI_Bcast failed", mpi_code);

#ifdef USE_GPFS_HINTS
    {
        /* Free all byte range tokens. This is a good thing to do if raw data is aligned on 256kB boundaries (a GPFS page is
         * 256kB). Care should be taken that there aren't too many sub-page writes, or the mmfsd may become overwhelmed.  This
         * should probably eventually be passed down here as a property. The gpfs_fcntl() will most likely fail if `fd' isn't
         * on a GPFS file system. */
        struct {
            gpfsFcntlHeader_t   hdr;
            gpfsFreeRange_t     fr;
        } hint;
        memset(&hint, 0, sizeof hint);
        hint.hdr.totalLength = sizeof hint;
        hint.hdr.fcntlVersion = GPFS_FCNTL_CURRENT_VERSION;
        hint.fr.structLen = sizeof hint.fr;
        hint.fr.structType = GPFS_FREE_RANGE;
        hint.fr.start = 0;
        hint.fr.length = 0;
        
        if (gpfs_fcntl(fd, &hint)<0)
            HGOTO_ERROR(H5E_FILE, H5E_FCNTL, NULL, "failed to send hints to GPFS");

        if (0==mpi_rank)
            fprintf(stderr, "HDF5: using GPFS hint mechanism...\n");
    }
#endif

    /* Build the file struct and initialize it */
    if (NULL==(file=H5MM_calloc(sizeof(H5FD_mpiposix_t))))
        HGOTO_ERROR(H5E_RESOURCE, H5E_NOSPACE, NULL, "memory allocation failed");

#ifdef REPORT_IO
    fprintf(stderr, "open:  rank=%d name=%s file=0x%08lx\n", mpi_rank, name, (unsigned long)file);
#endif

    /* Set the general file information */
    file->fd = fd;
    file->eof = sb.st_size;
    file->blksize = sb.st_blksize;

    /* Set the MPI information */
    file->comm = fa->comm;
    file->mpi_rank = mpi_rank;
    file->mpi_size = mpi_size;
    file->mpi_round = 0;        /* Start metadata writes with process 0 */

    /* Reset the last file I/O operation */
    file->pos = HADDR_UNDEF;
    file->op = OP_UNKNOWN;

    /* Set the information for the file's device and inode */
#ifdef WIN32
    filehandle = _get_osfhandle(fd);
    results = GetFileInformationByHandle((HANDLE)filehandle, &fileinfo);
    file->fileindexhi = fileinfo.nFileIndexHigh;
    file->fileindexlo = fileinfo.nFileIndexLow;
#else
    file->device = sb.st_dev;
    file->inode = sb.st_ino;
#endif

    /* Indicate success */
    ret_value=(H5FD_t *)file;

done:
    /* Error cleanup */
    if(ret_value==NULL) {
        /* Close the file if it was left open */
        if(fd!=(-1))
            HDclose(fd);
    } /* end if */

    FUNC_LEAVE(ret_value);
} /* end H5FD_mpiposix_open() */


/*-------------------------------------------------------------------------
 * Function:    H5FD_mpiposix_close
 *
 * Purpose:     Closes a file.
 *
 * Return:      Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:  Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_mpiposix_close(H5FD_t *_file)
{
    H5FD_mpiposix_t	*file = (H5FD_mpiposix_t*)_file;

    FUNC_ENTER(H5FD_mpiposix_close, FAIL);

    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);

    /* Close the unix file */
    if (HDclose(file->fd)<0)
        HRETURN_ERROR(H5E_IO, H5E_CANTCLOSEFILE, FAIL, "unable to close file");

    H5MM_xfree(file);

    FUNC_LEAVE(SUCCEED);
} /* end H5FD_mpiposix_close() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_cmp
 *
 * Purpose:	Compares two files belonging to this driver using an
 *		arbitrary (but consistent) ordering.
 *
 * Return:	Success:	A value like strcmp()
 *		Failure:	never fails (arguments were checked by the
 *				caller).
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static int
H5FD_mpiposix_cmp(const H5FD_t *_f1, const H5FD_t *_f2)
{
    const H5FD_mpiposix_t	*f1 = (const H5FD_mpiposix_t*)_f1;
    const H5FD_mpiposix_t	*f2 = (const H5FD_mpiposix_t*)_f2;
    int ret_value=0;

    FUNC_ENTER(H5FD_mpiposix_cmp, H5FD_VFD_DEFAULT);

#ifdef WIN32
    if (f1->fileindexhi < f2->fileindexhi) ret_value= -1;
    if (f1->fileindexhi > f2->fileindexhi) ret_value= 1;

    if (f1->fileindexlo < f2->fileindexlo) ret_value= -1;
    if (f1->fileindexlo > f2->fileindexlo) ret_value= 1;

#else
    if (f1->device < f2->device) ret_value= -1;
    if (f1->device > f2->device) ret_value= 1;

    if (f1->inode < f2->inode) ret_value= -1;
    if (f1->inode > f2->inode) ret_value= 1;
#endif

    FUNC_LEAVE(ret_value);
} /* end H5FD_mpiposix_cmp() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_query
 *
 * Purpose:	Set the flags that this VFL driver is capable of supporting.
 *              (listed in H5FDpublic.h)
 *
 * Return:	Success:	non-negative
 *		Failure:	negative
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_mpiposix_query(const H5FD_t UNUSED *_file, unsigned long *flags /* out */)
{
    herr_t ret_value=SUCCEED;

    FUNC_ENTER(H5FD_mpiposix_query, FAIL);

    /* Set the VFL feature flags that this driver supports */
    if(flags) {
        *flags=0;
        *flags|=H5FD_FEAT_AGGREGATE_METADATA; /* OK to aggregate metadata allocations */

        /* Distinguish between updating the metadata accumulator on writes and
         * reads.  This is particularly (perhaps only, even) important for MPI-I/O
         * where we guarantee that writes are collective, but reads may not be.
         * If we were to allow the metadata accumulator to be written during a
         * read operation, the application would hang.
         */
        *flags|=H5FD_FEAT_ACCUMULATE_METADATA_WRITE; /* OK to accumulate metadata for faster writes */

        *flags|=H5FD_FEAT_AGGREGATE_SMALLDATA; /* OK to aggregate "small" raw data allocations */
    } /* end if */

    FUNC_LEAVE(ret_value);
} /* end H5FD_mpiposix_query() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_get_eoa
 *
 * Purpose:	Gets the end-of-address marker for the file. The EOA marker
 *		is the first address past the last byte allocated in the
 *		format address space.
 *
 * Return:	Success:	The end-of-address marker.
 *		Failure:	HADDR_UNDEF
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static haddr_t
H5FD_mpiposix_get_eoa(H5FD_t *_file)
{
    H5FD_mpiposix_t *file = (H5FD_mpiposix_t*)_file;

    FUNC_ENTER(H5FD_mpiposix_get_eoa, HADDR_UNDEF);
    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);

    FUNC_LEAVE(file->eoa);
} /* end H5FD_mpiposix_get_eoa() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_set_eoa
 *
 * Purpose:	Set the end-of-address marker for the file. This function is
 *		called shortly after an existing HDF5 file is opened in order
 *		to tell the driver where the end of the HDF5 data is located.
 *
 * Return:	Success:	non-negative
 *		Failure:	negative
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_mpiposix_set_eoa(H5FD_t *_file, haddr_t addr)
{
    H5FD_mpiposix_t	*file = (H5FD_mpiposix_t*)_file;

    FUNC_ENTER(H5FD_mpiposix_set_eoa, FAIL);
    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);

    file->eoa = addr;

    FUNC_LEAVE(SUCCEED);
} /* end H5FD_mpi_posix_set_eoa() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_get_eof
 *
 * Purpose:	Gets the end-of-file marker for the file. The EOF marker
 *		is the real size of the file.
 *
 *		The MPIPOSIX driver doesn't bother keeping this field updated
 *		since that's a relatively expensive operation. Fortunately
 *		the library only needs the EOF just after the file is opened
 *		in order to determine whether the file is empty, truncated,
 *		or okay.
 *
 * Return:	Success:	The end-of-address marker.
 *		Failure:	HADDR_UNDEF
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static haddr_t
H5FD_mpiposix_get_eof(H5FD_t *_file)
{
    H5FD_mpiposix_t	*file = (H5FD_mpiposix_t*)_file;

    FUNC_ENTER(H5FD_mpiposix_get_eof, HADDR_UNDEF);
    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);

    FUNC_LEAVE(MAX(file->eof,file->eoa));
} /* end H5FD_mpiposix_get_eof() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_read
 *
 * Purpose:	Reads SIZE bytes of data from FILE beginning at address ADDR
 *		into buffer BUF according to data transfer properties in
 *		DXPL_ID using potentially complex file and buffer types to
 *		effect the transfer.
 *
 *		Reading past the end of the file returns zeros instead of
 *		failing.
 *
 * Return:	Success:	Non-negative. Result is stored in caller-supplied
 *				buffer BUF.
 *		Failure:	Negative, Contents of buffer BUF are undefined.
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_mpiposix_read(H5FD_t *_file, H5FD_mem_t UNUSED type, hid_t UNUSED dxpl_id, haddr_t addr, hsize_t size,
	       void *buf/*out*/)
{
    H5FD_mpiposix_t	*file = (H5FD_mpiposix_t*)_file;
    ssize_t	        nbytes;         /* Number of bytes read each I/O call */
    herr_t             	ret_value=SUCCEED;

    FUNC_ENTER(H5FD_mpiposix_read, FAIL);

    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);
    assert(buf);

    /* Check for overflow conditions */
    if (HADDR_UNDEF==addr)
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "addr undefined");
    if (REGION_OVERFLOW(addr, size))
        HGOTO_ERROR(H5E_ARGS, H5E_OVERFLOW, FAIL, "addr overflow");
    if (addr+size>file->eoa)
        HGOTO_ERROR(H5E_ARGS, H5E_OVERFLOW, FAIL, "addr overflow");

#ifdef REPORT_IO
    {
        int commrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &commrank);
        fprintf(stderr, "read:  rank=%d file=0x%08lx type=%d, addr=%lu size=%lu\n",
                commrank, (unsigned long)file, (int)type, (unsigned long)addr, (unsigned long)size);
    }
#endif

    /* Seek to the correct location */
    if ((addr!=file->pos || OP_READ!=file->op) &&
            file_seek(file->fd, (file_offset_t)addr, SEEK_SET)<0)
        HGOTO_ERROR(H5E_IO, H5E_SEEKERROR, FAIL, "unable to seek to proper position");

    /*
     * Read data, being careful of interrupted system calls, partial results,
     * and the end of the file.
     */
    while (size>0) {
        do {
            H5_CHECK_OVERFLOW(size,hsize_t,size_t);
            nbytes = HDread(file->fd, buf, (size_t)size);
        } while (-1==nbytes && EINTR==errno);
        if (-1==nbytes)
            HGOTO_ERROR(H5E_IO, H5E_READERROR, FAIL, "file read failed");
        if (0==nbytes) {
            /* end of file but not end of format address space */
            H5_CHECK_OVERFLOW(size,hsize_t,size_t);
            HDmemset(buf, 0, (size_t)size);
            size = 0;
        } /* end if */
        assert(nbytes>=0);
        assert((hsize_t)nbytes<=size);
        size -= (hsize_t)nbytes;
        addr += (haddr_t)nbytes;
        buf = (char*)buf + nbytes;
    }
    
    /* Update current position */
    file->pos = addr;
    file->op = OP_READ;

done:
    /* Check for error */
    if(ret_value<0) {
        /* Reset last file I/O information */
        file->pos = HADDR_UNDEF;
        file->op = OP_UNKNOWN;
    } /* end if */

    FUNC_LEAVE(ret_value);
} /* end H5FD_mpiposix_read() */


/*-------------------------------------------------------------------------
 * Function:	H5FD_mpiposix_write
 *
 * Purpose:	Writes SIZE bytes of data to FILE beginning at address ADDR
 *		from buffer BUF according to data transfer properties in
 *		DXPL_ID using potentially complex file and buffer types to
 *		effect the transfer.
 *
 * Return:	Success:	non-negative
 *		Failure:	negative
 *
 * Programmer:	Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *              Quincey Koziol - 2002/07/18
 *              Added "block_before_meta_write" dataset transfer flag, which
 *              is set during writes from a metadata cache flush and indicates
 *              that all the processes must sync up before (one of them)
 *              writing metadata.
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_mpiposix_write(H5FD_t *_file, H5FD_mem_t type, hid_t dxpl_id, haddr_t addr,
		hsize_t size, const void *buf)
{
    H5FD_mpiposix_t	*file = (H5FD_mpiposix_t*)_file;
    const H5D_xfer_t	*xfer_parms = NULL;     /* Dataset transfer plist */
    int			mpi_code;	/* MPI return code */
    ssize_t	        nbytes;         /* Number of bytes written each I/O call */
    herr_t             	ret_value=SUCCEED;      /* Return value */

    FUNC_ENTER(H5FD_mpiposix_write, FAIL);

    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);
    assert(buf);

    /* Check for overflow conditions */
    if (HADDR_UNDEF==addr) 
        HGOTO_ERROR(H5E_ARGS, H5E_BADVALUE, FAIL, "addr undefined");
    if (REGION_OVERFLOW(addr, size))
        HGOTO_ERROR(H5E_ARGS, H5E_OVERFLOW, FAIL, "addr overflow");
    if (addr+size>file->eoa)
        HGOTO_ERROR(H5E_ARGS, H5E_OVERFLOW, FAIL, "addr overflow");
    
    /* Get the dataset transfer property list */
    if (H5P_DEFAULT == dxpl_id) {
        xfer_parms = &H5D_xfer_dflt;
    } /* end if */
    else if (H5P_DATASET_XFER != H5P_get_class(dxpl_id) ||
	       NULL == (xfer_parms = H5I_object(dxpl_id))) {
        HGOTO_ERROR(H5E_ARGS, H5E_BADTYPE, FAIL, "not xfer parms");
    } /* end if */

    /* Metadata specific actions */
    if(type!=H5FD_MEM_DRAW) {
        /* Check if we need to syncronize all processes before attempting metadata write
         * (Prevents race condition where the process writing the metadata goes ahead
         * and writes the metadata to the file before all the processes have
         * read the data, "transmitting" data from the "future" to the reading
         * process. -QAK )
         */
        if(xfer_parms->block_before_meta_write) {
            if (MPI_SUCCESS!= (mpi_code=MPI_Barrier(file->comm)))
                HMPI_GOTO_ERROR(FAIL, "MPI_Barrier failed", mpi_code);
        } /* end if */

#if 1
        /* Only p<round> will do the actual write if all procs in comm write same data */
        if (H5_mpiposix_1_metawrite_g) {
            if (file->mpi_rank != file->mpi_round)
                HGOTO_DONE(SUCCEED) /* skip the actual write */
        } /* end if */
#else
        /* Only task zero will do the actual write if all tasks in comm write same data */
        if (H5_mpiposix_1_metawrite_g && 0!=file->mpi_rank)
            HGOTO_DONE(SUCCEED);
#endif
    } /* end if */

#ifdef REPORT_IO
    {
        int commrank;
        MPI_Comm_rank(MPI_COMM_WORLD, &commrank);
        fprintf(stderr, "write: rank=%d file=0x%08lx type=%d, addr=%lu size=%lu %s\n",
                commrank, (unsigned long)file, (int)type, (unsigned long)addr, (unsigned long)size,
                0==file->naccess?"(FIRST ACCESS)":"");
    }
#endif

    if (0==file->naccess++) {
        /* First write access to this file */
#ifdef USE_GPFS_HINTS
        struct {
            gpfsFcntlHeader_t           hdr;
            gpfsMultipleAccessRange_t   mar;
        } hint;
        memset(&hint, 0, sizeof hint);
        hint.hdr.totalLength = sizeof hint;
        hint.hdr.fcntlVersion = GPFS_FCNTL_CURRENT_VERSION;
        hint.mar.structLen = sizeof hint.mar;
        hint.mar.structType = GPFS_MULTIPLE_ACCESS_RANGE;
        hint.mar.accRangeCnt = 1;
        hint.mar.accRangeArray[0].blockNumber = addr / file->blksize;
        hint.mar.accRangeArray[0].start = addr % file->blksize;
        hint.mar.accRangeArray[0].length = MIN(file->blksize-hint.mar.accRangeArray[0].start, size);
        hint.mar.accRangeArray[0].isWrite = 1;
        if (gpfs_fcntl(file->fd, &hint)<0)
            HGOTO_ERROR(H5E_FILE, H5E_FCNTL, NULL, "failed to send hints to GPFS");
#endif
    }
        
    /* Seek to the correct location */
    if ((addr!=file->pos || OP_WRITE!=file->op) &&
            file_seek(file->fd, (file_offset_t)addr, SEEK_SET)<0)
        HGOTO_ERROR(H5E_IO, H5E_SEEKERROR, FAIL,
		      "unable to seek to proper position");

    /*
     * Write the data, being careful of interrupted system calls and partial
     * results
     */
    while (size>0) {
        do {
            H5_CHECK_OVERFLOW(size,hsize_t,size_t);
            nbytes = HDwrite(file->fd, buf, (size_t)size);
        } while (-1==nbytes && EINTR==errno);
        if (-1==nbytes)
            HGOTO_ERROR(H5E_IO, H5E_WRITEERROR, FAIL, "file write failed");
        assert(nbytes>0);
        assert((hsize_t)nbytes<=size);
        size -= (hsize_t)nbytes;
        addr += (haddr_t)nbytes;
        buf = (const char*)buf + nbytes;
    } /* end while */

    /* Update current last file I/O information */
    file->pos = addr;
    file->op = OP_WRITE;

done:
    /* Check for error */
    if(ret_value<0) {
        /* Reset last file I/O information */
        file->pos = HADDR_UNDEF;
        file->op = OP_UNKNOWN;
    } /* end if */
    /* Guard against getting into metadata broadcast in failure cases */
    else {
        /* if only p<round> writes, need to broadcast the ret_value to other processes */
        if ((type!=H5FD_MEM_DRAW) && H5_mpiposix_1_metawrite_g) {
            if (MPI_SUCCESS != (mpi_code= MPI_Bcast(&ret_value, sizeof(ret_value), MPI_BYTE, file->mpi_round, file->comm)))
                HMPI_GOTO_ERROR(FAIL, "MPI_Bcast failed", mpi_code);

            /* Round-robin rotate to the next process */
            file->mpi_round = (file->mpi_round+1)%file->mpi_size;
        } /* end if */
    } /* end else */

    FUNC_LEAVE(ret_value);
} /* end H5FD_mpiposix_write() */


/*-------------------------------------------------------------------------
 * Function:    H5FD_mpiposix_flush
 *
 * Purpose:     Makes sure that all data is on disk.  This is collective.
 *
 * Return:      Success:	Non-negative
 * 		Failure:	Negative
 *
 * Programmer:  Quincey Koziol
 *              Thursday, July 11, 2002
 *
 * Modifications:
 *
 *-------------------------------------------------------------------------
 */
static herr_t
H5FD_mpiposix_flush(H5FD_t *_file, hid_t UNUSED dxpl_id)
{
    H5FD_mpiposix_t	*file = (H5FD_mpiposix_t*)_file;
#ifdef WIN32
    HFILE filehandle;   /* Windows file handle */
    LARGE_INTEGER li;   /* 64-bit integer for SetFilePointer() call */
#endif /* WIN32 */
    int			mpi_code;	/* MPI return code */
    herr_t              ret_value=SUCCEED;

    FUNC_ENTER(H5FD_mpiposix_flush, FAIL);

    assert(file);
    assert(H5FD_MPIPOSIX==file->pub.driver_id);

    /* Extend the file to make sure it's large enough */
    if(file->eoa>file->last_eoa) {
        /* Use the round-robin process to truncate (extend) the file */
        if(file->mpi_rank == file->mpi_round) {
#ifdef WIN32
            /* Map the posix file handle to a Windows file handle */
            filehandle = _get_osfhandle(fd);

            /* Translate 64-bit integers into form Windows wants */
            /* [This algorithm is from the Windows documentation for SetFilePointer()] */
            li.QuadPart = file->eoa;
            SetFilePointer((HANDLE)filehandle,li.LowPart,&li.HighPart,FILE_BEGIN);
            if(SetEndOfFile((HANDLE)filehandle)==0)
                HGOTO_ERROR(H5E_IO, H5E_SEEKERROR, FAIL, "unable to extend file properly");
#else /* WIN32 */
            if(-1==file_truncate(file->fd, (file_offset_t)file->eoa))
                HGOTO_ERROR(H5E_IO, H5E_SEEKERROR, FAIL, "unable to extend file properly");
#endif /* WIN32 */
        } /* end if */

        /* Don't let any proc return until all have extended the file.
         * (Prevents race condition where some processes go ahead and write
         * more data to the file before all the processes have finished making
         * it the shorter length, potentially truncating the file and dropping
         * the new data written)
         */
        if (MPI_SUCCESS!= (mpi_code=MPI_Barrier(file->comm)))
            HMPI_GOTO_ERROR(FAIL, "MPI_Barrier failed", mpi_code);

        /* Update the 'last' eoa and eof values */
        file->last_eoa=file->eoa;
        file->eof = file->eoa;

        /* Reset last file I/O information */
        file->pos = HADDR_UNDEF;
        file->op = OP_UNKNOWN;
    } /* end if */

done:
    FUNC_LEAVE(ret_value);
} /* end H5FD_mpiposix_flush() */

#endif /*H5_HAVE_PARALLEL*/

