! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!   Copyright by The HDF Group.                                               *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the files COPYING and Copyright.html.  COPYING can be found at the root   *
!   of the source code distribution tree; Copyright.html can be found at the  *
!   root level of an installed copy of the electronic HDF5 document set and   *
!   is linked from the top-level documents page.  It can also be found at     *
!   http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have          *
!   access to either file, you may request a copy from help@hdfgroup.org.     *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
!
!
! 
!    This file contains subroutines which are used in
!    all the hdf5 fortran tests
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: check
!DEC$endif

         SUBROUTINE check(string,error,total_error)
            CHARACTER(LEN=*) :: string
            INTEGER :: error, total_error
            if (error .lt. 0) then
                total_error=total_error+1
                write(*,*) string, " failed"
            endif
            RETURN
         END SUBROUTINE check   


!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: verify
!DEC$endif
         SUBROUTINE verify(string,value,correct_value,total_error)
            CHARACTER(LEN=*) :: string
            INTEGER :: value, correct_value, total_error
            if (value .ne. correct_value) then
                total_error=total_error+1
                write(*,*) string
            endif
            RETURN
         END SUBROUTINE verify

!----------------------------------------------------------------------
! Name:		h5_fixname_f 
!
! Purpose:	Create a file name from the a file base name.
!               It is a fortran counterpart for the h5_fixname in ../../test/h5test.c
!
! Inputs:  
!		base_name	- base name of the file 
!               fapl		- file access property list 
! Outputs:  
!		full_name	- full file name
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
!
! Programmer:	Elena Pourmal
!		September 13, 2002
!
!
!----------------------------------------------------------------------
          SUBROUTINE h5_fixname_f(base_name, full_name, fapl, hdferr)
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5_fixname_f
!DEC$endif
            USE H5GLOBAL           
            IMPLICIT NONE
            CHARACTER(LEN=*), INTENT(IN) :: base_name   ! base name 
            CHARACTER(LEN=*), INTENT(IN) :: full_name   ! full name 
            INTEGER, INTENT(OUT) :: hdferr         ! Error code 
            INTEGER(HID_T), INTENT(IN) :: fapl ! file access property list

            INTEGER(SIZE_T) :: base_namelen ! Length of the base name character string
            INTEGER(SIZE_T) :: full_namelen ! Length of the full name character string
!            INTEGER(HID_T) :: fapl_default

            INTERFACE
              INTEGER FUNCTION h5_fixname_c(base_name, base_namelen, fapl, &
                               full_name, full_namelen)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !DEC$ ATTRIBUTES C,reference,decorate,alias:'H5_FIXNAME_C':: h5_fixname_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: base_name 
              !DEC$ATTRIBUTES reference :: full_name 
              CHARACTER(LEN=*), INTENT(IN) :: base_name
              INTEGER(SIZE_T) :: base_namelen
              INTEGER(HID_T), INTENT(IN) :: fapl
              CHARACTER(LEN=*), INTENT(IN) :: full_name
              INTEGER(SIZE_T) :: full_namelen
              END FUNCTION h5_fixname_c
            END INTERFACE

            base_namelen = LEN(base_name)
            full_namelen = LEN(full_name)
            hdferr = h5_fixname_c(base_name, base_namelen, fapl, &
                     full_name, full_namelen) 

          END SUBROUTINE h5_fixname_f
          
!----------------------------------------------------------------------
! Name:		h5_cleanup_f 
!
! Purpose:	Cleanups tests files
!               It is a fortran counterpart for the h5_cleanup in ../../test/h5test.c
!
! Inputs:  
!		base_name	- base name of the file 
!               fapl		- file access property list 
! Outputs:  
!		hdferr:		- error code		
!				 	Success:  0
!				 	Failure: -1   
!
! Programmer:	Elena Pourmal
!		September 19, 2002
!
!
!----------------------------------------------------------------------
          SUBROUTINE h5_cleanup_f(base_name, fapl, hdferr)
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5_cleanup_f
!DEC$endif
            USE H5GLOBAL           
            IMPLICIT NONE
            CHARACTER(LEN=*), INTENT(IN) :: base_name   ! base name 
            INTEGER, INTENT(OUT) :: hdferr         ! Error code 
            INTEGER(HID_T), INTENT(IN) :: fapl ! file access property list

            INTEGER(SIZE_T) :: base_namelen ! Length of the base name character string

            INTERFACE
              INTEGER FUNCTION h5_cleanup_c(base_name, base_namelen, fapl)
              USE H5GLOBAL
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !DEC$ ATTRIBUTES C,reference,decorate,alias:'H5_CLEANUP_C':: h5_cleanup_c
              !DEC$ ENDIF
              !DEC$ATTRIBUTES reference :: base_name 
              CHARACTER(LEN=*), INTENT(IN) :: base_name
              INTEGER(SIZE_T) :: base_namelen
              INTEGER(HID_T), INTENT(IN) :: fapl
              END FUNCTION h5_cleanup_c
            END INTERFACE

            base_namelen = LEN(base_name)
            hdferr = h5_cleanup_c(base_name, base_namelen, fapl)

          END SUBROUTINE h5_cleanup_f

!----------------------------------------------------------------------
! Name:		h5_exit_f 
!
! Purpose:	Exit application
!               It is a fortran counterpart for the standard C 'exit()' routine
! 		Be careful not to overflow the exit value range since
! 		UNIX supports a very small range such as 1 byte.
!		Therefore, exit(256) may end up as exit(0).
!
! Inputs:  
!		status	- Status to return from application
!
! Outputs:  
!		none
!
! Programmer:	Quincey Koziol
!		December 14, 2004
!
!
!----------------------------------------------------------------------
          SUBROUTINE h5_exit_f(status)
!
!This definition is needed for Windows DLLs
!DEC$if defined(BUILD_HDF5_DLL)
!DEC$attributes dllexport :: h5_exit_f
!DEC$endif
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: status         ! Return code

            INTERFACE
              SUBROUTINE h5_exit_c(status)
              !DEC$ IF DEFINED(HDF5F90_WINDOWS)
              !DEC$ ATTRIBUTES C,reference,decorate,alias:'H5_EXIT_C':: h5_exit_c
              !DEC$ ENDIF
              INTEGER, INTENT(IN) :: status
              END SUBROUTINE h5_exit_c
            END INTERFACE

            CALL h5_exit_c(status)

          END SUBROUTINE h5_exit_f

