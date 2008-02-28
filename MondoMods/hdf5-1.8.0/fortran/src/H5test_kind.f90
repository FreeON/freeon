! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!   Copyright by the Board of Trustees of the University of Illinois.         *
!   All rights reserved.                                                      *
!                                                                             *
!   This file is part of HDF5.  The full HDF5 copyright notice, including     *
!   terms governing use, modification, and redistribution, is contained in    *
!   the files COPYING and Copyright.html.  COPYING can be found at the root   *
!   of the source code distribution tree; Copyright.html can be found at the  *
!   root level of an installed copy of the electronic HDF5 document set and   *
!   is linked from the top-level documents page.  It can also be found at     *
!   http://hdf.ncsa.uiuc.edu/HDF5/doc/Copyright.html.  If you do not have     *
!   access to either file, you may request a copy from hdfhelp@ncsa.uiuc.edu. *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
! H5test_kind.f90
!
! This fortran program generates H5fortran_detect.f90
! 
!
         program test_kind
         integer :: i, j, ii, last, kind_numbers(10)
         integer :: jr, jd
         last = -1
         ii = 0
         j = selected_int_kind(18)
!         write(*,*) j
         do i = 1,100
            j = selected_int_kind(i)
            if(j .ne. last) then
              if(last .ne. -1) then
                  ii = ii + 1
                  kind_numbers(ii) = last
              endif
            last = j
            if(j .eq. -1) exit
            endif
          enddo
!          write(*,*) kind_numbers(1:ii)
! Generate a program
          write(*,*) "program int_kind"
          write(*,*) "write(*,*) "" /*generating header file*/ """
             j = 0
             write(*, "("" call i"", i2.2,""()"")") j
             jr = 0
             write(*, "("" call r"", i2.2,""()"")") jr
             jd = 0
             write(*, "("" call d"", i2.2,""()"")") jd
          do i = 1, ii
             j = kind_numbers(i)
             write(*, "("" call i"", i2.2,""()"")") j
          enddo
          write(*,*) "end program int_kind"
              j = 0
             write(*, "("" subroutine i"" i2.2,""()"")") j
             write(*,*)"   implicit none"
             write(*,*)"   integer :: a = 0"
             write(*,*)"   integer :: a_size"
             write(*,*)"   a_size = bit_size(a)"
             write(*,*)"   if (a_size .eq. 8) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_NATIVE_1"" "
             write(*,*)"   endif"
             write(*,*)"   if (a_size .eq. 16) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_NATIVE_2"" "
             write(*,*)"   endif"
             write(*,*)"   if (a_size .eq. 32) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_NATIVE_4"" "
             write(*,*)"   endif"
             write(*,*)"   if (a_size .eq. 64) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_NATIVE_8"" "
             write(*,*)"   endif"
             write(*,*)"   if (a_size .eq. 128) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_NATIVE_16"" "
             write(*,*)"   endif"
             write(*,*)"   return"
             write(*,*)"   end subroutine"    
              jr = 0
             write(*, "("" subroutine r"" i2.2,""()"")") j
             write(*,*)"   implicit none"
             write(*,*)"   real :: b(1) = 0"
             write(*,*)"   integer :: a(1) = 0"
             write(*,*)"   integer :: a_size"
             write(*,*)"   integer :: real_size"
             write(*,*)"   integer :: ab_size ! How many integers needed to hold a real"
             write(*,*)"   integer :: ba_size ! How many reals needed to hold an integer"
             write(*,*)"   a_size = bit_size(a(1)) ! Size in bits for integer"
             write(*,*)"   ab_size = size(transfer(b,a))"
             write(*,*)"   ba_size = size(transfer(a,b))"
             write(*,*)"   if (ab_size .eq. ba_size) real_size=a_size"
             write(*,*)"   if (ab_size .gt. ba_size) real_size=a_size*ba_size"
             write(*,*)"   if (ab_size .lt. ba_size) real_size=a_size/ba_size"
             write(*,*)"   if (real_size .eq. 32) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_REAL_NATIVE_4"" "
             write(*,*)"   endif"
             write(*,*)"   if (real_size .eq. 64) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_REAL_NATIVE_8"" "
             write(*,*)"   endif"
             write(*,*)"   if (real_size .eq. 128) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_REAL_NATIVE_16"" "
             write(*,*)"   endif"
             write(*,*)"   return"
             write(*,*)"   end subroutine"    
              jd = 0
             write(*, "("" subroutine d"" i2.2,""()"")") jd
             write(*,*)"   implicit none"
             write(*,*)"   double precision :: b = 0"
             write(*,*)"   integer :: a(8) = 0"
             write(*,*)"   integer :: a_size"
             write(*,*)"   integer :: b_size"
             write(*,*)"   a_size = bit_size(a(1))"
             write(*,*)"   b_size = size(transfer(b,a))*a_size"
             write(*,*)"   if (b_size .eq. 64) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_DOUBLE_NATIVE_8"" "
             write(*,*)"   endif"
             write(*,*)"   if (b_size .eq. 128) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_DOUBLE_NATIVE_16"" "
             write(*,*)"   endif"
             write(*,*)"   return"
             write(*,*)"   end subroutine"    
          do i = 1, ii
              j = kind_numbers(i)
             write(*, "("" subroutine i"" i2.2,""()"")") j
             write(*,*)"   implicit none"
             write(*,*)"   integer(",j,") :: a = 0"
             write(*,*)"   integer :: a_size"
             write(*,*)"   a_size = bit_size(a)"
             write(*,*)"   if (a_size .eq. 8) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_INTEGER_1"" "
             write(*,*)"   endif"
             write(*,*)"   if (a_size .eq. 16) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_INTEGER_2"" "
             write(*,*)"   endif"
             write(*,*)"   if (a_size .eq. 32) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_INTEGER_4"" "
             write(*,*)"   endif"
             write(*,*)"   if (a_size .eq. 64) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_INTEGER_8"" "
             write(*,*)"   endif"
             write(*,*)"   if (a_size .eq. 128) then"
             write(*,*)"       write(*,*) ""#define H5_FORTRAN_HAS_INTEGER_16"" "
             write(*,*)"   endif"
             write(*,*)"   return"
             write(*,*)"   end subroutine"    
          enddo
          end program
              
            

