$!#
$!# Copyright by The HDF Group.
$!# Copyright by the Board of Trustees of the University of Illinois.
$!# All rights reserved.
$!#
$!# This file is part of HDF5.  The full HDF5 copyright notice, including
$!# terms governing use, modification, and redistribution, is contained in
$!# the files COPYING and Copyright.html.  COPYING can be found at the root
$!# of the source code distribution tree; Copyright.html can be found at the
$!# root level of an installed copy of the electronic HDF5 document set and
$!# is linked from the top-level documents page.  It can also be found at
$!# http://hdfgroup.org/HDF5/doc/Copyright.html.  If you do not have
$!# access to either file, you may request a copy from help@hdfgroup.org.
$!#
$!
$!
$! This file builds C, Frtran, C++ HDF5 librraies and runs the tests
$! Specify location of the top HDF5 source directory
$
$ hdf5top == "sys$sysusers:[pourmale.hdf5]"
$ len = F$LENGTH(hdf5top)
$ tmp = F$EXTRACT(0, len-1, hdf5top)
$ hdf5vms     = tmp + ".VMS]"
$ hdf5ctest   = tmp + ".TEST]"
$ hdf5f90test = tmp + ".FORTRAN.TEST]"
$ hdf5cxxtest = tmp  + ".C__.TEST]"
$ hdf5toolstest = tmp  + ".TOOLS.TESTFILES]"
$ h5importtest  = tmp  + ".TOOLS.H5IMPORT.TESTFILES]"
$ set def 'hdf5vms'
$@make
$ set def 'hdf5ctest'
$@check
$ set def 'hdf5f90test'
$@check
$ set def 'hdf5cxxtest'
$@check
$ set def 'hdf5toolstest'
$ copy [-.h5dump]check_h5dump.com     check_h5dump.com 
$ copy [-.h5ls]check_h5ls.com         check_h5ls.com
$ copy [-.h5diff]check_h5diff.com     check_h5diff.com
$ copy [-.h5repack]check_h5repack.com check_h5repack.com
$@check_h5dump.com
$@check_h5ls.com
$@check_h5diff.com
$@check_h5repack.com
$!
$ set def 'h5importtest'
$ copy [.-]check_h5import.com check_h5import.com
$@check_h5import.com
$
$ exit
