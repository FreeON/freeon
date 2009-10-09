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
$! Makefile for VMS systems.
$!
$ type sys$input
       Running examples
$ run  h5_write
$ run  h5_read
$ run  h5_extend_write
$ run  h5_chunk_read
$ run  h5_compound
$ run  h5_group
$ run  h5_select
$ run  h5_attribute
$ run  h5_mount
$ run  h5_reference
$ run  h5_ref2reg
$ run  h5_drivers
$ exit
