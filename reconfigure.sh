#!/bin/bash
#
# Use this script to reconfigure the sources after checking them out from git.

autoreconf --no-recursive --install --symlink --warnings=no-portability . LinearAlgebraLibrary MondoMods/hdf5-1.8.3
