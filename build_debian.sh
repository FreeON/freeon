#!/bin/bash

make dist
mv freeon-1.0.3.tar.bz2 freeon_1.0.3.orig.tar.bz2
rm -rf freeon-1.0.3
tar xf freeon_1.0.3.orig.tar.bz2
rsync -av debian freeon-1.0.3/
cd freeon-1.0.3
debuild -us -uc
