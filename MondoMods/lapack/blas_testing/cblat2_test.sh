#!/bin/bash
./xblat2c < ${srcdir}/cblat2.in
if test `grep -i fail cblat2.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
