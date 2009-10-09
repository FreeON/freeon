#!/bin/bash
./xblat2z < ${srcdir}/zblat2.in
if test `grep -i fail zblat2.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
