#!/bin/bash
./xblat3z < ${srcdir}/zblat3.in
if test `grep -i fail zblat3.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
