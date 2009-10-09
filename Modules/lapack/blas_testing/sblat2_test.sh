#!/bin/bash
./xblat2s < ${srcdir}/sblat2.in
if test `grep -i fail sblat2.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
