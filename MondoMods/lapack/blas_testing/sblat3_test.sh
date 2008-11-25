#!/bin/bash
./xblat3s < ${srcdir}/sblat3.in
if test `grep -i fail sblat3.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
