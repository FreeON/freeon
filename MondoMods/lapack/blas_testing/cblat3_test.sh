#!/bin/bash
./xblat3c < ${srcdir}/cblat3.in
if test `grep -i fail cblat3.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
