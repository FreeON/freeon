#!/bin/bash
./xblat3d < ${srcdir}/dblat3.in
if test `grep -i fail dblat3.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
