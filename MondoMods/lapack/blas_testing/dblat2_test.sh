#!/bin/bash
./xblat2d < ${srcdir}/dblat2.in
if test `grep -i fail dblat2.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
