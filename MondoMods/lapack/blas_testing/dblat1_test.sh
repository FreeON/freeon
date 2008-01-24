#!/bin/bash
./xblat1d > dblat1.out
if test `grep -i fail dblat1.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
