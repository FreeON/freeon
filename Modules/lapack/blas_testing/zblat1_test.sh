#!/bin/bash
./xblat1z > zblat1.out
if test `grep -i fail zblat1.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
