#!/bin/bash
./xblat1s > sblat1.out
if test `grep -i fail sblat1.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
