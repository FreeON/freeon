#!/bin/bash
./xblat1c > cblat1.out
if test `grep -i fail cblat1.out | wc -l` -ne 0; then
  exit 1
fi

exit 0
