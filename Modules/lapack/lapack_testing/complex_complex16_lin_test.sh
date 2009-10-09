#!/bin/bash
. ${srcdir}/test_functions.sh

run_test lin/xlintstzc ${srcdir}/zctest.in zctest.out || exit 1
