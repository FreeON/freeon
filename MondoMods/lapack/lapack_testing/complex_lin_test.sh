#!/bin/bash
. ${srcdir}/test_functions.sh

run_test lin/xlintstc ${srcdir}/ctest.in ctest.out || exit 1
