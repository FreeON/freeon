#!/bin/bash
. ${srcdir}/test_functions.sh

run_test lin/xlintstds ${srcdir}/dstest.in dstest.out || exit 1
