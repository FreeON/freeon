#!/bin/bash
. ${srcdir}/test_functions.sh

run_test lin/xlintstd ${srcdir}/dtest.in dtest.out || exit 1
