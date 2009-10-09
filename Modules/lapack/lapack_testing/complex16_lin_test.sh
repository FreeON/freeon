#!/bin/bash
. ${srcdir}/test_functions.sh

run_test lin/xlintstz ${srcdir}/ztest.in ztest.out || exit 1
