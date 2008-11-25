#!/bin/bash
. ${srcdir}/test_functions.sh

run_test lin/xlintsts ${srcdir}/stest.in stest.out || exit 1
