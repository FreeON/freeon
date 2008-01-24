#!/bin/bash
. test_functions.sh

run_test lin/xlintstzc zctest.in zctest.out || exit 1
