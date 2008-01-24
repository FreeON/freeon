#!/bin/bash
. test_functions.sh

run_test lin/xlintstz ztest.in ztest.out || exit 1
