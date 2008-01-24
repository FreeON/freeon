#!/bin/bash
. test_functions.sh

run_test lin/xlintstd dtest.in dtest.out || exit 1
