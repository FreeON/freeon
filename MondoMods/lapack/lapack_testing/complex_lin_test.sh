#!/bin/bash
. test_functions.sh

run_test lin/xlintstc ctest.in ctest.out || exit 1
