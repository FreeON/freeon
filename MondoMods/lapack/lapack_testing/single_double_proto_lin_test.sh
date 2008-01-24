#!/bin/bash
. test_functions.sh

run_test lin/xlintstds dstest.in dstest.out || exit 1
