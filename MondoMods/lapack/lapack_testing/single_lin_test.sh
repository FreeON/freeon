#!/bin/bash
. test_functions.sh

run_test lin/xlintsts stest.in stest.out || exit 1
