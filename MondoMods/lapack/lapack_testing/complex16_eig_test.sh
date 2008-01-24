#!/bin/bash
. test_functions.sh

run_test eig/xeigtstz nep.in   znep.out  || exit 1
run_test eig/xeigtstz sep.in   zsep.out  || exit 1
run_test eig/xeigtstz svd.in   zsvd.out  || exit 1
run_test eig/xeigtstz zec.in   zec.out   || exit 1
run_test eig/xeigtstz zed.in   zed.out   || exit 1
run_test eig/xeigtstz zgg.in   zgg.out   || exit 1
run_test eig/xeigtstz zgd.in   zgd.out   || exit 1
run_test eig/xeigtstz zsb.in   zsb.out   || exit 1
run_test eig/xeigtstz zsg.in   zsg.out   || exit 1
run_test eig/xeigtstz zbal.in  zbal.out  || exit 1
run_test eig/xeigtstz zbak.in  zbak.out  || exit 1
run_test eig/xeigtstz zgbal.in zgbal.out || exit 1
run_test eig/xeigtstz zgbak.in zgbak.out || exit 1
run_test eig/xeigtstz zbb.in   zbb.out   || exit 1
run_test eig/xeigtstz glm.in   zglm.out  || exit 1
run_test eig/xeigtstz gqr.in   zgqr.out  || exit 1
run_test eig/xeigtstz gsv.in   zgsv.out  || exit 1
run_test eig/xeigtstz lse.in   zlse.out  || exit 1
