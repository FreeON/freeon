#!/bin/bash
. test_functions.sh

run_test eig/xeigtstd nep.in   dnep.out  || exit 1
run_test eig/xeigtstd sep.in   dsep.out  || exit 1
run_test eig/xeigtstd svd.in   dsvd.out  || exit 1
run_test eig/xeigtstd dec.in   dec.out   || exit 1
run_test eig/xeigtstd ded.in   ded.out   || exit 1
run_test eig/xeigtstd dgg.in   dgg.out   || exit 1
run_test eig/xeigtstd dgd.in   dgd.out   || exit 1
run_test eig/xeigtstd dsb.in   dsb.out   || exit 1
run_test eig/xeigtstd dsg.in   dsg.out   || exit 1
run_test eig/xeigtstd dbal.in  dbal.out  || exit 1
run_test eig/xeigtstd dbak.in  dbak.out  || exit 1
run_test eig/xeigtstd dgbal.in dgbal.out || exit 1
run_test eig/xeigtstd dgbak.in dgbak.out || exit 1
run_test eig/xeigtstd dbb.in   dbb.out   || exit 1
run_test eig/xeigtstd glm.in   dglm.out  || exit 1
run_test eig/xeigtstd gqr.in   dgqr.out  || exit 1
run_test eig/xeigtstd gsv.in   dgsv.out  || exit 1
run_test eig/xeigtstd lse.in   dlse.out  || exit 1
