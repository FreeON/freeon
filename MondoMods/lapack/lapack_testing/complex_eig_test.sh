#!/bin/bash
. test_functions.sh

run_test eig/xeigtstc nep.in   cnep.out  || exit 1
run_test eig/xeigtstc sep.in   csep.out  || exit 1
run_test eig/xeigtstc svd.in   csvd.out  || exit 1
run_test eig/xeigtstc cec.in   cec.out   || exit 1
run_test eig/xeigtstc ced.in   ced.out   || exit 1
run_test eig/xeigtstc cgg.in   cgg.out   || exit 1
run_test eig/xeigtstc cgd.in   cgd.out   || exit 1
run_test eig/xeigtstc csb.in   csb.out   || exit 1
run_test eig/xeigtstc csg.in   csg.out   || exit 1
run_test eig/xeigtstc cbal.in  cbal.out  || exit 1
run_test eig/xeigtstc cbak.in  cbak.out  || exit 1
run_test eig/xeigtstc cgbal.in cgbal.out || exit 1
run_test eig/xeigtstc cgbak.in cgbak.out || exit 1
run_test eig/xeigtstc cbb.in   cbb.out   || exit 1
run_test eig/xeigtstc glm.in   cglm.out  || exit 1
run_test eig/xeigtstc gqr.in   cgqr.out  || exit 1
run_test eig/xeigtstc gsv.in   cgsv.out  || exit 1
run_test eig/xeigtstc lse.in   clse.out  || exit 1
