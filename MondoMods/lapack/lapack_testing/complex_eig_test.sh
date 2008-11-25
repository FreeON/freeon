#!/bin/bash
. ${srcdir}/test_functions.sh

run_test eig/xeigtstc ${srcdir}/nep.in   cnep.out  || exit 1
run_test eig/xeigtstc ${srcdir}/sep.in   csep.out  || exit 1
run_test eig/xeigtstc ${srcdir}/svd.in   csvd.out  || exit 1
run_test eig/xeigtstc ${srcdir}/cec.in   cec.out   || exit 1
run_test eig/xeigtstc ${srcdir}/ced.in   ced.out   || exit 1
run_test eig/xeigtstc ${srcdir}/cgg.in   cgg.out   || exit 1
run_test eig/xeigtstc ${srcdir}/cgd.in   cgd.out   || exit 1
run_test eig/xeigtstc ${srcdir}/csb.in   csb.out   || exit 1
run_test eig/xeigtstc ${srcdir}/csg.in   csg.out   || exit 1
run_test eig/xeigtstc ${srcdir}/cbal.in  cbal.out  || exit 1
run_test eig/xeigtstc ${srcdir}/cbak.in  cbak.out  || exit 1
run_test eig/xeigtstc ${srcdir}/cgbal.in cgbal.out || exit 1
run_test eig/xeigtstc ${srcdir}/cgbak.in cgbak.out || exit 1
run_test eig/xeigtstc ${srcdir}/cbb.in   cbb.out   || exit 1
run_test eig/xeigtstc ${srcdir}/glm.in   cglm.out  || exit 1
run_test eig/xeigtstc ${srcdir}/gqr.in   cgqr.out  || exit 1
run_test eig/xeigtstc ${srcdir}/gsv.in   cgsv.out  || exit 1
run_test eig/xeigtstc ${srcdir}/lse.in   clse.out  || exit 1
