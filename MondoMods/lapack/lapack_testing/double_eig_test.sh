#!/bin/bash
. ${srcdir}/test_functions.sh

run_test eig/xeigtstd ${srcdir}/nep.in   dnep.out  || exit 1
run_test eig/xeigtstd ${srcdir}/sep.in   dsep.out  || exit 1
run_test eig/xeigtstd ${srcdir}/svd.in   dsvd.out  || exit 1
run_test eig/xeigtstd ${srcdir}/dec.in   dec.out   || exit 1
run_test eig/xeigtstd ${srcdir}/ded.in   ded.out   || exit 1
run_test eig/xeigtstd ${srcdir}/dgg.in   dgg.out   || exit 1
run_test eig/xeigtstd ${srcdir}/dgd.in   dgd.out   || exit 1
run_test eig/xeigtstd ${srcdir}/dsb.in   dsb.out   || exit 1
run_test eig/xeigtstd ${srcdir}/dsg.in   dsg.out   || exit 1
run_test eig/xeigtstd ${srcdir}/dbal.in  dbal.out  || exit 1
run_test eig/xeigtstd ${srcdir}/dbak.in  dbak.out  || exit 1
run_test eig/xeigtstd ${srcdir}/dgbal.in dgbal.out || exit 1
run_test eig/xeigtstd ${srcdir}/dgbak.in dgbak.out || exit 1
run_test eig/xeigtstd ${srcdir}/dbb.in   dbb.out   || exit 1
run_test eig/xeigtstd ${srcdir}/glm.in   dglm.out  || exit 1
run_test eig/xeigtstd ${srcdir}/gqr.in   dgqr.out  || exit 1
run_test eig/xeigtstd ${srcdir}/gsv.in   dgsv.out  || exit 1
run_test eig/xeigtstd ${srcdir}/lse.in   dlse.out  || exit 1
