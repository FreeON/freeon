#!/bin/bash
. ${srcdir}/test_functions.sh

run_test eig/xeigtsts ${srcdir}/nep.in   snep.out  || exit 1
run_test eig/xeigtsts ${srcdir}/sep.in   ssep.out  || exit 1
run_test eig/xeigtsts ${srcdir}/svd.in   ssvd.out  || exit 1
run_test eig/xeigtsts ${srcdir}/sec.in   sec.out   || exit 1
run_test eig/xeigtsts ${srcdir}/sed.in   sed.out   || exit 1
run_test eig/xeigtsts ${srcdir}/sgg.in   sgg.out   || exit 1
run_test eig/xeigtsts ${srcdir}/sgd.in   sgd.out   || exit 1
run_test eig/xeigtsts ${srcdir}/ssb.in   ssb.out   || exit 1
run_test eig/xeigtsts ${srcdir}/ssg.in   ssg.out   || exit 1
run_test eig/xeigtsts ${srcdir}/sbal.in  sbal.out  || exit 1
run_test eig/xeigtsts ${srcdir}/sbak.in  sbak.out  || exit 1
run_test eig/xeigtsts ${srcdir}/sgbal.in sgbal.out || exit 1
run_test eig/xeigtsts ${srcdir}/sgbak.in sgbak.out || exit 1
run_test eig/xeigtsts ${srcdir}/sbb.in   sbb.out   || exit 1
run_test eig/xeigtsts ${srcdir}/glm.in   sglm.out  || exit 1
run_test eig/xeigtsts ${srcdir}/gqr.in   sgqr.out  || exit 1
run_test eig/xeigtsts ${srcdir}/gsv.in   sgsv.out  || exit 1
run_test eig/xeigtsts ${srcdir}/lse.in   slse.out  || exit 1
