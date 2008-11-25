#!/bin/bash
. ${srcdir}/test_functions.sh

run_test eig/xeigtstz ${srcdir}/nep.in   znep.out  || exit 1
run_test eig/xeigtstz ${srcdir}/sep.in   zsep.out  || exit 1
run_test eig/xeigtstz ${srcdir}/svd.in   zsvd.out  || exit 1
run_test eig/xeigtstz ${srcdir}/zec.in   zec.out   || exit 1
run_test eig/xeigtstz ${srcdir}/zed.in   zed.out   || exit 1
run_test eig/xeigtstz ${srcdir}/zgg.in   zgg.out   || exit 1
run_test eig/xeigtstz ${srcdir}/zgd.in   zgd.out   || exit 1
run_test eig/xeigtstz ${srcdir}/zsb.in   zsb.out   || exit 1
run_test eig/xeigtstz ${srcdir}/zsg.in   zsg.out   || exit 1
run_test eig/xeigtstz ${srcdir}/zbal.in  zbal.out  || exit 1
run_test eig/xeigtstz ${srcdir}/zbak.in  zbak.out  || exit 1
run_test eig/xeigtstz ${srcdir}/zgbal.in zgbal.out || exit 1
run_test eig/xeigtstz ${srcdir}/zgbak.in zgbak.out || exit 1
run_test eig/xeigtstz ${srcdir}/zbb.in   zbb.out   || exit 1
run_test eig/xeigtstz ${srcdir}/glm.in   zglm.out  || exit 1
run_test eig/xeigtstz ${srcdir}/gqr.in   zgqr.out  || exit 1
run_test eig/xeigtstz ${srcdir}/gsv.in   zgsv.out  || exit 1
run_test eig/xeigtstz ${srcdir}/lse.in   zlse.out  || exit 1
