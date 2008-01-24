#!/bin/bash
. test_functions.sh

run_test eig/xeigtsts nep.in   snep.out  || exit 1
run_test eig/xeigtsts sep.in   ssep.out  || exit 1
run_test eig/xeigtsts svd.in   ssvd.out  || exit 1
run_test eig/xeigtsts sec.in   sec.out   || exit 1
run_test eig/xeigtsts sed.in   sed.out   || exit 1
run_test eig/xeigtsts sgg.in   sgg.out   || exit 1
run_test eig/xeigtsts sgd.in   sgd.out   || exit 1
run_test eig/xeigtsts ssb.in   ssb.out   || exit 1
run_test eig/xeigtsts ssg.in   ssg.out   || exit 1
run_test eig/xeigtsts sbal.in  sbal.out  || exit 1
run_test eig/xeigtsts sbak.in  sbak.out  || exit 1
run_test eig/xeigtsts sgbal.in sgbal.out || exit 1
run_test eig/xeigtsts sgbak.in sgbak.out || exit 1
run_test eig/xeigtsts sbb.in   sbb.out   || exit 1
run_test eig/xeigtsts glm.in   sglm.out  || exit 1
run_test eig/xeigtsts gqr.in   sgqr.out  || exit 1
run_test eig/xeigtsts gsv.in   sgsv.out  || exit 1
run_test eig/xeigtsts lse.in   slse.out  || exit 1
