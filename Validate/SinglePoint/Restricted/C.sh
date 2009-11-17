#!/bin/bash

FreeON C.inp | ${srcdir}/../../parseTestResults.py --reference C.reference
