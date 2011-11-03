#!/usr/bin/python
#
# This code is part of the MondoSCF suite of programs for linear scaling
# electronic structure theory and ab initio molecular dynamics.
#
# Copyright (2004). The Regents of the University of California. This
# material was produced under U.S. Government contract W-7405-ENG-36
# for Los Alamos National Laboratory, which is operated by the University
# of California for the U.S. Department of Energy. The U.S. Government has
# rights to use, reproduce, and distribute this software.  NEITHER THE
# GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS OR IMPLIED,
# OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version. Accordingly, this program is distributed in
# the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
# the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU General Public License at www.gnu.org for details.
#
# While you may do as you like with this software, the GNU license requires
# that you clearly mark derivative software.  In addition, you are encouraged
# to return derivative works to the MondoSCF group for review, and possible
# disemination in future releases.
#
# Regression tests for FreeON.
#
# Nicolas Bock <nbock@lanl.gov>

import logging
import optparse
import os
import re
import sys

parser = optparse.OptionParser()

parser.add_option("--reference",
    metavar = "FILE",
    help = "load reference values from FILE",
    dest = "reference")

parser.add_option("--output",
    metavar = "FILE",
    help = "load FreeON output from FILE instead of from standard input",
    dest = "output")

parser.add_option("--verbose", "-v",
    help = "print lots of stuff.",
    dest = "verbose",
    action = "store_true")

( options, arguments ) = parser.parse_args()

# Set up logging.
log = logging.getLogger("regressionTest")
log.setLevel(logging.DEBUG)

# Set console logger.
logHandler = logging.StreamHandler()
logFormatter = logging.Formatter("[%(name)s] %(message)s")
logHandler.setFormatter(logFormatter)
logHandler.setLevel(logging.INFO)

if options.verbose:
  logHandler.setLevel(logging.DEBUG)

log.addHandler(logHandler)

# Set file logger.
logHandler = logging.FileHandler("regressionTest.log")
logFormatter = logging.Formatter("%(asctime)s [%(name)s %(levelname)s] %(message)s", "%y-%m-%d %H:%M:%S")
logHandler.setFormatter(logFormatter)
logHandler.setLevel(logging.DEBUG)

log.addHandler(logHandler)

log.debug("test")

if not options.reference:
  log.error("missing reference file")
  parser.print_help()
  sys.exit(1)

if not options.output:
  # Load output from standard input.
  log.debug("loading output from standard input")
  output = sys.stdin.readlines()

else:
  # Load output from file.
  log.debug("loading output from file: " + options.output)
  fd = open(options.output)
  output = fd.readlines()
  fd.close()

# Parse output against reference.
log.debug("loading reference tags")
reference = {}

fd = open(options.reference)
lines = fd.readlines()

# numberErrors:
#   A matched tag, i.e. a line that matches the regular expression given, is
#   checked against the value and error given in the reference file. If the
#   value is not within the error, numberErrors is incremented by one.
numberErrors = 0

# numberUnmatched:
#   If the regular expression matches a tag, but all lines corresponding to
#   this tag have already been found in the output, then numberUnmatched is
#   incremented by one.
numberUnmatched = 0

# numberMatched:
#   If the regular expression matches a tag, numberMatched is incremented by
#   one, regardless whether the value was within the error or not.
numberMatched = 0

# numberMissing:
#   A line in the reference file that can not be matched with the output
#   causes numberMissing to be > 0.
numberMissing = 0

linenumber = 0
for line in lines:
  line = line.strip()
  linenumber += 1

  check = re.compile("(^[^#]*)(#.+$)").search(line)
  if check:
    line = check.group(1).strip()

  if len(line) == 0:
    continue

  check = re.compile("\"(.*)\" *= *([0-9.DdEe+\-]*) *\((([0-9.DdEe+\-]*))\)").search(line)
  if check:
    log.debug("read tag from reference \"" + check.group(1) + "\", value = " + check.group(2) + ", error = " + check.group(3))

    # Convert fortran exponential "D" to "E".
    value = check.group(2)
    error = check.group(3)
    value = value.lower().replace("d", "e")
    error = error.lower().replace("d", "e")
    new_tag = { "value": float(value), "error": float(error) }
    if check.group(1) in reference:
      reference[check.group(1)]["values"].append(new_tag)
    else:
      reference[check.group(1)] = { "index": 0, "values": [ new_tag ] }
    numberMissing += 1

log.debug("checking tags: " + str(reference))

successfullyTerminated = False
scratchDirectory = None

linenumber = 0
for line in output:
  linenumber += 1
  line = line.rstrip()

  check = re.compile("Successful FreeON run").search(line)
  if check:
    successfullyTerminated = True

  check = re.compile("scratch directory at (.*)").search(line)
  if check:
    scratchDirectory = check.group(1)

  for tag in reference.keys():
    check = re.compile(tag).search(line)
    if check:
      if reference[tag]["index"] >= len(reference[tag]["values"]):
        log.error("unmatched key on line " + str(linenumber) + ": " + line)
        numberUnmatched += 1
        continue

      # Convert number from f90 format.
      value = float(check.group(1).lower().replace("d", "e"))

      # Compare to reference.
      ref_value = reference[tag]["values"][reference[tag]["index"]]["value"]
      ref_error = reference[tag]["values"][reference[tag]["index"]]["error"]
      if abs(value-ref_value) > ref_error:
        numberErrors += 1
        log.error("line " + str(linenumber) + ", " + line)
        log.error("--> wrong value " + str(value))
        log.error("--> expected " + str(ref_value) + " +- " + str(ref_error))
        log.error("--> observed difference = " + str(abs(value-ref_value)))
        log.error("--> tag \"" + tag + "\", " + "index " + str(reference[tag]["index"]))
      else:
        log.debug("line " + str(linenumber) + ", found tag " + \
            tag + ": " + line + " <--> value verified")

      reference[tag]["index"] += 1
      numberMatched += 1
      numberMissing -= 1

if options.output:
  log.info("output file: " + options.output)

log.info("reference file: " + options.reference)
log.info("scratch directory: " + scratchDirectory)
log.info("log file: " + os.path.join(os.getcwd(), "regressionTest.log"))

log.info("found " + str(numberErrors) + " value errors, " \
    + str(numberMatched)   + " matched entries, " \
    + str(numberUnmatched) + " unmatched entries, i.e found tag in output but did not have reference for it, " \
    + str(numberMissing)   + " entries in reference that were missing from the output")

if successfullyTerminated:
  log.info("FreeON successfully terminated")
else:
  log.info("FreeON did not successfully terminate")

# Exit with a proper exit code.
if numberErrors == 0 and numberMissing == 0 and numberUnmatched == 0:
  sys.exit(0)
else:
  sys.exit(99)
