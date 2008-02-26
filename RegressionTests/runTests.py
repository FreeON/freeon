#!/usr/bin/python
#
# Regression tests for MondoSCF.
#
# Nicolas Bock <nbock@lanl.gov>

import optparse, sys, logging
import lib

parse = optparse.OptionParser(usage = "%prog [options] input.test", \
    description = """This scripts reads input from the file given on the command
line and executes the tests in there.""")

parse.add_option("--tar", \
    help = "Build Mondo from tarfile.", \
    metavar = "tarfile", \
    dest = "tar", \
    type = "string")

option, argument = parse.parse_args()

logging.basicConfig(level = logging.DEBUG, \
    filename = "regression_test.log", \
    format = "%(asctime)s [%(name)s %(levelname)s] %(message)s", \
    datefmt='%y-%m-%d %H:%M:%S')

console_handler = logging.StreamHandler()
console_formatter = logging.Formatter("[%(name)s] %(message)s")
console_handler.setFormatter(console_formatter)

console_handler.setLevel(logging.INFO)

log = logging.getLogger("main")

log.info("starting new regression test")

if len(argument) != 1:
  log.error("I need exactly 1 input file")
  sys.exit(1)

# Parse the input file.
testfile = lib.parse_input(argument[0])

# Run the tests.
