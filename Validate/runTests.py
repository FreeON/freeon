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
# Regression tests for MondoSCF.
#
# Nicolas Bock <nbock@lanl.gov>

import logging
import optparse
import os.path
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

# Some variables.
#dirprefix = "/tmp/MondoRegressionTest"
dirprefix = tempfile.mkdtemp()
builddir = os.path.join(dirprefix, "build")
installdir = os.path.join(dirprefix, "install")
rundirbase = tempfile.mkdtemp(dir=dirprefix, prefix="run")
rundir = ""
Mondo_inputbase = ""
inputfield = {}

# Start of the main program.
parse = optparse.OptionParser(usage = "%prog [options] input.test", \
    description = """This scripts reads input from the file given on the command
line and executes the tests in there.""")

parse.add_option("--tar", \
    help = "Build Mondo from tarfile.", \
    metavar = "tarfile", \
    dest = "tar", \
    type = "string")

parse.add_option("--verbose", "-v", \
    help = "print lots of stuff.", \
    dest = "verbose", \
    action = "store_true")

parse.add_option("--clean-run", \
    help = "clean the run directory", \
    dest = "clean_run", \
    action = "store_true")

parse.add_option("--clean-build", \
    help = "clean the build directory", \
    dest = "clean_build", \
    action = "store_true")

parse.add_option("--pretend", \
    help = "do not run the test, simply anaylize the output file", \
    dest = "pretend", \
    action = "store_true")

parse.add_option("--FreeON",
    help = "specify explicitly where the FreeON executable is",
    dest = "executable")

option, argument = parse.parse_args()

# Set up logging.
logging.basicConfig(level = logging.INFO, \
    filename = "regression_test.log", \
    format = "%(asctime)s [%(name)s %(levelname)s] %(message)s", \
    datefmt='%y-%m-%d %H:%M:%S')

console_handler = logging.StreamHandler()
console_formatter = logging.Formatter("[%(name)s] %(message)s")
console_handler.setFormatter(console_formatter)

if option.verbose:
  console_handler.setLevel(logging.DEBUG)
else:
  console_handler.setLevel(logging.INFO)

logging.getLogger().addHandler(console_handler)

log = logging.getLogger("main")

# Start...
log.info("starting new regression test")

if len(argument) != 1:
  log.error("I need exactly 1 input file")
  sys.exit(1)

# Parse the input file.
inputfile = argument[0]
log.debug("reading from file " + inputfile)

if not os.path.exists(inputfile):
  log.error("input file " + inputfile + " does not exist")
  sys.exit(1)

fd = open(inputfile, 'r')
lines = fd.readlines()
fd.close()

# Parse through input file.
linenumber = 1
for line in lines:
  # Remove comments.
  line = line.strip()
  check = re.compile("(^[^#]*)(#.+$)").search(line)
  if check:
    line = check.group(1)

  # Search for known keywords.
  for fieldname in [
      "Mondo_input",
      "Mondo_reference",
      "Mondo_executable",
      "Mondo_tar",
      "configure_options",
      "make_options",
      "logfile",
      "suppress_stdout" ]:
    check = re.compile("(^\s*" + fieldname + "\s*=)(.*$)").search(line)
    if check:
      if not check.group(2):
        log.error("syntax error on line " + str(linenumber))
        sys.exit(1)
      else:
        inputfield[fieldname] = check.group(2).strip()
        log.debug("found input tag " + fieldname + " with value " + inputfield[fieldname])

  linenumber += 1

# Set some values straight.
if "suppress_stdout" not in inputfield:
  inputfield["suppress_stdout"] = False

else:
  if inputfield["suppress_stdout"].lower() in [ "no", "false" ]:
    inputfield["suppress_stdout"] = False
  elif inputfield["suppress_stdout"].lower() in [ "yes", "true" ]:
    inputfield["suppress_stdout"] = True
  else:
    log.error("I do not understand your input for suppress_stdout")
    sys.exit(1)

# Run the tests.
if "Mondo_tar" in inputfield:
  if not os.path.exists(inputfield["Mondo_tar"]):
    log.error("tar file does not exist: " + inputfield["Mondo_tar"])
    sys.exit(1)

  log.debug("checking tar")

  # Extract tar file version name.
  tarnamebase = os.path.basename(inputfield["Mondo_tar"])
  log.debug("stripped leading directories: " + tarnamebase)

  tarext = []
  tarname = tarnamebase
  while(True):
    tarext.append("")
    tarname, tarext[-1] = os.path.splitext(tarname)
    if not tarext[-1]:
      del tarext[-1]
      break

    if not (tarext[-1] in [ ".tar", ".bz2" ]):
      tarname = tarname + tarext[-1]
      del tarext[-1]
      break

  log.debug("tar file version name is " + tarname)
  log.debug("tar extensions are " + str(tarext))

  # Check for clean_build option.
  if option.clean_build and os.path.exists(builddir):
    log.info("wiping builddir " + builddir)
    shutil.rmtree(builddir)

  # Check whether we already built this source.
  if os.path.exists(os.path.join(builddir, tarname)):
    log.info("sources are apparently already built.")

  else:
    log.info("building sources")
    try:
      os.makedirs(builddir)
    except:
      log.error("error making builddir")

    log.debug("copying tar file to builddir " + builddir)
    shutil.copy(inputfield["Mondo_tar"], builddir)

    log.debug("unpacking sources")
    src = tarfile.open(inputfield["Mondo_tar"], "r:bz2")
    for member in src:
      src.extract(member, builddir)
    src.close()

    log.debug("configuring")
    configure_arguments = [ "./configure" ]

    if "configure_options" in inputfield:
      configure_arguments += inputfield["configure_options"].split()

    configure_arguments.append("--prefix=" + installdir)

    log.debug("running " + str(configure_arguments))
    configure = subprocess.Popen(configure_arguments, \
        stdin = subprocess.PIPE, \
        stdout = subprocess.PIPE, \
        stderr = subprocess.PIPE, \
        cwd = os.path.join(builddir, tarname))

    configure_stdout = configure.stdout.readlines()
    configure_stderr = configure.stderr.readlines()

    configure.wait()

    if configure.returncode != 0:
      log.error("configure failed with return code " \
          + str(configure.returncode))

      log.error("standard output:")
      for line in configure_stdout:
        log.error(line.strip())

      log.error("standard error:")
      for line in configure_stderr:
        log.error(line.strip())

      sys.exit(1)

    else:
      log.debug("configure done")

    log.debug("building")
    make_arguments = [ "make" ]

    if "make_options" in inputfield:
      make_arguments += inputfield["make_options"].split()

    log.debug("running " + str(make_arguments))
    make = subprocess.Popen(make_arguments, \
        stdin = subprocess.PIPE, \
        stdout = subprocess.PIPE, \
        stderr = subprocess.PIPE, \
        cwd = os.path.join(builddir, tarname))

    make_stdout = make.stdout.readlines()
    make_stderr = make.stderr.readlines()

    make.wait()

    if make.returncode != 0:
      log.error("make failed with return code " + str(make.returncode))

      log.error("standard output:")
      for line in make_stdout:
        log.error(line.strip())

      log.error("standard error:")
      for line in make_stderr:
        log.error(line.strip())

      sys.exit(1)

    else:
      log.debug("building done")

    log.debug("installing binaries")
    install = subprocess.Popen([ "make", "install" ], \
        stdin = subprocess.PIPE, \
        stdout = subprocess.PIPE, \
        stderr = subprocess.PIPE, \
        cwd = os.path.join(builddir, tarname))

    install_stdout = install.stdout.readlines()
    install_stderr = install.stderr.readlines()

    install.wait()

    if install.returncode != 0:
      log.error("install failed with return code " + str(install.returncode))

      log.error("standard output:")
      for line in install_stdout:
        log.error(line.strip())

      log.error("standard error:")
      for line in install_stderr:
        log.error(line.strip())

      sys.exit(1)

    else:
      log.debug("installing done")

    log.info("sources build and install correctly")

if "Mondo_executable" in inputfield:
  log.debug("running " + inputfield["Mondo_executable"])
elif "Mondo_tar" in inputfield:
  inputfield["Mondo_executable"] = os.path.join(installdir, "bin", "MondoSCF")
  log.debug("constructed Mondo_executable " + inputfield["Mondo_executable"])
else:
  log.info("no explicit FreeON executable given, using simply FreeON")
  inputfield["Mondo_executable"] = "FreeON"

if not "Mondo_input" in inputfield:
  log.info("No Mondo_input given, we are done")
  sys.exit(0)

# Construct input filename.
Mondo_inputbase = os.path.basename(inputfield["Mondo_input"])

# Run the test.
log.info("creating run directory " + rundirbase)
if not os.path.exists(rundirbase):
  try:
    os.makedirs(rundirbase)
  except:
    log.error("error making " + rundirbase)

# Create run directory for this test.
rundir = os.path.join(rundirbase, Mondo_inputbase)
while(True):
  rundir, ext = os.path.splitext(rundir)
  if not ext:
    break

if not option.pretend:
  log.debug("creating directory " + rundir)
  if os.path.exists(rundir) and not option.clean_run:
    log.error("rundir already exists " + rundir)
    sys.exit(1)
  if os.path.exists(rundir) and option.clean_run:
    log.info("rundir " + rundir + " already exists but I am told to wipe it")
    shutil.rmtree(rundir)
  os.mkdir(rundir)

  log.debug("copying input file to rundir " + rundir)
  shutil.copy(inputfield["Mondo_input"], rundir)

  Mondo_arguments = [ inputfield["Mondo_executable"], Mondo_inputbase ]

  log.info("running " + str(Mondo_arguments))
  Mondo = subprocess.Popen(Mondo_arguments, \
      stdin = subprocess.PIPE, \
      stdout = subprocess.PIPE, \
      stderr = subprocess.PIPE, \
      cwd = rundir)

  Mondo_stdout = []
  while not Mondo.stdout.closed:
    Mondo_stdout.append(Mondo.stdout.readline())
    if not inputfield["suppress_stdout"]:
      print Mondo_stdout[-1].strip()

  Mondo_stderr = Mondo.stderr.readlines()

  Mondo.wait()

  if Mondo.returncode != 0:
    log.error("Mondo failed with return code " + str(Mondo.returncode))

    log.error("standard output:")
    for line in Mondo_stdout:
      log.error(line.strip())

    log.error("standard error:")
    for line in Mondo_stderr:
      log.error(line.strip())

    sys.exit(1)

  else:
    log.debug("Mondo done")

# Verify reference result.
if not "Mondo_reference" in inputfield:
  log.info("no reference given, we are done")
  sys.exit(0)

# Read in the reference file.
log.debug("reading references from " + inputfield["Mondo_reference"])
fd = open(inputfield["Mondo_reference"])
lines = fd.readlines()
fd.close()

# The references.
last_tag = None
ref = {}
output_file = None

log.info("anaylyzing result")
linenumber = 0
for line in lines:
  line = line.strip()
  linenumber += 1

  check = re.compile("(^[^#]*)(#.+$)").search(line)
  if check:
    line = check.group(1).strip()

  if len(line) == 0:
    continue

  check = re.compile("output_file *= *(.*)$").search(line)
  if check:
    log.debug("found output_file " + check.group(1).strip())
    output_file = check.group(1).strip()
    continue

  check = re.compile("\"(.*)\" *= *([0-9.DdEe+\-]*) *\((([0-9.DdEe+\-]*))\)").search(line)
  if check:
    log.debug("found tag \"" + check.group(1) + "\", value = " + \
        check.group(2) + ", error = " + check.group(3))
    # Convert f90 exponential "D" to "E".
    value = check.group(2)
    error = check.group(3)
    value = value.lower().replace("d", "e")
    error = error.lower().replace("d", "e")
    new_tag = { "value": float(value), "error": float(error) }
    if check.group(1) in ref:
      ref[check.group(1)]["values"].append(new_tag)
    else:
      ref[check.group(1)] = { "index": 0, "values": [ new_tag ] }

log.debug("checking tags: " + str(ref))

# Find the output file.
files = os.listdir(rundir)
output_files = []
for file in files:
  check = re.compile(output_file).search(file)
  if check:
    log.debug("found possible output file: " + file)
    output_files.append(file)

if len(output_files) > 1:
  log.error("considered " + str(files) + " for output file")
  log.error("found more than 1 possible output file matching " + output_file)
  log.error(str(output_files))
  sys.exit(1)

if len(output_files) < 1:
  log.error("considered " + str(files) + " for output file")
  log.error("could not find output file matching " + output_file)
  sys.exit(1)

log.info("analyzing output file " + os.path.join(rundir, output_files[0]))
fd = open(os.path.join(rundir, output_files[0]))
lines = fd.readlines()
fd.close()

number_errors = 0
number_unmatched = 0

linenumber = 0
for line in lines:
  linenumber += 1
  line = line.strip()

  for tag in ref.keys():
    check = re.compile(tag).search(line)
    if check:
      if ref[tag]["index"] >= len(ref[tag]["values"]):
        log.error("unmatched key on line " + str(linenumber) + ": " + line)
        number_unmatched += 1
        continue

      # Convert number from f90 format.
      value = float(check.group(1).lower().replace("d", "e"))

      # Compare to reference.
      ref_value = ref[tag]["values"][ref[tag]["index"]]["value"]
      ref_error = ref[tag]["values"][ref[tag]["index"]]["error"]
      if abs(value-ref_value) > ref_error:
        number_errors += 1
        log.error("line " + str(linenumber) + ", " + line + \
            " <--> wrong value " + str(value) + ", " + \
            "expected " + str(ref_value) + ", " + \
            "difference = " + str(abs(value-ref_value)) + ", " + \
            "tag \"" + tag + "\", " + \
            "index " + str(ref[tag]["index"]) + ", " + \
            " +- " + str(ref_error))
      else:
        log.debug("line " + str(linenumber) + ", found tag " + \
            tag + ": " + line + " <--> value verified")

      ref[tag]["index"] += 1

log.info("done analyzing, found " + str(number_errors) + " errors and " + \
    str(number_unmatched) + " unmatched tags")

# Exit with a proper exit code.
if number_errors == 0:
  sys.exit(0)
else:
  sys.exit(1)
