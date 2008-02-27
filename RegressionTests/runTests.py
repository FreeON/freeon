#!/usr/bin/python
#
# Regression tests for MondoSCF.
#
# Nicolas Bock <nbock@lanl.gov>

import optparse
import sys
import logging
import re
import os.path
import shutil
import tarfile
import subprocess

# Some variables.
dirprefix = "/tmp/MondoRegressionTest"
builddir = os.path.join(dirprefix, "build")
installdir = os.path.join(dirprefix, "install")
rundirbase = os.path.join(dirprefix, "run")
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

option, argument = parse.parse_args()

# Set up logging.
logging.basicConfig(level = logging.DEBUG, \
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
  for fieldname in [ "Mondo_input", \
      "Mondo_executable", \
      "Mondo_tar", \
      "configure_options", \
      "make_options", \
      "logfile" ]:
    check = re.compile("(^\s*" + fieldname + "\s*=)(.*$)").search(line)
    if check:
      if not check.group(2):
        log.error("syntax error on line " + str(linenumber))
        sys.exit(1)
      else:
        inputfield[fieldname] = check.group(2).strip()
        log.debug("found input tag " + fieldname + " with value " + inputfield[fieldname])

  linenumber += 1

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
  log.error("no tarfile and no exectuable, I need something here....")
  sys.exit(1)

if not "Mondo_input" in inputfield:
  log.info("No Mondo_input given, we are done")
  sys.exit(0)
else:
  Mondo_inputbase = os.path.basename(inputfield["Mondo_input"])

# Run the test.
log.debug("creating run directory")
try:
  os.makedirs(rundirbase)
except:
  log.error("error making rundir")

# Create run directory for this test.
rundir = os.path.join(rundirbase, Mondo_inputbase)
while(True):
  rundir, ext = os.path.splitext(rundir)
  if not ext:
    break

log.debug("creating directory " + rundir)
if os.path.exists(rundir):
  log.error("rundir already exists " + rundir)
  sys.exit(1)
os.mkdir(rundir)

log.debug("copying input file to rundir " + rundir)
shutil.copy(inputfield["Mondo_input"], rundir)

Mondo_arguments = [ inputfield["Mondo_executable"], Mondo_inputbase ]

log.debug("running " + str(Mondo_arguments))
Mondo = subprocess.Popen(Mondo_arguments, \
    stdin = subprocess.PIPE, \
    stdout = subprocess.PIPE, \
    stderr = subprocess.PIPE, \
    cwd = rundir)

Mondo_stdout = Mondo.stdout.readlines()
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
