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
import regressionlib

dirprefix = "/tmp/MondoRegressionTest"

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
console_handler.setLevel(logging.DEBUG)

logging.getLogger().addHandler(console_handler)

log = logging.getLogger("main")

log.info("starting new regression test")

if len(argument) != 1:
  log.error("I need exactly 1 input file")
  sys.exit(1)

# Parse the input file.
testfile = regressionlib.parse_input(argument[0])

# Run the tests.
if "Mondo_tar" in testfile.field:
  if not os.path.exists(testfile.field["Mondo_tar"]):
    log.error("tar file does not exist: " + testfile.field["Mondo_tar"])
    sys.exit(1)

  log.debug("checking tar")

  # Come up with paths.
  builddir = os.path.join(dirprefix, "build")
  installdir = os.path.join(dirprefix, "install")

  # Extract tar file version name.
  tarnamebase = os.path.basename(testfile.field["Mondo_tar"])
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
    shutil.copy(testfile.field["Mondo_tar"], builddir)

    log.debug("unpacking sources")
    src = tarfile.open(testfile.field["Mondo_tar"], "r:bz2")
    for member in src:
      src.extract(member, builddir)
    src.close()

    log.debug("configuring")
    configure_arguments = [ "./configure" ]

    if "configure_options" in testfile.field:
      configure_arguments += testfile.field["configure_options"].split()

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

    if "make_options" in testfile.field:
      make_arguments += testfile.field["make_options"].split()

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
      log.error("make stage failed with return code " + str(make.returncode))

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
      log.error("install stage failed with return code " + str(install.returncode))

      log.error("standard output:")
      for line in install_stdout:
        log.error(line.strip())

      log.error("standard error:")
      for line in install_stderr:
        log.error(line.strip())

      sys.exit(1)

    else:
      log.debug("installing done")

    log.info("sources built and install correctly")
