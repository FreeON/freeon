#!/usr/bin/python
#
# Script to extract the diffusion coefficient from an MD calculation.
#
# Nicolas Bock <nbock@lanl.gov>

import sys, optparse, parseMDO, logging

parser = optparse.OptionParser(usage = "%prog [options] file.MDO", \
    description = """This script parses an MD output file (file.MDO) and
calculates the diffusion coefficient from the positions and velocities.""")

option, argument = parser.parse_args()

if len(argument) != 1:
  print("I need exactly one MDO file")
  sys.exit(1)

# Set up the logger.
logging.basicConfig(level = logging.INFO)
log = logging.getLogger("diffusionCoefficient")

# Parse positions and velocities.
MDO = parseMDO.parseMDO(argument[0])
time, position, velocity = MDO.parse()
log.info("analyzing " + str(len(time)) + " MD steps")

# Calculate something.

