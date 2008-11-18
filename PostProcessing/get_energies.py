#!/usr/bin/python

import sys, re, optparse

parser = optparse.OptionParser(usage = "%prog [options] outputfile ...", \
    description = """This script parses the output file(s) given on the command
line for the final SCF energy. The special filename '-' means standard input.""")

option, argument = parser.parse_args()

if len(argument) == 0:
  print >> sys.stderr, "reading from standard input"
  argument.append("-")

for file in argument:
  if file == "-":
    lines = sys.stdin.readlines()
  else:
    fd = open(file)
    lines = fd.readlines()
    fd.close()

  lastGeometry = 1
  lastSCFCycle = 0
  Geometry = 1
  SCFCycle = 0
  SCFEnergy = 0.0

  for line in lines:
    result = re.compile("SCFCycle #([0-9]+),.+Geometry #([0-9]+)").search(line)
    if result:
      lastGeometry = Geometry
      lastSCFCycle = SCFCycle
      SCFCycle = int(result.group(1))
      Geometry = int(result.group(2))

      if Geometry > lastGeometry:
        print lastGeometry, SCFEnergy, lastSCFCycle

      continue

    result = re.compile("<SCF> += (.+)$").search(line)
    if result:
      SCFEnergy = float(result.group(1))
