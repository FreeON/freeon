#!/usr/bin/python
#
# vim: tw=0

import sys, re, optparse

parser = optparse.OptionParser(usage = "%prog [options] outputfile ...",
    description = """This script parses the output file(s) given on the command
line for the final SCF energy. The special filename '-' means standard input.
The output is written to standard output. There are 4 columns: (1) geometry,
(2) SCF Energy, (3) SCFCycle, (4) dD, the difference in electronic density""")

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

  lastSCFCycle = 0
  lastGeometry = 1
  lastBasis = 1

  SCFCycle = 0
  Basis = 1
  Geometry = 1
  SCFEnergy = 0.0
  dD = 0.0

  lineNumber = 0

  for line in lines:
    lineNumber += 1
    result = re.compile("SCF.+\[([0-9]+),([0-9]+),([0-9]+)\].+<SCF> = ([-+.0-9dDeE]+).*dD = ([-+.0-9dDeE]+).*").search(line)
    if result:
      lastSCFCycle = SCFCycle
      lastBasis    = Basis
      lastGeometry = Geometry
      lastdD       = dD
      SCFCycle = int(result.group(1))
      Basis    = int(result.group(2))
      Geometry = int(result.group(3))

      if Geometry > lastGeometry:
        print lastGeometry, SCFEnergy, lastSCFCycle, dD

      try:
        SCFEnergy = float(result.group(4).lower().replace("d", "e"))
      except ValueError, message:
        print "anaylzing line " + str(lineNumber) + ": " + line.strip()
        print result.group(1) + " " + result.group(2) + " " + result.group(3) + " " + result.group(4) + " " + result.group(5)
        print "ValueError (SCF) in input file " + file + " on line " + str(lineNumber) + ": \"" + result.group(4) + "\" -> " + str(message)
        sys.exit(1)

      try:
        dD = float(result.group(5).lower().replace("d", "e"))
      except ValueError, message:
        print "anaylzing line " + str(lineNumber) + ": " + line.strip()
        print result.group(1) + " " + result.group(2) + " " + result.group(3) + " " + result.group(4) + " " + result.group(5)
        print "ValueError (dD) in input file " + file + " on line " + str(lineNumber) + ": \"" + result.group(5) + "\" -> " + str(message)
        sys.exit(1)
