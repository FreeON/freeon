#!/usr/bin/python
#
# vim: tw=0

import sys, re, optparse

parser = optparse.OptionParser(usage = "%prog [options] outputfile ...",
    description = """This script parses the output file(s) given on the command
line for the final SCF energy. The special filename '-' means standard input.
The output is written to standard output. There these columns:
(1) MDTime,
(4) MDEtot,
(3) SCF Energy,
(5) MDEkin,
(6) MDEpot,
(7) MDTemp,
(2) geometry,
(8) SCFCycle,
(9) dD, the difference in electronic density""")

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

  SCFEnergy = "0.0"
  dD        = "0.0"
  MDTime    = "0.0"
  MDEtot    = "0.0"
  MDEkin    = "0.0"
  MDEpot    = "0.0"
  MDTemp    = "0.0"

  lineNumber = 0

  for line in lines:
    lineNumber += 1

    result = re.compile("MDTime += ([-+.0-9dDeE]+)").search(line)
    if result:
      MDTime = result.group(1)

    result = re.compile("MDEtot[^=]+= ([-+.0-9dDeE]+)").search(line)
    if result:
      MDEtot = result.group(1)

    result = re.compile("MDEkin += ([-+.0-9dDeE]+)").search(line)
    if result:
      MDEkin = result.group(1)

    result = re.compile("MDEpot += ([-+.0-9dDeE]+)").search(line)
    if result:
      MDEpot = result.group(1)

    result = re.compile("T += ([-+.0-9dDeE]+)").search(line)
    if result:
      MDTemp = result.group(1)

    # Get SCF energy and other stuff.
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
        print MDTime, MDEtot, SCFEnergy, MDEkin, MDEpot, MDTemp, lastGeometry, lastSCFCycle, dD

      SCFEnergy = result.group(4)
      dD = result.group(5)
