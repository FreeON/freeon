#!/usr/bin/python

import sys, optparse, re

parser = optparse.OptionParser()

options, arguments = parser.parse_args()

if len(arguments) != 1:
  print >> sys.stderr, "what file?"

fd = open(arguments[0])
line = fd.readline()

geometry = 1
clone = -1
lastClone = -1
fMax = []

while line:
  result = re.compile("Clone ([0-9]+).+fMax = ([^ ]+)").search(line)
  if result:
    clone = int(result.group(1))
    if clone < lastClone:
      print geometry,
      for i in range(len(fMax)):
        print fMax[i],
        if i == len(fMax)-1:
          print
      fMax = []
      geometry += 1

    fMax.append(float(result.group(2).lower().replace("d", "e")))

  lastClone = clone
  line = fd.readline()

fd.close()
