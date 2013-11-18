#!/usr/bin/python

import sys, optparse, re

parser = optparse.OptionParser()

options, arguments = parser.parse_args()

if len(arguments) != 1:
  print("what file?", file = sys.stderr)

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
      print(geometry, end = ' ')
      for i in range(len(fMax)):
        print(fMax[i], end = ' ')
        if i == len(fMax)-1:
          print()
      fMax = []
      geometry += 1

    fMax.append(float(result.group(2).lower().replace("d", "e")))

  lastClone = clone
  line = fd.readline()

fd.close()
