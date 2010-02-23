#!/usr/bin/python
#
# Get particular frames from an xyz file.

import optparse, sys

parser = optparse.OptionParser(usage = "%prog [options] FILE",
    description = """Loads coordinates from FILE and extracts select frames from
it. The extracted frames are written standard output in xyz format.""")

parser.add_option("-f", "--frame",
    metavar = "N[...]",
    help = """Load frame N from xyz file. The following argument,
--frames=1,3,4-10,44 would load frames 1, 3, 4, 5, 6, 7, 8, 9, 10, 44.""",
    dest = "frameString",
    type = "string")

options, arguments = parser.parse_args()

if not options.frameString:
  options.frameString = "1"

# Parse frame argument.
frameIDs = []
commas = options.frameString.split(",")
for comma in commas:
  dash = comma.split("-")
  if len(dash) == 1:
    frameIDs.append(int(dash[0]))
  elif len(dash) == 2:
    if int(dash[0]) >= int(dash[1]):
      print "syntax error in --frame argument, " + comma
      sys.exit(1)
    for i in range(int(dash[0]), int(dash[1])+1):
      frameIDs.append(i)
  else:
    print "syntax error in --frame argument, " + comma
    sys.exit(1)

print >> sys.stderr, "extracting frames " + str(frameIDs)

frames = {}

if len(arguments) == 0:
  print "missing xyz file..."
  print
  parser.print_help()
  sys.exit(1)

if len(arguments) > 1:
  print >> sys.stderr, "will only read from first file"

fd = open(arguments[0], "r")
lines = fd.readlines()
fd.close()

frameNumber = 0
lineNumber = 0
while lineNumber < len(lines):
  # Update frame number.
  frameNumber += 1

  # Read in size of frame.
  N = int(lines[lineNumber])
  lineNumber += 1
  comment = lines[lineNumber].rstrip()
  lineNumber += 1

  # Read frame.
  if frameNumber in frameIDs:
    frames[frameNumber] = {}
    frames[frameNumber]["comment"] = comment
    frames[frameNumber]["coordinates"] = []
    print >> sys.stderr, "reading frame " + str(frameNumber)
    for i in range(N):
      frames[frameNumber]["coordinates"].append(lines[lineNumber].rstrip())
      lineNumber += 1
  else:
    lineNumber += N

# Output frames.
for frameID in frames:
  print len(frames[frameID])
  print frames[frameID]["comment"]
  for line in frames[frameID]["coordinates"]:
    print line
