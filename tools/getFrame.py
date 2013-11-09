#!/usr/bin/python
#
# Get particular frames from an xyz file.

import optparse, re, sys

parser = optparse.OptionParser(usage = "%prog [options] FILE",
    description = """Loads coordinates from FILE and extracts select frames from
it. The extracted frames are written standard output in xyz format.""")

parser.add_option("-f", "--frame",
    metavar = "N[...]",
    help = """Load frame N from xyz file. Here are some examples of arguments:
(1) --frame=1,3,4-10,44 loads frames 1, 3, 4, 5, 6, 7, 8, 9, 10, 44. (2)
--frame=2/5 loads every 5th frame starting with frame 2.""",
    dest = "frameString",
    type = "string")

parser.add_option("--last",
    metavar = "N",
    help = "Load the last N frames",
    dest = "last",
    type = "int",
    default = -1)

options, arguments = parser.parse_args()

if len(arguments) == 0:
  print("missing xyz file...")
  print()
  parser.print_help()
  sys.exit(1)

if len(arguments) > 1:
  print("will only read from first file", file = sys.stderr)

# Parse frame argument.
frameIDs = []
if options.frameString:
  commas = options.frameString.split(",")
  for comma in commas:
    # Check for a "-" or a "/".
    if re.compile("-").search(comma) and re.compile("/").search(comma):
      print("syntax error in --frame argument, " + comma)
      sys.exit(1)

    if re.compile("-").search(comma):
      token = comma.split("-")
      if len(token) == 2:
        if int(token[0]) >= int(token[1]):
          print("syntax error in --frame argument, " + comma)
          sys.exit(1)
        for i in range(int(token[0]), int(token[1])+1):
          frameIDs.append(i)
      else:
        print("syntax error in --frame argument, " + comma)
        sys.exit(1)

    elif re.compile("/").search(comma):
      token = comma.split("/")
      if len(token) == 2:
        for i in range(int(token[0]), 1000, int(token[1])):
          frameIDs.append(i)
      else:
        print("syntax error in --frame argument, " + comma)
        sys.exit(1)

    else:
      frameIDs.append(int(comma))

if len(frameIDs) > 0:
  print("extracting frames " + str(frameIDs), file = sys.stderr)
if options.last >= 0:
  print("extracting last " + str(options.last) + " frames", file = sys.stderr)

# Read the xyz file.
fd = open(arguments[0], "r")
lines = fd.readlines()
fd.close()

# Count the total number of frames.
totalFrames = 0
lineNumber = 0
while lineNumber < len(lines):
  # Update total number of frames.
  totalFrames += 1

  # Read in size of frame.
  N = int(lines[lineNumber])
  lineNumber += 1
  comment = lines[lineNumber].rstrip()
  lineNumber += 1

  # Advance in file to next frame.
  lineNumber += N

print("counted " + str(totalFrames) + " total frames", file = sys.stderr)

frames = []
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
  if (frameNumber in frameIDs) or (frameNumber > totalFrames-options.last):
    frames.append({})
    frames[-1]["comment"] = comment
    frames[-1]["coordinates"] = []
    print("reading frame " + str(frameNumber), file = sys.stderr)
    for i in range(N):
      frames[-1]["coordinates"].append(lines[lineNumber].rstrip())
      lineNumber += 1
  else:
    lineNumber += N

print("parsed through " + str(frameNumber) + " frames", file = sys.stderr)

# Output frames.
for frame in frames:
  print(len(frame["coordinates"]))
  print(frame["comment"])
  for line in frame["coordinates"]:
    print(line)
