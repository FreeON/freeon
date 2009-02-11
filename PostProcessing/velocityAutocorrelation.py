#!/usr/bin/python
#
# Nicolas Bock <nbock@lanl.gov>

import optparse
import re
import sys

parser = optparse.OptionParser()

parser.add_option("--t0",
    metavar = "t0",
    help = "Start calculating velocity autocorrelation at time t0",
    type = "float",
    dest = "t0")

parser.add_option("--tf",
    metavar = "tf",
    help = "Stop calculating velocity autocorrelation at time tf",
    type = "float",
    dest = "tf")

parser.add_option("--dt",
    metavar = "dt",
    help = "break velocity autocorrelation function into pieces dt long",
    type = "float",
    dest = "dt")

parser.set_defaults(t0 = 0.0)

(options, arguments) = parser.parse_args()

if len(arguments) == 0:
  print("what xyz file should I read from?")
  print("")
  parser.print_help()
  sys.exit(1)

time = []
position = []
velocity = []

for file in arguments:
  print >> sys.stderr, "reading " + file

  try:
    fd = open(file)
    lines = fd.readlines()
    fd.close()

  except IOError, e:
    print >> sys.stderr, "error opening file: " + e.strerror
    sys.exit(1)

  numberOfAtoms = int(lines[0].strip())
  print >> sys.stderr, "file contains information on " + str(numberOfAtoms) + " atoms"

  lineNumber = 0
  while(lineNumber < len(lines)):
    position.append([])
    velocity.append([])

    N = int(lines[lineNumber].strip())
    if N != numberOfAtoms:
      print >> sys.stderr, "line " + str(lineNumber+1) + ": number of atoms is different from previous values (" + str(N) + ")"
      sys.exit(1)

    lineNumber += 1
    comment = lines[lineNumber].strip()

    result = re.compile("t = ([-+.0-9eEdD]+)").search(comment)
    if result:
      timestring = result.group(1)
      timestring = timestring.lower().replace("d", "e")
      try:
        time.append(float(timestring))

      except ValueError, e:
        print >> sys.stderr, "can not convert time value on line " + str(lineNumber+1) + ": " + comment
        print >> sys.stderr, e.message
        sys.exit(1)

    else:
      time.append(float(len(time)))

    lineNumber += 1

    for i in range(numberOfAtoms):
      tokens = lines[lineNumber+i].split()
      if len(tokens) < 7:
        print >> sys.stderr, str(lineNumber+i+1) + ": I need positions and velocities on each line"
        sys.exit(1)

      try:
        pos = [ float(tokens[1]), float(tokens[2]), float(tokens[3]) ]
        vel = [ float(tokens[4]), float(tokens[5]), float(tokens[6]) ]

      except ValueError, e:
        print >> sys.stderr, "can not convert line " + str(lineNumber+i+1) + ": " + lines[lineNumber+i+1].strip()
        print(e.message)
        sys.exit(1)

      position[-1].append(pos)
      velocity[-1].append(vel)

    lineNumber += numberOfAtoms

print >> sys.stderr, "read " + str(len(time)) + " configurations"

# Find index of t0.
t0Index = len(time)
for i in range(len(time)):
  if time[i] > options.t0:
    t0Index = i-1
    break

if t0Index == len(time):
  print >> sys.stderr, "t0 too large"
  sys.exit(1)

# Find index of tf.
tfIndex = 0
if options.tf:
  for i in range(len(time)):
    if time[i] > options.tf:
      tfIndex = i
      break

  if tfIndex == 0:
    print >> sys.stderr, "tf too small"
    sys.exit(1)

else:
  tfIndex = len(time)

if t0Index >= tfIndex:
  print >> sys.stderr, "syntax error: tf <= t0, " + str(options.tf) + " <= " + str(options.t0)
  sys.exit(1)

# Calculate correlation function.
C = [[]]
for i in range(t0Index, tfIndex):
  C[-1].append(0)
  for n in range(len(velocity[i])):
    try:
      C[-1][-1] += velocity[t0Index][n][0]*velocity[i][n][0] \
                +  velocity[t0Index][n][1]*velocity[i][n][1] \
                +  velocity[t0Index][n][2]*velocity[i][n][2]

    except IndexError, e:
      print >> sys.stderr, "index error: " + e.message
      print >> sys.stderr, "i       = " + str(i)
      print >> sys.stderr, "n       = " + str(n)
      print >> sys.stderr, "t0Index = " + str(t0Index)
      print >> sys.stderr, "tfIndex = " + str(tfIndex)
      sys.exit(1)

  C[-1][-1] /= float(len(velocity[i]))

  print time[i], C[-1][-1]
