#!/usr/bin/python
#
# Convert an MDO output file (what you get when running an MD simulation) into
# an xyz file.

import optparse, sys, re

parser = optparse.OptionParser(usage = "%prog [options] file.MDO [file.xyz]")

option, arguments = parser.parse_args()

if len(arguments) not in [1, 2]:
  print >> sys.stderr, "wrong number of arguments"
  parser.print_help()
  sys.exit(1)

input_filename = arguments[0]

if len(arguments) < 2:
  output_filename = input_filename + ".xyz"
else:
  output_filename = arguments[1]

print >> sys.stderr, "opening MDO file " + input_filename \
  + " and write output into " + output_filename

file = open(input_filename)
lines = file.readlines()
file.close()

line_number = 0
block_started = False
frames = []

for line in lines:
  line_number += 1
  result_start = re.compile("^Start: Atom Position").search(line)
  result_end = re.compile("^End: Atom Position").search(line)

  if result_start != None:
    if not block_started:
      block_started = True
      frames.append([])
      continue
    else:
      print >> sys.stderr, "line " + str(line_number) + ": block already started"
      sys.exit(1)

  if result_end != None:
    if block_started:
      block_started = False

  if block_started:
    token = line.split()
    frames[-1].append([token[0], float(token[1]), float(token[2]), float(token[3])])

print >> sys.stderr, "read " + str(len(frames)) + " frames"

file = open(output_filename, 'w')
frame_counter = 0
for frame in frames:
  frame_counter += 1
  print >> file, len(frame)
  print >> file, "MDO frame " + str(frame_counter)
  for atom in frame:
    print >> file, "%s %f %f %f" % (atom[0], atom[1], atom[2], atom[3])

file.close()
