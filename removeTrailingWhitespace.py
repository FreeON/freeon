#!/usr/bin/python
#
# Call this with a list of files to consider.

import optparse, sys, re, tempfile

parser = optparse.OptionParser()

options, arguments = parser.parse_args()

if len(arguments) == 0:
  print("what files should I fix?")
  sys.exit(1)

for file in arguments:
  fd = open(file)
  lines = fd.readlines()
  fd.close()

  lineNumber = 0
  fileNeedsFixing = False
  for line in lines:
    lineNumber += 1
    if re.compile("[ \t]+$").search(line):
      fileNeedsFixing = True
      print("trailing whitespace in file " + file + " on line " + str(lineNumber) + ": " + line.strip())

  if fileNeedsFixing:
    print("fixing file " + file)
    fd = open(file, "w")
    for line in lines:
      line = re.sub("[ \t\n]+$", "", line)
      print >> fd, line
    fd.close()
