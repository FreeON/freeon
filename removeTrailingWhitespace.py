#!/usr/bin/python
#
# Call this with a list of files to consider.

import optparse, sys, re, tempfile

parser = optparse.OptionParser()

parser.add_option("--pretend",
    help = "Just check, but not change any files",
    action = "store_true",
    default = False,
    dest = "pretend")

options, arguments = parser.parse_args()

if len(arguments) == 0:
  print("what files should I fix?")
  sys.exit(1)

for file in arguments:
  try:
    fd = open(file)
    lines = fd.readlines()
    fd.close()

  except:
    print("error opening file " + file)
    sys.exit(1)

  lineNumber = 0
  fileNeedsFixing = False
  for line in lines:
    lineNumber += 1
    if re.compile("[ \t]+$").search(line):
      fileNeedsFixing = True
      print("trailing whitespace in file " + file + " on line " + str(lineNumber))

  if not options.pretend and fileNeedsFixing:
    print("fixing file " + file)
    fd = open(file, "w")
    for line in lines:
      line = re.sub("\s+$", "", line)
      print >> fd, line
    fd.close()
