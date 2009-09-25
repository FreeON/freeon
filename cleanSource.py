#!/usr/bin/python
#
# Call this with a list of files to consider.

import optparse, sys, re, tempfile

parser = optparse.OptionParser(description = """This script cleans a source
file. It does (1) remove trailing white space characters, and (2) replace a
TAB character with space characters in cases where the TAB character follows a
leading space character.""")

parser.add_option("--pretend",
    help = "Just check, but not change any files",
    action = "store_true",
    default = False,
    dest = "pretend")

parser.add_option("--tab",
    help = "Replace a TAB character with N spaces (default N = 2)",
    default = 2,
    type = "int",
    dest = "tab",
    metavar = "N")

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
