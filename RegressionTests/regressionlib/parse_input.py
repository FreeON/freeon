# Parse the input file.

import logging, re, sys

class parse_input:

  # The logger.
  log = None

  # The field values read from the input file.
  field = {}

  def __init__(self, testfile):
    self.log = logging.getLogger("parse_input")
    self.log.debug("reading from file " + testfile)

    # Read input file.
    fd = open(testfile, 'r')
    lines = fd.readlines()
    fd.close()

    # Parse through input file.
    linenumber = 1
    for line in lines:
      # Remove comments.
      line = line.strip()
      check = re.compile("(^[^#]*)(#.+$)").search(line)
      if check:
        line = check.group(1)

      # Search for known keywords.
      for fieldname in [ "input", "Mondo_executable", "Mondo_tar", \
          "configure_options" ]:
        check = re.compile("(^\s*" + fieldname + "\s*=)(.*$)").search(line)
        if check:
          if not check.group(2):
            self.log.error("syntax error on line " + str(linenumber))
            sys.exit(1)
          else:
            self.field[fieldname] = check.group(2).strip()
            self.log.debug("found input tag " + fieldname + " with value " + self.field[fieldname])

      linenumber += 1
