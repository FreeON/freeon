# Parse the input file.

import logging

class parse_input:

  # The logger.
  log = None

  def __init__(self, testfile):
    self.log = logging.getLogger("parse_input")
    self.log.debug("reading from file " + testfile)
