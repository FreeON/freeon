#!/bin/csh

#
#  This script is an example of how to run the PHiPAC search scripts
#  with the "-custom" parameter. The output is simultaneously sent
#  to stderr and a file called "LOG" using the UNIX tee command.
#

set echo

perl ../../search-2.2/search.pl -custom \
	-machine machine_specs -ccopt compiler_options \
	-prec double -level 0 \
	|& tee -a LOG

# eof
