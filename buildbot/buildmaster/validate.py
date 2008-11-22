import re
from buildbot.steps.shell import Test
from buildbot.status.builder import SUCCESS, FAILURE

class Validate(Test):

  name = "validate"
  description = ["validating"]
  descriptionDone = ["validate"]
  command = ["make", "validate"]

  def evaluateCommand(self, cmd):

    # Reset counters.
    totalCount = 0
    failedCount = 0
    passedCount = 0
    warningsCount = 0

    # Read input.
    lines = self.getLog('stdio').readlines()

    warnings = []
    #warnings.append("analyzing " + str(len(lines)) + " lines of output")

    for line in lines:
      # Find errors.
      #warnings.append("anaylzing line " + line.strip())
      result = re.compile("done analyzing.+([0-9]+) errors.+([0-9]+) unmatched tags").search(line)
      if result != None:
        warnings.append("errors: " + str(int(result.group(1))))
        warnings.append("unmatched tags: " + str(int(result.group(2))))
        totalCount += 1
        if int(result.group(1)) > 0:
          failedCount += 1
        else:
          passedCount += 1

        if int(result.group(2)) > 0:
          warningsCount += int(result.group(2))

    # Set the results.
    self.setTestResults(total = totalCount,
                        failed = failedCount,
                        passed = passedCount,
                        warnings = warningsCount)
    self.addCompleteLog("warnings", "\n".join(warnings) + "\n")

    # Return
    if failedCount > 0:
      return FAILURE
    else:
      return SUCCESS
