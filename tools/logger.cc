/** @file
 *
 * The implementation of the logging functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#include "config.h"
#include "logger.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/** @internal The time the logger was initialized (hopefully the same as the
 * time the main program started).
 */
struct timespec loggerStartTime;

/** Initialize the logger. This function should be called very early in the
 * main method, for example inside an initnode entry method.
 */
void initializeLogger (void)
{
  if(clock_gettime(CLOCKTYPE, &loggerStartTime) < 0)
  {
    printf("[%s:%d] can not get time\n", __FILE__, __LINE__);
    exit(1);
  }
}

/** The main logging function. It is usually more convenient to use the macros
 * DEBUG, INFO, and ABORT instead.
 *
 * @param filename The name of the file where the message occurred.
 * @param linenumber The line number in that file.
 * @param function_name The name of the function.
 * @param tag A tag.
 * @param format The format string. See printf() for details.
 */
void logger (const char *const filename,
    const int linenumber,
    const char *const function_name,
    const char *const tag,
    const char *const format, ...)
{
  const int format_length = 2000;
  const int string_length = 5000;

  va_list ap;
  char new_format[format_length];
  char output_string[string_length];

  struct timespec nowTime;

  if(clock_gettime(CLOCKTYPE, &nowTime) < 0)
  {
    printf("[%s:%d] can not get time\n", __FILE__, __LINE__);
    exit(1);
  }
  double walltime = (nowTime.tv_sec+nowTime.tv_nsec/1.0e9)
    -(loggerStartTime.tv_sec+loggerStartTime.tv_nsec/1.0e9);

  if(strlen(tag) > 0)
  {
    snprintf(new_format, format_length, "[%s:%d (%s) %1.2fs %s] %s", filename,
        linenumber, function_name, walltime, tag, format);
  }

  else
  {
    snprintf(new_format, format_length, "[%s:%d (%s) %1.2fs] %s", filename,
        linenumber, function_name, walltime, format);
  }

  va_start(ap, format);
  vsnprintf(output_string, string_length, new_format, ap);
  va_end(ap);
  printf("%s", output_string);
}
