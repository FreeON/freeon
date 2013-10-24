/** @file
 *
 * The header file for the logging functions.
 *
 * @author Nicolas Bock <nicolas.bock@freeon.org>
 * @author Matt Challacombe <matt.challacombe@freeon.org>
 */

#ifndef __LOGGER_H
#define __LOGGER_H

#include "config.h"

#include <string>

void initializeLogger (void);

void logger (const char *const filename,
    const int linenumber,
    const char *const function_name,
    const char *const tag,
    const char *const format, ...);

/** A convenience macro for debugging messages. */
#ifdef DEBUG_OUTPUT
#define DEBUG(format, ...) logger(__FILE__, __LINE__, __func__, "", format, ##__VA_ARGS__)
#else
#define DEBUG(format, ...)
#endif

/** A convenience macro for info messages. */
#define INFO(format, ...) logger(__FILE__, __LINE__, __func__, "INFO", format, ##__VA_ARGS__)

/** A convenience macro for aborting messages. This macro will also terminate
 * the program. */
#define ABORT(format, ...) logger(__FILE__, __LINE__, __func__, "ERROR", format, ##__VA_ARGS__); printf("this is fatal\n"); exit(-1)

#endif
