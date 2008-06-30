#if ! defined (__MONDOLOGGER_H)
#define __MONDOLOGGER_H 1

#include "MondoConfig.h"

#if MAX_LINE_LENGTH

#define LOG_PLAIN(msg)   CALL MondoLog(DEBUG_NONE, "", msg)
#define LOG_NONE(msg)    CALL MondoLog(DEBUG_NONE, __FILE__//":"//TRIM(IntToChar(__LINE__)), msg)
#define LOG_MINIMUM(msg) CALL MondoLog(DEBUG_MINIMUM, __FILE__//":"//TRIM(IntToChar(__LINE__)), msg)
#define LOG_MEDIUM(msg)  CALL MondoLog(DEBUG_MEDIUM, __FILE__//":"//TRIM(IntToChar(__LINE__)), msg)
#define LOG_CHKSUMS(msg) CALL MondoLog(DEBUG_CHKSUMS, __FILE__//":"//TRIM(IntToChar(__LINE__)), msg)
#define LOG_MAXIMUM(msg) CALL MondoLog(DEBUG_MAXIMUM, __FILE__//":"//TRIM(IntToChar(__LINE__)), msg)

#else

#define LOG_PLAIN(msg)   CALL MondoLogPlain(msg)
#define LOG_NONE(msg)    CALL MondoLogPlain(msg)
#define LOG_MINIMUM(msg) CALL MondoLogPlain(msg)
#define LOG_MEDIUM(msg)  CALL MondoLogPlain(msg)
#define LOG_CHKSUMS(msg) CALL MondoLogPlain(msg)
#define LOG_MAXIUM(msg)  CALL MondoLogPlain(msg)

#endif

#endif
