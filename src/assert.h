#include <R.h>

/* Modified from assert.h included with max os x */

#undef assert

#ifdef NDEBUG
#define assert(e) ((void)0)
#else
#define macroevolution_assert(e, file, line) \
    error("%s:%u: failed assertion `%s'\n", file, line, e)
#define assert(e) \
    ((void) ((e) ? ((void)0) : macroevolution_assert(#e, __FILE__, __LINE__)))
#endif
