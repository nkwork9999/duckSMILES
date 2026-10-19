#pragma once

// DuckDB 2 also force-includes this shim in generated C entry points.
#if defined(_MSC_VER) && defined(__cplusplus)
#include <iterator>
#include <string>

// DuckDB v1.5.3 vendors an older fmt release that uses
// stdext::checked_array_iterator whenever _SECURE_SCL is defined. VS 18
// still defines the compatibility macro but no longer provides stdext. Include
// core STL headers first so they cannot define _SECURE_SCL again while fmt is
// being parsed.
#ifdef _SECURE_SCL
#undef _SECURE_SCL
#endif
#endif
