if(MSVC)
    # DuckDB 1.5.x defaults core targets to C++11; the vendored fmt needs C++17
    # under recent MSVC.
    if(NOT DEFINED CMAKE_CXX_STANDARD OR CMAKE_CXX_STANDARD LESS 17)
        set(CMAKE_CXX_STANDARD 17 CACHE STRING "C++ standard to enforce" FORCE)
        set(CMAKE_CXX_STANDARD_REQUIRED ON CACHE BOOL "Require the configured C++ standard" FORCE)
    endif()
    # MSVC 14.51+ (VS 18) removed stdext::checked_array_iterator but still
    # defines _SECURE_SCL, breaking DuckDB v1.5.3's vendored fmt
    # (duckdb/duckdb#22704). Force-include a shim that undefines the macro
    # before fmt headers are parsed. Drop once the pinned DuckDB has the fix.
    if(MSVC_VERSION GREATER_EQUAL 1950)
        add_compile_options("/FI${CMAKE_CURRENT_LIST_DIR}/cmake/msvc_fmt_compat.hpp")
    endif()
endif()

duckdb_extension_load(ducksmiles
    SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}
    LOAD_TESTS
)
