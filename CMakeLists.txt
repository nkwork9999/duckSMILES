cmake_minimum_required(VERSION 3.5...3.29)

set(TARGET_NAME ducksmiles)
set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)

project(${TARGET_NAME})

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ---------------------------------------------------------------------------
# Build Rust static library via cargo
# ---------------------------------------------------------------------------
set(RUST_DIR ${CMAKE_CURRENT_SOURCE_DIR})

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(CARGO_PROFILE release)
    set(CARGO_FLAG --release)
else()
    set(CARGO_PROFILE debug)
    set(CARGO_FLAG "")
endif()

if(DEFINED ENV{CARGO_TARGET})
    set(CARGO_TARGET_DIR "target/$ENV{CARGO_TARGET}/${CARGO_PROFILE}")
    set(CARGO_TARGET_FLAG "--target" "$ENV{CARGO_TARGET}")
else()
    set(CARGO_TARGET_DIR "target/${CARGO_PROFILE}")
    set(CARGO_TARGET_FLAG "")
endif()

set(RUST_STATIC_LIBS
    "${RUST_DIR}/${CARGO_TARGET_DIR}/libducksmiles_smiles.a"
    "${RUST_DIR}/${CARGO_TARGET_DIR}/libducksmiles_inchi.a"
    "${RUST_DIR}/${CARGO_TARGET_DIR}/libducksmiles_mol.a"
    "${RUST_DIR}/${CARGO_TARGET_DIR}/libducksmiles_pdb.a"
    "${RUST_DIR}/${CARGO_TARGET_DIR}/libducksmiles_selfies.a"
)

add_custom_command(
    OUTPUT ${RUST_STATIC_LIBS}
    COMMAND cargo build ${CARGO_FLAG} ${CARGO_TARGET_FLAG} --workspace
    WORKING_DIRECTORY ${RUST_DIR}
    COMMENT "Building Rust workspace (duckSmiles crates)"
    VERBATIM
)

add_custom_target(rust_lib DEPENDS ${RUST_STATIC_LIBS})

# ---------------------------------------------------------------------------
# C++ extension sources
# ---------------------------------------------------------------------------
set(EXTENSION_SOURCES
    src/ducksmiles_extension.cpp
)

include_directories(src/include)

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})

add_dependencies(${EXTENSION_NAME} rust_lib)
add_dependencies(${LOADABLE_EXTENSION_NAME} rust_lib)

foreach(RUST_LIB ${RUST_STATIC_LIBS})
    target_link_libraries(${EXTENSION_NAME} ${RUST_LIB})
    target_link_libraries(${LOADABLE_EXTENSION_NAME} ${RUST_LIB})
endforeach()

if(APPLE)
    target_link_libraries(${EXTENSION_NAME} "-framework Security" "-framework CoreFoundation")
    target_link_libraries(${LOADABLE_EXTENSION_NAME} "-framework Security" "-framework CoreFoundation")
endif()

install(
  TARGETS ${EXTENSION_NAME}
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")
