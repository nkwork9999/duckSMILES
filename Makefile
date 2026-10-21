PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))
EXT_NAME=ducksmiles
EXT_CONFIG=${PROJ_DIR}extension_config.cmake
include extension-ci-tools/makefiles/duckdb_extension.Makefile

# Native prebuilt libraries can carry a staging version for the pinned commit.
# Append after the shared include so its defaults cannot erase this flag.
ifneq ($(strip $(DUCKDB_PREBUILT_LIBRARY)),)
ducksmiles_shell_quote = '$(subst ','"'"',$(1))'
DUCKSMILES_PREBUILT_VERSION := $(strip $(shell python3 $(call ducksmiles_shell_quote,$(PROJ_DIR)scripts/prebuilt_duckdb_version.py) --library $(call ducksmiles_shell_quote,$(value DUCKDB_PREBUILT_LIBRARY)) --source-dir $(call ducksmiles_shell_quote,$(DUCKDB_SRCDIR)) --requested-version $(call ducksmiles_shell_quote,$(DUCKDB_VERSION))))
ifeq ($(DUCKSMILES_PREBUILT_VERSION),)
$(error Unable to verify the native prebuilt DuckDB version)
endif
override EXT_RELEASE_FLAGS += -DDUCKDB_EXPLICIT_VERSION=$(DUCKSMILES_PREBUILT_VERSION)
endif
