#!/bin/sh
set -eu
root=$(CDPATH= cd -- "$(dirname -- "$0")/../../.." && pwd)
platform=${PLATFORM_NAME:-iphonesimulator}
case "$platform" in
  iphoneos) target=aarch64-apple-ios ;;
  iphonesimulator) target=aarch64-apple-ios-sim ;;
  *) echo "Unsupported platform: $platform" >&2; exit 1 ;;
esac
if [ "${CHEMOB_USE_PREBUILT_ENGINE:-0}" = "1" ]; then
  test -f "$root/apps/Chemob/Generated/$platform/libducksmiles_smiles.a" || { echo "Build the engine before reusing it" >&2; exit 1; }
  echo "Using the previously built DuckSMILES library for UI iteration"
  exit 0
fi
export PATH="$HOME/.cargo/bin:/opt/homebrew/bin:$PATH"
export CARGO_TARGET_DIR="${CHEMOB_CARGO_TARGET_DIR:-/private/tmp/chemob-apple-engine-target}"
export SDKROOT
SDKROOT=$(xcrun --sdk "$platform" --show-sdk-path)
cargo build --manifest-path "$root/Cargo.toml" -p ducksmiles_smiles --release --target "$target"
mkdir -p "$root/apps/Chemob/Generated/$platform"
cp "$CARGO_TARGET_DIR/$target/release/libducksmiles_smiles.a" "$root/apps/Chemob/Generated/$platform/"
