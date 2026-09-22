#!/usr/bin/env bash
# Fill the recipe's sha256 from the published GitHub release tarball (run after pushing tag vX.Y.Z).
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
recipe="$here/rnaseq-pipeline/meta.yaml"
version="$(sed -n 's/^{% set version = "\(.*\)" %}$/\1/p' "$recipe")"
url="https://github.com/vipin839/rnaseq_pipeline/archive/refs/tags/v${version}.tar.gz"
sha="$(curl -fsSL "$url" | sha256sum | cut -d' ' -f1)"
sed -i "s/^  sha256: .*/  sha256: ${sha}/" "$recipe"
echo "v${version}: sha256 ${sha} written to ${recipe}"
