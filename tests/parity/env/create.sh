#!/usr/bin/env bash
# Create the pinned conda env used to run the Panaroo parity reference.
#   ./create.sh            create (or update) the env
#   conda env remove -n panaroo-parity     to undo
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
conda env create -f "$here/environment.yml" --yes 2>/dev/null \
  || conda env update -f "$here/environment.yml" --prune
echo
echo "created. record the resolved versions:"
conda list -n panaroo-parity --explicit > "$here/environment.lock.txt"
echo "  $here/environment.lock.txt"
