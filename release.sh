#!/usr/bin/env bash
set -e  # Exit on error
ORIG_DIR="$(pwd)"
source activate geno

cd ../prstools
nbdev_bump_version
cd "$ORIG_DIR"
rsync -auv --exclude='.git/' --exclude-from='.gitignore' --existing ../prstools/ ./


git add -u
git commit -m "automatic release"
git push



