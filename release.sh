#!/bin/bash -i
source ~/.bashrc
set -e  # Exit on error 
ORIG_DIR="$(pwd)"
source ~/cmd/_geno
#conda activate geno
echo 'CREATING NEW PRSTOOLS RELEASE'

## Some notes:
# This release.sh code pulls code from private dev-repo and pushes it to release github and pypi
# and immediately tests it.
# Later a dynamic with a bleeding edge version needs to be incorporated.
# Also a small benchmark experiment should be part of the tests

cd ../prstools
nbdev_bump_version
cd "$ORIG_DIR"
rsync -auv --exclude='.git/' --exclude-from='.gitignore' --existing ../prstools/ ./

git add -u
git commit -m "automatic release"
git push
nbdev_pypi
pip install -U prstools
pytest -v -s --pyargs prstools



