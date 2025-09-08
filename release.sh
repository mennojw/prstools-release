#!/bin/bash -i
source ~/.bashrc
set -e  # Exit on error 
ORIG_DIR="$(pwd)"
source ~/cmd/_geno
#conda activate geno
echo 'THEREHEHEHTHTHE'

cd ../prstools
nbdev_bump_version
cd "$ORIG_DIR"
rsync -auv --exclude='.git/' --exclude-from='.gitignore' --existing ../prstools/ ./

git add -u
git commit -m "automatic release"
git push
nbdev_pypi
pip install -U prstools
pytest --pyargs prstools



