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

pip install -e ./
cd ../prstools
nbdev_bump_version
TODAY=$(date +%d-%m-%Y)
sed -i '' "s/^_date = .*/_date = \"$TODAY\"/" prstools/__init__.py
cd "$ORIG_DIR"
rsync -auv --exclude='.git/' --exclude-from='.gitignore' --existing ../prstools/ ./
cd ./prstools
python ./_cmd.py --dev
cd "$ORIG_DIR"

# Commiting and push to github & pypi
git add -u
git commit -m "automatic release"
git push
nbdev_pypi

# Test in a clean environment
echo ">>> Testing in clean conda env"
mamba create -y -n prstools_test python=3.11
source activate prstools_test
pip install -U prstools
pytest -v -s --pyargs prstools
conda deactivate
mamba env remove -n prstools_test -y

#pip uninstall -y prstools
#pip install -U prstools
#pytest -v -s --pyargs prstools



