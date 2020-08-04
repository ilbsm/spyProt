#!/bin/bash
# Resolve local directory of the script
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" && pwd )"


# Clean
rm -rf $DIR/venv
rm -rf $DIR/dist
rm -rf $DIR/_skbuild

# Setup environment
python3.8 -m venv $DIR/venv
source $DIR/venv/bin/activate
pip install wheel pytest scikit_build
pip install -r $DIR/requirements.txt
# Run tests
cd $DIR
pytest
# Build package
python setup.py bdist_wheel
