#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

if [ $# -ne 1 ]; then
  echo "usage: $0 BOOST_DIR"
  echo "BOOST_DIR should be a directory containing the full unpacked boost source."
  echo "example: $0 ~/Downloads/boost_1_55_0"
  exit 1
fi

boost=$1

# Boostrap boost and build the bcp tool.
pushd $boost
./bootstrap.sh
./bjam tools/bcp
popd

# Extract boost.
tmp=$(mktemp -d boost.XXXXXX)
$boost/dist/bin/bcp --boost=$boost --scan $(ls *.cpp *.hh) $tmp
$boost/dist/bin/bcp --boost=$boost --scan $(find deweylab/ -name '*.hh') $tmp
$boost/dist/bin/bcp --boost=$boost --scan $(find deweylab/ -name '*.cc') $tmp
$boost/dist/bin/bcp --boost=$boost build $tmp
rm -rf $tmp/tools/build/v2/engine/bootstrap/
rm -rf $tmp/tools/build/v2/engine/bin.linuxx86_64/
rm -rf $tmp/tools/build/v2/engine/bin.macosxx86_64/

# Move newly extracted files to boost directory, and move the old boost
# directory elsewhere.
if [ -e boost ]; then
  old=$(mktemp -d old_boost.XXXXXX)
  mv boost $old
  echo "Moved old ./boost to $old/boost"
fi
mv $tmp boost
echo "Moved newly extracted boost to ./boost"

# Restored our own files.
git checkout -- boost/.gitignore
git checkout -- boost/Makefile
echo "Restored boost/.gitignore and boost/Makefile"
