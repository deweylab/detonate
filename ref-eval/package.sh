#!/bin/bash
set -o nounset
set -o pipefail
set -o errexit

if [ ! -f ref-eval.cpp ]; then
  echo "This script should be run from within a check out of the detonate repo."
  exit 1;
fi

version=`date "+%Y%m%d"`
pkg=ref-eval-$version

rm -rf $pkg
git clone . $pkg
rm -rf $pkg/.git
tar cvfz $pkg.tar.gz $pkg
rm -rf $pkg
mv $pkg.tar.gz ../

echo "REF-EVAL version $version has been saved to ../$pkg.tar.gz"
