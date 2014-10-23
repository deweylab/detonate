#!/bin/bash
set -o nounset
set -o pipefail
set -o errexit

if [ ! -f config.sh ]; then
  echo "This script should be run from within a check out of the detonate repo."
  exit 1;
fi

# Read the RSEM-EVAL and REF-EVAL versions, and set the DETONATE version.
. config.sh
pkg=detonate-$version

# Init the package directory.
rm -rf $pkg
mkdir $pkg

# Copy RSEM-EVAL into the package.
cp -rp rsem-eval $pkg/
rm -rf $pkg/rsem-eval/.git

# Copy REF-EVAL into the package.
cp -rp ref-eval $pkg/
rm -rf $pkg/ref-eval/.git

# Record the versions used.
echo "DETONATE version is $pkg" > $pkg/VERSION

# Make/copy the documentation.
bash make_doc.sh $pkg

# Copy the makefile.
cp -p Makefile $pkg/

# Make the tarball.
tar cvfz $pkg.tar.gz $pkg/

# Clean up.
#rm -rf $pkg
echo "DETONATE version $version has been saved to $pkg.tar.gz"
