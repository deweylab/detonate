#!/bin/bash
set -o nounset
set -o pipefail
set -o errexit

if [[ "$#" -eq 1 ]]; then
  pkg=$1
else
  pkg=.
fi

if [[ ! -d rsem-eval && ! -d ref-eval ]]; then
  echo "This script should be run from within a check out of the detonate repo."
  exit 1;
fi

mkdir -p $pkg/html
python3 $pkg/ref-eval/make_doc.py \
  --template index.template.html \
  --html     $pkg/html/index.html \
  --text     $pkg/README \
  --cxx      /dev/null

python3 $pkg/ref-eval/make_doc.py \
  --template vignette.template.html \
  --html     $pkg/html/vignette.html \
  --text     $pkg/VIGNETTE \
  --cxx      /dev/null

pushd $pkg/ref-eval; make doc; popd
cp $pkg/ref-eval/ref-eval.html $pkg/html/
cp $pkg/ref-eval/ref-eval-estimate-true-assembly.html $pkg/html/

cat > rsem-eval.template.html <<-EOF
<?xml version="1.0"?>
<html>
<head>
<title>RSEM-EVAL: A novel reference-free transcriptome assembly evaluation measure</title>
</head>
<body>
<h1>RSEM-EVAL: A novel reference-free transcriptome assembly evaluation measure</h1>
EOF
perl Markdown_1.0.1/Markdown.pl < $pkg/rsem-eval/README.md | \
  awk '/h2.*Introduction/ { f=1 } f' >> rsem-eval.template.html
cat >> rsem-eval.template.html <<-EOF
</body>
</html>
EOF
python3 $pkg/ref-eval/make_doc.py \
  --template rsem-eval.template.html \
  --html     $pkg/html/rsem-eval.html \
  --text     /dev/null \
  --cxx      /dev/null
