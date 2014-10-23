set -o errexit
set -o nounset
set -o pipefail

out=$1
cmd=$2
shift 2
$cmd "$@" | tee $out
