if [ $# -ne 2 ]; then
  echo "usage: bash $0 server port"
  echo "example: bash $0 mammoth-1 9123"
  exit 1
fi

server=$1
port=$2

for i in `seq 1 40`; do
  nfwq start_condor_worker --server $server --port $port --client pulsar-$i --num_processes 6 --frac_mem 0.90 --timeout 60s
done

# These have 64G memory each.
for i in 35 36 37 `seq 39 45` 47 48; do
  nfwq start_condor_worker --server $server --port $port --client quasar-$i --num_processes 6 --frac_mem 0.90 --timeout 60s
done

# These have 128G memory each.
for i in 46 49 50 51 52; do
  nfwq start_condor_worker --server $server --port $port --client quasar-$i --num_processes 10 --frac_mem 0.90 --timeout 60s
done

nfwq start_condor_worker --server $server --port $port --client babel --num_processes 8 --frac_mem 0.90 --timeout 60s

for client in broman-10 broman-9; do
  nfwq start_condor_worker --server $server --port $port --client $client --num_processes 20 --frac_mem 0.85 --frac_cpu 0.75 --timeout 60s
done
