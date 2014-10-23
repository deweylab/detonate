import os, sys, nfwq, work_queue, argparse, assemblies_meta
from collections import namedtuple

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True)
args = p.parse_args()

# Figure out where the input and output files should live.
indir = os.path.join(args.workdir, "detonate_assemblies")
outdir = os.path.join(args.workdir, "ref_eval_scores")

# Figure out where the executables live.
script_dir = os.path.dirname(os.path.abspath(__file__))
wrapper = script_dir + "/wrapper.sh"
ref_eval = os.path.dirname(script_dir) + "/ref-eval"

# Initialize the DAG. We'll add our jobs here.
dag = nfwq.Dag(args.workdir + "/dag.db")
dag.init()

Reads = namedtuple("Reads", ("id", "ss", "readlen", "num_reads"))
reads = [
  Reads(id="ensembl_sim_oases",  ss=False, readlen=76, num_reads=24130103),
  Reads(id="oases_mouse_real",   ss=False, readlen=76, num_reads=41612875),
  Reads(id="real",               ss=True,  readlen=76, num_reads=52402378),
  Reads(id="trinity_yeast_real", ss=True,  readlen=68, num_reads=49999986)
]

for r in reads:
  if r.ss:
    assemblies = assemblies_meta.ss_assemblies_meta(r.id)
  else:
    assemblies = assemblies_meta.all_assemblies_meta(r.id)
  for m in assemblies:

    # Figure out where the input and output files live in the workdir.
    A_seqs = os.path.join(indir, r.id, m.assembler_id, m.assembly_id + ".fa")
    B_seqs = os.path.join(indir, r.id, "cc", "cc_0.fa")
    A_expr = os.path.join(indir, r.id, m.assembler_id, m.assembly_id + "_expression.isoforms.results")
    B_expr = os.path.join(indir, r.id, "cc", "cc_0_expression.isoforms.results")
    A_to_B = os.path.join(indir, r.id, m.assembler_id, m.assembly_id + "_to_cc_0.psl")
    B_to_A = os.path.join(indir, r.id, m.assembler_id, "cc_0_to_" + m.assembly_id + ".psl")
    scores = os.path.join(outdir, r.id, m.assembler_id, m.assembly_id + "_to_cc_0.scores")

    # Make sure the output dir exists.
    os.makedirs(os.path.dirname(scores), exist_ok=True)

    # Specify the inputs and outputs for this job.
    inputs = [
      nfwq.input_file(A_seqs, "A_seqs", cache=False),
      nfwq.input_file(B_seqs, "B_seqs", cache=False),
      nfwq.input_file(A_expr, "A_expr", cache=False),
      nfwq.input_file(B_expr, "B_expr", cache=False),
      nfwq.input_file(A_to_B, "A_to_B", cache=False),
      nfwq.input_file(B_to_A, "B_to_A", cache=False)
    ]
    outputs = [
      nfwq.output_file(scores, "stdout", cache=False)
    ]
    bins = [
      nfwq.input_file(wrapper,  "wrapper.sh", cache=False),
      nfwq.input_file(ref_eval, "ref-eval",   cache=False)
    ]

    # Specify the ref-eval command for this job. The wrapper.sh script
    # redirects the stdout of the command specified by all but the first
    # argument to the file specified by the first argument.
    cmd = ["bash", "wrapper.sh", "stdout",
      "./ref-eval",
      "--scores=nucl,pair,contig,kmer,kc",
      "--weighted=both",
      "--A-seqs", "A_seqs",
      "--B-seqs", "B_seqs",
      "--A-expr", "A_expr",
      "--B-expr", "B_expr",
      "--A-to-B", "A_to_B",
      "--B-to-A", "B_to_A",
      "--num-reads", str(r.num_reads),
      "--readlen", str(r.readlen),
      "--kmerlen", str(r.readlen),
      "--min-segment-len", str(r.readlen)]
    if r.ss:
      cmd += ["--strand-specific"]

    # Actually add this job to the dag.
    tag = r.id + "_" + m.assembly_id
    dag.add(tag=tag, cmd=cmd, files=inputs + outputs + bins)
