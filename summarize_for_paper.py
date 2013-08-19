import subprocess, argparse, sys, os

class ADFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawDescriptionHelpFormatter):
  pass

p = argparse.ArgumentParser(description="""

This script computes the reference-based scores used in the current paper
draft. These scores are: 

- Nucleotide F1
- Transcript F1 
- WKR-ICR ("weighted kmer recall minus inverse compression ratio"), aka "Novel score"

Before running this script, you need to run the following blat commands, where
A.fa is the assembly, and B.fa is the set of reference transcripts:

  blat -minIdentity=80 B.fa A.fa A_to_B.psl
  blat -minIdentity=80 A.fa B.fa B_to_A.psl

Then run:

  mkdir tmpdir
  python3 summarize_for_paper.py \\
    --A-seqs A.fa \\
    --B-seqs B.fa \\
    --A-expr A_expression/expression.isoforms.results \\
    --B-expr B_expression/expression.isoforms.results \\
    --A-to-B A_to_B.psl \\
    --B-to-A B_to_A.psl \\
    --readlen 76 \\
    --num_reads 52402378 \\
    --num_nucls 14625717 \\
    --output tmpdir

where appropriate values are filled in for readlen, num_reads, and num_nucls.
The argument "--strand-specific" should be included if appropriate; in this
case, it is assumed that everything is in the forward orientation.

Output useful for debugging, etc., is written to several files in tmpdir. The
main scores are written to standard output. Status messages are written to
standard error.

""", formatter_class=ADFormatter)
p.add_argument("--A-seqs", required=True, help="The assembly sequences, in FASTA format.")
p.add_argument("--B-seqs", required=True, help="The oracleset sequences, in FASTA format.")
p.add_argument("--A-expr", help="The assembly expression, as produced by RSEM in a file called *.isoforms.results.")
p.add_argument("--B-expr", help="The oracleset expression, as produced by RSEM in a file called *.isoforms.results.")
p.add_argument("--no-expr", action="store_true", help="Do not use expression at all. No weighted scores will be produced.")
p.add_argument("--A-to-B", required=True, help="The alignments of A to B.")
p.add_argument("--B-to-A", required=True, help="The alignments of B to A.")
p.add_argument("--strand-specific", action="store_true", help="Ignore alignments that are to the reverse strand.")
p.add_argument("--readlen", required=True, type=int, help="The read length.")
p.add_argument("--num_reads", required=True, type=int, help="The number of reads in the dataset used to make the assembly A.")
p.add_argument("--num_nucls", required=True, type=int, help="The number of nucleotides in the assembly A.")
p.add_argument("--frac-identity-thresh", default="0.99", help="The threshold for frac_identity (wrt both a and b) below which an alignment is not counted.")
p.add_argument("--frac-indel-thresh", default="0.01", help="The threshold for frac_indel (wrt both a and b) above which an alignment is not counted.")
p.add_argument("--bindir", default=os.path.dirname(os.path.abspath(__file__)), help="The directory where the summarize_* executables are located.")
p.add_argument("--output", required=True, help="Output useful for debugging will be written here.")
args = p.parse_args()

if args.no_expr:
  if args.A_expr or args.B_expr:
    print("Error: If --no-expr is given, then --A-expr and --B-expr cannot be given.", file=sys.stderr)
    sys.exit(1)
else:
  if not args.A_expr or not args.B_expr:
    print("Error: If --no-expr is not given, then --A-expr and --B-expr must be given.", file=sys.stderr)
    sys.exit(1)

if args.strand_specific:
  strand_specific = "--strand_specific"
else:
  strand_specific = ""

if args.no_expr:
  args_expr = ["--no-expr"]
else:
  args_expr = ["--A-expr", args.A_expr, "--B-expr", args.B_expr]

print("Running summarize_matched.", file=sys.stderr)
matched_output = subprocess.check_output([
  args.bindir + "/summarize_matched",
  "--A-seqs", args.A_seqs,
  "--B-seqs", args.B_seqs,
  ] + args_expr + [
  "--A-to-B", args.A_to_B,
  "--B-to-A", args.B_to_A,
  "--alignment-type", "psl",
  strand_specific,
  "--readlen", str(args.readlen),
  "--plot-output", args.output + "/summarize_matched_plot_output"]).decode("utf-8").strip("\n")
print(matched_output, end="", file=open(args.output + "/summarize_matched_output", "w"))

print("Running summarize_oomatched.", file=sys.stderr)
oomatched_output = subprocess.check_output([
  args.bindir + "/summarize_oomatched",
  "--A-seqs", args.A_seqs,
  "--B-seqs", args.B_seqs,
  ] + args_expr + [
  "--A-to-B", args.A_to_B,
  "--B-to-A", args.B_to_A,
  "--alignment-type", "psl",
  strand_specific,
  "--frac-identity-thresh", args.frac_identity_thresh,
  "--frac-indel-thresh", args.frac_indel_thresh,
  "--output", args.output + "/summarize_oomatched_matching"]).decode("utf-8").strip("\n")
print(oomatched_output, end="", file=open(args.output + "/summarize_oomatched_output", "w"))

if not args.no_expr:
  print("Running summarize_kmer.", file=sys.stderr)
  kmer_output = subprocess.check_output([
    args.bindir + "/summarize_kmer",
    "--A-seqs", args.A_seqs,
    "--B-seqs", args.B_seqs,
    ] + args_expr + [
    strand_specific,
    "--readlen", str(args.readlen)]).decode("utf-8").strip("\n")
  print(kmer_output, end="", file=open(args.output + "/summarize_kmer_output", "w"))

matched_dict   = dict(l.strip("\n").split("\t") for l in matched_output.split("\n"))
oomatched_dict = dict(l.strip("\n").split("\t") for l in oomatched_output.split("\n"))

if not args.no_expr:
  kmer_dict    = dict(l.strip("\n").split("\t") for l in kmer_output.split("\n"))
  WKR = float(kmer_dict["weighted_kmer_frac_present_at_one"])
  ICR = 1.0 * args.num_nucls / (args.num_reads * args.readlen)

print("nucl_precision\t{}".format(matched_dict["unweighted_matched_nucl_precision"]))
print("nucl_recall\t{}"   .format(matched_dict["unweighted_matched_nucl_recall"]))
print("nucl_F1\t{}"       .format(matched_dict["unweighted_matched_nucl_F1"]))

print("tran_precision\t{}".format(oomatched_dict["unweighted_oomatched_tran_precision"]))
print("tran_recall\t{}"   .format(oomatched_dict["unweighted_oomatched_tran_recall"]))
print("tran_F1\t{}"       .format(oomatched_dict["unweighted_oomatched_tran_F1"]))

if not args.no_expr:
  print("WKR\t{}"    .format(WKR))
  print("ICR\t{}"    .format(ICR))
  print("WKR-ICR\t{}".format(WKR-ICR))
