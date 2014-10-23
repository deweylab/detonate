import sys, csv, assemblies_meta
from collections import namedtuple

# This script is just a quick and dirty tool to extract a subset of the full
# (large) CSV file containing the results described in our paper, for use by
# this regression test suite.
#
# This script is not intended for end-user use.

print("Do you really want to run this?")
sys.exit(1)

columns_map = [
  ("reads_id", "reads_id"),
  ("assembler_id", "assembler_id"),
  ("assembly_id", "assembly_id"),
  ("weighted_matched_nucl_precision", "weighted_nucl_precision"),
  ("weighted_matched_nucl_recall",    "weighted_nucl_recall"),
  ("weighted_matched_nucl_F1",        "weighted_nucl_F1"),
  ("unweighted_matched_nucl_precision", "unweighted_nucl_precision"),
  ("unweighted_matched_nucl_recall",    "unweighted_nucl_recall"),
  ("unweighted_matched_nucl_F1",        "unweighted_nucl_F1"),
  ("weighted_matched_pair_precision", "weighted_pair_precision"),
  ("weighted_matched_pair_recall",    "weighted_pair_recall"),
  ("weighted_matched_pair_F1",        "weighted_pair_F1"),
  ("unweighted_matched_pair_precision", "unweighted_pair_precision"),
  ("unweighted_matched_pair_recall",    "unweighted_pair_recall"),
  ("unweighted_matched_pair_F1",        "unweighted_pair_F1"),
  ("weighted_oomatched_tran_precision", "weighted_contig_precision"),
  ("weighted_oomatched_tran_recall",    "weighted_contig_recall"),
  ("weighted_oomatched_tran_F1",        "weighted_contig_F1"),
  ("unweighted_oomatched_tran_precision", "unweighted_contig_precision"),
  ("unweighted_oomatched_tran_recall",    "unweighted_contig_recall"),
  ("unweighted_oomatched_tran_F1",        "unweighted_contig_F1"),
  ("weighted_kmer_simple_recall_at_one", "weighted_kmer_recall"),
  #( "", "inverse_compression_rate"), # no equivalent
  ("WKSR-NN/(L*NR)", "kmer_compression_score"),
  ("weighted_kmer_probs_KL_A_to_M_at_one", "weighted_kmer_KL_A_to_M"),
  ("weighted_kmer_probs_KL_B_to_M_at_one", "weighted_kmer_KL_B_to_M"),
  ("weighted_kmer_probs_JS_at_one", "weighted_kmer_jensen_shannon"),
  ("weighted_kmer_probs_hellinger_at_one", "weighted_kmer_hellinger"),
  #("", "weighted_kmer_total_variation"), # no equivalent
  ("unweighted_kmer_probs_KL_A_to_M_at_one", "unweighted_kmer_KL_A_to_M"),
  ("unweighted_kmer_probs_KL_B_to_M_at_one", "unweighted_kmer_KL_B_to_M"),
  ("unweighted_kmer_probs_JS_at_one", "unweighted_kmer_jensen_shannon"),
  ("unweighted_kmer_probs_hellinger_at_one", "unweighted_kmer_hellinger"),
  #("", "unweighted_kmer_total_variation"), # no equivalent
]

fp = open("/tier2/deweylab/scratch/nathanae/summary_nov13.csv")
reader = csv.DictReader(fp, delimiter="\t")
ls = [l for l in reader]
df = {l["reads_id"] + "_" + l["assembly_id"]: l for l in ls}

Reads = namedtuple("Reads", ("id", "ss"))
reads = [
  Reads(id="ensembl_sim_oases",  ss=False),
  Reads(id="oases_mouse_real",   ss=False),
  Reads(id="real",               ss=True),
  Reads(id="trinity_yeast_real", ss=True),
]

fo = open("summary_nov13_subset.csv", "w")
print("\t".join(new for _, new in columns_map), file=fo)
output = []
for r in reads:
  if r.ss:
    assemblies = assemblies_meta.ss_assemblies_meta(r.id)
  else:
    assemblies = assemblies_meta.all_assemblies_meta(r.id)
  for m in assemblies:
    rec = df[r.id + "_" + m.assembly_id]
    subrec = [rec[old] for old, _ in columns_map]
    print("\t".join(subrec), file=fo)
