import sys, os, argparse, csv, assemblies_meta
from collections import namedtuple

p = argparse.ArgumentParser()
p.add_argument("--workdir", required=True)
args = p.parse_args()

# Read old results. They are stored as a map from (reads_id,
# assembly_id) to a record.
script_dir = os.path.dirname(os.path.abspath(__file__))
fp = open(script_dir + "/summary_nov13_subset.csv")
dr = csv.DictReader(fp, delimiter="\t")
df = {(l["reads_id"], l["assembly_id"]): l for l in dr}

Reads = namedtuple("Reads", ("id", "ss"))
reads = [
  Reads(id="ensembl_sim_oases",  ss=False),
  Reads(id="oases_mouse_real",   ss=False),
  Reads(id="real",               ss=True),
  Reads(id="trinity_yeast_real", ss=True),
]

# Read the new results, and for each one, compare it to the
# corresponding old results.
num_failed, num_passed = 0, 0
d = os.path.join(args.workdir, "ref_eval_scores")
for r in reads:
  if r.ss:
    assemblies = assemblies_meta.ss_assemblies_meta(r.id)
  else:
    assemblies = assemblies_meta.all_assemblies_meta(r.id)
  for m in assemblies:
    fn = os.path.join(d, r.id, m.assembler_id, m.assembly_id + "_to_cc_0.scores")
    if not os.path.exists(fn):
      print("warning: {} doesn't exist, so skipping its tests".format(fn))
      continue
    fp = open(fn)
    ls = [l.strip("\n").split("\t") for l in fp]
    new_row = dict(ls)
    old_row = df[(r.id, m.assembly_id)]
    for k in old_row.keys() - ["reads_id", "assembler_id", "assembly_id"]:
      if k not in new_row:
        print("failure: {}, {}: {}: not found in new results".format(
          r.id, m.assembly_id, k))
        num_failed += 1
      elif abs(float(new_row[k]) - float(old_row[k])) > 1e-5:
        print("failure: {}, {}: {}: {} does not equal old {}".format(
          r.id, m.assembly_id, k, new_row[k], old_row[k]))
        num_failed += 1
      else:
        num_passed += 1
print("{} passed (i.e., {} fields had the same value as in the DETONATE paper)".format(num_passed, num_passed))
print("{} failed (i.e., {} fields had a different value as in the DETONATE paper)".format(num_failed, num_failed))
