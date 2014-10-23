import os, sys, subprocess, assemblies_meta

# This script is just a quick and dirty tool to extract the assemblies
# evaluated in the DETONATE paper, for use by this regression test suite.
#
# This script is not intended for end-user use.

print("Do you really want to run this?")
sys.exit(1)

def vcc(cmd):
  print(cmd)
  subprocess.check_call(cmd)

indir  = "/tier2/deweylab/scratch/nathanae"
outdir = "/scratch/nathanae/detonate_assemblies"

for reads_id in ["oases_mouse_real", "ensembl_sim_oases", "real", "trinity_yeast_real"]:
  vcc(["mkdir", "-p", os.path.join(outdir, reads_id)])
  vcc(["mkdir", "-p", os.path.join(outdir, reads_id, "cc")])
  vcc(["cp", "-p",
      os.path.join(indir,  reads_id, "cc", "cc_0.fa"),
      os.path.join(outdir, reads_id, "cc", "cc_0.fa")])
  vcc(["cp", "-p",
      os.path.join(indir,  reads_id, "cc", "cc_0_expression/expression.isoforms.results"),
      os.path.join(outdir, reads_id, "cc", "cc_0_expression.isoforms.results")])
  if reads_id in ["real", "trinity_yeast_real"]:
    assemblies = assemblies_meta.ss_assemblies_meta(reads_id)
  else:
    assemblies = assemblies_meta.all_assemblies_meta(reads_id)
  for m in assemblies:
    vcc(["mkdir", "-p", os.path.join(outdir, reads_id, m.assembler_id)])
    vcc(["cp", "-p",
        os.path.join(indir,  reads_id, m.assembler_id, m.assembly_id + ".fa"),
        os.path.join(outdir, reads_id, m.assembler_id, m.assembly_id + ".fa")])
    vcc(["cp", "-p",
        os.path.join(indir,  reads_id, m.assembler_id, m.assembly_id + "_to_cc_0.psl"),
        os.path.join(outdir, reads_id, m.assembler_id, m.assembly_id + "_to_cc_0.psl")])
    vcc(["cp", "-p",
        os.path.join(indir,  reads_id, m.assembler_id, "cc_0_to_" + m.assembly_id + ".psl"),
        os.path.join(outdir, reads_id, m.assembler_id, "cc_0_to_" + m.assembly_id + ".psl")])
    vcc(["cp", "-p",
        os.path.join(indir,  reads_id, m.assembler_id, m.assembly_id + "_expression/expression.isoforms.results"),
        os.path.join(outdir, reads_id, m.assembler_id, m.assembly_id + "_expression.isoforms.results")])
