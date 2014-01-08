from collections import namedtuple

trinity_assembly_meta    = namedtuple("trinity_assembly_meta",    ("reads_id", "assembler_id", "assembly_id", "paramset"))
oases_assembly_meta      = namedtuple("oases_assembly_meta",      ("reads_id", "assembler_id", "assembly_id", "velvet_paramset", "oases_paramset"))
soap_assembly_meta       = namedtuple("soap_assembly_meta",       ("reads_id", "assembler_id", "assembly_id", "paramset", "assembler_variant"))
transabyss_assembly_meta = namedtuple("transabyss_assembly_meta", ("reads_id", "assembler_id", "assembly_id", "paramset"))

def trinity_assemblies_meta(reads_id):
  single_paramsets = [
    # (paramset, paramset_name)
    ("--min_glue 1", "min_glue_1"),
    ("--min_iso_ratio 0.01", "min_iso_ratio_0.01"),
    ("--min_iso_ratio 0.1", "min_iso_ratio_0.1"),
    ("--glue_factor 0.01", "glue_factor_0.01"),
    ("--glue_factor 0.1", "glue_factor_0.1"),
    ("--min_pct_read_mapping 1", "min_pct_read_mapping_1"),
    ("--min_pct_read_mapping 5", "min_pct_read_mapping_5"),
    ("--max_reads_per_graph 10000000", "max_reads_per_graph_10000000"),
    ("--max_reads_per_graph 50000000", "max_reads_per_graph_50000000"),
    ("--max_number_of_paths_per_node 1", "max_number_of_paths_per_node_1"),
    ("--max_number_of_paths_per_node 100", "max_number_of_paths_per_node_100"),
    ("--path_reinforcement_distance 1", "path_reinforcement_distance_1"),
    ("--path_reinforcement_distance 37", "path_reinforcement_distance_37")
  ]
  default_paramset = ("", "default")
  paramsets = [default_paramset] + list(single_paramsets)
  for p1, n1 in single_paramsets:
    for p2, n2 in single_paramsets:
      if p1 == p2:
        break
      if p1.split()[0] == p2.split()[0]: # same parameter, different value
        continue
      p = p1 + " " + p2
      n = n1 + "_" + n2
      paramsets.append((p, n))
  return [trinity_assembly_meta(
            reads_id = reads_id,
            assembler_id = "trinity",
            assembly_id = "trinity_" + paramset_name,
            paramset = paramset)
          for (paramset, paramset_name) in paramsets]

def oases_assemblies_meta(reads_id, only_single=False):
  single_paramsets = [
    # (velvet_paramset, oases_paramset, paramset_name)
    ("-cov_cutoff auto", "", "velvet_cov_cutoff_auto"),
    ("-cov_cutoff 5", "", "velvet_cov_cutoff_5"),
    ("-exp_cov auto", "", "velvet_exp_cov_auto"),
    ("-exp_cov 20", "", "velvet_exp_cov_20"),
    ("", "-cov_cutoff 0", "cov_cutoff_0"),
    ("", "-cov_cutoff 10", "cov_cutoff_10"),
    ("", "-edgeFractionCutoff 0.1", "edgeFractionCutoff_0.1"),
    ("", "-edgeFractionCutoff 0.5", "edgeFractionCutoff_0.5"),
    ("", "-degree_cutoff 1", "degree_cutoff_1"),
    ("", "-degree_cutoff 10", "degree_cutoff_10"),
    ("", "-min_pair_count 1", "min_pair_count_1"),
    ("", "-min_pair_count 10", "min_pair_count_10"),
  ]
  default_paramset = ("", "", "default")
  paramsets = [default_paramset] + list(single_paramsets)
  if not only_single:
    for vp1, op1, n1 in single_paramsets:
      for vp2, op2, n2 in single_paramsets:
        if vp1 == vp2 and op1 == op2: # don't have same param twice
          break
        if ((vp1 != "" and vp2 != "" and vp1.split()[0] == vp2.split()[0]) or # same parameter, different value
            (op1 != "" and op2 != "" and op1.split()[0] == op2.split()[0])):
          continue
        vp = vp1 + " " + vp2
        op = op1 + " " + op2
        n = n1 + "_" + n2
        paramsets.append((vp, op, n))
  return [oases_assembly_meta(
            reads_id = reads_id,
            assembler_id = "oases",
            assembly_id = "oases_" + paramset_name,
            velvet_paramset = velvet_paramset,
            oases_paramset = oases_paramset)
          for (velvet_paramset, oases_paramset, paramset_name) in paramsets]

def soap_assemblies_meta(reads_id):
  single_paramsets = [
    # (paramset, paramset_name, assembler_variant
    ("-K 19", "K_19", "31k"),
    #("-K 35", "K_35", "127"),
    ("-K 31", "K_31", "31k"),
    ("-M 0", "M_0", "31k"),
    ("-M 3", "M_3", "31k"),
    ("-d 1", "d_1", "31k"),
    ("-d 10", "d_10", "31k"),
    ("-D 0", "D_0", "31k"),
    ("-D 10", "D_10", "31k"),
    ("-e 0", "e_0", "31k"),
    ("-e 10", "e_10", "31k")
  ]
  default_paramset = ("", "default", "31k")
  paramsets = [default_paramset] + list(single_paramsets)
  for p1, n1, m1 in single_paramsets:
    for p2, n2, m2 in single_paramsets:
      if p1 == p2:
        break
      if p1.split()[0] == p2.split()[0]: # same parameter, different value
        continue
      p = p1 + " " + p2
      n = n1 + "_" + n2
      if m1 == "127" or m2 == "127":
        m = "127"
      else:
        m = "31k"
      paramsets.append((p, n, m))
  return [soap_assembly_meta(
            reads_id = reads_id,
            assembler_id = "soap",
            assembly_id = "soap_" + paramset_name,
            paramset = paramset,
            assembler_variant = assembler_variant)
          for (paramset, paramset_name, assembler_variant) in paramsets]

def transabyss_assemblies_meta(reads_id):
  return [transabyss_assembly_meta(
            reads_id = reads_id,
            assembler_id = "transabyss",
            assembly_id = "transabyss_k26_to_k64",
            paramset = "")]

def all_assemblies_meta(reads_id):
  return (trinity_assemblies_meta(reads_id) +
          oases_assemblies_meta(reads_id) +
          soap_assemblies_meta(reads_id) +
          transabyss_assemblies_meta(reads_id))

def ss_assemblies_meta(reads_id):
  return (trinity_assemblies_meta(reads_id) +
          oases_assemblies_meta(reads_id))

if __name__ == "__main__":
  
  import argparse, sys
  p = argparse.ArgumentParser()
  p.add_argument("--reads_ids", required=True)
  p.add_argument("--assembler_ids", required=True)
  args = p.parse_args(sys.argv[1:])

  for reads_id in args.reads_ids.split(","):
    for assembler_id in args.assembler_ids.split(","):
      for m in eval("{}_assemblies_meta(reads_id)".format(assembler_id)):
        print("{}\t{}\t{}".format(m.reads_id, m.assembler_id, m.assembly_id))
