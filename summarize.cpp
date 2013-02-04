#include "summarize_meat.hh"

int main(int argc, const char **argv)
{
  try {

    std::ios::sync_with_stdio(false);

    boost::program_options::variables_map vm;
    if (0) {
      int fake_argc = 14;
      const char *fake_argv[] = {
        "summarize",
        "--A-seqs", "../summarize_nucl/test_3/axolotl_trinity_default/trinity_default.fa",
        "--B-seqs", "../summarize_nucl/test_3/axolotl_trinity_default/xtrop.fa",
        "--A-expr", "../summarize_nucl/test_3/axolotl_trinity_default/trinity_default_expression/expression.isoforms.results",
        "--induce-B-expr",
        "--A-to-B", "../summarize_nucl/test_3/axolotl_trinity_default/test6_extra.out",
        "--alignment-type", "blast",
        "--plot-output", "plot.tmp"
      };
      parse_options(vm, fake_argc, fake_argv);
      argc = fake_argc;
      argv = fake_argv;
    } else if (0) {
      int fake_argc = 15;
      const char *fake_argv[] = {
        "summarize",
        "--A-seqs", "../summarize_nucl/ensembl_sim_oases_default/oases_default.fa",
        "--B-seqs", "../summarize_nucl/ensembl_sim_oases_default/cc_0.fa",
        "--A-expr", "../summarize_nucl/ensembl_sim_oases_default/oases_default_expression/expression.isoforms.results",
        "--B-expr", "../summarize_nucl/ensembl_sim_oases_default/cc_0_expression/expression.isoforms.results",
        "--A-to-B", "../summarize_nucl/ensembl_sim_oases_default/oases_default_to_cc_0.psl",
        "--alignment-type", "psl",
        "--plot-output", "plot.tmp"
      };
      parse_options(vm, fake_argc, fake_argv);
      argc = fake_argc;
      argv = fake_argv;
    } else if (0) {
      int fake_argc = 15;
      const char *fake_argv[] = {
        "summarize",
        "--A-seqs", "../summarize_nucl/test_3/A.fa",
        "--B-seqs", "../summarize_nucl/test_3/B.fa",
        "--A-expr", "../summarize_nucl/test_3/A_expression/expression.isoforms.results",
        "--B-expr", "../summarize_nucl/test_3/B_expression/expression.isoforms.results",
        "--A-to-B", "../summarize_nucl/test_3/A_to_B.psl",
        "--alignment-type", "psl",
        "--plot-output", "plot.tmp"
      };
      parse_options(vm, fake_argc, fake_argv);
      argc = fake_argc;
      argv = fake_argv;
    } else if (0) {
      int fake_argc = 15;
      const char *fake_argv[] = {
        "summarize",
        "--A-seqs", "../summarize_nucl/test_3/A.fa",
        "--B-seqs", "../summarize_nucl/test_3/A.fa",
        "--A-expr", "../summarize_nucl/test_3/A_expression/expression.isoforms.results",
        "--B-expr", "../summarize_nucl/test_3/A_expression/expression.isoforms.results",
        "--A-to-B", "../summarize_nucl/test_3/A_to_A.psl",
        "--alignment-type", "psl",
        "--plot-output", "plot.tmp"
      };
      parse_options(vm, fake_argc, fake_argv);
      argc = fake_argc;
      argv = fake_argv;
    }
    parse_options(vm, argc, argv);

    std::string alignment_type = vm["alignment-type"].as<std::string>();
    if (alignment_type == "blast")
      main_1<blast_alignment>(vm);
    else if (alignment_type == "psl")
      main_1<psl_alignment>(vm);

  } catch (const std::exception& x) {
    std::cerr << "Exception: " << x.what() << std::endl;
    return 1;
  }

  return 0;
}
