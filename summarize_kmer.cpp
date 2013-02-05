#include "summarize_kmer_meat.hh"

int main(int argc, const char **argv)
{
  try {

    std::ios::sync_with_stdio(false);

    boost::program_options::variables_map vm;
    if (0) {
      int fake_argc = 11;
      const char *fake_argv[] = {
        "summarize",
        "--A-seqs", "../summarize_nucl/ensembl_sim_oases_default/oases_default.fa",
        "--B-seqs", "../summarize_nucl/ensembl_sim_oases_default/cc_0.fa",
        "--A-expr", "../summarize_nucl/ensembl_sim_oases_default/oases_default_expression/expression.isoforms.results",
        "--B-expr", "../summarize_nucl/ensembl_sim_oases_default/cc_0_expression/expression.isoforms.results",
        "--readlen", "76"
      };
      parse_options(vm, fake_argc, fake_argv);
      argc = fake_argc;
      argv = fake_argv;
    } else if (0) {
      int fake_argc = 11;
      const char *fake_argv[] = {
        "summarize",
        "--A-seqs", "../summarize_nucl/test_3/A.fa",
        "--B-seqs", "../summarize_nucl/test_3/B.fa",
        "--A-expr", "../summarize_nucl/test_3/A_expression/expression.isoforms.results",
        "--B-expr", "../summarize_nucl/test_3/B_expression/expression.isoforms.results",
        "--readlen", "76"
      };
      parse_options(vm, fake_argc, fake_argv);
      argc = fake_argc;
      argv = fake_argv;
    } else if (0) {
      int fake_argc = 11;
      const char *fake_argv[] = {
        "summarize",
        "--A-seqs", "../summarize_nucl/test_3/A.fa",
        "--B-seqs", "../summarize_nucl/test_3/A.fa",
        "--A-expr", "../summarize_nucl/test_3/A_expression/expression.isoforms.results",
        "--B-expr", "../summarize_nucl/test_3/A_expression/expression.isoforms.results",
        "--readlen", "76"
      };
      parse_options(vm, fake_argc, fake_argv);
      argc = fake_argc;
      argv = fake_argv;
    }
    parse_options(vm, argc, argv);

    main_1(vm);

  } catch (const std::exception& x) {
    std::cerr << "Exception: " << x.what() << std::endl;
    return 1;
  }

  return 0;
}
