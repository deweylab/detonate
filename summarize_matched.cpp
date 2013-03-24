#include "summarize_matched_meat.hh"

int main(int argc, const char **argv)
{
  try {

    std::ios::sync_with_stdio(false);

    boost::program_options::variables_map vm;
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
