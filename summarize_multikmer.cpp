#include "summarize_multikmer_meat.hh"

int main(int argc, const char **argv)
{
  try {

    std::ios::sync_with_stdio(false);

    boost::program_options::variables_map vm;
    parse_options(vm, argc, argv);

    main_1(vm);

  } catch (const std::exception& x) {
    std::cerr << "Exception: " << x.what() << std::endl;
    return 1;
  }

  return 0;
}
