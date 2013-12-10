// Copyright (c) 2013
// Nathanael Fillmore (University of Wisconsin-Madison)
// nathanae@cs.wisc.edu
//
// This file is part of REF-EVAL.
//
// REF-EVAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// REF-EVAL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with REF-EVAL.  If not, see <http://www.gnu.org/licenses/>.

#include "summarize_meat.hh"

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
