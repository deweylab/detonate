#include <string>
#include <iostream>
#define BOOST_TEST_MODULE test_read_cluster_filter_alignments
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include "summarize_meat.hh"
#include "fake_alignment.hh"

void check_equal_sets(std::set<std::string> set1, std::set<std::string> set2)
{
  BOOST_CHECK(set1 == set2);
}

BOOST_AUTO_TEST_CASE(basic)
{
  std::vector<alignment_segment> segs = { {0, 100-1, 0, 100-1, {}, {}} };
  std::vector<fake_alignment> als {

    // read_alignments_and_filter_by_best_from_A

    {"a0", "b0",  1.0,  0.0, segs }, // ok

    {"a1", "b1",  0.9,  0.0, segs }, // frac_identity too low

    {"a2", "b20", 0.97, 0.2, segs }, // not chosen
    {"a2", "b21", 0.97, 0.0, segs }, // chosen

    {"a3", "b30", 0.98, 0.2, segs }, // chosen
    {"a3", "b31", 0.97, 0.0, segs }, // not chosen

    // cluster_best_alignments_to_B and filter_by_best_alignment_to_B

    {"A0",  "B0", 1.0,  0.0, segs }, // in cluster and filter

    {"A10", "B1", 0.99, 0.0, segs }, // in cluster
    {"A11", "B1", 1.0,  0.0, segs }, // in cluster and filter

    {"A20", "B2", 1.0,  0.0, segs }, // in cluster and filter
    {"A21", "B2", 1.0,  0.1, segs }, // in cluster

    {"A30", "B3", 1.0,  0.1, segs }, // in cluster and filter
    {"A31", "B3", 0.99, 0.0, segs }, // in cluster

  };
  fake_alignment::input_stream_type input_stream(als);

  std::set<std::string> A_names, B_names;
  std::map<std::string, size_t> A_names_to_idxs, B_names_to_idxs;
  for (const fake_alignment& al : als) { A_names.insert(al.a_name()); B_names.insert(al.b_name()); }
  size_t i = 0; for (const std::string& a_name : A_names) A_names_to_idxs.insert(std::make_pair(a_name, i++));
  size_t j = 0; for (const std::string& b_name : B_names) B_names_to_idxs.insert(std::make_pair(b_name, j++));

  std::vector<BestTuple<fake_alignment>> best_from_A(A_names_to_idxs.size());
  std::vector<std::vector<const fake_alignment *> > clustered_best_to_B(B_names_to_idxs.size());
  std::vector<std::vector<const fake_alignment *> > filtered_best_to_B(B_names_to_idxs.size());

  read_alignments_and_filter_by_best_from_A(best_from_A, input_stream, A_names_to_idxs);
  cluster_best_alignments_to_B(clustered_best_to_B, best_from_A, B_names_to_idxs);
  filter_by_best_alignment_to_B(filtered_best_to_B, best_from_A, B_names_to_idxs);

  // read_alignments_and_filter_by_best_from_A

  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a0"]].al.a_name(), "a0");
  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a0"]].al.b_name(), "b0");
  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a0"]].frac_identity, 1.0);
  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a0"]].frac_indel, 0.0);

  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a1"]].frac_identity, -1.0);

  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a2"]].al.a_name(), "a2");
  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a2"]].al.b_name(), "b21");
  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a2"]].frac_identity, 0.97);
  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a2"]].frac_indel, 0.0);

  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a3"]].al.a_name(), "a3");
  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a3"]].al.b_name(), "b30");
  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a3"]].frac_identity, 0.98);
  BOOST_CHECK_EQUAL(best_from_A[A_names_to_idxs["a3"]].frac_indel, 0.2);

  // cluster_best_alignments_to_B and filter_by_best_alignment_to_B

  std::vector<const fake_alignment *> cur;

  cur = clustered_best_to_B[B_names_to_idxs["B0"]];
  BOOST_CHECK_EQUAL(cur.size(), 1ul);
  BOOST_CHECK_EQUAL(cur[0]->a_name(), "A0");
  BOOST_CHECK_EQUAL(cur[0]->b_name(), "B0");

  cur = filtered_best_to_B[B_names_to_idxs["B0"]];
  BOOST_CHECK_EQUAL(cur.size(), 1ul);
  BOOST_CHECK_EQUAL(cur[0]->a_name(), "A0");
  BOOST_CHECK_EQUAL(cur[0]->b_name(), "B0");

  cur = clustered_best_to_B[B_names_to_idxs["B1"]];
  BOOST_CHECK_EQUAL(cur.size(), 2ul);
  check_equal_sets({cur[0]->a_name(), cur[1]->a_name()}, {"A10", "A11"});
  check_equal_sets({cur[0]->b_name(), cur[1]->b_name()}, {"B1"});

  cur = filtered_best_to_B[B_names_to_idxs["B1"]];
  BOOST_CHECK_EQUAL(cur.size(), 1ul);
  BOOST_CHECK_EQUAL(cur[0]->a_name(), "A11");
  BOOST_CHECK_EQUAL(cur[0]->b_name(), "B1");

  cur = clustered_best_to_B[B_names_to_idxs["B2"]];
  BOOST_CHECK_EQUAL(cur.size(), 2ul);
  check_equal_sets({cur[0]->a_name(), cur[1]->a_name()}, {"A20", "A21"});
  check_equal_sets({cur[0]->b_name(), cur[1]->b_name()}, {"B2"});

  cur = filtered_best_to_B[B_names_to_idxs["B2"]];
  BOOST_CHECK_EQUAL(cur.size(), 1ul);
  BOOST_CHECK_EQUAL(cur[0]->a_name(), "A20");
  BOOST_CHECK_EQUAL(cur[0]->b_name(), "B2");

  cur = clustered_best_to_B[B_names_to_idxs["B3"]];
  BOOST_CHECK_EQUAL(cur.size(), 2ul);
  check_equal_sets({cur[0]->a_name(), cur[1]->a_name()}, {"A30", "A31"});
  check_equal_sets({cur[0]->b_name(), cur[1]->b_name()}, {"B3"});

  cur = filtered_best_to_B[B_names_to_idxs["B3"]];
  BOOST_CHECK_EQUAL(cur.size(), 1ul);
  BOOST_CHECK_EQUAL(cur[0]->a_name(), "A30");
  BOOST_CHECK_EQUAL(cur[0]->b_name(), "B3");
}
