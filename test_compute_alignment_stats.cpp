#include <string>
#include <iostream>
#define BOOST_TEST_MODULE test_compute_alignment_stats
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include "summarize_meat.hh"
#include "fake_alignment.hh"

#if 1
#define CHECK_CLOSE(x, y) \
  do { \
    double w = x; \
    double v = y; \
    double z = fabs(w - v); \
    if (z > 1e-5) { \
      std::string sx = boost::lexical_cast<std::string>(w); \
      std::string sy = boost::lexical_cast<std::string>(v); \
      std::string sz = boost::lexical_cast<std::string>(z); \
      BOOST_ERROR(sx + " != " + sy + " [fabs(x - y) = " + sz + " > 1e-5]"); \
    } \
  } while (0)
#else
// segfaults on mac
#define CHECK_CLOSE(x, y) BOOST_CHECK_CLOSE(x, y, 0.01)
#endif

std::vector<std::string> make_seqs(std::vector<size_t> sizes)
{
  std::vector<std::string> v;
  for (size_t s : sizes)
    v.push_back(std::string("N", s));
  return v;
}

std::map<std::string, size_t> invert(std::vector<std::string> v)
{
  std::map<std::string, size_t> w;
  size_t i = 0;
  for (auto x : v)
    w.insert({x, i++});
  return w;
}

std::vector<std::vector<const fake_alignment *> > pointerize(const std::vector<std::vector<fake_alignment> >& v)
{
  std::vector<std::vector<const fake_alignment *> > w(v.size());
  for (size_t i = 0; i < v.size(); ++i) {
    for (size_t j = 0; j < v[i].size(); ++j)
      w[i].push_back(&(v[i][j]));
    assert(w[i].size() == v[i].size());
  }
  return w;
}

size_t choose_2(size_t n) { return n*(n-1)/2; }

// one perfect alignment from a0 -> b0, where these are the only seqs present
BOOST_AUTO_TEST_CASE(sanity)
{
  std::vector<fake_alignment> best_to_b0 = { fake_alignment{"a0", "b0", -1, -1, { alignment_segment{0, 100-1, 0, 100-1, {}, {}} }} };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({100});
  std::vector<double>      nu_A  = {1.0};
  std::vector<double>      tau_A = {1.0};

  std::vector<std::string> B     = make_seqs({100});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = {1.0};

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  for (auto x : {pair, nucl, tran}) {
    CHECK_CLOSE(x.precis, 1.0);
    CHECK_CLOSE(x.recall, 1.0);
    CHECK_CLOSE(x.F1, 1.0);
  }
  CHECK_CLOSE(b_frac_ones[0], 1.0);
}

// no alignments from a0 -> b0, where these are the only seqs present
BOOST_AUTO_TEST_CASE(sanity_2)
{
  std::vector<fake_alignment> best_to_b0 = { };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({100});
  std::vector<double>      nu_A  = {1.0};
  std::vector<double>      tau_A = {1.0};

  std::vector<std::string> B     = make_seqs({100});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = {1.0};

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  for (auto x : {pair, nucl, tran}) {
    CHECK_CLOSE(x.precis, 0.0);
    CHECK_CLOSE(x.recall, 0.0);
    CHECK_CLOSE(x.F1, 0.0);
  }
  CHECK_CLOSE(b_frac_ones[0], 0.0);
}

// a0 and a1 together cover b0 perfectly
//  a: 0000011111
//  b: ----------
BOOST_AUTO_TEST_CASE(perfect_two_to_one)
{
  std::vector<fake_alignment> best_to_b0 = {
    fake_alignment{"a0", "b0", -1, -1, { alignment_segment{0, 5000-1,    0,  5000-1, {}, {}} }},
    fake_alignment{"a1", "b0", -1, -1, { alignment_segment{0, 5000-1, 5000, 10000-1, {}, {}} }},
  };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({5000, 5000});
  std::vector<double>      nu_A  = {0.5, 0.5};
  std::vector<double>      tau_A = {0.5, 0.5};

  std::vector<std::string> B     = make_seqs({10000});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = {1.0};

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0", "a1"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  for (auto x : {nucl, tran}) {
    CHECK_CLOSE(x.precis, 1.0);
    CHECK_CLOSE(x.recall, 1.0);
    CHECK_CLOSE(x.F1, 1.0);
  }

  CHECK_CLOSE(pair.precis, 1.0);
  CHECK_CLOSE(pair.recall, 1.0*(choose_2(5000+1) + choose_2(5000+1))/choose_2(10000+1));

  CHECK_CLOSE(b_frac_ones[0], 1.0);
}

// a0 and a1 together cover b0 and overlap each other
//  a: 00000
//         111111
//  b: ----------
BOOST_AUTO_TEST_CASE(overlapping_two_to_one)
{
  std::vector<fake_alignment> best_to_b0 = {
    fake_alignment{"a0", "b0", -1, -1, { alignment_segment{0, 50-1,  0,  50-1, {}, {}} }},
    fake_alignment{"a1", "b0", -1, -1, { alignment_segment{0, 60-1, 40, 100-1, {}, {}} }},
  };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({50, 60});
  std::vector<double>      nu_A  = {0.3, 0.7};
  std::vector<double>      tau_A = {0.3, 0.7};

  std::vector<std::string> B     = make_seqs({100});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = {1.0};

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0", "a1"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  for (auto x : {nucl, tran}) {
    CHECK_CLOSE(x.precis, 1.0);
    CHECK_CLOSE(x.recall, 1.0);
    CHECK_CLOSE(x.F1, 1.0);
  }

  CHECK_CLOSE(pair.precis, 1.0);
  CHECK_CLOSE(pair.recall, 1.0*(choose_2(50+1) + choose_2(60+1) - choose_2(10+1))/choose_2(100+1));

  CHECK_CLOSE(b_frac_ones[0], 1.0);
}

// a0 and a1 together cover 97/100 of b0
// the alignment from a0->b0 covers 49/50 of a0
// the alignment from a1->b0 covers 48/50 of a1
//  a: 00000
//          11111
//  b: ----------
BOOST_AUTO_TEST_CASE(two_to_one_partial_coverage)
{
  std::vector<fake_alignment> best_to_b0 = {
    fake_alignment{"a0", "b0", -1, -1, { alignment_segment{0, 49-1,  0, 49-1, {}, {}} }},
    fake_alignment{"a1", "b0", -1, -1, { alignment_segment{0, 48-1, 50, 98-1, {}, {}} }},
  };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({50, 50});
  std::vector<double>      nu_A  = {0.3, 0.7};
  std::vector<double>      tau_A = {0.3, 0.7};

  std::vector<std::string> B     = make_seqs({100});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = {1.0};

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0", "a1"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  CHECK_CLOSE(pair.precis, 0.3*choose_2(49+1)/choose_2(50+1) + 0.7*choose_2(48+1)/choose_2(50+1));
  CHECK_CLOSE(pair.recall, 1.0*(choose_2(49+1) + choose_2(48+1))/choose_2(100+1));

  CHECK_CLOSE(nucl.precis, 0.3*49/50 + 0.7*48/50);
  CHECK_CLOSE(nucl.recall, 1.0*97/100);

  CHECK_CLOSE(tran.precis, 1.0);
  CHECK_CLOSE(tran.recall, 1.0);

  CHECK_CLOSE(b_frac_ones[0], 1.0*97/100);
}

// a0 and a1 together cover 97/1000 of b0
// the alignment from a0->b0 covers 49/500 of a0
// the alignment from a1->b0 covers 48/50 of a1
//  a: 00000
//          11111
//  b: ----------
//
// This is the same as two_to_one_partial_coverage, except that two of the
// sequences are longer, so they are ignored by the transcript variants.
BOOST_AUTO_TEST_CASE(two_to_one_partial_coverage_long)
{
  std::vector<fake_alignment> best_to_b0 = {
    fake_alignment{"a0", "b0", -1, -1, { alignment_segment{0, 49-1,  0, 49-1, {}, {}} }},
    fake_alignment{"a1", "b0", -1, -1, { alignment_segment{0, 48-1, 50, 98-1, {}, {}} }},
  };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({500, 50});
  std::vector<double>      nu_A  = {0.3, 0.7};
  std::vector<double>      tau_A = {0.3, 0.7};

  std::vector<std::string> B     = make_seqs({1000});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = {1.0};

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0", "a1"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  CHECK_CLOSE(pair.precis, 0.3*choose_2(49+1)/choose_2(500+1) + 0.7*choose_2(48+1)/choose_2(50+1));
  CHECK_CLOSE(pair.recall, 1.0*(choose_2(49+1) + choose_2(48+1))/choose_2(1000+1));

  CHECK_CLOSE(nucl.precis, 0.3*49/500 + 0.7*48/50);
  CHECK_CLOSE(nucl.recall, 1.0*97/1000);

  CHECK_CLOSE(tran.precis, 0.7);
  CHECK_CLOSE(tran.recall, 0.0);

  CHECK_CLOSE(b_frac_ones[0], 1.0*97/1000);
}

//  a0: ---x---
//  b0: -------
BOOST_AUTO_TEST_CASE(one_mismatch)
{
  std::vector<fake_alignment> best_to_b0 = {
    fake_alignment{"a0", "b0", -1, -1, { alignment_segment{0, 50-1,  0, 50-1, {25}, {25}} }},
  };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({50});
  std::vector<double>      nu_A  = {1.0};
  std::vector<double>      tau_A = nu_A;

  std::vector<std::string> B     = make_seqs({50});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = nu_B;

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  CHECK_CLOSE(pair.precis, 1.0*choose_2(49+1)/choose_2(50+1));
  CHECK_CLOSE(pair.recall, 1.0*choose_2(49+1)/choose_2(50+1));

  CHECK_CLOSE(nucl.precis, 1.0*49/50);
  CHECK_CLOSE(nucl.recall, 1.0*49/50);

  CHECK_CLOSE(tran.precis, 1.0);
  CHECK_CLOSE(tran.recall, 1.0);

  CHECK_CLOSE(b_frac_ones[0], 1.0*49/50);
}

//  a0: -x-x---
//  b0: -------
BOOST_AUTO_TEST_CASE(two_mismatches)
{
  std::vector<fake_alignment> best_to_b0 = {
    fake_alignment{"a0", "b0", -1, -1, { alignment_segment{0, 50-1,  0, 50-1, {15, 25}, {15, 25}} }},
  };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({50});
  std::vector<double>      nu_A  = {1.0};
  std::vector<double>      tau_A = nu_A;

  std::vector<std::string> B     = make_seqs({50});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = nu_B;

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  CHECK_CLOSE(pair.precis, 1.0*choose_2(48+1)/choose_2(50+1));
  CHECK_CLOSE(pair.recall, 1.0*choose_2(48+1)/choose_2(50+1));

  CHECK_CLOSE(nucl.precis, 1.0*48/50);
  CHECK_CLOSE(nucl.recall, 1.0*48/50);

  CHECK_CLOSE(tran.precis, 1.0);
  CHECK_CLOSE(tran.recall, 1.0);

  CHECK_CLOSE(b_frac_ones[0], 1.0*48/50);
}

//  a0: --xx---
//  a1: -----x-
//  b0: -------
BOOST_AUTO_TEST_CASE(disjoint_mismatches)
{
  std::vector<fake_alignment> best_to_b0 = {
    fake_alignment{"a0", "b0", -1, -1, { alignment_segment{0, 50-1,  0, 50-1, {24, 25}, {24, 25}} }},
    fake_alignment{"a1", "b0", -1, -1, { alignment_segment{0, 50-1,  0, 50-1, {30    }, {30    }} }},
  };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({50, 50});
  std::vector<double>      nu_A  = {0.5, 0.5};
  std::vector<double>      tau_A = nu_A;

  std::vector<std::string> B     = make_seqs({50});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = nu_B;

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0", "a1"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  CHECK_CLOSE(pair.precis, 0.5*choose_2(48+1)/choose_2(50+1) + 0.5*choose_2(49+1)/choose_2(50+1));
  CHECK_CLOSE(pair.recall, 1.0*(choose_2(48+1) + choose_2(49+1) - choose_2(47+1))/choose_2(50+1));

  CHECK_CLOSE(nucl.precis, 0.5*48/50 + 0.5*49/50);
  CHECK_CLOSE(nucl.recall, 1.0);

  CHECK_CLOSE(tran.precis, 1.0);
  CHECK_CLOSE(tran.recall, 1.0);

  CHECK_CLOSE(b_frac_ones[0], 1.0);
}

//  a0: ---x-x-
//  a1:    --x--x-
//  b0: ----------
BOOST_AUTO_TEST_CASE(more_complicated_mismatches)
{
  std::vector<fake_alignment> best_to_b0 = {
    fake_alignment{"a0", "b0", -1, -1, { alignment_segment{0, 50-1,   0,  50-1, {25, 30}, {25, 30}} }},
    fake_alignment{"a1", "b0", -1, -1, { alignment_segment{5, 80-1,  25, 100-1, {10, 50}, {30, 70}} }},
  };
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0 };

  std::vector<std::string> A     = make_seqs({50, 80});
  std::vector<double>      nu_A  = {0.3, 0.7};
  std::vector<double>      tau_A = nu_A;

  std::vector<std::string> B     = make_seqs({100});
  std::vector<double>      nu_B  = {1.0};
  std::vector<double>      tau_B = nu_B;

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0", "a1"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  CHECK_CLOSE(pair.precis, 0.3*choose_2(48+1)/choose_2(50+1) + 0.7*choose_2(73+1)/choose_2(80+1));
  CHECK_CLOSE(pair.recall, 1.0*(choose_2(48+1) + choose_2(73+1) - choose_2(23+1))/choose_2(100+1));
  // overlap: 25, -1 missing in both, -1 missing in  a0 ==> 23

  CHECK_CLOSE(nucl.precis, 0.3*48/50 + 0.7*73/80);
  CHECK_CLOSE(nucl.recall, 1.0*98/100);

  CHECK_CLOSE(tran.precis, 0.3);
  CHECK_CLOSE(tran.recall, 1.0);

  CHECK_CLOSE(b_frac_ones[0], 1.0*98/100);
}

// A:  00000           3333                        6666
//      --1111111            5555
//      4444            
//
// B:  000000000000   111111111111  22222222222
// 
//  a0: ----
//  a1:    --x--x-
//  b0: ----------
BOOST_AUTO_TEST_CASE(several_B_elements)
{
  std::vector<fake_alignment> best_to_b0 = {
    fake_alignment{"a0", "b0", -1, -1, { alignment_segment{ 0, 30-1,   0,  30-1, {}, {}} }},
    fake_alignment{"a1", "b0", -1, -1, { alignment_segment{10, 60-1,  20,  70-1, {}, {}} }},
    fake_alignment{"a4", "b0", -1, -1, { alignment_segment{ 0, 20-1,  10,  20-1, {}, {}} }},
  };
  std::vector<fake_alignment> best_to_b1 = {
    fake_alignment{"a3", "b1", -1, -1, { alignment_segment{ 0, 15-1,   5,  20-1, {}, {}} }},
    fake_alignment{"a5", "b1", -1, -1, { alignment_segment{ 0, 25-1,  55,  80-1, {}, {}} }},
  };
  std::vector<fake_alignment> best_to_b2 = {};
  std::vector<std::vector<fake_alignment> > best_to_B = { best_to_b0, best_to_b1, best_to_b2 };

  std::vector<std::string> A     = make_seqs({30, 60, 300, 15, 20, 25, 10});
  std::vector<double>      nu_A  = {0.3, 0.25, 0.2, 0.15, 0.03, 0.07};
  std::vector<double>      tau_A = nu_A;

  std::vector<std::string> B     = make_seqs({100, 90, 40});
  std::vector<double>      nu_B  = {0.1, 0.4, 0.5};
  std::vector<double>      tau_B = nu_B;

  std::map<std::string,size_t> A_names_to_idxs = invert({"a0", "a1", "a2", "a3", "a4", "a5", "a6"});

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(best_to_B.size());
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones,
                                         pointerize(best_to_B),
                                         A, nu_A, tau_A,
                                         B, nu_B, tau_B,
                                         A_names_to_idxs);

  CHECK_CLOSE(pair.precis, nu_A[0] + nu_A[1]*choose_2(50+1)/choose_2(60+1) + nu_A[4] +
                           nu_A[3] + nu_A[5]);
  CHECK_CLOSE(pair.recall, nu_B[0]*(choose_2(30+1) + choose_2(50+1) - choose_2(10+1))/choose_2(100+1) +
                           nu_B[1]*(choose_2(15+1) + choose_2(25+1))/choose_2(90+1));

  CHECK_CLOSE(nucl.precis, nu_A[0] + nu_A[1]*50/60 + nu_A[4] + nu_A[3] + nu_A[5]);
  CHECK_CLOSE(nucl.recall, nu_B[0]*70/100 + nu_B[1]*(15+25)/90);

  CHECK_CLOSE(tran.precis, tau_A[0] + tau_A[4] + tau_A[3] + tau_A[5]);
  CHECK_CLOSE(tran.recall, 0.0);

  CHECK_CLOSE(b_frac_ones[0], 1.0*70/100);
  CHECK_CLOSE(b_frac_ones[1], 1.0*(15+25)/90);
  CHECK_CLOSE(b_frac_ones[2], 0.0);
}
