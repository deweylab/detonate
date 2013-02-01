#include <string>
#include <iostream>
#define BOOST_TEST_MODULE test_pslx
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include "pslx.hh"

BOOST_AUTO_TEST_CASE(alignment_concept)
{
  pslx_alignment al;
  std::string line = "79	0	0	0	3	10	2	7	+	Locus_1_Transcript_1/48_Confidence_0.136_Length_1158	1158	268	357	Locus_1_Transcript_25/41_Confidence_0.222_Length_1952	1952	1866	1952	4	17,47,6,9,	268,286,335,348,	1866,1883,1931,1943,	catgtgcatcaacatct,ggacacagacctcggaatactggtccacgtaggttcttgcatatgtc,gacacg,ccccaagcg,	catgtgcatcaacatct,ggacacagacctcggaatactggtccacgtaggttcttgcatatgtc,gacacg,ccccaagcg,";
  al.parse_line(line);

  BOOST_CHECK_EQUAL(al.a_name(), "Locus_1_Transcript_1/48_Confidence_0.136_Length_1158");
  BOOST_CHECK_EQUAL(al.b_name(), "Locus_1_Transcript_25/41_Confidence_0.222_Length_1952");
  BOOST_CHECK_EQUAL(al.frac_identity(), 1.0*(79+0)/1158);
  BOOST_CHECK_EQUAL(al.frac_indel(), 1.0*(10+7)/1158);
}

BOOST_AUTO_TEST_CASE(single_interval_realistic_usage)
{
  pslx_alignment al;
  std::string line = "32	1	0	0	0	0	0	0	+	Locus_54_Transcript_14/14_Confidence_0.000_Length_639	639	460	493	Locus_51_Transcript_1/5_Confidence_0.800_Length_1299	1299	1235	1268	1	33,	460,	1235,	aaaaaagaaaaaaaaaaaaaaaaaaaaaaaaaa,	aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa,";
  al.parse_line(line);
  size_t i = 0;
  BOOST_FOREACH(const alignment_segment& seg, al.segments("", "")) {
    BOOST_CHECK_EQUAL(seg.a_start, 460ul);
    BOOST_CHECK_EQUAL(seg.a_end,   460ul + 33 - 1);
    BOOST_CHECK_EQUAL(seg.b_start, 1235ul);
    BOOST_CHECK_EQUAL(seg.b_end,   1235ul + 33 - 1);
    BOOST_CHECK_EQUAL(seg.a_mismatches.size(), 1ul);
    BOOST_CHECK_EQUAL(seg.b_mismatches.size(), 1ul);
    BOOST_CHECK_EQUAL(seg.a_mismatches[0], 460ul + 6);
    BOOST_CHECK_EQUAL(seg.b_mismatches[0], 1235ul + 6);
    ++i;
  }
  BOOST_CHECK_EQUAL(i, 1ul);
}

BOOST_AUTO_TEST_CASE(several_intervals_forward_strand)
{
  pslx_alignment al;
  std::string line = "79	0	0	0	3	10	2	7	+	Locus_1_Transcript_1/48_Confidence_0.136_Length_1158	1158	268	357	Locus_1_Transcript_25/41_Confidence_0.222_Length_1952	1952	1866	1952	4	"
    "17,47,6,9,\t"             // block sizes
    "268,286,335,348,\t"       // a starts
    "1866,1883,1931,1943,\t"   // b starts
    "ccccccccccccccccc,ccccccccccccccccccccccccccccccccccccccccccccccc,cgcgcc,ccccccccc,\t" // a segment seqs
    "ccccccccccccccccc,gcccccccccccccccccccccccccccccccccccccccccccccc,cccccc,ccccccccc,";  // b segment seqs
  al.parse_line(line);
  pslx_alignment::segments_type segs = al.segments("", "");
  pslx_alignment::segments_type::const_iterator it = segs.begin(), end = segs.end();
  size_t block_sizes[] = {17,47,6,9};
  size_t a_starts[] = {268,286,335,348};
  size_t b_starts[] = {1866,1883,1931,1943};
  for (size_t i = 0; i < 4; ++it, ++i) {
    BOOST_CHECK(it != end);
    BOOST_CHECK_EQUAL(it->a_start, a_starts[i]);
    BOOST_CHECK_EQUAL(it->a_end,   a_starts[i] + block_sizes[i] - 1);
    BOOST_CHECK_EQUAL(it->b_start, b_starts[i]);
    BOOST_CHECK_EQUAL(it->b_end,   b_starts[i] + block_sizes[i] - 1);
    if (i == 0 || i == 3) {
      BOOST_CHECK_EQUAL(it->a_mismatches.size(), 0ul);
      BOOST_CHECK_EQUAL(it->b_mismatches.size(), 0ul);
    } else if (i == 1) {
      BOOST_CHECK_EQUAL(it->a_mismatches.size(), 1ul);
      BOOST_CHECK_EQUAL(it->b_mismatches.size(), 1ul);
      BOOST_CHECK_EQUAL(it->a_mismatches[0], a_starts[i]);
      BOOST_CHECK_EQUAL(it->b_mismatches[0], b_starts[i]);
    } else if (i == 2) {
      BOOST_CHECK_EQUAL(it->a_mismatches.size(), 2ul);
      BOOST_CHECK_EQUAL(it->b_mismatches.size(), 2ul);
      BOOST_CHECK_EQUAL(it->a_mismatches[0], a_starts[i]+1);
      BOOST_CHECK_EQUAL(it->b_mismatches[0], b_starts[i]+1);
      BOOST_CHECK_EQUAL(it->a_mismatches[1], a_starts[i]+3);
      BOOST_CHECK_EQUAL(it->b_mismatches[1], b_starts[i]+3);
    }
  }
  BOOST_CHECK(it == end);
}


BOOST_AUTO_TEST_CASE(several_intervals_reverse_strand)
{
  pslx_alignment al;
  std::string line = "47	1	0	0	3	7	0	0	-	Locus_54_Transcript_14/14_Confidence_0.000_Length_639	639	450	505	Locus_651_Transcript_1/4_Confidence_0.833_Length_280	280	229	277	4	"
    "5,4,33,6,\t"         // block sizes
    "134,141,146,183,\t"  // a starts
    "229,234,238,271,\t"  // b starts
    "ggggt,gggg,ggggggggggggggggggggggggggggggggg,gtggtg,\t" // a segment seqs
    "ggggg,tggg,ggggggggggggggggggggggggggggggggg,gggtgg,";  // b segment seqs
  al.parse_line(line);
  pslx_alignment::segments_type segs = al.segments("", "");
  pslx_alignment::segments_type::const_iterator it = segs.begin(), end = segs.end();
  size_t block_sizes[] = {5,4,33,6};
  size_t a_starts[] = {639-(134+1),639-(141+1),639-(146+1),639-(183+1)};
  size_t b_starts[] = {229,234,238,271};
  for (size_t i = 0; i < 4; ++it, ++i) {
    BOOST_CHECK(it != end);
    BOOST_CHECK_EQUAL(it->a_start, a_starts[i]);
    BOOST_CHECK_EQUAL(it->a_end,   a_starts[i] - (block_sizes[i] - 1));
    BOOST_CHECK_EQUAL(it->b_start, b_starts[i]);
    BOOST_CHECK_EQUAL(it->b_end,   b_starts[i] + (block_sizes[i] - 1));
    if (i == 0) {
      BOOST_CHECK_EQUAL(it->a_mismatches.size(), 1ul);
      BOOST_CHECK_EQUAL(it->b_mismatches.size(), 1ul);
      BOOST_CHECK_EQUAL(it->a_mismatches[0], a_starts[i]-4);
      BOOST_CHECK_EQUAL(it->b_mismatches[0], b_starts[i]+4);
    } else if (i == 1) {
      BOOST_CHECK_EQUAL(it->a_mismatches.size(), 1ul);
      BOOST_CHECK_EQUAL(it->b_mismatches.size(), 1ul);
      BOOST_CHECK_EQUAL(it->a_mismatches[0], a_starts[i]);
      BOOST_CHECK_EQUAL(it->b_mismatches[0], b_starts[i]);
    } else if (i == 2) {
      BOOST_CHECK_EQUAL(it->a_mismatches.size(), 0ul);
      BOOST_CHECK_EQUAL(it->b_mismatches.size(), 0ul);
    } else if (i == 3) {
      BOOST_CHECK_EQUAL(it->a_mismatches.size(), 3ul);
      BOOST_CHECK_EQUAL(it->b_mismatches.size(), 3ul);
      BOOST_CHECK_EQUAL(it->a_mismatches[0], a_starts[i]-1);
      BOOST_CHECK_EQUAL(it->b_mismatches[0], b_starts[i]+1);
      BOOST_CHECK_EQUAL(it->a_mismatches[1], a_starts[i]-3);
      BOOST_CHECK_EQUAL(it->b_mismatches[1], b_starts[i]+3);
      BOOST_CHECK_EQUAL(it->a_mismatches[2], a_starts[i]-4);
      BOOST_CHECK_EQUAL(it->b_mismatches[2], b_starts[i]+4);
    }
  }
  BOOST_CHECK(it == end);
}
