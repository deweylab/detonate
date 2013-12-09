#include <string>
#include <iostream>
#define BOOST_TEST_MODULE test_blast
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include "blast.hh"

BOOST_AUTO_TEST_CASE(alignment_concept)
{
  blast_alignment al;
  std::string line = "comp12_c1_seq1	3219	gi|301621625|gb|XP_002940146|	1125	1155	2213	459	899	HLPSPVTAQKYRCELLYEGPPDDEAAMGIKNCDPKGPLMMYISKMVPTTDKGRFYAFGRVFSGVVSTGLKVRIMGPNYTPGKKEDLYLKPIQRT-------------------ILMMGRYVEPIEDVPCGNIVGLVGVDQFLVKTGTITTFEHAHNLRVMKFSVSPVVRVAVEAKNPADLPKLVEGLKRLAKSDPMVQCIIEESGEHIIAGAGELHLEICLKDLEEDHACIPIKKSDPVVSYRET-------------------------VAEESSV-----------LCLSKSPNKHNRLYMKARPFQDGLAEDIDKGEVAAR--------------------------QELKTRARYLAEKYEWDVTEARKIWCFGPDGTGPNMLTDITKGV----------------QYLNEIKDSVVAGFQWATKEGALCEENMRAVRFDVHDVTL	QLPAEQVPSSNECEM----PANKENTKGDKQTSHRDPEVVKPEKQ-EAENKSYFIAFARVFSGIVRRGQKIFVLGPKYDPA--ETLLKLPLNCTPACDLPGIPHMACCSLDNLYLLMGRELENLEEVPVGNVLGIGGLEEFVLKSATLSTSPACPPFIPLNFEATPIVRVAVEPKHPSEMPQLVKGMKLLNQADPCVEVLIQETGEHVLITAGEVHLQRCLDDLRERFAKIQVSSSAPIIPFRETIIRPPKVDMVNEEIGKQQKIAVIHQVKEEQSKCPEGVQVDPDGLVTLTTPNKLATLSVRAMPLPEEVTQLLEKNSDLIRTMEQINMALNEGSYTIHFNESAIERITAFKSNLQQLLQGRRWR-NAVDQIWSFGPRRYGPNILLNRIEGYDRPSVWQCLEKSIREGKYRN-FDNSIVSGFQLATLAGPMCEEPLMGVCFIVEKLDL	8e-46	 180	456	450	28.89	130	214	203	10	106	45.11	3	0";
  al.parse_line(line);

  BOOST_CHECK_EQUAL(al.a_name(), "comp12_c1_seq1");
  BOOST_CHECK_EQUAL(al.b_name(), "gi|301621625|gb|XP_002940146|");
  BOOST_CHECK_EQUAL(al.frac_identity_wrt_a(), 1.0*130*3/3219);
  BOOST_CHECK_EQUAL(al.frac_indel_wrt_a(), 1.0*106*3/3219);
  BOOST_CHECK_EQUAL(al.frac_identity_wrt_b(), 1.0*130/1125);
  BOOST_CHECK_EQUAL(al.frac_indel_wrt_b(), 1.0*106/1125);
}

BOOST_AUTO_TEST_CASE(one_interval)
{
  blast_alignment al;
  std::string line = "comp3_c0_seq1	932	gi|89268196|gb|CAJ82654|	140	479	877	8	140	GGSGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG	5e-89	 266	680	133	99.25	132	1	133	0	0	100.00	2	0";
  al.parse_line(line);

  // Check non-segment parts of alignment concept
  // todo

  blast_alignment::segments_type segs = al.segments("", "");
  blast_alignment::segments_type::const_iterator it = segs.begin(), end = segs.end();
  BOOST_CHECK_EQUAL(it != end, true);
  BOOST_CHECK_EQUAL(it->a_start, 478ul);
  BOOST_CHECK_EQUAL(it->a_end, 876ul);
  BOOST_CHECK_EQUAL(it->b_start, 7ul);
  BOOST_CHECK_EQUAL(it->b_end, 139ul);
  BOOST_CHECK_EQUAL(it->a_mismatches.size(), 3ul);
  BOOST_CHECK_EQUAL(it->a_mismatches[0], 478ul + 2*3ul);
  BOOST_CHECK_EQUAL(it->a_mismatches[1], 478ul + 2*3ul + 1ul);
  BOOST_CHECK_EQUAL(it->a_mismatches[2], 478ul + 2*3ul + 2ul);
  BOOST_CHECK_EQUAL(it->b_mismatches.size(), 1ul);
  BOOST_CHECK_EQUAL(it->b_mismatches[0], 7ul + 2ul);
  ++it;
  BOOST_CHECK_EQUAL(it == end, true);
}

BOOST_AUTO_TEST_CASE(several_intervals)
{
  blast_alignment al;
  std::string line = "comp12_c0_seq1	3219	gi|301610695|gb|XP_002934886|	636	3079	2756	23	116\t"
  // 0         1         2           3         4         5         6         7         8         9         0
  // 012345678901234567890123456789  012345678901234567890123456789012345678901234567890123456789012345678901234567
    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII--SSSSSSSSSSSSSSSSSSSSSSSSSSSSSFFFFFFFFFFFFFFFFDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\t"
    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIILLSSSSSSSSSSSSSSSSSSSSSSSSSSSSS----------------DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\t"
  // 0123456789012345678901234567890123456789012345678901234567890                123456789012345678901234567890123
  // 0         1         2         3         4         5         6                         7         8         9
    "2e-11	68.9	167	110	39.09	43	49	61	2	18	55.45	-3	0";
  al.parse_line(line);
  blast_alignment::segments_type segs = al.segments("", "");
  blast_alignment::segments_type::const_iterator it = segs.begin(), end = segs.end();
  BOOST_CHECK_EQUAL(it != end, true);
  BOOST_CHECK_EQUAL(it->a_start, 3078ul);
  BOOST_CHECK_EQUAL(it->a_end,   3078ul - 3*29ul - 2);
  BOOST_CHECK_EQUAL(it->b_start, 22ul);
  BOOST_CHECK_EQUAL(it->b_end,   22ul + 29ul);
  ++it;
  BOOST_CHECK_EQUAL(it != end, true);
  BOOST_CHECK_EQUAL(it->a_start, 3078ul - 3*30ul);
  BOOST_CHECK_EQUAL(it->a_end,   3078ul - 3*58ul - 2);
  BOOST_CHECK_EQUAL(it->b_start, 22ul + 32ul);
  BOOST_CHECK_EQUAL(it->b_end,   22ul + 60ul);
  ++it;
  BOOST_CHECK_EQUAL(it != end, true);
  BOOST_CHECK_EQUAL(it->a_start, 3078ul - 3*75ul);
  BOOST_CHECK_EQUAL(it->a_end,   3078ul - 3*107ul - 2);
  BOOST_CHECK_EQUAL(it->b_start, 22ul + 61ul);
  BOOST_CHECK_EQUAL(it->b_end,   22ul + 93ul);
  ++it;
  BOOST_CHECK_EQUAL(it == end, true);
}

BOOST_AUTO_TEST_CASE(iterator_equality)
{
  // The details of this alignment are irrelevant except that there are 3
  // segments.
  blast_alignment al;
  std::string line = "comp12_c0_seq1	3219	gi|301610695|gb|XP_002934886|	636	3079	2756	23	116\t"
  // 0         1         2           3         4         5         6         7         8         9         0
  // 012345678901234567890123456789  012345678901234567890123456789012345678901234567890123456789012345678901234567
    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIII--SSSSSSSSSSSSSSSSSSSSSSSSSSSSSFFFFFFFFFFFFFFFFDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\t"
    "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIILLSSSSSSSSSSSSSSSSSSSSSSSSSSSSS----------------DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD\t"
  // 0123456789012345678901234567890123456789012345678901234567890                123456789012345678901234567890123
  // 0         1         2         3         4         5         6                         7         8         9
    "2e-11	68.9	167	110	39.09	43	49	61	2	18	55.45	-3	0";
  al.parse_line(line);
  blast_alignment::segments_type segs = al.segments("", "");

  // it1, it2 = seg 0, seg 0
  blast_alignment::segments_type::const_iterator it1 = segs.begin(), it2 = segs.begin(), end = segs.end();
  BOOST_CHECK(it1 == it2);
  BOOST_CHECK(it1 != end);
  BOOST_CHECK(it2 != end);

  // it1, it2 = seg 1, seg 0
  ++it1;
  BOOST_CHECK(it1 != it2);
  BOOST_CHECK(it1 != end);

  // it1, it2 = seg 1, seg 1
  ++it2;
  BOOST_CHECK(it1 == it2);
  BOOST_CHECK(it1 != end);

  // it1, it2 = seg 2, seg 1
  ++it1;
  BOOST_CHECK(it1 != it2);
  BOOST_CHECK(it1 != end);

  // it1, it2 = seg 2, seg 2
  ++it2;
  BOOST_CHECK(it1 == it2);
  BOOST_CHECK(it2 != end);

  // it1, it2 = end, seg 2
  ++it1;
  BOOST_CHECK(it1 != it2);
  BOOST_CHECK(it1 == end);
}
