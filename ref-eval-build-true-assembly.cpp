// Copyright (c) 2013
// Nathanael Fillmore and Bo Li (University of Wisconsin-Madison)
// nathanae@cs.wisc.edu
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <vector>
#include <list>
#include <string>
#include <stdio.h>
#include <time.h>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include "sam/bam.h"
#include "fasta.hh"
#include "util.hh"
#include "re_bta_help.hh"

struct probs_and_lidxs {
  std::vector<size_t> lidxs;
  std::vector<double> probs;
};

struct ReadStruct
{
  int rid, tid, pos, readlen;
  int cid, rpos; // id in contigs and right most end's contig position

  ReadStruct(int rid, int tid, int pos, int readlen) {
    this->rid = rid;
    this->tid = tid;
    this->pos = pos;
    this->readlen = readlen;
    cid = rpos = -1;
  }

  bool operator< (const ReadStruct& o) const {
    if (tid != o.tid) return tid < o.tid;
    return pos < o.pos;
  }
};

struct ContigStruct
{
  int tid, pos, len;

  ContigStruct(int tid, int pos, int len) {
    this->tid = tid;
    this->pos = pos;
    this->len = len;
  }
};

struct BreakPointStruct
{
  int cid, pos, overlap;

  BreakPointStruct(int cid, int pos, int overlap) {
    this->cid = cid;
    this->pos = pos;
    this->overlap = overlap;
  }
};

bool is_ok(const bam1_t *rec)
{
  // // Exclude ``right'' mates.
  // // we only need to use one of the pairs of reads; let's choose the
  // // left/upstream ones
  // // XXX should we really do this?
  // if (rec->core.flag & BAM_FREAD2) // this is read2
  //   return false;

  // Exclude noise reads.
  if (rec->core.qual == 0) // this is a noise read
    return false;
  if (rec->core.tid < 0)
    throw std::runtime_error("Found an alignment with qual > 0 "
        "but tid < 0. Was your BAM file produced by RSEM?");

  return true;
}

// Calculate the posterior probability that this read follows this alignment,
// as calculated by RSEM.
double compute_prob(const bam1_t *rec)
{
  // According to rsem-calculate-expression -?, ``The MAPQ field of each
  // alignment is set to min(100, floor(-10 * log10(1.0 - w) + 0.5)), where w
  // is the posterior probability of that alignment being the true mapping of
  // a read.''
  //
  // To go backwards, we solve:
  //        min(100, floor(-10 * log10(1.0 - w) + 0.5) = q
  // i.e.,           round(-10 * log10(1.0 - w))       = q
  // i.e.,                 -10 * log10(1.0 - w)        = q
  // i.e.,                       log10(1.0 - w)        = -q/10
  // i.e.,                             1.0 - w         = 10^(-q/10)
  // i.e.,                                 - w         = 10^(-q/10) - 1.0
  // i.e.,                                   w         = -( 10^(-q/10) - 1.0 )
  //
  // As a doublecheck: In sam_rsem_cvt.h, there is the following function:
  //
  //   uint8_t getMAPQ(double val) {
  //       double err = 1.0 - val;
  //       if (err <= 1e-10) return 100;
  //       return (uint8_t)(-10 * log10(err) + .5); // round it
  //   }
  //
  // If we plug in val = -(10^(-q/10) - 1.0), we get
  //
  //    -10 * log10(err)
  // =  -10 * log10(1.0 - val)
  // =  -10 * log10(1.0 - (-(10^(-q/10) - 1.0)))
  // =  -10 * log10(1.0 + 10^(-q/10) - 1.0)
  // =  -10 * log10(      10^(-q/10)      )
  // =  -10 * (               -q/10       )
  // =  -     (               -q          )
  // =                         q           
  //
  // as expected if the above is correct.
  double prob = -(pow(10.0, -rec->core.qual/10.0) - 1.0);
  if (prob < 0.0 || prob > 1.0)
    throw std::runtime_error("Posterior probability of an alignment is < 0 or > 1.");
  return prob;
}

size_t extract_alignment_probs(
    std::map<std::string, probs_and_lidxs>& read_to_probs_and_lidxs,
    const std::string& bam_fname,
    double min_alignment_prob)
{
  // Open BAM file and read header.
  bamFile fi = bam_open(bam_fname.c_str(), "r");
  if (!fi)
    throw std::runtime_error("Cannot open BAM file: " + bam_fname);
  bam_header_t *header = bam_header_read(fi);

  // Iterate over lines of the bamfile.
  bam1_t *rec = bam_init1();
  size_t lidx;
  for (lidx = 0;; ++lidx) { // alignment idx
    
    // Read next record. If no more records, we're done.
    if (bam_read1(fi, rec) < 0)
      break;

    // Check that this isn't an alignment to the noise transcript, etc.
    if (!is_ok(rec))
      continue;

    // Calculate the posterior probability that this read follows this
    // alignment, as calculated by RSEM.
    double prob = compute_prob(rec);

    // Check that the posterior probability is big enough.
    if (prob < min_alignment_prob)
      continue;

    // Actually add the alignment.
    const char *read_name = bam1_qname(rec);
    probs_and_lidxs& pal = read_to_probs_and_lidxs[read_name];
    pal.probs.push_back(prob);
    pal.lidxs.push_back(lidx);

  }

  // Clean up.
  bam_destroy1(rec);
  bam_header_destroy(header);

  // Return the number of alignments.
  return lidx;
}

void choose_alignments(
    std::vector<bool>& is_lidx_chosen,
    boost::random::mt19937& rng,
    const std::map<std::string, probs_and_lidxs>& read_to_probs_and_lidxs,
    const std::string& alignment_policy)
{
  std::map<std::string, probs_and_lidxs>::const_iterator it;

  if (alignment_policy == "sample") {
    for (it = read_to_probs_and_lidxs.begin();
         it != read_to_probs_and_lidxs.end();
         ++it) {
      boost::random::discrete_distribution<> dd(it->second.probs);
      size_t i = dd(rng);
      size_t chosen_lidx = it->second.lidxs[i];
      is_lidx_chosen[chosen_lidx] = true;
    }
  }

  else if (alignment_policy == "best") {
    for (it = read_to_probs_and_lidxs.begin();
         it != read_to_probs_and_lidxs.end();
         ++it) {
      size_t num_alignments = it->second.probs.size();
      size_t best = 0;
      if (num_alignments > 1) {
        for (size_t i = 1; i < num_alignments; ++i)
          if (it->second.probs[i] > it->second.probs[best])
            best = i;
      }
      size_t chosen_lidx = it->second.lidxs[best];
      is_lidx_chosen[chosen_lidx] = true;
    }
  }

  else if (alignment_policy == "all") {
    for (it = read_to_probs_and_lidxs.begin();
         it != read_to_probs_and_lidxs.end();
         ++it) {
      size_t num_alignments = it->second.probs.size();
      for (size_t i = 0; i < num_alignments; ++i)
        is_lidx_chosen[it->second.lidxs[i]] = true;
    }
  }
}

void extract_reads(
    std::vector<ReadStruct>& reads,
    const std::vector<bool>& is_lidx_chosen,
    const std::string& bam_fname)
{
  // Open BAM file and read header.
  bamFile fi = bam_open(bam_fname.c_str(), "r");
  if (!fi)
    throw std::runtime_error("Cannot open BAM file: " + bam_fname);
  bam_header_t *header = bam_header_read(fi);

  // Iterate over lines of the bamfile.
  bam1_t *rec = bam_init1();
  for (size_t lidx = 0;; ++lidx) { // alignment idx

    // Read next record. If no more records, we're done.
    if (bam_read1(fi, rec) < 0)
      break;

    // Check that this alignment has been chosen.
    if (!is_lidx_chosen[lidx])
      continue;

    // Extract the transcript idx, the leftmost position of the read, and
    // the read length.
    int tid = rec->core.tid;
    int pos = rec->core.pos;
    int readlen = rec->core.l_qseq;

    // Record the above info about this read. This is the main info that
    // will be relevant for making the contigs later on.
    reads.push_back(ReadStruct(-1, tid, pos, readlen));

  }

  // Clean up.
  bam_destroy1(rec);
  bam_header_destroy(header);
}

void assemble(
  std::vector<ContigStruct>& contigs,
  std::vector<ReadStruct>& reads,
  int min_overlap)
{
  int s = -1;
  size_t num_reads = reads.size();

  for (size_t i = 0; i < num_reads; i++) {
    if (contigs.size() > 0 && contigs[s].tid == reads[i].tid) {
      assert(contigs[s].pos <= reads[i].pos);
      assert(contigs[s].pos + contigs[s].len <= reads[i].pos + reads[i].readlen);
    }
    if (contigs.size() == 0 ||
        contigs[s].tid != reads[i].tid ||
        contigs[s].pos + contigs[s].len - reads[i].pos < min_overlap) {
      ++s;
      contigs.push_back(ContigStruct(reads[i].tid, reads[i].pos, reads[i].readlen));
    }
    else {
      contigs[s].len = reads[i].pos + reads[i].readlen - contigs[s].pos;
    }

    reads[i].cid = s;
    reads[i].rpos = contigs[s].len - 1;
  }

  ++s; // there are s contigs

  // if (breakpoints != NULL) {
  //   for (int i = 0; i < num_reads - 1; i++) {
  //     // read is the end of a contig
  //     if (reads[i].rpos + 1 == contigs[reads[i].cid].len)
  //       continue;
  //     // there is another read coincide with this read
  //     if (reads[i].rpos == reads[i + 1].rpos)
  //       continue;
  //     assert(reads[i].cid == reads[i + 1].cid);
  //     breakpoints->push_back(BreakPointStruct(reads[i].cid, reads[i].rpos, readLen - (reads[i + 1].rpos - reads[i].rpos)));
  //   }
  // }
}

void output(
  const std::string& output_fname,
  const std::vector<ContigStruct>& contigs,
  const fasta& ref_fa,
  const std::string& bam_fname)
{
  // Open BAM file and read header.
  bamFile fi = bam_open(bam_fname.c_str(), "r");
  if (!fi)
    throw std::runtime_error("Cannot open BAM file: " + bam_fname);
  bam_header_t *header = bam_header_read(fi);

  // Open file to write output.
  std::ofstream out(output_fname.c_str());
  if (!out)
    throw std::runtime_error("Cannot open file to write assembly: " + output_fname);

  // Assemble each contig.
  for (size_t i = 0; i < contigs.size(); i++) {

    // We have recorded, in contigs[i].tid, the idx of the transcript that this
    // contig comes from, relative to the ordering of transcripts in the BAM
    // file. However, the transcript sequences are stored in the FASTA file,
    // and the order here could be different. Thus, here, we will figure out
    // the idx of the transcript relative to the order in the FASTA file.
    std::string tname(header->target_name[contigs[i].tid]);
    std::map<std::string, size_t>::const_iterator it;
    it = ref_fa.names_to_idxs.find(tname);
    if (it == ref_fa.names_to_idxs.end())
      throw std::runtime_error("The following transcript name does not "
          "correspond to a reference sequence: " + tname);
    size_t tidx = it->second;

    // Get the transcript sequence and extract the contig.
    const std::string& refseq = ref_fa.seqs[contigs[i].tid];
    if (contigs[i].pos + contigs[i].len > static_cast<int>(refseq.length())) {
      std::ostringstream ss;
      ss << "Contig extends past end of the transcript sequence. "
         << "Info for debugging: "
         << "tname="  << tname << ", " 
         << "pos="    << contigs[i].pos << ", " 
         << "length=" << contigs[i].len << ", " 
         << "transcript_length=" << refseq.length();
      throw std::runtime_error(ss.str());
    }
    std::string seq = refseq.substr(contigs[i].pos, contigs[i].len);

    // Output the contig in FASTA format.
    out << ">contig_" << i
        << " tname="  << tname
        << " pos="    << contigs[i].pos 
        << " length=" << contigs[i].len << "\n";
    out << seq << "\n";

  }

  // Clean up.
  bam_header_destroy(header);

  // if (outputBreakPoints) {
  //   sprintf(breakPointsF, "%s_%d.bp", breakPointsFN, windowSize);
  //   fo = fopen(breakPointsF, "w");
  //   if (!fo)
  //     throw std::runtime_error(std::string("Cannot open file to write breakpoints: ") + breakPointsF);
  //   int nbp = (int)breakPoints.size();
  //   fprintf(fo, "%d\n", nbp);
  //   for (int i = 0; i < nbp; i++) 
  //     fprintf(fo, "%d %d %d\n", breakPoints[i].cid, breakPoints[i].pos, breakPoints[i].overlap);
  //   fclose(fo);
  //   printf("All break points are written to %s!\n", breakPointsF);
  // }
}

void print_help()
{
  std::cout << get_help_string() << std::endl;
}

boost::program_options::options_description describe_options()
{
  namespace po = boost::program_options;
  po::options_description desc;
  desc.add_options()
    ("help,?", "Display this information.")
    ("reference", po::value<std::string>())
    ("expression", po::value<std::string>())
    ("assembly", po::value<std::string>())
    ("min-overlap", po::value<std::string>())
    ("min-alignment-prob", po::value<double>())
    ("alignment-policy", po::value<std::string>())
    ;
  return desc;
}

struct opts
{
  std::string reference;
  std::string expression;
  std::string assembly;
  size_t min_overlap_begin;
  size_t min_overlap_end;
  double min_alignment_prob;
  std::string alignment_policy;
};

void parse_options(opts& o, const boost::program_options::variables_map& vm)
{
  namespace po = boost::program_options;

  if (!vm.count("reference"))
    throw po::error("--reference is required.");
  o.reference = vm["reference"].as<std::string>();

  if (!vm.count("expression"))
    throw po::error("--expression is required.");
  o.expression = vm["expression"].as<std::string>();

  if (!vm.count("assembly"))
    throw po::error("--assembly is required.");
  o.assembly = vm["assembly"].as<std::string>();

  std::string mo;
  if (vm.count("min-overlap"))
    mo = vm["min-overlap"].as<std::string>();
  else
    mo = "0";

  size_t tpos = mo.find_first_of(',');
  if (tpos == std::string::npos) {
    o.min_overlap_begin = boost::lexical_cast<size_t>(mo);
    o.min_overlap_end   = o.min_overlap_begin;
  }
  else {
    o.min_overlap_begin = boost::lexical_cast<size_t>(mo.substr(0, tpos));
    o.min_overlap_end   = boost::lexical_cast<size_t>(mo.substr(tpos + 1));
  }

  if (vm.count("alignment-policy")) {
    o.alignment_policy = vm["alignment-policy"].as<std::string>();
    if (o.alignment_policy != "random" &&
        o.alignment_policy != "best" &&
        o.alignment_policy != "all")
      throw po::error("Invalid value for --alignment-policy: " + o.alignment_policy);
  } else {
    o.alignment_policy = "random";
  }

  if (vm.count("min-alignment-prob")) {
    o.min_alignment_prob = vm["min-alignment-prob"].as<double>();
    if (o.min_alignment_prob < 0.0 || o.min_alignment_prob > 1.0)
      throw po::error("Invalid value for --min-alignment-prob: " + vm["min-alignment-prob"].as<std::string>());
  } else {
    o.min_alignment_prob = 0.0;
  }

}

int main(int argc, const char **argv)
{
  namespace po = boost::program_options;

  try {

    po::options_description desc = describe_options();
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (argc == 1 || vm.count("help")) {
      print_help();
      exit(0);
    }

    opts o;
    parse_options(o, vm);
    notify(vm);

    // Read transcripts and extract a map from transcript names to transcript
    // idxs. The idxs start at 1, with 0 being the noise transcript.
    std::cerr << curtime() << "Reading the transcripts..." << std::flush;
    fasta ref_fa;
    std::string ref_fname = o.reference + ".transcripts.fa";
    read_fasta(ref_fa, ref_fname);
    std::cerr << "done." << std::endl;

    // Read alignments and extract the posterior probability of each alignment.
    std::cerr << curtime() << "Reading the alignments and extracting the posterior probability of each alignment..." << std::flush;
    std::map<std::string, probs_and_lidxs> read_to_probs_and_lidxs;
    std::string bam_fname = o.expression + ".transcript.sorted.bam";
    size_t num_alignments = extract_alignment_probs(read_to_probs_and_lidxs, bam_fname, o.min_alignment_prob);
    std::cerr << "done." << std::endl;

    // Select a subset of alignments (currently: by sampling one alignment for
    // each read, according to the posterior probability of each alignment).
    std::cerr << curtime() << "Selecting a subset of alignments..." << std::flush;
    std::vector<bool> is_lidx_chosen(num_alignments);
    boost::random::mt19937 rng(time(0));
    choose_alignments(is_lidx_chosen, rng, read_to_probs_and_lidxs, o.alignment_policy);
    std::cerr << "done." << std::endl;

    // Read the alignments again and, for the reads chosen above, extract info
    // about the read locations relative to the reference transcripts.
    std::cerr << curtime() << "Reading the alignments and extracting info about read locations..." << std::flush;
    std::vector<ReadStruct> reads;
    extract_reads(reads, is_lidx_chosen, bam_fname);
    std::cerr << "done." << std::endl;

    for (size_t mo = o.min_overlap_begin; mo <= o.min_overlap_end; ++mo) {
      // Assemble the reads into contigs.
      std::cerr << curtime() << "Assembling the reads into contigs at minimum overlap " << mo << "..." << std::flush;
      std::vector<ContigStruct> contigs;
      assemble(contigs, reads, mo);

      // Output the contigs.
      std::ostringstream ss;
      ss << o.assembly << "_" << mo << ".fa";
      std::cerr << "writing " << contigs.size() << " contigs..." << std::flush;
      output(ss.str(), contigs, ref_fa, bam_fname);
      std::cerr << "done." << std::endl;
    }

    std::cerr << curtime() << "Done with everything." << std::endl;
    return 0;

  } catch (const boost::program_options::error& x) {

    std::cerr << std::endl;
    std::cerr << argv[0] << ": Error: " << x.what() << std::endl;
    std::cerr << "Check " << argv[0] << " --help for more information." << std::endl;
    return 1;

  } catch (const std::exception& x) {

    std::cerr << std::endl;
    std::cerr << argv[0] << ": Error: " << x.what() << std::endl;
    return 1;

  }
}
