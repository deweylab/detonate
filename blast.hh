#pragma once
#include <vector>
#include <string>
#include <boost/iterator.hpp>
#include "lazycsv.hh"
#include "line_stream.hh"
#include "alignment_segment.hh"

namespace detail
{
  class blast_alignment_input_stream;
  class blast_alignment_segments;
}

////////////////////////////////////////////////////////////////////////////
// blast_alignment represents a single blast alignment. It realizes the
// Alignment concept.
////////////////////////////////////////////////////////////////////////////

class blast_alignment
{
public:
  // Realization of Alignment concept
  void parse_line(const std::string& line) { lazy_csv.parse_line(line); }
  std::string a_name()        const { return lazy_csv.at<std::string>(0); }
  std::string b_name()        const { return lazy_csv.at<std::string>(2); }
  double      frac_identity() const { return 1.0 * nident() * (qframe() == 0 ? 1 : 3) / qlen(); }
  double      frac_indel()    const { return 1.0 * gaps()   * (qframe() == 0 ? 1 : 3) / qlen(); }
  typedef detail::blast_alignment_input_stream input_stream_type;
  typedef detail::blast_alignment_segments     segments_type;
  segments_type segments(const std::string& a, const std::string& b) const; // defined below

  double      frac_identity_reverse() const { return 1.0 * nident() / slen(); }

  // For use by blast-specific algorithms:
  inline  int         qlen  ()        const { return lazy_csv.at<int>        ( 1); }
  inline  int         slen  ()        const { return lazy_csv.at<int>        ( 3); }
  inline  int         qstart()        const { return lazy_csv.at<int>        ( 4) - 1; } // 1-based -> 0-based
  inline  int         qend  ()        const { return lazy_csv.at<int>        ( 5) - 1; }
  inline  int         sstart()        const { return lazy_csv.at<int>        ( 6) - 1; }
  inline  int         send  ()        const { return lazy_csv.at<int>        ( 7) - 1; }
  inline  std::string qseq  ()        const { return lazy_csv.at<std::string>( 8); }
  inline  std::string sseq  ()        const { return lazy_csv.at<std::string>( 9); }
  inline  double      evalue()        const { return lazy_csv.at<double>     (10); }
  inline  int         nident()        const { return lazy_csv.at<int>        (15); }
  inline  int         gaps  ()        const { return lazy_csv.at<int>        (19); }
  inline  int         qframe()        const { return lazy_csv.at<int>        (21); }
  inline  int         sframe()        const { return lazy_csv.at<int>        (22); }

private:
  lazycsv<23, '\t'> lazy_csv;
};

namespace detail
{
  ////////////////////////////////////////////////////////////////////////////
  // blast_alignment_input_stream is an input stream that produces
  // blast_alignment objects, one per line of the input stream.
  ////////////////////////////////////////////////////////////////////////////

  class blast_alignment_input_stream
  {
  public:
    blast_alignment_input_stream(std::istream& is) : ls(is) {}

    blast_alignment_input_stream& operator>>(blast_alignment& al)
    {
      if (ls >> line)
        al.parse_line(line);
      return *this;
    }

    inline operator bool() { return ls; }
    inline bool operator!() { return !ls; }

  private:
    line_stream ls;
    std::string line;
  };

  ////////////////////////////////////////////////////////////////////////////
  // blast_alignment_segment_iterator is a forward iterator over alignment
  // segments of a blast alignment. This class should never be instantiated
  // directly by the end user. Instead, use the segments method of the
  // Alignment concept to get a blast_alignment_segments proxy, and use that
  // proxy's begin and end methods to get an iterator pair.
  ////////////////////////////////////////////////////////////////////////////

  class blast_alignment_segment_iterator
    : public boost::iterator_facade<blast_alignment_segment_iterator,
                                    const alignment_segment,
                                    std::forward_iterator_tag>
  {
  public:
    blast_alignment_segment_iterator() : at_end(true), al(NULL) { /* std::cout<< "in default constructor" << std::endl; */ }

    blast_alignment_segment_iterator(const blast_alignment& al_) : at_end(false), al(&al_)
    {
      a_seq        = al->qseq();                         b_seq        = al->sseq();
      a_start      = al->qstart();                       b_start      = al->sstart();
      a_end        = al->qend();                         b_end        = al->send();
      a_is_nucl    = (al->qframe() != 0);                b_is_nucl    = (al->sframe() != 0);
      a_incr_small = (a_start <= a_end) ? +1 : -1;       b_incr_small = (b_start <= b_end) ? +1 : -1; 
      a_incr       = a_incr_small * (a_is_nucl ? 3 : 1); b_incr       = b_incr_small * (b_is_nucl ? 3 : 1); 
      a_pos        = a_start;                            b_pos        = b_start;
      i = 0;

      if (a_seq.size() != b_seq.size())
        throw std::runtime_error("Failed assumption that a_seq and b_seq are the same length.");
      if (a_seq.size() == 0)
        throw std::runtime_error("Failed assumption that a_seq and b_seq of positive length.");
      if (a_seq[0] == '-' || a_seq[0] == '-')
        throw std::runtime_error("Failed assumption that alignment does not start with a gap.");
      if (a_is_nucl && b_is_nucl)
        throw std::runtime_error("Either a or b (or both) needs to be of type prot if you use blast alignments.");

      increment();
    }

  private:
    friend class boost::iterator_core_access;

    // ITERATOR FACADE PRIVATE INTERFACE STUFF

    void increment()
    {
      // Precondition: (1) (a_seq[i], b_seq[i], a_pos, b_pos) are at start of segment, or (2) at end of string.
      // Postcondition: (1) they are at the start of the next segment or (and seg is filled in properly), or (2) at_end == true.

      // check for end
      if (i == a_seq.size()) {
        at_end = true;
        return;
      }

      seg.a_mismatches.clear();
      seg.b_mismatches.clear();

      seg.a_start = a_pos;
      seg.b_start = b_pos;

      // look for end of segment and mismatches within it
      for (; i != a_seq.size(); ) {
        if (a_seq[i] == '-' || b_seq[i] == '-') // gap
          break;
        else {
          if (a_seq[i] != b_seq[i]) {
            seg.a_mismatches.push_back(a_pos);
            if (a_is_nucl) {
              seg.a_mismatches.push_back(a_pos + a_incr_small);
              seg.a_mismatches.push_back(a_pos + 2*a_incr_small);
            }
            seg.b_mismatches.push_back(b_pos);
            if (b_is_nucl) {
              seg.b_mismatches.push_back(b_pos + b_incr_small);
              seg.b_mismatches.push_back(b_pos + 2*b_incr_small);
            }
          }
        }
        ++i;
        a_pos += a_incr;
        b_pos += b_incr;
      }

      seg.a_end = a_pos - a_incr_small;
      seg.b_end = b_pos - b_incr_small;

      // skip past gap
      for (; i != a_seq.size(); ) {
        if      (a_seq[i] != '-' && b_seq[i] != '-') break;
        else if (a_seq[i] == '-' && b_seq[i] != '-') b_pos += b_incr;
        else if (a_seq[i] != '-' && b_seq[i] == '-') a_pos += a_incr;
        else throw std::runtime_error("Unexpected alignment of gap to gap.");
        ++i;
      }
    }

    bool equal(const blast_alignment_segment_iterator& other) const
    {
      if (!at_end && other.at_end)
        return false;
      else if (at_end && other.at_end)
        return true;
      else
        throw std::runtime_error("General equal operator is not implemented.");
    }

    const alignment_segment& dereference() const { return seg; }

    // BOOKKEEPING DATA

    bool at_end;
    const blast_alignment *al;
    alignment_segment seg;

    std::string a_seq,        b_seq;
    int         a_start,      b_start;
    int         a_end,        b_end;
    bool        a_is_nucl,    b_is_nucl;
    int         a_incr_small, b_incr_small;
    int         a_incr,       b_incr;
    int         a_pos,        b_pos;
    size_t i;
  };

  ////////////////////////////////////////////////////////////////////////////
  // blast_alignment_segments is a proxy class, which looks like a standard
  // container that can be iterated over.
  ////////////////////////////////////////////////////////////////////////////

  class blast_alignment_segments
  {
  public:
    blast_alignment_segments(const blast_alignment& al) : al(al) {}
    typedef const alignment_segment&         const_reference;
    typedef blast_alignment_segment_iterator const_iterator;
    const_iterator begin() const { return blast_alignment_segment_iterator(al); }
    const_iterator end()   const { return blast_alignment_segment_iterator(); }

  private:
    const blast_alignment& al;
  };

} // namespace detail

blast_alignment::segments_type blast_alignment::segments(const std::string& /*a*/, const std::string& /*b*/) const
{
  return blast_alignment::segments_type(*this);
}

// SUMMARY OF BLAST FORMAT:
//
// std::string qseqid;   // 0
// int         qlen;     // 1
// std::string sseqid;   // 2
// int         slen;     // 3
// int         qstart;   // 4
// int         qend;     // 5
// int         sstart;   // 6
// int         send;     // 7
// std::string qseq;     // 8
// std::string sseq;     // 9
// double      evalue;   // 10
// int         bitscore; // 11
// int         score;    // 12
// int         length;   // 13
// double      pident;   // 14
// int         nident;   // 15
// int         mismatch; // 16
// int         positive; // 17
// int         gapopen;  // 18
// int         gaps;     // 19
// double      ppos;     // 20
// int         qframe;   // 21 // This is -1, -2, -3, 1, 2, 3 if the query is nucleotide, and 0 if the query is protein.
// int         sframe;   // 22 // Similarly for the subject.
