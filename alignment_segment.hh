#pragma once
#include <vector>
#include <iostream>
#include <boost/foreach.hpp>

struct alignment_segment
{
  size_t a_start, a_end, b_start, b_end; // the segment is [a_start, a_end] -> [b_start, b_end]
  std::vector<size_t> a_mismatches, b_mismatches;
};

std::ostream& operator<<(std::ostream& out, const alignment_segment& seg)
{
  out << "(a_start=" << seg.a_start
      << ", a_end=" << seg.a_end
      << ", b_start=" << seg.b_start
      << ", b_end=" << seg.b_end
      << ", a_mismatches=( ";
  BOOST_FOREACH(size_t x, seg.a_mismatches)
    out << x << " ";
  out << "), b_mismatches=( ";
  BOOST_FOREACH(size_t x, seg.b_mismatches)
    out << x << " ";
  out << "))";
  return out;
}

bool operator==(const alignment_segment& seg1, const alignment_segment& seg2)
{
  return seg1.a_start      == seg2.a_start      &&
         seg1.a_end        == seg2.a_end        &&
         seg1.b_start      == seg2.b_start      &&
         seg1.b_end        == seg2.b_end        &&
         seg1.a_mismatches == seg2.a_mismatches &&
         seg1.b_mismatches == seg2.b_mismatches;
}

struct segment_ops_wrt_b
{
  static inline size_t l(const alignment_segment& seg) { return std::min(seg.b_start, seg.b_end); }
  static inline size_t r(const alignment_segment& seg) { return std::max(seg.b_start, seg.b_end); }

  static inline size_t& b_start(alignment_segment& seg) { return seg.b_start; }
  static inline size_t& b_end  (alignment_segment& seg) { return seg.b_end  ; }

  static inline size_t& a_start(alignment_segment& seg) { return seg.a_start; }
  static inline size_t& a_end  (alignment_segment& seg) { return seg.a_end  ; }
};

struct segment_ops_wrt_a
{
  static inline size_t l(const alignment_segment& seg) { return std::min(seg.a_start, seg.a_end); }
  static inline size_t r(const alignment_segment& seg) { return std::max(seg.a_start, seg.a_end); }

  static inline size_t& b_start(alignment_segment& seg) { return seg.a_start; }
  static inline size_t& b_end  (alignment_segment& seg) { return seg.a_end  ; }

  static inline size_t& a_start(alignment_segment& seg) { return seg.b_start; }
  static inline size_t& a_end  (alignment_segment& seg) { return seg.b_end  ; }
};

template<typename H>
void set_l_after_r(alignment_segment& seg2, const alignment_segment& seg1)
{
  bool b2_is_forward = (H::b_start(seg2) <= H::b_end(seg2));
  bool a2_is_forward = (H::a_start(seg2) <= H::a_end(seg2));

  // lb1---------rb1             -->     lb1---------rb1      
  //       sb2>>>>>>>>>>>eb2                            sb2>>eb2
  if (b2_is_forward) {
    int delta = H::r(seg1) + 1 - H::b_start(seg2);
    H::b_start(seg2) += delta;

    //     sb2>>>>>>>>>>>eb2     -->                    sb2>>eb2
    //     |               |                            |      |
    //     sa2>>>>>>>>>>>ea2                            sa2>>ea2
    if (a2_is_forward)
      H::a_start(seg2) += delta;

    //     sb2>>>>>>>>>>>eb2     -->                    sb2>>eb2
    //     |               /                            |    /
    //      \-------------\                     /-----------/
    //     /-------------/ |                   |       /    
    //     ea2<<<<<<<<<<<sa2                   ea2<<sa2       
    else
      H::a_start(seg2) -= delta;
  }

  // lb1---------rb1             -->     lb1---------rb1      
  //       eb2<<<<<<<<<<<sb2                            eb2<<sb2
  else {
    int delta = H::r(seg1) + 1 - H::b_end(seg2);
    H::b_end(seg2) += delta;

    //     eb2<<<<<<<<<<<sb2     -->                    eb2<<sb2
    //     |               /                           /     /
    //      \-------------\                           |     /
    //     /-------------/ |                   /-----------/
    //     sa2>>>>>>>>>>>ea2                   sa2>>ea2
    if (a2_is_forward)
      H::a_end(seg2) -= delta;

    //     eb2<<<<<<<<<<<sb2     -->                    eb2<<sb2
    //     |               |                            |      |
    //     ea2<<<<<<<<<<<sa2                            ea2<<sa2
    else
      H::a_end(seg2) += delta;
  }
}

template<typename H>
void set_r_before_l(alignment_segment& seg2, const alignment_segment& seg1)
{
  bool b2_is_forward = (H::b_start(seg2) <= H::b_end(seg2));
  bool a2_is_forward = (H::a_start(seg2) <= H::a_end(seg2));

  //           lb1-------rb1                     lb1-----rb1
  //   sb2>>>>>>>>>>>eb2         -->     sb2>>eb2
  if (b2_is_forward) {
    int delta = H::b_end(seg2) + 1 - H::l(seg1);
    H::b_end(seg2) -= delta;

    // sb2>>>>>>>>>>>eb2         -->     sb2>>eb2
    // |               |                 |      |
    // sa2>>>>>>>>>>>ea2                 sa2>>ea2
    if (a2_is_forward)
      H::a_end(seg2) -= delta;

    // sb2>>>>>>>>>>>eb2         -->     sb2>>eb2
    // |               /                 \       \          .
    //  \-------------\                   --------------\   .
    // /-------------/ |                          |      |
    // ea2<<<<<<<<<<<sa2                          ea2<<sa2
    else
      H::a_end(seg2) += delta;
  }

  //           lb1-------rb1                     lb1-----rb1
  //   eb2<<<<<<<<<<<sb2         -->     eb2<<sb2
  else {
    int delta = H::b_start(seg2) + 1 - H::l(seg1);
    H::b_start(seg2) -= delta;

    // eb2<<<<<<<<<<<sb2          -->     eb2<<sb2
    // |               /                  \       \          .
    //  \-------------\                    --------------\   .
    // /-------------/ |                           |      |
    // sa2>>>>>>>>>>>ea2                           sa2>>ea2
    if (a2_is_forward)
      H::a_start(seg2) += delta;

    // eb2<<<<<<<<<<<sb2          -->     eb2<<sb2
    // |               |                  |      |
    // ea2>>>>>>>>>>>sa2                  ea2<<sa2
    else
      H::a_start(seg2) -= delta;
  }
}

template<typename H>
bool intersects(const std::vector<alignment_segment>& segs2,
                const std::vector<alignment_segment>& segs1)
{
  BOOST_FOREACH(const alignment_segment& seg1, segs1) {
    BOOST_FOREACH(const alignment_segment& seg2, segs2) {

      // these conditions are copied from subtract below, q.v.

      if (H::l(seg1) <= H::l(seg2) && H::l(seg2) <= H::r(seg1) && H::r(seg1) <= H::r(seg2))
        return true;

      else if (H::l(seg2) <= H::l(seg1) && H::l(seg1) <= H::r(seg2) && H::r(seg2) <= H::r(seg1))
        return true;

      else if (H::l(seg1) <= H::l(seg2) && H::r(seg2) <= H::r(seg1))
        return true;

      else if (H::l(seg2) <= H::l(seg1) && H::r(seg1) <= H::r(seg2))
        return true;
    
    }
  }
  return false;
}

template<typename H>
void subtract_in_place(      std::vector<alignment_segment>& segs2,
                       const std::vector<alignment_segment>& segs1)
{
  std::vector<alignment_segment>::iterator it;
  std::vector<alignment_segment> new_segs;

  BOOST_FOREACH(const alignment_segment& seg1, segs1) {

    for (it = segs2.begin(); it != segs2.end(); ) {
      alignment_segment& seg2 = *it;

      // Possibilities:
      //
      // l1---------r1           -->     l1---------r1          |    if (l1 <= l2 && l2 <= r1 && r1 <= r2)
      //       l2---------r2                          l2--r2    |      l2 = r1 + 1
      //
      //
      //
      //       l1---------r1                   l1---------r1    |    if (l2 <= l1 && l1 <= r2 && r2 <= r1)
      // l2---------r2           -->     l2--r2                 |      r2 = l1 - 1
      //
      //
      //
      // l1--------->r1                  l1----------r1         |    if (l1 <= l2 && r2 <= r1)
      //     l2--r2              -->        null                |      delete segment 2
      //
      //
      //
      //       l1--r1            -->           l1--r1           |    if (l2 <= l1 && r1 <= r2)
      // l2--------------r2              l2--r2      l3--r3     |      [l3,r3] = copy([l2,r2])
      //                                                        |      r2 = l1 - 1
      //                                                        |      l3 = r1 + 1

      if (H::l(seg1) <= H::l(seg2) && H::l(seg2) <= H::r(seg1) && H::r(seg1) < H::r(seg2)) {
        set_l_after_r<H>(seg2, seg1);
        ++it;
      }

      else if (H::l(seg2) < H::l(seg1) && H::l(seg1) <= H::r(seg2) && H::r(seg2) <= H::r(seg1)) {
        set_r_before_l<H>(seg2, seg1);
        ++it;
      }

      else if (H::l(seg1) <= H::l(seg2) && H::r(seg2) <= H::r(seg1)) {
        std::cout << "erasing seg: " << seg2 << std::endl;
        it = segs2.erase(it);
      }

      else if (H::l(seg2) < H::l(seg1) && H::r(seg1) < H::r(seg2)) {
        alignment_segment seg3(seg2);
        set_r_before_l<H>(seg2, seg1);
        set_l_after_r <H>(seg3, seg1);
        new_segs.push_back(seg3);
        ++it;
      }

      else // no overlap
        ++it;
    
    }

    BOOST_FOREACH(const alignment_segment& seg3, new_segs)
      segs2.push_back(seg3);
    new_segs.clear();
  }
}

template<typename H>
std::vector<alignment_segment> subtract(      std::vector<alignment_segment>  segs2, // copy
                                        const std::vector<alignment_segment>& segs1)
{
  subtract_in_place<H>(segs2, segs1);
  return segs2;
}
