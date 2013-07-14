#pragma once
#include <set>
#include <vector>
#include <numeric>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

class kmerset
{
public:
  kmerset(size_t k, size_t top) : k(k), top(top), sz(0), data(top+1-k, false) {}

  size_t size() { return sz; }

  template<typename ConstIterator>
  void add_kmers_with_exceptions(size_t lo, size_t hi, ConstIterator exceptions_begin, ConstIterator exceptions_end)
  {
    if (lo > hi)
      std::swap(lo, hi);

    // If this segment is shorter than the kmer length k, then we should skip
    // adding it.
    // 
    // E.g., k=3, lo=4, hi=6: 456  -> ok,     6-4+1 = 3 >= 3
    //       k=3, lo=4, hi=5: 45   -> not ok, 5-4+1 = 2 < 3
    //       k=4, lo=4, hi=6: 456  -> not ok, 6-4+1 = 3 < 4
    //       k=3, lo=4, hi=7: 4567 -> ok,     7-4+1 = 4 >= 3
    // So the requirement is that hi-lo+1 >= k, i.e., hi+1-lo >= k.
    if (hi+1-lo < k)
      return;

    // full tran:   0123456789
    // lo=3, hi=7      .   .
    // exceptions:     01234
    // length = 7+1-3 = 5 = hi+1-lo
    // *x=2 => ignore
    // *x=5 => 5-3=2:    .
    // *x=3 => 3-3=0:  .
    // *x=7 => 7-3=4:      .
    // Thus, add *x-lo to exceptions.
    std::vector<bool> is_exception(hi+1-lo, false);
    for (ConstIterator x = exceptions_begin; x != exceptions_end; ++x)
      if (*x >= lo && *x <= hi)
        is_exception[*x-lo] = true;

    // E.g., k=3, lo=10, hi=15, then we can add (10,12), (11,13), (12,14), (13,15).
    // Note 13 = 15-3+1, i.e., hi-k+1, or (for size_t) hi+1-k.
    // Note 12 = 10+3-1, i.e, i+k-1.
    for (size_t i = lo; i <= hi+1-k; ++i) {
      if (i >= top+1-k)
        throw std::runtime_error("Out of bounds.");
      size_t num_exceptions = std::accumulate(is_exception.begin()+(i-lo),
                                              is_exception.begin()+(i+k-1-lo),
                                              static_cast<size_t>(0));
      if (num_exceptions <= 5) {
        if (data[i] == false) {
          ++sz;
          data[i] = true;
        }
      }
    }
  }

private:
  size_t k, top, sz;
  std::vector<bool> data;
};
