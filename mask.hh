#pragma once
#include <vector>
#include <set>
#include <boost/foreach.hpp>

class mask
{
public:
  mask(size_t length) : data(length, false), num_1s(0) {}

  // Adds the elements [start, end] - exceptions to the mask
  template<typename ConstIterator>
  void add_interval_with_exceptions(size_t start, size_t end, ConstIterator exceptions_begin, ConstIterator exceptions_end)
  {
    if (start > end)
      std::swap(start, end);
    std::set<size_t> exceptions(exceptions_begin, exceptions_end);
    for (size_t i = start; i <= end; ++i) {
      if (exceptions.count(i) == 0) {
        if (!data[i]) // excluded -> included => increment num ones
          ++num_1s;
        data[i] = true;
      }
    }
  }

  size_t num_ones() { return num_1s; }

private:
  std::vector<bool> data;
  size_t num_1s;
};
