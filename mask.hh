#pragma once
#include <vector>
#include <boost/foreach.hpp>

class mask
{
public:
  mask(size_t length) : data(length, false), num_1s(0) {}

  void add_interval(size_t start, size_t end)
  {
    for (size_t i = start; i <= end; ++i) {
      if (!data[i]) // excluded -> included => increment num ones
        ++num_1s;
      data[i] = true;
    }
  }

  void remove(const std::vector<size_t>& mismatches)
  {
    BOOST_FOREACH(size_t i, mismatches) {
      if (data[i]) // included -> excluded => decrement num ones
        --num_1s;
      data[i] = false;
    }
  }

  size_t num_ones() { return num_1s; }

private:
  std::vector<bool> data;
  size_t num_1s;
};
