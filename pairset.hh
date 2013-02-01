#pragma once
#include <set>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

/*
template<typename T>
struct pairset
{
  // Returns the number of pairs in the set.
  size_t size();

  // Adds the pairs { (x,y) : x,y in [lo, hi] - exceptions } to the set.
  void add_square_with_exceptions(T lo, T hi, const std::vector<T>& exceptions);

  // OR:

  // Adds the pairs { (x,y) : x,y in [lo, hi] } to the set.
  void add_square(T lo, T hi);

  // Removes the pair (x,y) from the set.
  void remove_pair(T x, T y);
};
*/

struct brute_force_pairset
{
  typedef size_t T;

  brute_force_pairset(size_t /*top*/) {}

  size_t size() { return set.size(); }

  void add_square(T lo, T hi)
  {
    if (lo > hi)
      std::swap(lo, hi);
    for (T x = lo; x <= hi; ++x)
      for (T y = x; y <= hi; ++y)
        set.insert(std::make_pair(x, y));
  }

  void remove_pair(T x, T y) { set.erase(std::make_pair(x, y)); }

  void remove_all_pairs(const std::vector<T>& xs)
  {
    BOOST_FOREACH(T x, xs)
      BOOST_FOREACH(T y, xs)
        if (x <= y)
          remove_pair(x, y);
  }

private:
  std::set<std::pair<T, T> > set;
};

struct big_matrix_pairset
{
  typedef size_t T;

  big_matrix_pairset(size_t top) : top(top), sz(0), data(top*top, false) {}

  size_t size() {
    /*
    size_t sz = 0;
    BOOST_FOREACH(bool b, data)
      if (b)
        ++sz;
    */
    return sz;
  }

  void add_square(T lo, T hi)
  {
    if (lo > hi)
      std::swap(lo, hi);
    for (T x = lo; x <= hi; ++x)
      for (T y = x; y <= hi; ++y)
        assign(x, y, true);
  }

  void remove_pair(T x, T y) { assign(x, y, false); }

  void remove_all_pairs(const std::vector<T>& xs)
  {
    BOOST_FOREACH(T x, xs)
      BOOST_FOREACH(T y, xs)
        if (x <= y)
          remove_pair(x, y);
  }

private:
  size_t top, sz;
  std::vector<bool> data;
  void assign(T x, T y, bool val)
  {
    if (x >= top)
      throw std::runtime_error("x = " + boost::lexical_cast<std::string>(x) + " is not in range [0, " + boost::lexical_cast<std::string>(top) + ").");
    if (y >= top)
      throw std::runtime_error("y = " + boost::lexical_cast<std::string>(y) + " is not in range [0, " + boost::lexical_cast<std::string>(top) + ").");

    bool old = data[x*top + y];
    if (!old && val)       // excluded -> included => increment size
      ++sz;
    else if (old && !val)  // included -> excluded => decrement size
      --sz;

    data[x*top + y] = val;
  }
};
