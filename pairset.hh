#pragma once
#include <set>
#include <vector>
#include <queue>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

// Pairset concept:
//
// General notes:
// - All pairs are of the form (i,j) where i <= j.
// - In other words, the pairs (i,j) and (j,i) are considered equivalent and
//       are not counted twice.
// - "Pairs between (l,l) and (h,h)" means all pairs (i,j) such that i and j
//       are both in the interval [l,h] and i <= j.
//
// Interface:
// - Pairset(size_t top) - constructs a new pairset with possible pairs (0,0),
//       (0,1), and (top-1,top-1), inclusive.
//
// - size_t size() - returns how many pairs are in the pairset
//
// - void add_square_with_exceptions(size_t lo, size_t hi, ConstIterator
//       exceptions_begin, ConstIterator exceptions_end) - adds all pairs
//       between (lo,lo) and (hi,hi), except all pairs (i,x) or (x,i) such that
//       x is an exception and i is between lo and hi. Note that the iterators
//       are over exceptional coordinates, not pairs.

class brute_force_pairset
{
public:
  brute_force_pairset(size_t /*top*/) {}

  size_t size() { return set.size(); }

  template<typename ConstIterator>
  void add_square_with_exceptions(size_t lo, size_t hi, ConstIterator exceptions_begin, ConstIterator exceptions_end)
  {
    if (lo > hi)
      std::swap(lo, hi);

    std::set<std::pair<size_t,size_t> > exceptions;
    for (ConstIterator x = exceptions_begin; x != exceptions_end; ++x) {
      for (size_t i = lo; i <= hi; ++i) {
        if (*x <= i)
          exceptions.insert(std::make_pair(*x, i));
        else
          exceptions.insert(std::make_pair(i, *x));
      }
    }

    for (size_t i = lo; i <= hi; ++i)
      for (size_t j = i; j <= hi; ++j)
        if (exceptions.count(std::make_pair(i,j)) == 0)
          set.insert(std::make_pair(i,j));
  }

private:
  std::set<std::pair<size_t,size_t> > set;
};


class big_matrix_pairset
{
public:
  big_matrix_pairset(size_t top) : top(top), sz(0), data(top*top, false) {}

  size_t size() { return sz; }

  template<typename ConstIterator>
  void add_square_with_exceptions(size_t lo, size_t hi, ConstIterator exceptions_begin, ConstIterator exceptions_end)
  {
    if (lo > hi)
      std::swap(lo, hi);

    std::set<size_t> exceptions;
    for (ConstIterator x = exceptions_begin; x != exceptions_end; ++x)
      exceptions.insert(*x);

    for (size_t i = lo; i <= hi; ++i)
      for (size_t j = i; j <= hi; ++j)
        if (exceptions.count(i) == 0 && exceptions.count(j) == 0)
          assign(i, j, true);
  }

private:
  size_t top, sz;
  std::vector<bool> data;

  void assign(size_t x, size_t y, bool val)
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

namespace detail
{
  inline size_t choose_2(size_t n) { return n*(n-1)/2; }

  struct coord
  {
    size_t x, y;
    coord(size_t x, size_t y) : x(x), y(y) {}
  };

  struct rect
  {
    coord lo, hi;
    rect(coord lo, coord hi) : lo(lo), hi(hi) {}
  };
}

namespace std
{
  template<>
  struct less<detail::coord>
  {
    bool operator()(const detail::coord& a,
                    const detail::coord& b) const
    {
      return (a.x < b.x) || (a.x == b.x && a.y < b.y);
    }
  };

  template<>
  struct less<detail::rect>
  {
    bool operator()(const detail::rect& a,
                    const detail::rect& b) const
    {
      return (a.lo.x <  b.lo.x)
          || (a.lo.x == b.lo.x && a.lo.y <  b.lo.y)
          || (a.lo.x == b.lo.x && a.lo.y == b.lo.y && a.hi.x <  b.hi.x)
          || (a.lo.x == b.lo.x && a.lo.y == b.lo.y && a.hi.x == b.hi.x && a.hi.y < b.hi.y);
    }
  };
}


class smart_pairset
{
public:
  smart_pairset(size_t /*top*/) {}

  size_t size() const
  {
    // Make a copy of the squares so that we can modify it.
    std::list<detail::rect> rects(squares.begin(), squares.end());

    // Remove any squares that are completely contained inside another one.
    std::list<detail::rect>::iterator A, B, end = rects.end();
    for (A = rects.begin(); A != end; ) {
      bool did_erase = false;
      for (B = rects.begin(); B != end; ++B) {
        if (A != B && B->lo.x <= A->lo.x && A->hi.x <= B->hi.x) { // B contains A
          assert(B->lo.y <= A->lo.y && A->hi.y <= B->hi.y);
          rects.erase(A++);
          did_erase = true;
          break;
        }
      }
      if (!did_erase)
        ++A;
    }

    // We now know that if square B starts before square A starts, then B ends
    // before A ends too. (This holds because all squares are on the diagonal
    // and none is contained inside another.)

    // Go from left to right, and if a square overlaps a predecessor, chop it
    // off on the left, so that they no longer overlap. (Remember that pair
    // (i,j) is equivalent to (j,i), so we don't have to worry about what
    // happens "above the diagonal".) We can go from left to right by just
    // iterating over rects, because the sorting criterion says that "B < A"
    // if, first of all, B starts to the left of where A starts.
    for (A = rects.begin(); A != end; ++A) {
      for (B = rects.begin(); B != A; ++B) { // iterate over all B such that B starts before A
        if (A->lo.x <= B->hi.x) {            // - if they overlap, i.e., B ends after (or =) A starts
          A->lo.x = B->hi.x + 1;             //   - truncate A so they don't overlap
          assert(A->lo.x <= A->hi.x);
        }
      }
    }

    // Count how many pairs each rectangle contributes. On the "diagonal
    // square" of the rectangle, it contributes (n+1) choose 2 pairs, where n
    // is the side length of the diagonal square. On the rest of the rectangle,
    // it contributes n*m pairs, where n and m are the side lengths of the rest
    // of the rectangle.
    size_t sz = 0;
    for (A = rects.begin(); A != end; ++A) {
      size_t diagonal_side_length = A->hi.x - A->lo.x + 1;
      size_t remainder_side_length_1 = A->lo.x - A->lo.y;
      size_t remainder_side_length_2 = diagonal_side_length;
      sz += detail::choose_2(diagonal_side_length+1) + remainder_side_length_1*remainder_side_length_2;
    }

    // Subtract the exceptions
    sz -= excepted_pairs.size();

    return sz;
  }

  template<typename ConstIterator>
  void add_square_with_exceptions(size_t lo, size_t hi, ConstIterator exceptions_begin, ConstIterator exceptions_end)
  {
    if (lo > hi)
      std::swap(lo, hi);

    // First make a set of the new exceptions for easy lookup.
    std::set<size_t> exceptions;
    for (ConstIterator x = exceptions_begin; x != exceptions_end; ++x)
      exceptions.insert(*x);

    // Remove any exceptions (x,y) from the main list (i.e., existing, "old"
    // exceptions) that are covered by the new square and are not new
    // exceptions.
    for (std::set<detail::coord>::const_iterator it = excepted_pairs.begin(), end = excepted_pairs.end(); it != end;) {
      if (lo <= it->x && it->x <= hi &&     // (x,y) is in the new square,
          lo <= it->y && it->y <= hi &&     // and
          exceptions.count(it->x) == 0 &&   // (x,y) is not a new exception
          exceptions.count(it->y) == 0) {
        excepted_pairs.erase(it++);
      } else {
        ++it;
      }
    }

    // Add any new exceptions that are not covered by an old square to the main
    // list. Skip old exceptions, since no change would be made by adding them
    // again.
    for (ConstIterator i = exceptions_begin; i != exceptions_end; ++i) {
      for (size_t j = lo; j <= hi; ++j) {
        size_t x, y;
        if (*i <= j) {
          x = *i;
          y =  j;
        } else {
          x =  j;
          y = *i;
        }
        if (excepted_pairs.count(detail::coord(x, y))) // already an exception
          continue;
        else {
          // check if the exception is covered by some square
          bool is_covered = false;
          BOOST_FOREACH(const detail::rect& sq, squares) {
            if (sq.lo.x <= x && x <= sq.hi.x &&
                sq.lo.y <= y && y <= sq.hi.y) {
              is_covered = true;
              break;
            }
          }
          // if not, then the exception is truly an exception, so add it
          if (!is_covered)
            excepted_pairs.insert(detail::coord(x, y));
        }
      }
    }

    // Add the square itself. The square may overlap with existing squares.
    squares.insert(detail::rect(detail::coord(lo, lo), detail::coord(hi, hi)));
  }

private:
  std::set<detail::rect> squares;
  std::set<detail::coord> excepted_pairs;
};
