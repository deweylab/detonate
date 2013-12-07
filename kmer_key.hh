#include "city.h"

struct kmer_key
{
  const char *begin;
  const char *end;
  kmer_key() {}
  kmer_key(
      const char *begin,
      const char *end)
  : begin(begin),
    end(end)
  {}
};

struct kmer_key_hash
{
  size_t operator()(const kmer_key& k) const
  {
    size_t len = k.end - k.begin;
    return CityHash64(k.begin, len);
  }
};

struct kmer_key_equal_to
{
  bool operator()(const kmer_key& lhs, const kmer_key& rhs) const
  {
    const char *i = lhs.begin;
    const char *j = rhs.begin;
    for (;;) {
      // If i == end and j == end, return true.
      // If i != end and j != end and *i == *j, continue;
      // In all other cases, return false.
      if (i == lhs.end) {
        if (j == rhs.end)
          return true;    // i == end and j == end              -> true
        else
          return false;   // i == end and j != end              -> false
      } else {
        if (j == rhs.end)
          return false;   // i != end and j == end              -> false
        else {
          if (*i != *j)
            return false; // i != end and j != end and *i != *j -> false
          else {
          }               // i != end and j != end and *i == *j -> continue
        }
      }
      ++i, ++j;
    }
  }
};
