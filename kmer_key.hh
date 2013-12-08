#pragma once
#include "city.h"

typedef const char * kmer_key;

struct kmer_key_hash
{
  size_t kmerlen;

  kmer_key_hash(size_t kmerlen)
  : kmerlen(kmerlen)
  {}

  size_t operator()(const char *k) const
  {
    return CityHash64(k, kmerlen);
  }
};

struct kmer_key_equal_to
{
  size_t kmerlen;

  kmer_key_equal_to(size_t kmerlen)
  : kmerlen(kmerlen)
  {}

  bool operator()(const char *lhs, const char *rhs) const
  {
    for (size_t i = 0; i < kmerlen; ++i, ++lhs, ++rhs)
      if (*lhs != *rhs)
        return false;
    return true;
  }
};
