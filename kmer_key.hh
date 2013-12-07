struct kmer_key
{
  string::const_iterator begin;
  string::const_iterator end;
  kmer_key() {}
  kmer_key(
      const string::const_iterator& begin,
      const string::const_iterator& end)
  : begin(begin),
    end(end)
  {}
};

struct kmer_key_hash
{
  hash<std::string> h;

  size_t operator()(const kmer_key& k) const
  {
    return h(std::string(k.begin, k.end));
  }

  kmer_key_hash()
  : h()
  {}
};

struct kmer_key_equal_to
{
  bool operator()(const kmer_key& lhs, const kmer_key& rhs) const
  {
    string::const_iterator i = lhs.begin;
    string::const_iterator j = rhs.begin;
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

