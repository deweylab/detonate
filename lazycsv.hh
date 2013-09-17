#pragma once
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>

template<size_t num_fields, char sep = ','>
class lazycsv {
public:
  lazycsv() {}
  lazycsv(const std::string& line) : line(line) { update_starts(); }

  // Get the value of the t'th field, converted to type T.
  template<typename T>
  inline T at(size_t t) const
  {
    size_t len = (starts[t+1] - 1) - starts[t]; // "- 1" so as not to include sep char
    return boost::lexical_cast<T>(line.substr(starts[t], len));
  }

  void parse_line(const std::string& line_)
  {
    line = line_;
    update_starts();
  }

  bool operator==(const lazycsv& other) const { return line == other.line; }
  bool operator!=(const lazycsv& other) const { return !(*this == other); }

private:
  std::string line;
  boost::array<size_t, num_fields + 1> starts;

  // Figure out where fields start (and end).
  //
  // Note: we record field start positions, not iterators pointing to those
  // positions, because it makes copying this object easier and more efficient.
  void update_starts()
  {
    #if 0
    // much slower for some reason
    starts[0] = 0;
    size_t f = 1;
    for (size_t i = 0; i < line.size(); ++i) {
      if (line[i] == sep) {
        starts[f] = i + 1;
        ++f;
      }
    }
    #else
    starts[0] = 0;
    size_t f = 1;
    std::string::size_type i = 0;
    for (i = line.find(sep, 0);
         i != std::string::npos;
         i = line.find(sep, i + 1)) {
      starts[f] = i+1;
      ++f;
    }
    #endif

    starts[num_fields] = line.size() + 1; // "+ 1" to imitate sep char
    if (f != num_fields)
      throw std::runtime_error("Invalid number of fields in line: '" + line + "'");
  }
};
