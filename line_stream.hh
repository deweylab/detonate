#pragma once
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <boost/shared_ptr.hpp>

class line_stream
{
public:
  line_stream(boost::shared_ptr<std::istream> istrm) : istrm(istrm) {}

  line_stream& operator>>(std::string& line)
  {
    getline(*istrm, line);
    if (line.size() == 0)
      istrm->peek(); // set EOF bit in istrm if next read will be EOF
    return *this;
  }

  inline operator bool() { return *istrm; }
  inline bool operator!() { return !*istrm; }

private:
  boost::shared_ptr<std::istream> istrm;
};
