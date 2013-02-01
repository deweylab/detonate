#pragma once
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

class line_stream
{
public:
  line_stream(std::istream& istrm) : istrm(istrm) {}

  line_stream& operator>>(std::string& line)
  {
    getline(istrm, line);
    if (line.size() == 0)
      istrm.peek(); // set EOF bit in istrm if next read will be EOF
    return *this;
  }

  inline operator bool() { return istrm; }
  inline bool operator!() { return !istrm; }

private:
  std::istream& istrm;
};
