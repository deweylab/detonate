// Copyright (c) 2013
// Nathanael Fillmore (University of Wisconsin-Madison)
// nathanae@cs.wisc.edu
//
// This file is part of REF-EVAL.
//
// REF-EVAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// REF-EVAL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with REF-EVAL.  If not, see <http://www.gnu.org/licenses/>.

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

  inline operator bool() { return static_cast<bool>(*istrm); }
  inline bool operator!() { return !bool(); }

private:
  boost::shared_ptr<std::istream> istrm;
};
