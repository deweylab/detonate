/* Copyright (c) 2007
   Colin Dewey (University of Wisconsin-Madison)
   cdewey@biostat.wisc.edu
   
   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "deweylab/util/string.hh"

#include <cstdio>
#include <vector>
#include <iostream>

namespace deweylab {
namespace util {
namespace string {

void reverse(std::string& s) {
	std::reverse(s.begin(), s.end());
}
		
std::string reverse_copy(const std::string& s) {
	std::string reversed(s);
	reverse(reversed);
	return reversed;
}
		
const std::string& StringStorer::store(const std::string& s) {
	StringSet::const_iterator pos = strings.find(s);
	if (pos == strings.end()) {
		pos = strings.insert(s).first;
	}
	return *pos;
}
	
std::string capitalize(const std::string& s) {
	std::string cs(s);
	cs[0] = std::toupper(cs[0]);
	return cs;
}

std::string center(const std::string& s, unsigned int width) {
	if (s.length() >= width) {
		return s;
	}
	std::string cs(width, ' ');
	cs.replace((width - s.length()) / 2, s.length(), s);
	return cs;
}

std::string ljust(const std::string& s,
				  const std::string::size_type width) {
	if (s.length() >= width) {
		return s;
	}
	std::string ljs(width, ' ');
	ljs.replace(0, s.length(), s);
	return ljs;
}
	
std::string rjust(const std::string& s,
				  const std::string::size_type width) {
	if (s.length() >= width) {
		return s;
	}
	std::string rjs(width, ' ');
	rjs.replace(width - s.length(), s.length(), s);
	return rjs;
}		
	
size_t count(const std::string& s,
			 const std::string& sub,
			 size_t start,
			 const size_t end) {
	unsigned int c = 0;
	while (true) {
		start = s.find(sub, start);
		if (start >= end) {
			break;
		}
		++c;
		++start;
	}
	return c;
}

struct UpperCaser {
	char operator()(const char c) const {
		return std::toupper(c);
	}
};

struct LowerCaser {
	char operator()(const char c) const {
		return std::tolower(c);
	}
};
	
void upperString(std::string& str) {
	std::transform(str.begin(), str.end(), str.begin(), UpperCaser());
}

void lowerString(std::string& str) {
	std::transform(str.begin(), str.end(), str.begin(), LowerCaser());
}

std::string upper(const std::string& s) {
	std::string us(s);
	upperString(us);
	return us;
}
	
std::string lower(const std::string& s) {
	std::string ls(s);
	lowerString(ls);
	return ls;
}
	
std::string stripLeft(const std::string& s,
					  const std::string& chars) {
    if (s.size() == 0)
      return "";
    std::string::size_type pos = s.find_first_not_of(chars);
    return (pos == std::string::npos ? "" : s.substr(pos));
}

std::string stripRight(const std::string& s,
					   const std::string& chars) {
    if (s.size() == 0)
      return "";
    std::string::size_type pos = s.find_last_not_of(chars);
    return (pos == std::string::npos ? "" : s.substr(0, pos + 1));
}
	
std::string strip(const std::string& s,
				  const std::string& chars) {
	return stripRight(stripLeft(s, chars), chars);
}
	
std::string validPosix(const std::string& s, char replacement) {
	std::string valid = s;
	std::string::size_type pos = 0;
	while ((pos = s.find_first_not_of(validPosixChars, pos))
		   != std::string::npos) {
		valid[pos] = replacement;
	}
	return valid;
}

std::string firstWord(const std::string& s,
					  std::string::size_type start) {
	start = s.find_first_not_of(whitespaceChars, start);
	if (start == std::string::npos) {
		return "";
	}
	std::string::size_type end = s.find_first_of(whitespaceChars, start);
	return (end == std::string::npos ?
			s.substr(start) :
			s.substr(start, end));
}

// Split the string S into its first word and the remaining string
std::pair<std::string, std::string>
split_first_word(const std::string& s) {
    std::pair<std::string, std::string> result;
    size_t first_end = s.find_first_of(whitespaceChars);
    if (first_end == std::string::npos) {
        result.first = s;
    } else {
        result.first = s.substr(0, first_end);
        size_t rest_begin = s.find_first_not_of(whitespaceChars, first_end);
        if (rest_begin != std::string::npos) {
            result.second = s.substr(rest_begin);
        }
    }
    return result;
}

bool startsWith(const std::string& str,
				const std::string& prefix,
				std::string::size_type start,
				std::string::size_type end) {
	std::string::const_iterator
		strBegin = str.begin() + start,
		strEnd = (end == std::string::npos ? str.end() : str.begin() + end),
		preBegin = prefix.begin(),
		preEnd = prefix.end();
	while (strBegin != strEnd &&
		   preBegin != preEnd &&
		   *strBegin == *preBegin) {
		++strBegin;
		++preBegin;
	}
	return (preBegin == preEnd);
}

bool endsWith(const std::string& str,
			  const std::string& suffix,
			  std::string::size_type start,
			  std::string::size_type end) {
	std::string::const_reverse_iterator
		strBegin(str.begin() + start),
		strEnd(end == std::string::npos ? str.end() : str.begin() + end),
		preBegin(suffix.begin()),
		preEnd(suffix.end());
	while (strBegin != strEnd &&
		   preBegin != preEnd &&
		   *strEnd == *preEnd) {
		++strEnd;
		++preEnd;
	}
	return (preBegin == preEnd);
}
		
struct StringJoiner {
	std::string sep;
	std::string& joinedStr;
	StringJoiner(std::string& joinedStr,
				 const std::string& sep)
		: sep(sep), joinedStr(joinedStr) {}
	void operator()(const std::string& s) {
		if (!joinedStr.empty()) {
			joinedStr.append(sep);
		}
		joinedStr.append(s);
	}
};
	
std::string fill(const std::string& s,
				 std::string::size_type width,
				 const std::string& initialIndent,
				 const std::string& subsequentIndent) {
	std::string joinedStr;
	StringJoiner joiner(joinedStr, "\n");
	wrap(s, boost::make_function_output_iterator(joiner),
		 width, initialIndent, subsequentIndent);
	return joinedStr;
}

// Return the percent encoding of a single character
std::string percent_encoded(char c) {
	char encoded[4]; // '%', two hex digits, and NULL character
	std::sprintf(encoded, "%%%02X", c);
	return std::string(encoded, 3);
}

// Returns a string with all characters percent-encoded
std::string percent_encoded(const std::string& s) {
	std::vector<char> encoded(s.length() * 3 + 1, 0);
	char* pos = &encoded[0];
	for (size_t i = 0; i < s.length(); ++i) {
		std::sprintf(pos, "%%%02X", s[i]);
		pos += 3;
	}
	return std::string(encoded.begin(), encoded.end() - 1);
}

// Return a string with all percent decoded characters in S decoded
std::string percent_decoded(const std::string& s) {
	std::string decoded;
	std::string::size_type pos = 0, next;
	while ((next = s.find_first_of("%", pos)) != std::string::npos) {
		if (s.length() - next < 3) {
			throw std::runtime_error("Invalid percent encoded string: " + s);
		}
		decoded.append(s, pos, next - pos);
		char hex[3] = {s[next + 1], s[next + 2], '\0'};
		char *endptr = NULL;
		char c = static_cast<char>(std::strtol(hex, &endptr, 16));
		decoded.push_back(c);
		pos = next + 3;
	}
	decoded.append(s, pos, s.length());
	
	return decoded;
}


} } }
