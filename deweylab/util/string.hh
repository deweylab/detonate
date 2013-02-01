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

#ifndef __DEWEYLAB_UTIL_STRING_HH__
#define __DEWEYLAB_UTIL_STRING_HH__

#include <string>
#include <sstream>

#include "boost/function_output_iterator.hpp"
#include "boost/tokenizer.hpp"
#include "boost/lexical_cast.hpp"
#include "boost/unordered_set.hpp"

namespace deweylab {
namespace util {
namespace string {

struct IsWhitespace {
	typedef const char argument_type;
	bool operator()(argument_type c) const {
		switch(c) {
		case ' ':
		case '\t':
		case '\n':
		case '\r':
		case '\f':
		case '\v':
			return true;
		default:
			return false;
		}
	}
};

static const std::string validPosixChars("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-" );
static const std::string whitespaceChars(" \t\n\r\f\v");
static const std::string digitChars("0123456789");

void reverse(std::string& s);
std::string reverse_copy(const std::string& s);

std::string capitalize(const std::string& s);
std::string center(const std::string& s,
				   const std::string::size_type width);
size_t count(const std::string& s,
			 const std::string& sub,
			 std::string::size_type start = 0,
			 const std::string::size_type end = std::string::npos);
std::string upper(const std::string& s);
std::string lower(const std::string& s);

std::string ljust(const std::string& s,
				  const std::string::size_type width);
std::string rjust(const std::string& s,
				  const std::string::size_type width);
	
std::string stripLeft(const std::string& s,
					  const std::string& chars = whitespaceChars);
std::string stripRight(const std::string& s,
					   const std::string& chars = whitespaceChars);
std::string strip(const std::string& s,
				  const std::string& chars = whitespaceChars);
std::string validPosix(const std::string& s,
					   char replacement = '_');
std::string firstWord(const std::string& s,
					  std::string::size_type start = 0);
	
bool startsWith(const std::string& str,
				const std::string& prefix,
				std::string::size_type start = 0,
				std::string::size_type end = std::string::npos);
bool endsWith(const std::string& str,
			  const std::string& suffix,
			  std::string::size_type start = 0,
			  std::string::size_type end = std::string::npos);
	
template<typename InputIterator>
std::string join(InputIterator begin,
				 InputIterator end,
				 const std::string& sep = "");

template<typename T>
std::string toString(const T& x);
	
void upperString(std::string& str);
void lowerString(std::string& str);


// // Percent encoding/decoding functions

std::string percent_encoded(char c);

// Returns a string with all characters percent-encoded
std::string percent_encoded(const std::string& s);

// Return a string with all characters evaluating to True with respect
// to PRED converted to percent encodings
template<typename ReplacePred>
std::string percent_encode(const std::string& s, const ReplacePred& pred);

// Return a string with all percent decoded characters in S decoded
std::string percent_decoded(const std::string& s);


// Split the string S into words separated by whitespace and
// output each word to OUT
template<typename OutputIterator>
void split(const std::string& s,
		   OutputIterator out);
	
// Split the string S by the characters in SEP and output each
// token to OUT
template<typename OutputIterator>
void split(const std::string& s,
		   OutputIterator out,
		   const std::string& sep);

// Split the string S into its first word and the remaining string
std::pair<std::string, std::string>
split_first_word(const std::string& s);

template<typename OutputIterator>
void wrap(const std::string& s,
		  OutputIterator out,
		  std::string::size_type width = 70,
		  const std::string& initialIndent = "",
		  const std::string& subsequentIndent = "");
		
std::string fill(const std::string& s,
				 std::string::size_type width = 70,
				 const std::string& initialIndent = "",
				 const std::string& subsequentIndent = "");
	
template<typename Target>
struct Converter {
	typedef Target result_type;
	Target operator()(const std::string& s) const {
		return boost::lexical_cast<Target, const std::string>(s);
	}
};
	
class StringStorer {
public:
	typedef boost::unordered_set<std::string> StringSet;
	typedef StringSet::const_iterator Iterator;
		
	StringStorer() {}
	const std::string& store(const std::string& s);
	Iterator begin() const { return strings.begin(); }
	Iterator end() const { return strings.end(); }
private:
	StringSet strings;
};

// Implementation

template<typename T>
std::string toString(const T& x) {
	std::ostringstream oss;
	oss << x;
	return oss.str();
}

template<typename InputIterator>
std::string join(InputIterator begin,
				 InputIterator end,
				 const std::string& sep) {
	std::string s;
	if (begin != end) {
		s = toString(*begin);
		++begin;
	}
	while (begin != end) {
		s += sep;
		s += toString(*begin);
		++begin;
	}
	return s;
}

// Split the string S into words separated by whitespace and
// output each word to OUT
template<typename OutputIterator>
void split(const std::string& s, OutputIterator out) {
	typedef boost::tokenizer<boost::char_separator<char> > CharTok;
	CharTok tabtok(s, boost::char_separator<char>(whitespaceChars.c_str(),
												  ""));
	std::copy(tabtok.begin(), tabtok.end(), out);
}
	
// Split the string S by the characters in SEP and output each
// token to OUT
template<typename OutputIterator>
void split(const std::string& s,
		   OutputIterator out,
		   const std::string& sep) {
	typedef boost::tokenizer<boost::char_separator<char> > CharTok;
	CharTok tabtok(s, boost::char_separator<char>(sep.c_str(),
												  "",
												  boost::keep_empty_tokens));
	std::copy(tabtok.begin(), tabtok.end(), out);
}
	
template<typename OutputIterator>
void wrap(const std::string& s,
		  OutputIterator out,
		  std::string::size_type width,
		  const std::string& initialIndent,
		  const std::string& subsequentIndent) {
	typedef boost::tokenizer<boost::char_separator<char> > CharTok;
	CharTok tabtok(s, boost::char_separator<char>(whitespaceChars.c_str(),
												  ""));

	std::string line;
	line.reserve(width);
	line.append(initialIndent);

	bool firstWord = true;
	CharTok::const_iterator word;
	for (word = tabtok.begin(); word != tabtok.end(); ++word) {
		if (line.length() + 1 + word->length() > width) {
			if (firstWord) {
				*out = *word;
				++out;
				line.clear();
				line.append(subsequentIndent);
			} else {
				*out = line;
				++out;
				line.clear();
				line.append(subsequentIndent);
				line.append(*word);
			}
		} else {
			if (firstWord) {
				line.append(*word);
				firstWord = false;
			} else {
				line.append(" ");
				line.append(*word);
			}
		}
	}

	if (!firstWord) {
		*out = line;
		++out;
	}
}

} } }

#include "deweylab/util/string.cc"

#endif // __DEWEYLAB_UTIL_STRING_HH__
