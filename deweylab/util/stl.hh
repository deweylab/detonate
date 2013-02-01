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

#ifndef __DEWEYLAB_UTIL_STL_HH__
#define __DEWEYLAB_UTIL_STL_HH__

#include <functional>
#include <numeric>
#include <iosfwd>
#include <string>
#include <algorithm>

namespace deweylab {
namespace util {
namespace stl {

template<typename InputIterator1, typename InputIterator2>
size_t matches(InputIterator1 first1, InputIterator1 last1,
			   InputIterator2 first2) {
	typedef typename std::iterator_traits<InputIterator1>::value_type
		value_type;
	return std::inner_product(first1, last1, first2, 0,
							  std::plus<size_t>(),
							  std::equal_to<value_type>());
}

template<typename InputIterator1, typename InputIterator2>
size_t mismatches(InputIterator1 first1, InputIterator1 last1,
				  InputIterator2 first2) {
	typedef typename std::iterator_traits<InputIterator1>::value_type
		value_type;
	return std::inner_product(first1, last1, first2, 0,
							  std::plus<size_t>(),
							  std::not_equal_to<value_type>());
}

template<typename Iterator>
void print_elements(std::ostream& strm,
					Iterator begin,
					Iterator end,
					const std::string& sep = " ") {
	Iterator pos = begin;
	while (pos != end) {
		if (pos != begin) {
			strm << sep;
		}
		strm << *pos;
		++pos;
	}
	strm << std::string("\n");
}

template<typename Container>
void print_elements(std::ostream& strm,
					const Container& c,
					const std::string& sep = " ") {
	print_elements(strm, c.begin(), c.end(), sep);
}
		
template<typename InputIterator>
typename InputIterator::value_type
sum(InputIterator begin, InputIterator end) {
	return std::accumulate(begin, end, 0);
}

template<typename InputIterator>
typename InputIterator::value_type
product(InputIterator begin, InputIterator end) {
	typedef typename InputIterator::value_type value_type;
	return std::accumulate(begin, end, 1, std::multiplies<value_type>());
}

template<typename InputIterator, typename UnaryFunc>
typename UnaryFunc::result_type
sum(InputIterator begin,
	InputIterator end,
	UnaryFunc f,
	typename UnaryFunc::result_type initValue = typename UnaryFunc::result_type())
{
	while (begin != end) {
		initValue += f(*begin);
		++begin;
	}
	return initValue;
}

template<typename InputIterator, typename UnaryFunc>
typename UnaryFunc::result_type
product(InputIterator begin,
		InputIterator end,
		UnaryFunc f,
		typename UnaryFunc::result_type initValue = typename UnaryFunc::result_type())
{
	while (begin != end) {
		initValue *= f(*begin);
		++begin;
	}
	return initValue;
}
		
template<class OP>
struct FunctorRef : public std::unary_function<typename OP::argument_type,
											   typename OP::result_type> {
	const OP& op;
	FunctorRef(const OP& op) : op(op) {}
	typename OP::result_type
	operator()(typename OP::argument_type x) {
		return op(x);
	}
};

template<class OP>
FunctorRef<OP> make_functor_ref(const OP& op) { return FunctorRef<OP>(op); }

template<typename Container>
void remove_duplicates(Container& c) {
	std::sort(c.begin(), c.end());
	c.erase(std::unique(c.begin(), c.end()), c.end());
}

// Implementation of copy_if from Effective STL
template<typename InputIterator,
		 typename OutputIterator,
		 typename Predicate>
OutputIterator copy_if(InputIterator begin,
					   InputIterator end,
					   OutputIterator destBegin,
					   Predicate p)
{
	while (begin != end) {
		if (p(*begin)) *destBegin++ = *begin;
		++begin;
	}
	return destBegin;
}

// Generic input iterator over an input stream
template<typename StreamType, typename ValueType>
class InputStreamIterator
	: public std::iterator<std::input_iterator_tag, ValueType>
{
private:
	StreamType* stream;
	ValueType value;
	bool ok;
		
public:      
	InputStreamIterator() : stream(NULL), ok(false) {
	}

	InputStreamIterator(StreamType& s) : stream(&s) {
		read();
	}

	InputStreamIterator(const InputStreamIterator& obj) 
		: stream(obj.stream), value(obj.value), ok(obj.ok) {
	}
		
	ValueType& operator*() {
		return value;
	}
		
	ValueType* operator->() {
		return &(operator*());
	}

	InputStreamIterator& operator++() {
		read(); 
		return *this; 
	}

	InputStreamIterator operator++(int)  {
		InputStreamIterator tmp = *this;
		read();
		return tmp;
	}

	bool operator==(const InputStreamIterator& other) const {
		return (ok == other.ok) and (not ok or stream == other.stream);
	}

	bool operator!=(const InputStreamIterator& other) const {
		return not (*this == other);
	}
	
private:      

	void read() {
		ok = (stream and *stream) ? true : false;
		if (ok) {
			*stream >> value;
			ok = *stream ? true : false;
		}
	}
};


// Generic output iterator to an output stream
template<typename StreamType, typename ValueType>
class OutputStreamIterator
	: public std::iterator<std::output_iterator_tag, ValueType>
{
private:
	StreamType* stream;
		
public:      
	OutputStreamIterator(StreamType& stream) : stream(&stream) {
	}

	OutputStreamIterator(const OutputStreamIterator& other) 
		: stream(other.stream) {
	}
		
	OutputStreamIterator& operator*() {
		return *this;
	}
		
	OutputStreamIterator& operator++() {
		return *this; 
	}

	OutputStreamIterator operator++(int)  {
		return *this;
	}
	
	OutputStreamIterator operator=(const ValueType& value) {
		return *stream << value;
	}
};

} } }

#endif // __DEWEYLAB_UTIL_STL_HH__
