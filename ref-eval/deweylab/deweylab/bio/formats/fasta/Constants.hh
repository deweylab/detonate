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

#ifndef __DEWEYLAB_BIO_FORMATS_FASTA_CONSTANTS_HH__
#define __DEWEYLAB_BIO_FORMATS_FASTA_CONSTANTS_HH__

#include <cstddef>

namespace deweylab {
namespace bio {
namespace formats {
namespace fasta {

struct Constants {
	static const size_t DEFAULT_LINE_WIDTH;
	static const char TITLE_LINE_PREFIX;
};

} } }

};

#include "deweylab/bio/formats/fasta/Constants.cc"

#endif // __DEWEYLAB_BIO_FORMATS_FASTA_CONSTANTS_HH__
