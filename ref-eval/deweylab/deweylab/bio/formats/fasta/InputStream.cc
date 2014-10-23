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

#include <istream>

#include "deweylab/bio/formats/fasta/InputStream.hh"
#include "deweylab/util/string.hh"

namespace deweylab {
namespace bio {
namespace formats {
namespace fasta {

InputStream::operator bool() const {
    return not atEnd;
}

bool InputStream::operator!() const {
    return atEnd;
}

void InputStream::setKeepWhitespace(bool b) {
    keepWhitespace = b;
}

InputStream::InputStream(std::istream& stream,
                         const size_t bufferSize)
    : stream(stream),
      buffer(bufferSize),
      pos(buffer.end()),
      atEnd(false),
      keepWhitespace(false)
{
}

void InputStream::fillBuffer() {
    // Fill buffer
    stream.read(&buffer[0], buffer.size());

    // If buffer is not completely filled, resize to fit
    if (static_cast<size_t>(stream.gcount()) < buffer.size()) {
        buffer.resize(stream.gcount());
    }
}

InputStream& InputStream::operator>>(Record& rec) {
    using util::string::IsWhitespace;
    using util::string::split_first_word;

    // Fill buffers until first '>' character is found or EOF
    while (pos == buffer.end() and stream) {
        fillBuffer();
        pos = std::find(buffer.begin(), buffer.end(), TITLE_LINE_PREFIX);
    }

    // If no > was found, mark stream as processed
    if (pos == buffer.end()) {
        atEnd = true;
        return *this;
    }

    std::string title;

    // Skip over '>' character
    std::vector<char>::iterator start = pos + 1;

    // Find newline char
    pos = std::find(start, buffer.end(), '\n');

    // Fill buffers until end of title is found
    while (stream && pos == buffer.end()) {
        title.append(start, pos);
        fillBuffer();
        start = buffer.begin();
        pos = std::find(start, buffer.end(), '\n');
    }
    title.append(start, pos);

    // Split title line into id and description
    std::pair<std::string, std::string> title_pair = split_first_word(title);
    rec.id = title_pair.first;
    rec.description = util::string::stripRight(title_pair.second);

    rec.sequence.clear();

    // Skip over newline
    start = pos + 1;

    while (true) {
        // Find the next title line, or end of buffer
        pos = std::find(start, buffer.end(), TITLE_LINE_PREFIX);

        // Append to sequence
        if (keepWhitespace) {
            rec.sequence.append(start, pos);
        } else {
            rec.sequence.append(start,
                                std::remove_if(start, pos, IsWhitespace()));
        }

        // Stop if we did not run into the end of the buffer or the
        // stream is exhausted
        if (pos != buffer.end() or not stream) { break; }
        fillBuffer();
        start = buffer.begin();
    }

    return *this;
}

} } } }
