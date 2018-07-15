//
// Strelka - Small Variant Caller
// Copyright (c) 2009-2017 Illumina, Inc.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//

/// \file
/// \author Chris Saunders
///

#pragma once

#include "blt_util/seq_util.hh"

#include <iosfwd>
#include <string>
#include <vector>


struct bed_record
{
    bed_record()
    {
        clear();
    }

    // set record from record string s, return false on error
    bool set(const char* s);

    void clear()
    {
        chrom.clear();
        begin=0;
        end=0;
        line=nullptr;
    }

    bool
    is_valid() const
    {
        return (begin <= end);
    }

    std::string chrom;
    int begin = 0;
    int end = 0;
    const char* line = nullptr;
};


std::ostream& operator<<(std::ostream& os, const bed_record& bedr);

