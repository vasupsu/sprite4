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

///
/// \author Chris Saunders
///

#pragma once

#include "BaseErrorCounts.hh"
#include "IndelErrorCounts.hh"


/// Stores sequencing error patterns observed from data
///
/// Used by downstream operations to estimate useful error
/// parameters for modeling, but no parameters are directly
/// stored in this object, this object is only a compressed
/// abstraction of the input data
///
struct SequenceErrorCounts
{
    BaseErrorCounts&
    getBaseCounts()
    {
        return _bases;
    }

    const BaseErrorCounts&
    getBaseCounts() const
    {
        return _bases;
    }

    IndelErrorCounts&
    getIndelCounts()
    {
        return _indels;
    }

    const IndelErrorCounts&
    getIndelCounts() const
    {
        return _indels;
    }

    const std::string&
    getSampleName() const
    {

        return _sampleName;
    }

    void
    merge(const SequenceErrorCounts& in);

    void
    clear()
    {
        _bases.clear();
        _indels.clear();
    }

    void
    setSampleName(const std::string& sampleName)
    {
        _sampleName = sampleName;
    }
    void
    save(const char* filename) const;

    void
    load(const char* filename);

    /// debug output
    void
    dump(std::ostream& os) const;

private:
    std::string _sampleName;
    BaseErrorCounts _bases;
    IndelErrorCounts _indels;
};
