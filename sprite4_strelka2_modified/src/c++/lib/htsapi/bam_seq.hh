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

#include "blt_util/blt_types.hh"
#include "blt_util/PolymorphicObject.hh"
#include "blt_util/reference_contig_segment.hh"
#include "blt_util/seq_util.hh"

#include <cassert>

#include <algorithm>
#include <iosfwd>
#include <string>


namespace BAM_BASE
{
enum index_t
{
    REF = 0x0,
    A = 0x1,
    C = 0x2,
    G = 0x4,
    T = 0x8,
    ANY = 0xF
};
}


inline
char
get_bam_seq_char(const uint8_t a)
{
    using namespace BAM_BASE;

    switch (a)
    {
    case REF:
        return '=';
    case A:
        return 'A';
    case C:
        return 'C';
    case G:
        return 'G';
    case T:
        return 'T';
    default:
        return 'N';
    }
}


inline
char
get_bam_seq_complement_char(const uint8_t a)
{
    using namespace BAM_BASE;

    switch (a)
    {
    case REF:
        return '=';
    case A:
        return 'T';
    case C:
        return 'G';
    case G:
        return 'C';
    case T:
        return 'A';
    default:
        return 'N';
    }
}


inline
uint8_t
get_bam_seq_code(const char c)
{
    using namespace BAM_BASE;

    switch (c)
    {
    case '=':
        return REF;
    case 'A':
        return A;
    case 'C':
        return C;
    case 'G':
        return G;
    case 'T':
        return T;
    default:
        return ANY;
    }
}


inline
uint8_t
bam_seq_code_to_id(const uint8_t a,
                   const uint8_t ref = BAM_BASE::ANY)
{
    using namespace BAM_BASE;

    switch (a)
    {
    case REF:
        return bam_seq_code_to_id(ref);
    case A:
        return 0;
    case C:
        return 1;
    case G:
        return 2;
    case T:
        return 3;
    case ANY:
        return 4;
    default:
        base_error("bam_seq_code_to_id",a);
        return 4;
    }
}


// interface to bam_seq -- allows us to pass either compressed
// sequences from bam files and regular strings using the same
// object:
//
struct bam_seq_base : public PolymorphicObject
{
    virtual uint8_t get_code(pos_t i) const = 0;

    virtual char get_char(const pos_t i) const = 0;

    virtual unsigned size() const = 0;

protected:
    bool
    is_in_range(const pos_t i) const
    {
        return ((i>=0) && (i<static_cast<pos_t>(size())));
    }
};

std::ostream& operator<<(std::ostream& os, const bam_seq_base& bs);


//
//
struct bam_seq : public bam_seq_base
{
    bam_seq(const uint8_t* s,
            const uint16_t init_size,
            const uint16_t offset=0)
        : _s(s), _size(init_size), _offset(offset) {}

#if 0
    bam_seq(const bam_seq bs,
            const uint16_t size,
            const uint16_t offset=0)
        : _s(bs.s), _size(size), _offset(bs.offset+offset)
    {
        assert((offset+size)<=bs.size);
    }
#endif

    uint8_t
    get_code(pos_t i) const override
    {
        if (! is_in_range(i)) return BAM_BASE::ANY;
        i += static_cast<pos_t>(_offset);
        return _s[(i/2)] >> 4*(1-(i%2)) & 0xf;
    }

    char
    get_char(const pos_t i) const override
    {
        return get_bam_seq_char(get_code(i));
    }

    char
    get_complement_char(const pos_t i) const
    {
        return get_bam_seq_complement_char(get_code(i));
    }

    std::string
    get_string() const
    {
        std::string s(_size,'N');
        for (unsigned i(0); i<_size; ++i)
        {
            s[i] = get_char(i);
        }
        return s;
    }

    // returns the reverse complement
    std::string
    get_rc_string() const
    {
        std::string s(_size,'N');
        for (unsigned i(0); i<_size; ++i)
        {
            s[i] = get_complement_char(i);
        }
        std::reverse(s.begin(),s.end());
        return s;
    }

    unsigned size() const override
    {
        return _size;
    }

private:
    const uint8_t* _s;
    uint16_t _size;
    uint16_t _offset;
};


//
//
struct string_bam_seq : public bam_seq_base
{
    explicit
    string_bam_seq(const std::string& s)
        : _s(s.c_str()), _size(s.size()) {}

    string_bam_seq(const char* s,
                   const unsigned init_size)
        : _s(s), _size(init_size) {}

    uint8_t
    get_code(pos_t i) const override
    {
        return get_bam_seq_code(get_char(i));
    }

    char
    get_char(const pos_t i) const override
    {
        if (! is_in_range(i)) return 'N';
        return _s[i];
    }

    unsigned size() const override
    {
        return _size;
    }

private:
    const char* _s;
    unsigned _size;
};


/// Coerce a reference_contig_segment to the bam_seq_base interface without copying
///
struct rc_segment_bam_seq : public bam_seq_base
{
    explicit
    rc_segment_bam_seq(const reference_contig_segment& r)
        : _r(r)
    {}

    uint8_t
    get_code(pos_t i) const override
    {
        return get_bam_seq_code(get_char(i));
    }

    char
    get_char(const pos_t i) const override
    {
        return _r.get_base(i);
    }

    unsigned size() const override
    {
        return _r.end();
    }

private:
    const reference_contig_segment& _r;
};

