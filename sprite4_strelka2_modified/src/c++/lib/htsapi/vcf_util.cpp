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

#include "htsapi/vcf_util.hh"

#include "blt_util/blt_exception.hh"
#include <cassert>
#include <cctype>
#include <ctime>

#include <iostream>
#include <sstream>



std::ostream&
vcf_fileDate(std::ostream& os)
{
    const time_t t(time(nullptr));
    struct tm* ct(localtime(&t));
    assert(nullptr != ct);

    static const unsigned dsize(64);
    char datebuf[dsize];
    const size_t ret(strftime(datebuf,dsize,"%Y%m%d",ct));
    assert(ret!=0);
    return os << datebuf;
}



void
write_vcf_filter(
    std::ostream& os,
    const char* id,
    const char* desc)
{
    os << "##FILTER=<ID=" << id << ",Description=\"" << desc << "\">\n";
}



struct gt_parse_helper
{
    // return is_valid_genotype
    static
    bool
    start(const char* gt,
          std::vector<int>& gti,
          const bool is_badend)
    {
        gti.clear();
        if (isdigit(*gt)) return digit(gt,gti,is_badend);

        switch (*gt)
        {
        case '.' :
            return unknown(gt,gti,is_badend);
        default:
            return false;
        }
    }

private:

    static
    bool
    unknown(const char* gt,
            std::vector<int>& gti,
            const bool is_badend)
    {
        gt++;
        gti.push_back(-1);
        switch (*gt)
        {
        case '\0' :
            return true;
        case '|' :
        case '/' :
            return sep(gt,gti,is_badend);
        default :
            return is_badend;
        }
    }

    static
    bool
    sep(const char* gt,
        std::vector<int>& gti,
        const bool is_badend)
    {
        gt++;
        if (isdigit(*gt)) return digit(gt,gti,is_badend);
        switch (*gt)
        {
        case '.' :
            return unknown(gt,gti,is_badend);
        default :
            return false;
        }
    }

    static
    bool
    digit(const char* gt,
          std::vector<int>& gti,
          const bool is_badend)
    {
        int val(0);
        while (isdigit(*gt))
        {
            val = val*10 + static_cast<int>(*gt-'0');
            gt++;
        }
        gti.push_back(val);

        switch (*gt)
        {
        case '\0' :
            return true;
        case '|' :
        case '/' :
            return sep(gt,gti,is_badend);
        default :
            return is_badend;
        }
    }
};



void
parse_gt(const char* gt,
         std::vector<int>& gti,
         const bool is_allow_bad_end_char)
{
    if (! gt_parse_helper::start(gt,gti,is_allow_bad_end_char))
    {
        std::ostringstream oss;
        oss << "ERROR: can't parse genotype string: '" << gt << "'\n";
        throw blt_exception(oss.str().c_str());
    }
}


std::ostream&
operator<<(std::ostream& os, const VcfGenotype& vcfGt)
{
    VcfGenotypeUtil::writeGenotype(vcfGt, os);
    return os;
}



void
VcfGenotypeUtil::
writeGenotype(
    const uint8_t allele0Index,
    std::ostream& os)
{
    os << static_cast<int>(allele0Index);
}



void
VcfGenotypeUtil::
writeGenotype(
    const uint8_t allele0Index,
    const uint8_t allele1Index,
    std::ostream& os)
{
    os << static_cast<int>(allele0Index) << '/' << static_cast<int>(allele1Index);
}



void
VcfGenotypeUtil::
writeGenotype(
    const VcfGenotype& vcfGt,
    std::ostream& os)
{
    if (vcfGt.isUnknown())
    {
        os << '.';
    }
    else
    {
        if (vcfGt.getPloidy() > 0)
        {
            os << static_cast<int>(vcfGt.getAllele0Index());
        }
        if (vcfGt.getPloidy() > 1)
        {
            os << (vcfGt.getIsPhased() ? '|' : '/');
            os << static_cast<int>(vcfGt.getAllele1Index());
        }
        assert (vcfGt.getPloidy() <= 2);
    }
}

