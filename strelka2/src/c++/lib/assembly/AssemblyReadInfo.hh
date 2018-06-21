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
/// \author Ole Schulz-Trieglaff
///

#pragma once

#include <vector>
#include <string>


/// Information added to each read in the process of assembly
///
struct AssemblyReadInfo
{
    bool isUsed = false;
    bool isFiltered = false; ///< if true, the read was 'used' but filtered out, so there is no meaningful contig id association
    bool isPseudo = false; ///< if true, the read was an assembled contig
    std::vector<unsigned> contigIds; ///< index of the contigs that this read is used in
};


typedef std::vector<std::string> AssemblyReadInput;
typedef std::vector<AssemblyReadInfo> AssemblyReadOutput;
