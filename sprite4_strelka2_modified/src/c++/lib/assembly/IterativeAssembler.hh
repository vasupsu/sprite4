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
/// \author Xiaoyu Chen
///

#pragma once

#include "assembly/AssembledContig.hh"
#include "assembly/AssemblyReadInfo.hh"
#include "options/IterativeAssemblerOptions.hh"


/// \brief run a de-bruijn graph assembler intended for small-scale allele discovery
///
/// the assembler iteratively builds multiple contigs through a range of word sizes
///
/// \param[in] opt assembly parameters
/// \param[in] reads the set of reads to use for the assembly
/// \param[out] assembledReadInfo for each read in 'reads', provide information on if and how it was assembled into a contig
/// \param[out] contigs zero to many assembled contigs
///
void
runIterativeAssembler(
    const IterativeAssemblerOptions& opt,
    AssemblyReadInput& reads,
    AssemblyReadOutput& assembledReadInfo,
    Assembly& contigs);


