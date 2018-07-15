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

#pragma once


template <typename ScoreType>
struct AlignmentScores
{
    AlignmentScores(
        ScoreType initMatch,
        ScoreType initMismatch,
        ScoreType initOpen,
        ScoreType initExtend,
        ScoreType initOffEdge,
        ScoreType initInsertDelete = 0,
        bool initIsAllowEdgeInsertion = false,
        bool initIsRequireEdgeDeletion = false) :
        match(initMatch),
        mismatch(initMismatch),
        open(initOpen),
        extend(initExtend),
        offEdge(initOffEdge),
        insertDelete(initInsertDelete),
        isAllowEdgeInsertion(initIsAllowEdgeInsertion),
        isRequireEdgeDeletion(initIsRequireEdgeDeletion)
    {}

    const ScoreType match; ///< match score
    const ScoreType mismatch; ///< mismatch score (should be negative)
    const ScoreType open; ///< gap open, gap of length N is scored (open + N * extend) (should be negative)
    const ScoreType extend; ///< gap extend, gap of length N is scored (open + N * extend) (should be negative or zero)
    const ScoreType offEdge; ///< score applied when query goes off the end of an edge (should be negative)
    const ScoreType insertDelete; ///< additional cost (beyond extend) applied to insert->delete transitions
    const bool isAllowEdgeInsertion; ///< are insertions allowed directly on the leading and trailing edges?
    const bool isRequireEdgeDeletion; ///< requires that any remaining reference sequence at the edge is counted as a deletion
};
