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

/// derived from ELAND implementation by Tony Cox


#include <cassert>

#ifdef DEBUG_ALN
#include "blt_util/log.hh"
#include <iostream>
#endif



template <typename ScoreType>
template <typename SymIter>
void
GlobalAligner<ScoreType>::
align(
    const SymIter queryBegin, const SymIter queryEnd,
    const SymIter refBegin, const SymIter refEnd,
    AlignmentResult<ScoreType>& result) const
{
    result.clear();

    const AlignmentScores<ScoreType>& scores(this->getScores());

    const size_t querySize(std::distance(queryBegin, queryEnd));
    const size_t refSize(std::distance(refBegin, refEnd));

    assert(0 != querySize);
    assert(0 != refSize);

    _score1.resize(querySize+1);
    _score2.resize(querySize+1);
    _ptrMat.resize(querySize+1, refSize+1);

    ScoreVec* thisSV(&_score1);
    ScoreVec* prevSV(&_score2);

    static const ScoreType badVal(-10000);

#if 0
    auto safeScore = [&] (const int score) -> ScoreType
    {
        return std::max(score, static_cast<int>(badVal));
    };
#endif

    // global alignment of query
    //
    // disallow start from the delete state, control start from insert state with flag
    //
    // query can 'fall-off' the end of a short reference, in which case it will
    // be soft-clipped and each base off the end will be scored as offEdge
    //
    for (unsigned queryIndex(0); queryIndex<=querySize; queryIndex++)
    {
        PtrVal& headPtr(_ptrMat.val(queryIndex,0));
        ScoreVal& val((*thisSV)[queryIndex]);
        headPtr.match = AlignState::MATCH;
        val.match = queryIndex * scores.offEdge;
        headPtr.del = AlignState::MATCH;
        val.del = badVal;
        if (not scores.isAllowEdgeInsertion)
        {
            headPtr.ins = AlignState::MATCH;
            val.ins = badVal;
        }
        else
        {
            headPtr.ins = AlignState::INSERT;
            val.ins = scores.open + (queryIndex * scores.extend);
        }
    }

#ifdef DEBUG_ALN_MATRIX
    // store full matrix of scores to print out later, don't turn this debug option on for large references!
    std::vector<ScoreVec> storeScores;

    storeScores.push_back(*thisSV);
#endif

    BackTrace<ScoreType> btrace;

    {
        unsigned refIndex(0);
        for (SymIter refIter(refBegin); refIter != refEnd; ++refIter, ++refIndex)
        {
            std::swap(thisSV,prevSV);

            {
                // control start from delete state with flag
                PtrVal& headPtr(_ptrMat.val(0,refIndex+1));
                ScoreVal& val((*thisSV)[0]);
                if (not scores.isRequireEdgeDeletion)
                {
                    headPtr.match = AlignState::MATCH;
                    val.match = 0;
                    headPtr.del = AlignState::MATCH;
                    val.del = badVal;
                }
                else
                {
                    headPtr.match = AlignState::MATCH;
                    val.match = badVal;
                    headPtr.del = AlignState::DELETE;
                    val.del = scores.open + ((refIndex+1) * scores.extend);
                }
                headPtr.ins = AlignState::MATCH;
                val.ins = badVal;
            }

            unsigned queryIndex(0);
            for (SymIter queryIter(queryBegin); queryIter != queryEnd; ++queryIter, ++queryIndex)
            {
                // update match
                ScoreVal& headScore((*thisSV)[queryIndex+1]);
                PtrVal& headPtr(_ptrMat.val(queryIndex+1,refIndex+1));
                {
                    const ScoreVal& sval((*prevSV)[queryIndex]);
                    headPtr.match = this->max3(
                                        headScore.match,
                                        sval.match,
                                        sval.del,
                                        sval.ins);

                    headScore.match += ((*queryIter==*refIter) ? scores.match : scores.mismatch);
                }

                // update delete
                {
                    const ScoreVal& sval((*prevSV)[queryIndex+1]);
                    headPtr.del = this->max3(
                                      headScore.del,
                                      sval.match + scores.open,
                                      sval.del,
                                      sval.ins + scores.insertDelete);

                    headScore.del += scores.extend;
                    if (0==refIndex) headScore.del = badVal;
                }

                // update insert
                {
                    const ScoreVal& sval((*thisSV)[queryIndex]);
                    headPtr.ins = this->max3(
                                      headScore.ins,
                                      sval.match + scores.open,
                                      badVal,
                                      sval.ins);

                    headScore.ins += scores.extend;
                    if (0==queryIndex) headScore.ins = badVal;
                }

#ifdef DEBUG_ALN
                log_os << "i1i2: " << queryIndex+1 << " " << refIndex+1 << "\n";
                log_os << headScore.match << ":" << headScore.del << ":" << headScore.ins << "/"
                       << static_cast<int>(headPtr.match) << static_cast<int>(headPtr.del) << static_cast<int>(headPtr.ins) << "\n";
#endif
            }
#ifdef DEBUG_ALN
            log_os << "\n";
#endif

#ifdef DEBUG_ALN_MATRIX
            storeScores.push_back(*thisSV);
#endif

            // record potential backtrace start point, unless full reference sequence must be explained
            if (not scores.isRequireEdgeDeletion)
            {
                const ScoreVal& sval((*thisSV)[querySize]);
                updateBacktrace(sval.match,refIndex+1,querySize,btrace);
            }
        }
    }

    // optionally require that full reference sequence is explained
    if (scores.isRequireEdgeDeletion)
    {
        const ScoreVal& sval((*thisSV)[querySize]);
        updateBacktrace(sval.match,refSize,querySize,btrace, AlignState::MATCH);
        updateBacktrace(sval.del,refSize,querySize,btrace, AlignState::DELETE);
    }

    // optionally allow for trailing insertion
    if (scores.isAllowEdgeInsertion)
    {
        const ScoreVal& sval((*thisSV)[querySize]);
        updateBacktrace(sval.ins,refSize,querySize,btrace, AlignState::INSERT);
    }

    // also allow for the case where query falls-off the end of the reference:
    for (unsigned queryIndex(0); queryIndex<querySize; queryIndex++)
    {
        const ScoreVal& sval((*thisSV)[queryIndex]);
        const ScoreType thisMax(sval.match + (querySize-queryIndex) * scores.offEdge);
        updateBacktrace(thisMax,refSize,queryIndex,btrace);
    }

#ifdef DEBUG_ALN_MATRIX
    std::vector<AlignState::index_t> dumpStates {AlignState::MATCH, AlignState::DELETE, AlignState::INSERT};
    this->dumpTables(queryBegin, queryEnd,
                     refBegin, refEnd,
                     querySize, _ptrMat,
                     dumpStates, storeScores);
#endif

    this->backTraceAlignment(
        queryBegin, queryEnd,
        refBegin, refEnd,
        querySize, refSize,
        _ptrMat,
        btrace, result);
}

