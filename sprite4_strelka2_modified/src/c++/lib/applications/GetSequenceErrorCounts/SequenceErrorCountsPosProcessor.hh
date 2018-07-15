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

#include "SequenceErrorCountsOptions.hh"
#include "SequenceErrorCountsStreams.hh"
#include "errorAnalysis/SequenceErrorCounts.hh"
#include "starling_common/starling_pos_processor_base.hh"

#include "blt_util/RecordTracker.hh"

///
///
struct SequenceErrorCountsPosProcessor : public starling_pos_processor_base
{
    typedef starling_pos_processor_base base_t;

    SequenceErrorCountsPosProcessor(
        const SequenceErrorCountsOptions& opt,
        const SequenceErrorCountsDerivOptions& dopt,
        const reference_contig_segment& ref,
        const SequenceErrorCountsStreams& fileStreams,
        RunStatsManager& statsManager);

    /// Notify this object that the region is being reset
    void reset() override;

    /// Notify this object that the process is completing (without error)
    ///
    /// This is the point where all one-time statistical summaries which
    /// may span multiple regions are completed.
    void completeProcessing();

    void resetRegion(
        const std::string& chromName,
        const known_pos_range2& reportRegion);

    void
    insertExcludedRegion(
        const known_pos_range2& excludedRange);

    void
    addKnownVariant(
        const vcf_record& knownVariant);

private:

    void
    process_pos_variants_impl(
        const pos_t pos,
        const bool isPosPrecedingReportableRange) override
    {
        if (isPosPrecedingReportableRange) return;
        process_pos_error_counts(pos);
    }

    /// Enumerate all SNV and indel counts for this position
    void
    process_pos_error_counts(
        const pos_t pos);

    bool
    derived_empty() const override
    {
        return _excludedRegions.empty();
    }

    const SequenceErrorCountsOptions& _opt;
    const SequenceErrorCountsDerivOptions& _dopt;
    const SequenceErrorCountsStreams& _streams;

    double _normChromDepth = 0.;
    double _maxChromDepth = 0.;

    SequenceErrorCounts _counts;

    double _max_candidate_normal_sample_depth = -1;

    RegionTracker _excludedRegions;
    RecordTracker _knownVariants;

    /// Count of all non-empty sites eligible to contribute to sequencing error stats
    unsigned long _nonEmptySiteCount = 0;
};
