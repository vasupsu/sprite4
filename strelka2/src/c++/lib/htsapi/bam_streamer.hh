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

//#include "starling_common/starling_ref_seq.hh"
#include "starling_common/starling_base_shared.hh"
#include "htsapi/bam_record.hh"
#include "htsapi/sam_util.hh"

#include "boost/utility.hpp"

#include <string>

typedef struct
{
	uint32_t pos;
	uint8_t seq[60];
	uint8_t quals[120];
	uint16_t flag;
	uint8_t qual;
	uint8_t  matchLen;
}fullRec;

typedef struct
{
	uint32_t pos;
	uint16_t flag;
	uint8_t seq[60];
	uint8_t quals[120];
	uint8_t qual;
	uint8_t n_cigar;
	uint16_t cigar[10];
}otherRec;
/// Stream bam records from CRAM/BAM/SAM files. For CRAM/BAM
/// files you can run an indexed stream from a specific genome region.
///
//
// Example use:
// bam_streamer stream("sample1.bam","hg19.fasta");
// stream.resetRegion("chr1:1000000-2000000");
// while (stream.next()) {
//     const bam_record& read(*(stream.get_record_ptr()));
//     if(read.is_unmapped()) unmappedCount++;
// }
//
struct bam_streamer : public boost::noncopyable
{
    /// \param filename CRAM/BAM input file
    /// \param referenceFilename Corresponding reference file. nullptr can be given here to indicate that the
    ///            the reference is not being provided, but many CRAM files cannot be read in this case.
    /// \param region Restrict the stream to iterate through a specific region. The BAM/CRAM input file must be
    ///            indexed for this option to work. If 'region' is not provided, the stream is configured to
    ///            iterate through the entire alignment file.
    bam_streamer(
        const char* filename,
        const char* referenceFilename,
        const char* region = nullptr);

    bam_streamer(
        const char* filename,
        const char* referenceFilename,
        const char* aebaibPrefix,
        const int maxReferenceSegs,
        int rank, int numTasks, 
        const char* region = nullptr);
    ~bam_streamer();

    /// \brief Set new region to iterate over, this will fail if the alignment file is not indexed
    ///
    /// \param region htslib-style region string in format: "chromName:beginPos-endPos", cannot be nullptr
    void
    resetRegion(const char* region);

    /// \brief Set new region to iterate over, this will fail if the alignment file is not indexed
    ///
    /// \param referenceContigId htslib zero-indexed contig id
    /// \param beginPos start position (zero-indexed, closed)
    /// \param endPos end position (zero-indexed, closed)
    void
    resetRegion(
        int referenceContigId,
        int beginPos,
        int endPos);

    bool next();

    const bam_record* get_record_ptr() const
    {
        if (_is_record_set) return &_brec;
        else                return nullptr;
    }

    const char* name() const
    {
        return _stream_name.c_str();
    }

    unsigned record_no() const
    {
        return _record_no;
    }

    void report_state(std::ostream& os) const;
    void aeb2bam ();
    void aib2bam ();
    void split_fasta();
    void split_aebaib();
    void split_aebaib2();
    void split_aebaib3();
    size_t getFirstRecordIndexForRange (FILE *_aeb_fp, int isaibfile, size_t startrec, size_t endrec, int rangePos);
    bam1_t *curBamRec;
    int startFile, endFile;
    int curTid, nContigs;
    uint32_t curStart;
    uint32_t curEnd;
    int cur_aeb_rec, cur_aib_rec, total_aeb_rec, total_aib_rec;
    int rank;
    int numTasks;
//    std::vector<AnalysisRegionInfo> regionInfoList;
//    std::vector<std::string> regions;
    regions_t regions;

    const char*
    target_id_to_name(const int32_t tid) const;

    int32_t
    target_name_to_id(const char* seq_name) const;

    const bam_hdr_t&
    get_header() const
    {
        return *(_hdr);
    }

private:
    void _load_index();

    bool _is_record_set;
    htsFile* _hfp;
    bam_hdr_t* _hdr;
    hts_idx_t* _hidx;
    hts_itr_t* _hitr;
    bam_record _brec;

    // track for debug only:
    unsigned _record_no;
    std::string _stream_name;
    bool _is_region;
    std::string _region;

    std::string aebaib_prefix;
    int totalSegments;
    int numOutChunks, numMaxChunks, maxChunkSize;
    FILE* _aeb_fp;
    FILE* _aib_fp;
    fullRec *fR;
    otherRec *oR;
    uint8_t *curData;
};
