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

#include "starling_run.hh"
#include "starling_pos_processor.hh"
#include "starling_streams.hh"

#include "appstats/RunStatsManager.hh"
#include "blt_util/id_map.hh"
#include "blt_util/log.hh"
#include "common/Exceptions.hh"
#include "htsapi/bam_header_util.hh"
#include "htsapi/vcf_record_util.hh"
#include "starling_common/HtsMergeStreamerUtil.hh"
#include "starling_common/ploidy_util.hh"
#include "starling_common/starling_pos_processor_util.hh"
#include <iostream>
#include <sys/time.h>
#ifdef USE_MPI
	#include <mpi.h>
#endif
#ifdef USE_OMP
	#include <pthread.h>
#endif

namespace INPUT_TYPE
{
enum index_t
{
    CANDIDATE_INDELS,
    FORCED_GT_VARIANTS,
    PLOIDY_REGION,
    NOCOMPRESS_REGION,
    CALL_REGION
};
}


/// \brief Parse sample names from ploidy VCF and match these to expected sample names
///
/// 1. Validate sample count/sample name match to input alignment files
/// 2. Construct a sample index translation map
///
/// \param[out] sampleIndexToPloidyVcfSampleIndex translate from the primary sample index order
///                 (derived from input alignment file order) to the VCF sample index (index in VCF sample columns)
///
/// TODO better place to move this fn?
static
void
mapVcfSampleIndices(
    const vcf_streamer& vcfStream,
    const std::vector<std::string>& sampleNames,
    std::vector<unsigned>& sampleIndexToPloidyVcfSampleIndex)
{
    using namespace illumina::common;

    id_set<std::string> vcfSampleSet;
    const unsigned vcfSampleCount(vcfStream.getSampleCount());
    for (unsigned vcfSampleIndex(0); vcfSampleIndex < vcfSampleCount; ++vcfSampleIndex)
    {
        const std::string vcfSampleName(vcfStream.getSampleName(vcfSampleIndex));

        if (vcfSampleSet.test_key(vcfSampleName))
        {
            std::ostringstream oss;
            oss << "ERROR: repeated entry for sample name '" << vcfSampleName << "' in ploidy VCF file: '" << vcfStream.name() << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }

        vcfSampleSet.insert_key(vcfSampleName);
    }

    sampleIndexToPloidyVcfSampleIndex.clear();
    for (const auto& sampleName : sampleNames)
    {
        const auto maybeId(vcfSampleSet.get_optional_id(sampleName));
        if (not maybeId)
        {
            std::ostringstream oss;
            oss << "ERROR: no entry for sample name '" << sampleName << "' in ploidy VCF file: '" << vcfStream.name() << "'\n";
            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
        }
        else
        {
            sampleIndexToPloidyVcfSampleIndex.push_back(*maybeId);
        }
    }
}



static
void
callRegion(
    const starling_options& opt,
    const AnalysisRegionInfo& regionInfo,
    const starling_streams& fileStreams,
    const std::vector<unsigned>& sampleIndexToPloidyVcfSampleIndex,
    const unsigned ploidyVcfSampleCount,
    starling_read_counts& readCounts,
    reference_contig_segment& ref,
    HtsMergeStreamer& streamData,
    starling_pos_processor& posProcessor)
{
    using namespace illumina::common;

    const unsigned sampleCount(opt.alignFileOpt.alignmentFilenames.size());

    posProcessor.resetRegion(regionInfo.regionChrom, regionInfo.regionRange);
    streamData.resetRegion(regionInfo.streamerRegion.c_str());
    setRefSegment(opt, regionInfo.regionChrom, regionInfo.refRegionRange, ref);
    int print =1;
    while (streamData.next())
    {
        const pos_t currentPos(streamData.getCurrentPos());
        const HTS_TYPE::index_t currentHtsType(streamData.getCurrentType());
        const unsigned currentIndex(streamData.getCurrentIndex());

        // wind posProcessor forward to position behind buffer head:
        posProcessor.set_head_pos(currentPos-1);

        if       (HTS_TYPE::BAM == currentHtsType)
        {
            // Remove the filter below because it's not valid for
            // RNA-Seq case, reads should be selected for the report
            // range by the bam reading functions
            //
            // /// get potential bounds of the read based only on current_pos:
            // const known_pos_range any_read_bounds(current_pos-maxIndelSize,current_pos+MAX_READ_SIZE+maxIndelSize);
            // if( posProcessor.is_range_outside_report_influence_zone(any_read_bounds) ) continue;

            // Approximate begin range filter: (removed for RNA-Seq)
            //if((current_pos+MAX_READ_SIZE+maxIndelSize) <= rlimit.begin_pos) continue;
            if (print)
            {
                const bam_streamer& bstr = streamData.getCurrentBamStreamer();
//                std::cout << " Process region " << regionInfo.streamerRegion << "C:" << bstr.curTid << ":" << bstr.startFile << "(" << bstr.cur_aeb_rec << "/" << bstr.total_aeb_rec << "," << bstr.cur_aib_rec << "/" << bstr.total_aib_rec << ")" << "-" << bstr.endFile << std::endl;
                print = 0;
            }
            processInputReadAlignment(opt, ref, streamData.getCurrentBamStreamer(),
                                      streamData.getCurrentBam(), currentPos,
                                      readCounts, posProcessor, currentIndex);
        }
        else if (HTS_TYPE::VCF == currentHtsType)
        {
            assertExpectedVcfReference(ref, streamData.getCurrentVcfStreamer());
            const vcf_record& vcfRecord(streamData.getCurrentVcf());
            if     (INPUT_TYPE::CANDIDATE_INDELS == currentIndex)     // process candidate indels input from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    process_candidate_indel(opt.maxIndelSize, vcfRecord, posProcessor);
                }
                else
                {
                    log_os << "WARNING: candidate indel vcf variant record cannot be categorized as indel:\n";
                    streamData.getCurrentVcfStreamer().report_state(log_os);
                }
            }
            else if (INPUT_TYPE::FORCED_GT_VARIANTS == currentIndex)     // process forced genotype tests from vcf file(s)
            {
                if (vcfRecord.is_indel())
                {
                    static const unsigned sample_no(0);
                    static const bool is_forced_output(true);
                    process_candidate_indel(opt.maxIndelSize, vcfRecord, posProcessor, sample_no, is_forced_output);
                }
                else if (vcfRecord.is_snv() or vcfRecord.is_ref_site())
                {
                    posProcessor.insert_forced_output_pos(vcfRecord.pos - 1);
                }
                else
                {
                    std::ostringstream oss;
                    oss << "ERROR: forcedGT vcf variant record cannot be categorized as SNV or indel:\n";
                    streamData.getCurrentVcfStreamer().report_state(oss);
                    BOOST_THROW_EXCEPTION(LogicException(oss.str()));
                }
            }
            else if (INPUT_TYPE::PLOIDY_REGION == currentIndex)
            {
                std::vector<unsigned> samplePloidy;
                known_pos_range2 ploidyRange;
                try
                {
                    parsePloidyFromVcf(ploidyVcfSampleCount, vcfRecord.line, ploidyRange, samplePloidy);
                }
                catch (...)
                {
                    log_os << "ERROR: Exception caught while parsing vcf ploidy record\n";
                    streamData.getCurrentVcfStreamer().report_state(log_os);
                    throw;
                }

                for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
                {
                    const unsigned ploidy(samplePloidy[sampleIndexToPloidyVcfSampleIndex[sampleIndex]]);
                    if ((ploidy == 0) || (ploidy == 1))
                    {
                        const bool retval(posProcessor.insert_ploidy_region(sampleIndex, ploidyRange, ploidy));
                        if (!retval)
                        {
                            std::ostringstream oss;
                            const auto& sampleName(fileStreams.getSampleNames()[sampleIndex]);
                            oss << "ERROR: ploidy vcf FORMAT/CN values conflict. Conflict detected in sample '"
                                << sampleName << "' at:\n";
                            streamData.getCurrentVcfStreamer().report_state(oss);
                            BOOST_THROW_EXCEPTION(LogicException(oss.str()));
                        }
                    }
                }
            }
            else
            {
                assert(false && "Unexpected hts index");
            }
        }
        else if (HTS_TYPE::BED == currentHtsType)
        {
            const bed_record& bedRecord(streamData.getCurrentBed());
            if (INPUT_TYPE::NOCOMPRESS_REGION == currentIndex)
            {
                known_pos_range2 range(bedRecord.begin,bedRecord.end);
                posProcessor.insert_nocompress_region(range);
            }
            else if (INPUT_TYPE::CALL_REGION == currentIndex)
            {
                known_pos_range2 range(bedRecord.begin,bedRecord.end);
                posProcessor.insertCallRegion(range);
            }
            else
            {
                assert(false && "Unexpected hts index");
            }
        }
        else
        {
            assert(false && "Invalid input condition");
        }
    }
}

typedef struct workerData
{
    int *procFinished;
    uint32_t *endRegion, *regionNumber, *numFinishedTasks;
    pthread_mutex_t *lock;
    pthread_cond_t *cond;
    int numTasks, rank;
}workerData;
#define RANGE_MSG_OFS 3
#define MIN_RANGES 3
void *thread_worker (void *data)
{
        workerData *wdata = (workerData *)data;
        int buf;
        MPI_Request req;
        MPI_Status stat;
        int rank = wdata->rank;
        int firstMsg = 1;

        int remainingCountReceivedCount=0, maxCount=0, maxRank=0;
        while (*wdata->numFinishedTasks < wdata->numTasks)
        {
            if (firstMsg)
            {
                MPI_Irecv (&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &req);
                while (1)
                {
                    int flag=0;
                    MPI_Test (&req, &flag, &stat);
                    if (flag != 0)
                        break;
                    else
                        sleep (10);
                }
                firstMsg=0;
            }
            else
                MPI_Recv (&buf, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
            int tag = stat.MPI_TAG;
            int source = stat.MPI_SOURCE;
            if (tag == 0)//source finished
            {
                if (wdata->procFinished[source]==0)
                {
                    (*wdata->numFinishedTasks)++;
                    wdata->procFinished[source]=1;
                }
//                std::cout << "At [" << rank << "] Rank " << source << " Done. numFinishedTasks: "  << *wdata->numFinishedTasks << std::endl;
                buf = (int)(*wdata->endRegion-*wdata->regionNumber)/wdata->numTasks;
                MPI_Isend (&buf, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &req);
            }
            else if (tag == 1)//get number of remaining regions to be processed
            {
//                std::cout << "At [" << rank << "] Rank " << source << " Remaining regions " << buf << std::endl;
                assert (remainingCountReceivedCount<wdata->numTasks);
                remainingCountReceivedCount++;
                if (maxCount < buf)
                {
                    maxCount = buf;
                    maxRank = source;
                }
                if (remainingCountReceivedCount == wdata->numTasks)
                {
                    if (maxCount > MIN_RANGES)
                    {
                        std::cout << "At [" << rank << "] MaxRemaining count " << maxCount << " maxRemRank " << maxRank << std::endl;
                        MPI_Send (&buf, 1, MPI_INT, maxRank, 2, MPI_COMM_WORLD);
                    }
                    remainingCountReceivedCount = 0;
                    maxCount = 0;
                    maxRank = 0;
                }
            }
            else if (tag == 2)//get request that source intends to steal work
            {
                std::cout << "At [" << rank << "] Received request for work stealing from " << source << std::endl;
                pthread_mutex_lock (wdata->lock);
                int rem = (int)(*wdata->endRegion-*wdata->regionNumber)/wdata->numTasks;
                if (rem > MIN_RANGES)
                {
                    int numRegionsToSend = rem/2;
                    int newEnd = *wdata->endRegion - numRegionsToSend*wdata->numTasks;
                    for (buf = newEnd + 1; buf <= *wdata->endRegion; buf++)
                    {
                        if (buf % wdata->numTasks == rank) break;
                    }
                    assert (buf % wdata->numTasks == rank);
                    for (int i=0; i<numRegionsToSend; i++)
                    {
                        assert (buf+i*wdata->numTasks <= *wdata->endRegion);
                    }
                    *wdata->endRegion = newEnd;
                    numRegionsToSend += RANGE_MSG_OFS;
                    MPI_Isend (&buf, 1, MPI_INT, source, numRegionsToSend, MPI_COMM_WORLD, &req);
                }
                pthread_mutex_unlock(wdata->lock);
            }
            else if (tag > RANGE_MSG_OFS)
            {
                std::cout << "At [" << rank << "] Received " << tag-RANGE_MSG_OFS << " work starting at " << buf << " from " << source << std::endl;
                pthread_mutex_lock (wdata->lock);
                *wdata->regionNumber = buf;
                *wdata->endRegion = buf+(tag-RANGE_MSG_OFS-1)*wdata->numTasks;
                pthread_cond_signal(wdata->cond);
                pthread_mutex_unlock(wdata->lock);
            }
            else
            {
                 std::cout << "At [" << rank << "] Rank " << source << "TAG 3 not possible\n";

            }
        }
        std::cout << "Rank " << rank << " Comm Done." << std::endl;
        pthread_mutex_lock(wdata->lock);
        pthread_cond_signal(wdata->cond);
        pthread_mutex_unlock(wdata->lock);
}

void
starling_run(
    const prog_info& pinfo,
    const starling_options& opt)
{
    using namespace illumina::common;

    int rank=0, numTasks=1;
#ifdef USE_MPI
    MPI_Comm_size(MPI_COMM_WORLD,&numTasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#endif
//    numTasks = 1536;
//    rank = 957;
    // ensure that this object is created first for runtime benchmark
    RunStatsManager statsManager(opt.segmentStatsFilename + "." + std::to_string(rank) + ".xml");

    opt.validate();

    const starling_deriv_options dopt(opt);
    starling_read_counts readCounts;
    reference_contig_segment ref;

    const unsigned sampleCount(opt.getSampleCount());

    ////////////////////////////////////////
    // setup streamData:
    //
    HtsMergeStreamer streamData(opt.referenceFilename);

    // additional data structures required in the region loop below, which are filled in as a side effect of
    // streamData initialization:
    std::vector<std::reference_wrapper<const bam_hdr_t>> bamHeaders;
    std::vector<std::string> sampleNames;
    unsigned ploidyVcfSampleCount(0);
    std::vector<unsigned> sampleIndexToPloidyVcfSampleIndex;

    {
        std::vector<unsigned> registrationIndices;
        for (unsigned sampleIndex(0); sampleIndex < sampleCount; ++sampleIndex)
        {
            registrationIndices.push_back(sampleIndex);
        }
        bamHeaders = registerAlignments(opt.alignFileOpt.alignmentFilenames, registrationIndices, streamData, rank, numTasks);

        assert(not bamHeaders.empty());
        const bam_hdr_t& referenceHeader(bamHeaders.front());

        static const bool noRequireNormalized(false);
        registerVcfList(opt.input_candidate_indel_vcf, INPUT_TYPE::CANDIDATE_INDELS, referenceHeader, streamData,
                        noRequireNormalized);
        registerVcfList(opt.force_output_vcf, INPUT_TYPE::FORCED_GT_VARIANTS, referenceHeader, streamData);

        for (const bam_hdr_t& bamHeader : bamHeaders)
        {
            sampleNames.push_back(get_bam_header_sample_name(bamHeader));
        }

        if (!opt.ploidy_region_vcf.empty())
        {
            const vcf_streamer& vcfStream(streamData.registerVcf(opt.ploidy_region_vcf.c_str(), INPUT_TYPE::PLOIDY_REGION));
            vcfStream.validateBamHeaderChromSync(referenceHeader);

            mapVcfSampleIndices(vcfStream, sampleNames, sampleIndexToPloidyVcfSampleIndex);
            ploidyVcfSampleCount = vcfStream.getSampleCount();
        }

        if (!opt.gvcf.nocompress_region_bedfile.empty())
        {
            streamData.registerBed(opt.gvcf.nocompress_region_bedfile.c_str(), INPUT_TYPE::NOCOMPRESS_REGION);
        }

        if (! opt.callRegionsBedFilename.empty())
        {
            streamData.registerBed(opt.callRegionsBedFilename.c_str(), INPUT_TYPE::CALL_REGION);
        }
    }

    starling_streams fileStreams(opt, pinfo, bamHeaders, sampleNames);
    starling_pos_processor posProcessor(opt, dopt, ref, fileStreams, statsManager);

    const bam_hdr_t& referenceHeader(bamHeaders.front());
    const bam_header_info referenceHeaderInfo(referenceHeader);

    // parse and sanity check regions
    unsigned supplementalRegionBorderSize(opt.maxIndelSize);
    if (opt.is_short_haplotyping_enabled)
    {
        supplementalRegionBorderSize += ActiveRegion::MaxRefSpanToBypassAssembly;
    }
    const auto& referenceAlignmentFilename(opt.alignFileOpt.alignmentFilenames.front());
    std::vector<AnalysisRegionInfo> regionInfoList;
    std::cout << "Before getStrelkaAnalysisRegions" << std::endl;
    std::cout << "bam_streamer::next() returns " << streamData.next() << std::endl;
    const HTS_TYPE::index_t currentHtsType(streamData.getCurrentType());
    if (HTS_TYPE::BAM == currentHtsType)
    {
        std::cout << "Current HtsType: HTS_TYPE::BAM"  << std::endl;
        const bam_streamer &bstr = streamData.getCurrentBamStreamer();
//        bstr.rank = rank;
        std::cout << "NumRegions: " << bstr.regions.size() << " Original " << opt.regions.size () << std::endl;
        getStrelkaAnalysisRegions(bstr, referenceAlignmentFilename, referenceHeaderInfo, supplementalRegionBorderSize, regionInfoList);
    }
//    getStrelkaAnalysisRegions(opt, referenceAlignmentFilename, referenceHeaderInfo, supplementalRegionBorderSize, regionInfoList);

#ifdef USE_MPI
        std::cout << "[" << rank << "] Number of regions:" << regionInfoList.size() << std::endl;
#endif

    uint32_t regionNumber=/*rank*regionInfoList.size()/numTasks*/rank, endRegion=/*((rank+1)*regionInfoList.size()/numTasks)-1*/regionInfoList.size()-1, numFinishedTasks=0;
/*    if (endRegion > (regionInfoList.size()-1)) 
        endRegion = regionInfoList.size()-1;*/
    printf ("[%d] [%d,%d]\n", rank, regionNumber, endRegion);
    int *procFinished = new int[numTasks];
    for (int i=0; i<numTasks; i++) procFinished[i]=0;
    
    pthread_mutex_t lock;
    pthread_cond_t cond = PTHREAD_COND_INITIALIZER;
    assert (pthread_mutex_init(&lock, NULL) == 0);

    workerData wdata;
    wdata.lock = &lock;
    wdata.cond = &cond;
    wdata.numTasks = numTasks;
    wdata.rank = rank;
    wdata.numFinishedTasks = &numFinishedTasks;
    wdata.procFinished = procFinished;
    wdata.endRegion = &endRegion;
    wdata.regionNumber = &regionNumber;
//#define COMM
#ifdef COMM
    pthread_t thr; 
    assert (pthread_create(&thr, NULL, thread_worker, &wdata) ==0);
/*    if (rank == 0)
    {
        for (uint32_t i=0; i<regionInfoList.size(); i++)
        {
            std::cout << regionInfoList[i].streamerRegion << std::endl;
        }
    }*/
    MPI_Barrier (MPI_COMM_WORLD);
        MPI_Request req;
        MPI_Status stat;
#endif
    struct timeval stime, etime;
    gettimeofday (&stime, NULL);
    while (numFinishedTasks < numTasks)
    {
            while (1)
            {
                const AnalysisRegionInfo& regionInfo = regionInfoList[regionNumber];
//                if (regionNumber%100 == 0)
//                if (/*(rank == 0) ||*/ (rank == 1)/* || (rank == 9) || (rank == 15) || (rank == 28) || (rank == 35) || (rank == 37) || (rank == 46)*/)
                    std::cout << "Rank " << rank << " Process region (" << regionNumber << ") " << regionInfo.streamerRegion << "( End " << endRegion << " )" << std::endl;
                assert (not opt.isUseCallRegions());
                {
//                    if (regionNumber == 17233)
                        callRegion(opt, regionInfo, fileStreams, sampleIndexToPloidyVcfSampleIndex, ploidyVcfSampleCount,
                               readCounts, ref, streamData, posProcessor);
                }
//                    std::cout << "Rank " << rank << " End region (" << regionNumber << ") " << regionInfo.streamerRegion << "( End " << endRegion << " )" << std::endl;
                int continueLoop=1;
                pthread_mutex_lock (&lock);
                
                    regionNumber+=numTasks;
//                regionNumber++;
                    if (regionNumber > endRegion)
                    {
                       continueLoop=0;
                       for (int i=0; i<numTasks; i++)
                       {
#ifdef COMM
                           MPI_Isend (&procFinished[rank], 1, MPI_INT, i, 0, MPI_COMM_WORLD, &req);
#endif
                       }
                    }
                pthread_mutex_unlock (&lock);
                if (!continueLoop) break;
            }
        pthread_mutex_lock (&lock);
        gettimeofday (&etime, NULL);
        std::cout << "Rank " << rank << " Comp Done. " << std::endl;

#ifdef COMM
	std::cout << "Unfinished: ";
	for (int i=0; i<numTasks; i++)
	{
		if (!procFinished[i]) std::cout << i << ", ";
	}
	std::cout << std::endl;
        pthread_cond_wait(&cond, &lock);
        std::cout << "[" << rank << "] after COND wait numFinishedTasks:" << numFinishedTasks << std::endl;
#endif
        pthread_mutex_unlock (&lock);
#ifndef COMM
        break;
#endif
    }
#ifdef COMM
    pthread_join (thr, NULL);
#endif
    long elapsedtime = (etime.tv_sec * 1000000 + etime.tv_usec) - (stime.tv_sec * 1000000 + stime.tv_usec);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "[" << rank << "] Elapsed time " << (double)elapsedtime/1000000 << " sec" << std::endl;
    delete procFinished;

    pthread_mutex_destroy(&lock);
    posProcessor.reset();
}
