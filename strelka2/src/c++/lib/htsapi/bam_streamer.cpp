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

#include "blt_util/blt_exception.hh"
#include "blt_util/log.hh"
#include "htsapi/bam_header_util.hh"
#include "htsapi/bam_streamer.hh"

#include <cassert>
#include <cstdlib>
#include <sys/time.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <sys/stat.h>


bam_streamer::
bam_streamer(
    const char* filename,
    const char* referenceFilename,
    const char* region)
    : _is_record_set(false),
      _hfp(nullptr),
      _aeb_fp(NULL), _aib_fp(NULL), fR(NULL), oR(NULL),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/ERR091571/opFiles/tmp1N_0"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/ERR194147/opFiles/tmp1N_0"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1),
      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/Venter/opFiles/tmp1N_0"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/ERR194147/opFiles/parsnip_out"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/ERR091571/opFiles/parsnip_out"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/Venter/opFiles/parsnip_out"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1),
      _hdr(nullptr),
      _hidx(nullptr),
      _hitr(nullptr),
      _record_no(0),
      _stream_name(filename),
      _is_region(false)
{
    assert(nullptr != filename);
    if ('\0' == *filename)
    {
        throw blt_exception("Can't initialize bam_streamer with empty filename\n");
    }

    _hfp = hts_open(filename, "rb");

    if (nullptr == _hfp)
    {
        std::ostringstream oss;
        oss << "Failed to open SAM/BAM/CRAM file for reading: '" << name() << "'";
        throw blt_exception(oss.str().c_str());
    }

    if (nullptr != referenceFilename)
    {
        const std::string referenceFilenameIndex(std::string(referenceFilename) + ".fai");
        hts_set_fai_filename(_hfp, referenceFilenameIndex.c_str());
    }

    _hdr = sam_hdr_read(_hfp);

    if (nullptr == _hdr)
    {
        std::ostringstream oss;
        oss << "Failed to parse header from SAM/BAM/CRAM file: " << name();
        throw blt_exception(oss.str().c_str());
    }
    curData = (uint8_t *)malloc (20000 * sizeof(uint8_t));
    assert (curData != NULL);
    curBamRec = bam_init1();
    std::cout << "BAM streamer constructor" << std::endl;
    split_fasta ();
    split_aebaib2 ();
    if (nullptr == region)
    {
        // setup to read the whole BAM file by default if resetRegion() is not called:
        if (_hdr->n_targets)
        {
            // parse any contig name so that header->hash is created
            // ignore returned tid value, so doesn't matter if fake name
            // exists
            target_name_to_id("fake_name");
        }
    }
    else
    {
	std::cout << "Region " << region << std::endl;
        // read a specific region of the bam file:
        resetRegion(region);
    }
}

bam_streamer::
bam_streamer(
    const char* filename,
    const char* referenceFilename,
    int r, int numTasks,
    const char* region)
    : _is_record_set(false),
      _hfp(nullptr),
      _aeb_fp(NULL), _aib_fp(NULL), fR(NULL), oR(NULL),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/ERR091571/opFiles/tmp1N_0"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1), rank(r), numTasks(numTasks),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/ERR194147/opFiles/tmp1N_0"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1), rank(r), numTasks(numTasks),
      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/Venter/opFiles/tmp1N_0"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1), rank(r), numTasks(numTasks),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/ERR194147/opFiles/parsnip_out"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1), rank(r), numTasks(numTasks),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/ERR091571/opFiles/parsnip_out"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1), rank(r), numTasks(numTasks),
//      curTid(-1), curData(NULL), curBamRec (NULL), totalSegments(500), numOutChunks(0), numMaxChunks(0), maxChunkSize(0), aebaib_prefix("/scratch/03201/tg824998/Venter/opFiles/parsnip_out"), startFile (-1), cur_aeb_rec (-1), total_aeb_rec (-1), endFile (-1), cur_aib_rec (-1), total_aib_rec (-1), rank(r), numTasks(numTasks),
      _hdr(nullptr),
      _hidx(nullptr),
      _hitr(nullptr),
      _record_no(0),
      _stream_name(filename),
      _is_region(false)
{
    assert(nullptr != filename);
    if ('\0' == *filename)
    {
        throw blt_exception("Can't initialize bam_streamer with empty filename\n");
    }

    _hfp = hts_open(filename, "rb");

    if (nullptr == _hfp)
    {
        std::ostringstream oss;
        oss << "Failed to open SAM/BAM/CRAM file for reading: '" << name() << "'";
        throw blt_exception(oss.str().c_str());
    }

    if (nullptr != referenceFilename)
    {
        const std::string referenceFilenameIndex(std::string(referenceFilename) + ".fai");
        hts_set_fai_filename(_hfp, referenceFilenameIndex.c_str());
    }

    _hdr = sam_hdr_read(_hfp);

    if (nullptr == _hdr)
    {
        std::ostringstream oss;
        oss << "Failed to parse header from SAM/BAM/CRAM file: " << name();
        throw blt_exception(oss.str().c_str());
    }
    curData = (uint8_t *)malloc (20000 * sizeof(uint8_t));
    assert (curData != NULL);
    curBamRec = bam_init1();
    std::cout << "BAM streamer constructor" << std::endl;
    split_fasta ();
//    split_aebaib2 ();
    split_aebaib3 ();
    if (nullptr == region)
    {
        // setup to read the whole BAM file by default if resetRegion() is not called:
        if (_hdr->n_targets)
        {
            // parse any contig name so that header->hash is created
            // ignore returned tid value, so doesn't matter if fake name
            // exists
            target_name_to_id("fake_name");
        }
    }
    else
    {
	std::cout << "Region " << region << std::endl;
        // read a specific region of the bam file:
        resetRegion(region);
    }
}

bam_streamer::
~bam_streamer()
{
    if (nullptr != _hitr) hts_itr_destroy(_hitr);
    if (nullptr != _hidx) hts_idx_destroy(_hidx);
    if (nullptr != _hdr) bam_hdr_destroy(_hdr);
    if (nullptr != _hfp)
    {
        const int retval = hts_close(_hfp);
        if (retval != 0)
        {
            log_os << "ERROR: Failed to close SAM/BAM/CRAM file: '" << name() << "'\n";
            std::exit(EXIT_FAILURE);
        }
    }
    if (nullptr != fR) free (fR);
    if (nullptr != oR) free (oR);
    bam_destroy1 (curBamRec);
}

void
bam_streamer::
split_aebaib()
{
    const bam_hdr_t &bhdr = get_header();
    size_t numAebs=0, numAibs=0;
    struct stat sbuf, sbuf2;
    for (int i=0; i<nContigs; i++)
    {
        for (int j=0; j<numMaxChunks; j++)
        {
            char fname[10000];
            sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), i, j);
            int ret = stat (fname, &sbuf);
            if (ret == 0)
                numAebs += sbuf.st_size/sizeof(fullRec);
            sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), i, j);
            ret = stat (fname, &sbuf);
            if (ret == 0)
                numAibs += sbuf.st_size/sizeof(otherRec);
        }
    }
    std::cout << "NumAebs " << numAebs << " NumAibs " << numAibs << std::endl;
    size_t rSize = 100000;
    size_t curSize=0, nextSize=rSize;

    std::cout << "regions " << std::endl;

    fullRec frec;
    otherRec orec;
    for (int i=0; i<nContigs; i++)
    {
        int minPos = -1, minChunk =-1;
        int maxPos = 0, maxChunk = -1;
        int lastPos = 0;
        if (i > 24) break;
        for (int j=0; j<numMaxChunks; j++)
        {
            char fname[10000];
            sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), i, j);
            FILE *tmpfp = fopen (fname, "rb");
            if (tmpfp != NULL)
            {
                assert (fread (&frec, sizeof(fullRec), 1, tmpfp) == 1);
                if (minPos == -1)
                {
                    minPos = frec.pos;
                    minChunk = j;
                }
//                std::cout << "0a-C" << i << ":" << frec.pos << std::endl;
                fclose (tmpfp);
            }
            sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), i, j);
            tmpfp = fopen (fname, "rb");
            if (tmpfp != NULL)
            {
                assert (fread (&orec, sizeof(otherRec), 1, tmpfp) == 1);
                if ((minPos == -1) || (orec.pos < minPos)) 
                {
                    minPos = orec.pos;
                    minChunk = j;
                }

//                std::cout << "0b-C" << i << ":" << orec.pos << std::endl;
                fclose (tmpfp);
            }
            if (minChunk == j)
            {
//                std::cout << "0-C" << i << ":" << minPos << std::endl;
                lastPos = minPos;
            }
            sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), i, j);
            size_t nrecs=0;
            int ret = stat (fname, &sbuf);
            if (ret == 0)
            {
                 nrecs = sbuf.st_size/sizeof(fullRec);///vasu
                 while (curSize+nrecs  > nextSize)
                 {
                     tmpfp = fopen (fname, "rb");
                     assert (tmpfp != NULL);
                     fseek (tmpfp, (nextSize-curSize)*sizeof(fullRec), SEEK_SET);
                     assert (fread (&frec, sizeof(fullRec), 1, tmpfp) == 1);
                     std::cout << "1-" << bhdr.target_name[i] << ":" << lastPos << "-" << frec.pos-1 << std::endl;
                     regions.push_back (std::string (bhdr.target_name[i]) + std::string (":") + std::string (std::to_string (lastPos)) + std::string ("-") + std::string (std::to_string (frec.pos-1)));
                     assert ((lastPos != frec.pos) && (lastPos < frec.pos-1));
                     lastPos = frec.pos;
                     fclose (tmpfp);
                     nextSize += rSize;
                 }
            }
            tmpfp = fopen (fname, "rb");
            if (tmpfp != NULL)
            {
                fseek (tmpfp, (nrecs-1)*sizeof(fullRec), SEEK_SET);
                assert (fread (&frec, sizeof(fullRec), 1, tmpfp) == 1);
                if (frec.pos > maxPos)
                    maxPos = frec.pos;
//                std::cout << "2a-C" << i << ":" << frec.pos << std::endl;
                fclose (tmpfp);
            }
            sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), i, j);
            tmpfp = fopen (fname, "rb");
            if (tmpfp != NULL)
            {
                fseek (tmpfp, 0, SEEK_END);
                fseek (tmpfp, ftell(tmpfp)-sizeof(otherRec), SEEK_SET);
                assert (fread (&orec, sizeof(otherRec), 1, tmpfp) == 1);
                if (orec.pos > maxPos) maxPos = orec.pos;
//                std::cout << "2b-C" << i << ":" << orec.pos << std::endl;
                fclose (tmpfp);
            }
            curSize += nrecs;
            nextSize = curSize + rSize;
        }
        if ((maxPos > 0) && (lastPos != maxPos))
            std::cout << "2-" << bhdr.target_name[i] << ":" << lastPos << "-" << maxPos << std::endl;
    }
}

void
bam_streamer::
split_aebaib3()
{
    const bam_hdr_t &bhdr = get_header();
    size_t rSize = 100000;
/*    size_t *numAebRecs = (size_t *)calloc (nContigs*numMaxChunks, sizeof(size_t));
    size_t *numAibRecs = (size_t *)calloc (nContigs*numMaxChunks, sizeof(size_t));
    size_t *prefixSumRecs = (size_t *)calloc (nContigs*numMaxChunks+1, sizeof(size_t));
    assert ((numAebRecs != NULL) && (numAibRecs != NULL) && (prefixSumRecs != NULL));
    size_t numAebs=0, numAibs=0;
    struct stat sbuf, sbuf2;
    for (int i=0; i<nContigs; i++)
    {
        for (int j=0; j<numMaxChunks; j++)
        {
            char fname[10000];
            sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), i, j);
            int ret = stat (fname, &sbuf);
            if (ret == 0)
            {
                numAebRecs[i*numMaxChunks + j] = sbuf.st_size/sizeof(fullRec);
                numAebs += numAebRecs[i*numMaxChunks + j];
            }
            sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), i, j);
            ret = stat (fname, &sbuf);
            if (ret == 0)
            {
                numAibRecs[i*numMaxChunks + j] = sbuf.st_size/sizeof(otherRec);
                numAibs += numAibRecs[i*numMaxChunks + j];
            }
            prefixSumRecs[i*numMaxChunks + j + 1] = prefixSumRecs[i*numMaxChunks + j] + numAebRecs[i*numMaxChunks + j] + numAibRecs[i*numMaxChunks + j];
        }
    }
    std::cout << "Rank " << rank << "/" << numTasks << " split_aebaib3: NumAebs " << numAebs << " NumAibs " << numAibs << " Total " << prefixSumRecs[nContigs*numMaxChunks] << std::endl;
    size_t startAebRec = rank * prefixSumRecs[nContigs*numMaxChunks] / numTasks;
    assert (startAebRec < prefixSumRecs[nContigs*numMaxChunks]);
    size_t endAebRec = (rank+1) * prefixSumRecs[nContigs*numMaxChunks] / numTasks;
    if (endAebRec > prefixSumRecs[nContigs*numMaxChunks]) endAebRec = prefixSumRecs[nContigs*numMaxChunks];
    int startSegment=-1, endSegment=-1;
    size_t sumSoFar=0;
    for (int i=0; i<nContigs*numMaxChunks; i++)
    {
        if ((prefixSumRecs[i] >= startAebRec) && (startSegment==-1))
            startSegment = i;
        if ((prefixSumRecs[i] >= endAebRec) && (endSegment==-1))
            endSegment = i;
    }
    std::cout << "Rank " << rank << "/" << numTasks << " split_aebaib3: startSegment " << startSegment << " endSegment " << endSegment << std::endl; */
    assert (rSize%50000 == 0);//50000 - hist file range created by sampa
    int numSegsPerRegion = rSize/50000;
    assert (numSegsPerRegion > 0);
    for (int s=0; s<nContigs*numMaxChunks; s++)
    {
        int i = s/numMaxChunks;
        if (i>24) break;
        int j = s%numMaxChunks;
        char fname[10000];
        sprintf (fname, "%s/C%d_%d_sorted.hist", aebaib_prefix.c_str(), i, j);
        int prevOfs, curOfs, curSeg=numSegsPerRegion;
        FILE *fp_hist = fopen (fname, "r");
        if (fp_hist != NULL)
        {
	    fread (&prevOfs, sizeof(int), 1, fp_hist);
            while (!feof(fp_hist))
            {
                while ((!feof (fp_hist)) && (curSeg > 0))
                {
	            fread (&curOfs, sizeof(int), 1, fp_hist);
                    curSeg--;
                }
                curSeg = numSegsPerRegion;
                if (curOfs > prevOfs)
                {
                    if (rank == 0)
                        std::cout << "Rank " << rank << "/" << numTasks << " split_aebaib3 " << bhdr.target_name[i] + std::string (":") + std::string (std::to_string (prevOfs)) + std::string ("-") + std::string (std::to_string (curOfs-1)) << std::endl;
                    regions.push_back (std::string (bhdr.target_name[i] + std::string (":") + std::string (std::to_string (prevOfs)) + std::string ("-") + std::string (std::to_string (curOfs-1))));
                }
                prevOfs = curOfs;
            }
            fclose (fp_hist);
        }
    }
/*    free (numAebRecs);
    free (numAibRecs);
    free (prefixSumRecs);*/
}

void
bam_streamer::
split_aebaib2()
{
    const bam_hdr_t &bhdr = get_header();
    size_t numAebs=0, numAibs=0;
    struct stat sbuf, sbuf2;
    for (int i=0; i<nContigs; i++)
    {
        for (int j=0; j<numMaxChunks; j++)
        {
            char fname[10000];
            sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), i, j);
            int ret = stat (fname, &sbuf);
            if (ret == 0)
                numAebs += sbuf.st_size/sizeof(fullRec);
            sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), i, j);
            ret = stat (fname, &sbuf);
            if (ret == 0)
                numAibs += sbuf.st_size/sizeof(otherRec);
        }
    }
    std::cout << "NumAebs " << numAebs << " NumAibs " << numAibs << std::endl;
    size_t rSize = 100000;
    size_t curSize=0, nextSize=rSize;

    std::cout << "regions " << std::endl;

    fullRec frec;
    otherRec orec;
    for (int i=0; i<nContigs; i++)
    {
        int minPos = -1, minChunk =-1;
        int maxPos = 0, maxChunk = -1;
        int lastPos = 0;
        if (i > 24) break;
        for (int j=0; j<numMaxChunks; j++)
        {
            char fname[10000];
            sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), i, j);
            FILE *tmpfp = fopen (fname, "rb");
            if (tmpfp != NULL)
            {
                assert (fread (&frec, sizeof(fullRec), 1, tmpfp) == 1);
                if (minPos == -1)
                {
                    minPos = frec.pos;
                    minChunk = j;
                }
//                std::cout << "0a-C" << i << ":" << frec.pos << std::endl;
                fclose (tmpfp);
            }
            sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), i, j);
            tmpfp = fopen (fname, "rb");
            if (tmpfp != NULL)
            {
                assert (fread (&orec, sizeof(otherRec), 1, tmpfp) == 1);
                if ((minPos == -1) || (orec.pos < minPos)) 
                {
                    minPos = orec.pos;
                    minChunk = j;
                }

//                std::cout << "0b-C" << i << ":" << orec.pos << std::endl;
                fclose (tmpfp);
            }
            if (minChunk == j)
            {
//                std::cout << "0-C" << i << ":" << minPos << std::endl;
                lastPos = minPos;
            }
            size_t nrecs=0, nrecs_aeb=0, nrecs_aib=0;
            sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), i, j);
            int ret = stat (fname, &sbuf);
            if (ret == 0)
            {
                 nrecs = sbuf.st_size/sizeof(otherRec);
                 nrecs_aib = sbuf.st_size/sizeof(otherRec);
            }
            sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), i, j);
            ret = stat (fname, &sbuf);
            if (ret == 0)
            {
                 nrecs_aeb = sbuf.st_size/sizeof(fullRec);
                 if (nrecs_aeb > nrecs_aib)
                     nrecs = sbuf.st_size/sizeof(fullRec);///vasu
            }
            float aeb_frac  = (float)nrecs_aeb/nrecs;
            float aib_frac = (float)nrecs_aib/nrecs;
            if (nrecs_aeb > nrecs_aib)
            {
                 sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), i, j);
                 while (curSize+nrecs  > nextSize)
                 {
                     tmpfp = fopen (fname, "rb");
                     assert (tmpfp != NULL);
                     fseek (tmpfp, (nextSize-curSize)*sizeof(fullRec), SEEK_SET);
                     assert (fread (&frec, sizeof(fullRec), 1, tmpfp) == 1);
                     if (rank ==0)
                         std::cout << "f1-" << bhdr.target_name[i] << ":" << lastPos << "-" << frec.pos-1 << " segNo " << j <<" nrecs_aeb " << nrecs_aeb << " nrecs_aib " << nrecs_aib << " aebofs: " << nextSize-curSize << std::endl;
                     regions.push_back (std::string (bhdr.target_name[i]) + std::string (":") + std::string (std::to_string (lastPos)) + std::string ("-") + std::string (std::to_string (frec.pos-1)));
                     assert ((lastPos != frec.pos) && (lastPos < frec.pos-1));
                     lastPos = frec.pos;
                     fclose (tmpfp);
                     nextSize += rSize;
                 }
            }
            else
            {
                 sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), i, j);
                 while (curSize+nrecs  > nextSize)
                 {
                     tmpfp = fopen (fname, "rb");
                     assert (tmpfp != NULL);
                     fseek (tmpfp, (nextSize-curSize)*sizeof(otherRec), SEEK_SET);
                     assert (fread (&orec, sizeof(otherRec), 1, tmpfp) == 1);
                     if (rank == 0)
                         std::cout << "o1-" << bhdr.target_name[i] << ":" << lastPos << "-" << orec.pos-1 << std::endl;
                     regions.push_back (std::string (bhdr.target_name[i]) + std::string (":") + std::string (std::to_string (lastPos)) + std::string ("-") + std::string (std::to_string (orec.pos-1)));
                     assert ((lastPos != orec.pos) && (lastPos < orec.pos-1));
                     lastPos = orec.pos;
                     fclose (tmpfp);
                     nextSize += rSize;
                 }
            }
            
            tmpfp = fopen (fname, "rb");
            if (tmpfp != NULL)
            {
                if (nrecs_aeb > nrecs_aib)
                {
                    fseek (tmpfp, ftell(tmpfp)-sizeof(fullRec), SEEK_SET);
                    assert (fread (&frec, sizeof(fullRec), 1, tmpfp) == 1);
                    if (frec.pos > maxPos)
                        maxPos = frec.pos;
//                  std::cout << "2a-C" << i << ":" << frec.pos << std::endl;
                }
                else
                {
                    fseek (tmpfp, ftell(tmpfp)-sizeof(otherRec), SEEK_SET);
                    assert (fread (&orec, sizeof(otherRec), 1, tmpfp) == 1);
                    if (orec.pos > maxPos)
                        maxPos = orec.pos;
//                  std::cout << "2a-C" << i << ":" << orec.pos << std::endl;
                }
                fclose (tmpfp);
            }
            if (nrecs_aeb > nrecs_aib)
                sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), i, j);
            else
                sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), i, j);
            tmpfp = fopen (fname, "rb");
            if (tmpfp != NULL)
            {
                fseek (tmpfp, 0, SEEK_END);
                if (nrecs_aeb > nrecs_aib)
                {
                    fseek (tmpfp, ftell(tmpfp)-sizeof(otherRec), SEEK_SET);
                    assert (fread (&orec, sizeof(otherRec), 1, tmpfp) == 1);
                    if (orec.pos > maxPos) maxPos = orec.pos;
                }
                else
                {
                    fseek (tmpfp, ftell(tmpfp)-sizeof(fullRec), SEEK_SET);
                    assert (fread (&frec, sizeof(fullRec), 1, tmpfp) == 1);
                    if (frec.pos > maxPos) maxPos = frec.pos;
                }
//                std::cout << "2b-C" << i << ":" << orec.pos << std::endl;
                fclose (tmpfp);
            }
            curSize += nrecs;
            nextSize = curSize + rSize;
        }
//        if ((maxPos > 0) && (lastPos != maxPos))
//            std::cout << "2-" << bhdr.target_name[i] << ":" << lastPos << "-" << maxPos << std::endl;
    }
}

void
bam_streamer::
split_fasta()
{
	const bam_hdr_t &bhdr = get_header();
	int numContigs = bhdr.n_targets;
	nContigs = numContigs;
//	std::cout << "NumContigs " << numContigs << std::endl;
	if (numContigs < 93) totalSegments=numContigs;
	int32_t maxContigSize = 0;
	int i=0;
	for (i=0; i<numContigs; i++)
	{
//		std::cout << "Contig " << i << " len " << bhdr.target_len[i] << std::endl;
		if (bhdr.target_len[i] > maxContigSize) maxContigSize = bhdr.target_len[i];
	}
	int curChunks = 0, curContigSize=maxContigSize;
	while (curChunks < totalSegments)
	{
		curChunks = 0;
		for (i=0; i<numContigs; i++)
		{
			curChunks += (bhdr.target_len[i]/curContigSize)+(bhdr.target_len[i]%curContigSize > 0);
		}
//		printf ("curChunkSize %d --> numChunks %d\n", curContigSize, curChunks);
		if (curChunks < totalSegments)
			curContigSize = curContigSize / 2;
	}
	if (curContigSize < maxContigSize)
		curContigSize = curContigSize * 2;
	int maxChunksPerContig = 0;
	curChunks = 0;
	for (i=0; i<numContigs; i++)
	{
		int curChunksPerContig = (bhdr.target_len[i]/curContigSize)+(bhdr.target_len[i]%curContigSize > 0);
		curChunks += curChunksPerContig;
		if (curChunksPerContig > maxChunksPerContig) maxChunksPerContig = curChunksPerContig;
	}
	numOutChunks = curChunks;
	numMaxChunks = maxChunksPerContig;
	maxChunkSize = curContigSize;
}

static
bool
fexists(const char* filename)
{
    std::ifstream ifile(filename);
    return (! ifile.fail());
}



static
bool
hasEnding(
    const std::string& fullString,
    const std::string& ending)
{
    if (fullString.length() < ending.length()) return false;
    return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
}



// load index if it hasn't been set already:
void
bam_streamer::
_load_index()
{
    /// TODO: Find out whether _hidx can be destroyed after the HTS
    /// iterator is created, in which case this could be a local
    /// variable. Until we know, _hidx should persist for the lifetime
    /// of _hiter
    if (nullptr != _hidx) return;

    std::string index_base(name());

    // hack to allow GATK/Picard bai name convention:
    if ((! fexists((index_base+".bai").c_str())) &&
        (! fexists((index_base+".csi").c_str())) &&
        (! fexists((index_base+".crai").c_str())))
    {
        static const std::string bamext(".bam");
        if (hasEnding(index_base,bamext))
        {
            index_base=index_base.substr(0,index_base.length()-bamext.length());
        }
    }

    _hidx = sam_index_load(_hfp, index_base.c_str());
    if (nullptr == _hidx)
    {
        std::ostringstream oss;
        oss << "BAM/CRAM index is not available for file: " << name();
        throw blt_exception(oss.str().c_str());
    }
}

size_t bam_streamer::getFirstRecordIndexForRange (FILE *fp, int isaibfile, size_t l1, size_t h1, int key)
{
    fullRec frec;
    otherRec orec;
    assert (h1 < 2147483640);
    int h = (int)h1;
    int l = (int)l1;
    int lastrec = h;
    size_t sizePerRec = isaibfile? sizeof(otherRec) : sizeof(fullRec);
    int m=l, iter=0;
    int mpos=-1, prevpos=-1, nextpos=-1;
    while (l <= h)
    {
        m = (int)(l+h)/2;
        mpos=-1; prevpos=-1; nextpos=-1;
        if (!isaibfile)
        {
//            std::cout << "Read aeb file record no. " << m << std::endl;
            fseek (fp, sizeof(fullRec)*m, SEEK_SET);
            assert (fread (&frec, sizeof(fullRec), 1, fp) == 1);
            mpos = frec.pos;
            if (m > 0)
            {
                fseek (fp, sizeof(fullRec)*(m-1), SEEK_SET);
                assert (fread (&frec, sizeof(fullRec), 1, fp) == 1);
                prevpos = frec.pos;
            }
            if (m < lastrec)
            {
                fseek (fp, sizeof(fullRec)*(m+1), SEEK_SET);
                assert (fread (&frec, sizeof(fullRec), 1, fp) == 1);
                nextpos = frec.pos;
            }
        }
        else
        {
            fseek (fp, sizeof(otherRec)*m, SEEK_SET);
            assert (fread (&orec, sizeof(otherRec), 1, fp) == 1);
            mpos = orec.pos;
            if (m > 0)
            {
                fseek (fp, sizeof(otherRec)*(m-1), SEEK_SET);
                assert (fread (&orec, sizeof(otherRec), 1, fp) == 1);
                prevpos = orec.pos;
            }
            if (m < lastrec)
            {
                fseek (fp, sizeof(otherRec)*(m+1), SEEK_SET);
                assert (fread (&orec, sizeof(otherRec), 1, fp) == 1);
                nextpos = orec.pos;
            }
        }
        if (mpos >= key)
        {
            if ((m>0) && (prevpos < key)) return m;
            h = m-1;
        }
        else if (mpos < key)
        {
            if ((m<lastrec) && (nextpos >= key)) return m+1;
            l = m+1;
        }
        iter++;
        if (iter > 25)
            break;
    }
    if (l != 0)
    {
       std::cout << "l,m,h: " << l << "(" << prevpos << ")"  << "," << m << "(" << mpos << ")" << "," << h << "(" << nextpos << ")" << " lastElem " << lastrec << " key " << key << std::endl;
    }
    assert ((l==0) || (l==(lastrec+1)));
    return l;
}

void
bam_streamer::
resetRegion(const char* region)
{
    int32_t referenceContigId, beginPos, endPos;
    parse_bam_region_from_hdr(_hdr, region, referenceContigId, beginPos, endPos);

    try
    {
        resetRegion(referenceContigId, beginPos, endPos);
        _region=region;
        curTid = referenceContigId;
        curStart = beginPos;
        curEnd = endPos;
//        std::cout << "curTid " << curTid << std::endl;
	
	startFile = beginPos / maxChunkSize;
	endFile = endPos / maxChunkSize;
	if (endPos % maxChunkSize == 0) endFile--;
        if (fR == NULL)
            fR = (fullRec *)malloc (sizeof(fullRec) * 10000);
        if (oR == NULL)
            oR = (otherRec *)malloc (sizeof(otherRec) * 10000);
        if ((nullptr == fR) || (nullptr == oR))
        {
            std::ostringstream oss;
            oss << "Cannot allocate memory for storing 1M AEB AIB records\n";
            throw blt_exception(oss.str().c_str());
        }
        cur_aeb_rec = cur_aib_rec = 0;
        total_aeb_rec = total_aib_rec = 0;

	while (startFile <= endFile)
        {
            char fname[10000];
            sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), referenceContigId, startFile);
/*            std::string fname (aebaib_prefix);
            fname << "/C" << referenceContigId << "_" << startFile  << "_sorted.aeb";*/
//            std::cout << "AEB: " << fname << std::endl;
            _aeb_fp = fopen (fname, "rb");
            sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), referenceContigId, startFile);
//            std::cout << "AEB: " << fname << std::endl;
/*            std::string fname (aebaib_prefix);
            fname << "/C" << referenceContigId << "_" << startFile  << "_sorted.aib";*/
            _aib_fp = fopen (fname, "rb");
            if ((_aeb_fp == NULL) && (_aib_fp == NULL))
            {
		startFile++;
            }
            else
                break;
        }
        if (_aeb_fp != NULL)
        {
            fseek (_aeb_fp, 0, SEEK_END);
            size_t numRecs = ftell(_aeb_fp)/sizeof(fullRec);
//            std::cout << "NumAEBRecs " << numRecs << std::endl;
            size_t firstRecordOfs = getFirstRecordIndexForRange (_aeb_fp, 0, 0, numRecs-1, beginPos);
//            if (firstRecordOfs < numRecs)
            {
                fseek (_aeb_fp, firstRecordOfs*sizeof(fullRec), SEEK_SET);
                total_aeb_rec = (int)fread (fR, sizeof(fullRec), 10000, _aeb_fp);
            }
            if (firstRecordOfs == numRecs)
            {
                std::cout << feof(_aeb_fp) << total_aeb_rec << "****\n";
            }

/*            if (total_aeb_rec == 0)
                std::cout << "File C" << referenceContigId << "_" << startFile << "_sorted.aeb size=0" << std::endl;
            while ((total_aeb_rec>0) && (fR[total_aeb_rec-1].pos < beginPos))
                total_aeb_rec = (int)fread (fR, sizeof(fullRec), 1000000, _aeb_fp);
            if (total_aeb_rec > 0)
            {
                while (fR[cur_aeb_rec].pos < beginPos) cur_aeb_rec++;
            }*/
            if (nContigs == 93)
            {
//                std::cout << "File C" << referenceContigId << "_" << startFile << "_sorted.aeb first,last " << fR[0].pos << "," << fR[total_aeb_rec-1].pos << " < [" << beginPos << "," << endPos << "] curRec " << cur_aeb_rec << " (" << firstRecordOfs << ")" << std::endl;
                assert ((total_aeb_rec == 0) || ((fR[total_aeb_rec-1].pos >= beginPos)));
            }
        }
        if (_aib_fp != NULL)
        {
            fseek (_aib_fp, 0, SEEK_END);
            size_t numRecs = ftell(_aib_fp)/sizeof(otherRec);
//            std::cout << "NumAIBRecs " << numRecs << std::endl;
            size_t firstRecordOfs = getFirstRecordIndexForRange (_aib_fp, 1, 0, numRecs-1, beginPos);

//            if (firstRecordOfs < numRecs)
            {
                fseek (_aib_fp, firstRecordOfs*sizeof(otherRec), SEEK_SET);
                total_aib_rec = (int)fread (oR, sizeof(otherRec), 10000, _aib_fp);
            }
            if (firstRecordOfs == numRecs)
            {
//                fseek (_aib_fp, 0, SEEK_END);
                std::cout << feof(_aib_fp) << total_aib_rec << "++++\n";
            }

/*            if (total_aib_rec == 0)
                std::cout << "File C" << referenceContigId << "_" << startFile << "_sorted.aib size=0" << std::endl;
            while ((total_aib_rec>0) && (oR[total_aib_rec-1].pos < beginPos))
                total_aib_rec = (int)fread (oR, sizeof(otherRec), 1000000, _aib_fp);
            if (total_aib_rec > 0)
            {
                while (oR[cur_aib_rec].pos < beginPos) cur_aib_rec++;
            }*/
            if (nContigs == 93)
            {
//                std::cout << "File C" << referenceContigId << "_" << startFile << "_sorted.aib first,last " << oR[0].pos << "," << oR[total_aib_rec-1].pos << " < [" << beginPos << "," << endPos << "] curRec " << cur_aib_rec << " (" << firstRecordOfs << ")" << std::endl;
                assert ((total_aib_rec == 0) || ((oR[total_aib_rec-1].pos >= beginPos)));
            }
        }
    }
    catch (const std::exception& /*e*/)
    {
        log_os << "ERROR: exception while fetching BAM/CRAM region: '" << region
               << "' from file '" << name() << "'\n";
        throw;
    }
    std::cout << "bam_streamer::resetRegion " << region << "-" << std::endl;
}



void
bam_streamer::
resetRegion(
    int referenceContigId,
    int beginPos,
    int endPos)
{
    if (nullptr != _hitr) hts_itr_destroy(_hitr);

    _load_index();

    if (referenceContigId < 0)
    {
        std::ostringstream oss;
        oss << "Invalid region (contig index: " << referenceContigId << ") specified for BAM/CRAM file: " << name();
        throw blt_exception(oss.str().c_str());
    }

    _hitr = sam_itr_queryi(_hidx,referenceContigId,beginPos,endPos);
    if (_hitr == nullptr)
    {
        std::ostringstream oss;
        oss << "Failed to fetch region: #" << referenceContigId << ":" << beginPos << "-" << endPos << " specified for BAM/CRAM file: " << name();
        throw blt_exception(oss.str().c_str());
    }
    _is_region = true;
    _region.clear();

    _is_record_set = false;
    _record_no = 0;
}

void
bam_streamer::
aeb2bam () {
    bam1_core_t *c = &curBamRec->core;
    assert (curTid != -1);
    c->tid = curTid;
    c->pos = fR[cur_aeb_rec].pos-1;

    c->l_qname = 20;
    sprintf ((char *)curData, "%u", _record_no);
    curData[19]=0;

    c->n_cigar = 1;
    uint32_t cigar = fR[cur_aeb_rec].matchLen << 4;
    *((uint32_t *)(curData+20)) = cigar;

    c->l_qseq = fR[cur_aeb_rec].matchLen;
    memcpy (&curData[20+4], fR[cur_aeb_rec].seq, (c->l_qseq+1)>>1);
    memcpy (&curData[20+4+((c->l_qseq+1)>>1)], fR[cur_aeb_rec].quals, c->l_qseq);


//    curBamRec->l_aux=0;//1?
    curData[20+4+((c->l_qseq+1)>>1)+c->l_qseq]=0;

    curBamRec->l_data = 20+4+((c->l_qseq+1)>>1)+c->l_qseq/*+1*/;
    curBamRec->m_data = curBamRec->l_data;
    curBamRec->m_data = kroundup32(curBamRec->m_data);
    curBamRec->data = curData;


    uint32_t endPos = bam_calend(c, (uint32_t *)(curData+1));
    c->bin = bam_reg2bin(c->pos, endPos);
    c->qual = fR[cur_aeb_rec].qual;
    c->flag = fR[cur_aeb_rec].flag;
    c->mtid = -1;
    c->mpos = -1;
    c->isize = 0;
}

void
bam_streamer::
aib2bam () {
    bam1_core_t *c = &curBamRec->core;
    assert (curTid != -1);
    c->tid = curTid;
    c->pos = oR[cur_aib_rec].pos-1;

    c->l_qname = 20;
    sprintf ((char *)curData, "%u", _record_no);
    curData[19]=0;

    c->n_cigar = oR[cur_aib_rec].n_cigar;
    uint32_t *cigars = (uint32_t *)(curData+20);
    c->l_qseq = 0;
    for (int i=0; i<c->n_cigar; i++)
    {
        cigars[i]=oR[cur_aib_rec].cigar[i];
        char c1="MIDNSHP=XB"[oR[cur_aib_rec].cigar[i] & 0xF];
        uint32_t len = (oR[cur_aib_rec].cigar[i])>>4;
        if ((c1=='M') || (c1=='I') || (c1=='S')) c->l_qseq+=len;
    }

    memcpy (&curData[20+(c->n_cigar*4)], oR[cur_aib_rec].seq, (c->l_qseq+1)>>1);
    memcpy (&curData[20+(c->n_cigar*4)+((c->l_qseq+1)>>1)], oR[cur_aib_rec].quals, c->l_qseq);


//    curBamRec->l_aux=0;//1?
    curData[20+(c->n_cigar*4)+((c->l_qseq+1)>>1)+c->l_qseq]=0;

    curBamRec->l_data = 20+(c->n_cigar*4)+((c->l_qseq+1)>>1)+c->l_qseq/*+1*/;
    curBamRec->m_data = curBamRec->l_data;
    curBamRec->m_data = kroundup32(curBamRec->m_data);
    curBamRec->data = curData;


    uint32_t endPos = bam_calend(c, (uint32_t *)(curData+1));
    c->bin = bam_reg2bin(c->pos, endPos);
    c->qual = oR[cur_aib_rec].qual;
    c->flag = oR[cur_aib_rec].flag;
    c->mtid = -1;
    c->mpos = -1;
    c->isize = 0;
}


bool
bam_streamer::
next()
{
    if (nullptr == _hfp) return false;

    if (!_hitr)
        std::cout << "bam_streamer::next " << _hitr << std::endl;
    int ret=-1;
    if (nullptr == _hitr)
    {

        ret = sam_read1(_hfp,_hdr, _brec._bp);

        // Semi-documented sam_read1 API: -1 is expected read failure at end of stream, any other negative value
        // is an error
        if (ret < -1)
        {
            std::ostringstream oss;
            oss << "ERROR: Unknown htslib error value in sam_read1 '" << ret << "' while attempting to read BAM/CRAM file:\n";
            report_state(oss);
            throw blt_exception(oss.str().c_str());
        }
    }
    else
    {
/*        ret = sam_itr_next(_hfp, _hitr, _brec._bp);

        // Re sam_itr_next API: -1 is expected read failure at end of stream. As of v1.5 errors also give a return
        // value of -1. If PR #575 is accepted then errors should return a value less than -1.
        if (ret < -1)
        {
            std::ostringstream oss;
            oss << "ERROR: Unknown htslib error value in sam_read1 '" << ret << "' while attempting to read BAM/CRAM file:\n";
            report_state(oss);
            throw blt_exception(oss.str().c_str());
        }*/
        if (cur_aeb_rec == total_aeb_rec) 
            total_aeb_rec = cur_aeb_rec = 0;
        if (cur_aib_rec == total_aib_rec) 
            total_aib_rec = cur_aib_rec = 0;
    
	struct timeval stime, etime;
        int aebeof=0, aibeof=0;
        if (total_aeb_rec == 0)
        {
            if (_aeb_fp && (!feof (_aeb_fp)))
            {
		gettimeofday (&stime, NULL);
//                std::cout << "aebread start\n";
                total_aeb_rec = fread (fR, sizeof(fullRec), 1000, _aeb_fp);
                gettimeofday (&etime, NULL);
                long elapsed = (etime.tv_sec * 1000000 + etime.tv_usec) - (stime.tv_sec * 1000000 + stime.tv_usec);
//                std::cout << "aebread end" << (double)elapsed/1000000 << "\n";
                if (total_aeb_rec == 0) aebeof = 1;
            }
            else
                aebeof = 1;
//            std::cout << "numAeb " << total_aeb_rec << " numAib " << total_aib_rec << std::endl;
        }
    
        if (total_aib_rec == 0)
        {
            if (_aib_fp && (!feof (_aib_fp)))
            {
		gettimeofday (&stime, NULL);
//                std::cout << "aibread start\n";
                total_aib_rec = fread (oR, sizeof(otherRec), 1000, _aib_fp);
                gettimeofday (&etime, NULL);
                long elapsed = (etime.tv_sec * 1000000 + etime.tv_usec) - (stime.tv_sec * 1000000 + stime.tv_usec);
//                std::cout << "aibread end" << (double)elapsed/1000000 << "\n";
                if (total_aib_rec == 0) aibeof = 1;
            }
            else
                aibeof = 1;
//            std::cout << "numAeb " << total_aeb_rec << " numAib " << total_aib_rec << std::endl;
        }
    
        if (aebeof && aibeof)
        {
            assert ((cur_aeb_rec == 0) && (cur_aib_rec == 0) && (total_aeb_rec == 0) && (total_aib_rec == 0));
            startFile ++;
            if (_aeb_fp) fclose (_aeb_fp);
            if (_aib_fp) fclose (_aib_fp);
            
            while (startFile <= endFile)
            {
                char fname[10000];
                sprintf (fname, "%s/C%d_%d_sorted.aeb", aebaib_prefix.c_str(), curTid, startFile);
        //        std::cout << "AEB: " << fname << std::endl;
                _aeb_fp = fopen (fname, "rb");
                sprintf (fname, "%s/C%d_%d_sorted.aib", aebaib_prefix.c_str(), curTid, startFile);
        //        std::cout << "AIB: " << fname << std::endl;
                _aib_fp = fopen (fname, "rb");
                if ((_aeb_fp == NULL) && (_aib_fp == NULL))
               	    startFile++;
                else
                {
//                    std::cout << "aebaibread start\n";
                    if (_aeb_fp)
                        total_aeb_rec = (int)fread (fR, sizeof(fullRec), 10000, _aeb_fp);
                    if (_aib_fp)
                        total_aib_rec = (int)fread (oR, sizeof(otherRec), 10000, _aib_fp);
//                    std::cout << "aebaibread end\n";
        
                    break;
                }
            }
        }
        int aebValid = 0, aibValid = 0;
        if ((cur_aeb_rec < total_aeb_rec) && (fR[cur_aeb_rec].pos < curEnd)) aebValid = 1;
        if ((cur_aib_rec < total_aib_rec) && (oR[cur_aib_rec].pos < curEnd)) aibValid = 1;
        if (aebValid && (!aibValid || (fR[cur_aeb_rec].pos <= oR[cur_aib_rec].pos)))
        {
            aeb2bam();
            cur_aeb_rec ++;
        }
        else
        {
            aib2bam ();
            cur_aib_rec ++;
        }
    
        if (aebValid || aibValid)
        {
            bam1_t *bR = _brec._bp;
            bR->core = curBamRec->core;
            if (bR->m_data < curBamRec->m_data)
            {
               bR->data = (uint8_t *)realloc(bR->data, curBamRec->m_data);
               assert (bR->data != NULL);
            }
            memcpy (bR->data, curBamRec->data, curBamRec->l_data);
    //        bR->l_aux=0;
            bR->l_data = curBamRec->l_data;
            bR->m_data = curBamRec->m_data;
            ret = 0;
        }
        else
        {
            std::cout << "Region " << curTid << ":" << curStart << "-" << curEnd << " NumRecsProcessed: " << _record_no << " numAEBs " << cur_aeb_rec << " numAIBs " << cur_aib_rec << std::endl;
        }
    }

    _is_record_set=(ret >= 0);
    if (_is_record_set) _record_no++;

    return _is_record_set;
}



const char*
bam_streamer::
target_id_to_name(const int32_t tid) const
{
    // assert(tid < _bfp->header->n_targets);
    if (tid<0)
    {
        static const char unmapped[] = "*";
        return unmapped;
    }
    return _hdr->target_name[tid];
}



int32_t
bam_streamer::
target_name_to_id(const char* seq_name) const
{
    return bam_name2id(_hdr,seq_name);
}



void
bam_streamer::
report_state(std::ostream& os) const
{
    const bam_record* bamp(get_record_ptr());

    os << "\tbam_stream_label: " << name() << "\n";
    if (_is_region && (! _region.empty()))
    {
        os << "\tbam_stream_selected_region: " << _region << "\n";
    }
    if (nullptr != bamp)
    {
        os << "\tbam_stream_record_no: " << record_no() << "\n";
        os << "\tbam_record QNAME/read_number: " << bamp->qname() << "/" << bamp->read_no() << "\n";
        const char* chrom_name(target_id_to_name(bamp->target_id()));
        os << "\tbam record RNAME: " << chrom_name << "\n";
        os << "\tbam record POS: " << bamp->pos() << "\n";
    }
    else
    {
        os << "\tno bam record currently set\n";
    }
}
