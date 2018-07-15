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

#include "blt_util/log.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"
#include "indelModelVariantAndBinomialMixtureError.hh"

//#define CODEMIN_DEBUG
#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

#include <iomanip>
#include <iostream>

namespace
{

namespace MIN_PARAMS3
{
enum index_t
{
    LN_INSERT_ERROR_RATE,
    LN_DELETE_ERROR_RATE,
    LN_NOISY_LOCUS_RATE,
    LN_THETA,
    SIZE
};
}


//#define DEBUG_MODEL3



static
double
getObsLogLhood(
    const double logHomPrior,
    const double logHetPrior,
    const double logAltHetPrior,
    const double logNoIndelPrior,
    const double logInsertErrorRate,
    const double logDeleteErrorRate,
    const double logNoIndelRefRate,
    const ExportedIndelObservations& obs)
{
    static const double log0(-std::numeric_limits<double>::infinity());

    static const double homAltRate(0.99);
    static const double hetAltRate(0.5);

    static const double logHomAltRate(std::log(homAltRate));
    static const double logHomRefRate(std::log(1. - homAltRate));
    static const double logHetRate(std::log(hetAltRate));

    // get lhood of homref GT:
    double noindel(log0);
    {
        unsigned totalInsertObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::INSERT_1); altIndex < INDEL_SIGNAL_TYPE::DELETE_1; ++altIndex)
        {
            totalInsertObservations += obs.altObservations[altIndex];
        }

        unsigned totalDeleteObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::DELETE_1); altIndex < INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            totalDeleteObservations += obs.altObservations[altIndex];
        }

        noindel = (
                      logInsertErrorRate * totalInsertObservations +
                      logDeleteErrorRate * totalDeleteObservations +
                      logNoIndelRefRate * obs.refObservations);
    }

    unsigned maxIndex(0);
    for (unsigned altIndex(1); altIndex < INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
    {
        if (obs.altObservations[altIndex] > obs.altObservations[maxIndex]) maxIndex = altIndex;
    }

    // get lhood of het and hom GT:
    double het(log0);
    double hom(log0);
    {
        // approximate that the most frequent observations is the only potential variant allele:

        unsigned remainingInsertObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::INSERT_1); altIndex < INDEL_SIGNAL_TYPE::DELETE_1; ++altIndex)
        {
            if (altIndex == maxIndex) continue;
            remainingInsertObservations += obs.altObservations[altIndex];
        }

        unsigned remainingDeleteObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::DELETE_1); altIndex < INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            if (altIndex == maxIndex) continue;
            remainingDeleteObservations += obs.altObservations[altIndex];
        }

        // compute lhood of het/hom states given that maxIndex is the variant allele:
        het = (logHetRate * (obs.refObservations + obs.altObservations[maxIndex]) +
               logInsertErrorRate * remainingInsertObservations +
               logDeleteErrorRate * remainingDeleteObservations);

        hom = (logHomAltRate * obs.altObservations[maxIndex] +
               logHomRefRate * obs.refObservations +
               logInsertErrorRate * remainingInsertObservations +
               logDeleteErrorRate * remainingDeleteObservations);
    }

    // get lhood of althet GT:
    double althet(log0);
    {
        // approximate that the two most frequent observations are the only potential variant alleles:
        assert(INDEL_SIGNAL_TYPE::SIZE > 1);
        unsigned maxIndex2(maxIndex == 0 ? 1 : 0);
        for (unsigned altIndex(maxIndex2 + 1); altIndex < INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            if (altIndex == maxIndex) continue;
            if (obs.altObservations[altIndex] > obs.altObservations[maxIndex2]) maxIndex2 = altIndex;
        }

        unsigned remainingInsertObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::INSERT_1); altIndex < INDEL_SIGNAL_TYPE::DELETE_1; ++altIndex)
        {
            if (altIndex == maxIndex) continue;
            if (altIndex == maxIndex2) continue;
            remainingInsertObservations += obs.altObservations[altIndex];
        }

        unsigned remainingDeleteObservations(0);
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::DELETE_1); altIndex < INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            if (altIndex == maxIndex) continue;
            if (altIndex == maxIndex2) continue;
            remainingDeleteObservations += obs.altObservations[altIndex];
        }

        // compute lhood of het/hom states given that maxIndex is the variant allele:
        althet = (logHetRate * (obs.altObservations[maxIndex] + obs.altObservations[maxIndex2]) +
                  logHomRefRate * obs.refObservations +
                  logInsertErrorRate * remainingInsertObservations +
                  logDeleteErrorRate * remainingDeleteObservations);
    }

    return log_sum( log_sum(logHomPrior+hom,logHetPrior+het), log_sum(logNoIndelPrior+noindel,logAltHetPrior+althet) );
}



static
double
contextLogLhood(
    const std::vector<ExportedIndelObservations>& observations,
    const double logInsertErrorRate,
    const double logDeleteErrorRate,
    const double logNoisyLocusRate,
    const double logTheta)
{
#ifdef DEBUG_MODEL3
    log_os << "MODEL3: loghood input:"
           << " insert: " << std::exp(logInsertErrorRate)
           << " delete: " << std::exp(logDeleteErrorRate)
           << " noise: " << std::exp(logNoisyLocusRate)
           << " theta: " << std::exp(logTheta)
           << "\n";
#endif

    static const double log2(std::log(2));
    const double logHomPrior(logTheta - log2);
    const double logHetPrior(logTheta);
    const double logAltHetPrior(logTheta * 2);
    const double theta(std::exp(logTheta));
    const double logNoIndelPrior(std::log(1 - (theta * 3. / 2. + (theta * theta))));

    //const double logNoIndelRefRate(std::log(1-std::exp(logInsertErrorRate))+std::log(1-std::exp(logDeleteErrorRate)));
    const double logNoIndelRefRate(std::log(1-std::exp(logInsertErrorRate)-std::exp(logDeleteErrorRate)));

    const double logCleanLocusRate(std::log(1 - std::exp(logNoisyLocusRate)));

    static const double cleanLocusIndelRate(1e-8);
    static const double logCleanLocusIndelRate(std::log(cleanLocusIndelRate));
    static const double logCleanLocusRefRate(std::log(1 - cleanLocusIndelRate));

    double logLhood(0.);
    for (const auto& obs : observations)
    {
        const double noisyMix(getObsLogLhood(logHomPrior, logHetPrior, logAltHetPrior, logNoIndelPrior,
                                             logInsertErrorRate, logDeleteErrorRate, logNoIndelRefRate, obs));
        const double cleanMix(getObsLogLhood(logHomPrior, logHetPrior, logAltHetPrior, logNoIndelPrior,
                                             logCleanLocusIndelRate, logCleanLocusIndelRate, logCleanLocusRefRate, obs));

        const double mix(log_sum(logCleanLocusRate + cleanMix, logNoisyLocusRate + noisyMix));

#ifdef DEBUG_MODEL3
        log_os << "MODEL3: loghood obs: noisy/clean/mix/delta: " << noisyMix << " " << cleanMix << " " << mix << " " << (mix*obs.repeatCount) << "\n";
#endif

        logLhood += (mix * obs.observationCount);
    }

#ifdef DEBUG_MODEL3
    log_os << "MODEL3: loghood output:" << logLhood << "\n";
#endif

    return logLhood;
}


struct error_minfunc_model3 : public codemin::minfunc_interface<double>
{
    explicit
    error_minfunc_model3(
        const std::vector<ExportedIndelObservations>& observations,
        const bool isLockTheta = false)
        : _obs(observations), _isLockTheta(isLockTheta)
    {}

    unsigned dim() const override
    {
        return (_isLockTheta ? (MIN_PARAMS3::SIZE - 1) : MIN_PARAMS3::SIZE);
    }

    double val(const double* in) override
    {
        argToParameters(in, _params);
        return -contextLogLhood(_obs,
                                _params[MIN_PARAMS3::LN_INSERT_ERROR_RATE],
                                _params[MIN_PARAMS3::LN_DELETE_ERROR_RATE],
                                _params[MIN_PARAMS3::LN_NOISY_LOCUS_RATE],
                                (_isLockTheta ? defaultLogTheta : _params[MIN_PARAMS3::LN_THETA]));
    }

    /// normalize the minimization values back to usable parameters
    ///
    /// most values are not valid on [-inf,inf] -- the minimizer doesn't
    /// know this. here is where we fill in the gap:
    ///
    static
    void
    argToParameters(
        const double* in,
        double* out)
    {
        auto rateSmoother = [](double a) -> double
        {
            static const double triggerVal(1e-3);
            static const double logTriggerVal(std::log(triggerVal));
            if (a > logTriggerVal)
            {
                a = std::log1p(a - logTriggerVal) + logTriggerVal;
            }
            return (a > maxLogRate ? maxLogRate - std::abs(a - maxLogRate) : a);
        };

        auto locusRateSmoother = [](double a) -> double
        {
            static const double triggerVal(0.8);
            static const double logTriggerVal(std::log(triggerVal));
            if (a > logTriggerVal)
            {
                a = std::log1p(a - logTriggerVal) + logTriggerVal;
            }
            return (a > maxLogLocusRate ? maxLogLocusRate - std::abs(a - maxLogLocusRate) : a);
        };

        // A lot of conditioning is required to keep the model from winding
        // theta around zero and getting confused, here we start applying a
        // second log to the delta above triggerTheta, and finally put a hard stop
        // at logMaxTheta -- hard stops are obviously bad b/c the model can get lost
        // on the flat plane even if the ML value is well below this limit, but
        // in practice this is such a ridiculously high value for theta, that
        // I don't see the model getting trapped.
        auto thetaSmoother = [](double a) -> double
        {
            static const double triggerVal(1e-3);
            static const double logTriggerVal(std::log(triggerVal));

            if (a > logTriggerVal)
            {
                a = std::log1p(a - logTriggerVal) + logTriggerVal;
            }
            return (a > maxLogTheta ? maxLogTheta - std::abs(a - maxLogTheta) : a);
        };

        for (unsigned paramIndex(MIN_PARAMS3::LN_INSERT_ERROR_RATE); paramIndex<MIN_PARAMS3::LN_NOISY_LOCUS_RATE; ++paramIndex)
        {
            out[paramIndex] = rateSmoother(in[paramIndex]);
        }
        out[MIN_PARAMS3::LN_NOISY_LOCUS_RATE] = locusRateSmoother(in[MIN_PARAMS3::LN_NOISY_LOCUS_RATE]);
        out[MIN_PARAMS3::LN_THETA] = thetaSmoother(in[MIN_PARAMS3::LN_THETA]);
    }

#if 0
    // this should help in theory, but in practice the minimizer is more likely to get stuck
    bool
    is_val_computable(
        const double* in) override
    {
        if (in[MIN_PARAMS3::LN_INSERT_ERROR_RATE]>maxLogRate) return false;
        if (in[MIN_PARAMS3::LN_DELETE_ERROR_RATE]>maxLogRate) return false;
        if (in[MIN_PARAMS3::LN_NOISY_LOCUS_RATE]>maxLogLocusRate) return false;
        if (in[MIN_PARAMS3::LN_THETA]>maxLogTheta) return false;
        return true;
    }
#endif

    static const double defaultLogTheta;
    static const double maxLogTheta;
    static const double maxLogRate;
    static const double maxLogLocusRate;

private:
    const std::vector<ExportedIndelObservations>& _obs;
    bool _isLockTheta;
    double _params[MIN_PARAMS3::SIZE];
};

const double error_minfunc_model3::defaultLogTheta = std::log(1e-4);
const double error_minfunc_model3::maxLogTheta = std::log(0.4);
const double error_minfunc_model3::maxLogRate = std::log(0.5);
const double error_minfunc_model3::maxLogLocusRate = std::log(1.0);



struct SignalGroupTotal
{
    double ref = 0;
    double alt = 0;
    double locus = 0;
};



static
void
getAltSigTotal(
    const std::vector<ExportedIndelObservations>& observations,
    const unsigned altBeginIndex,
    const unsigned altEndIndex,
    SignalGroupTotal& sigTotal)
{
    for (const ExportedIndelObservations& obs : observations)
    {
        unsigned totalAltObservations(0);
        for (unsigned altIndex(altBeginIndex); altIndex < altEndIndex; ++altIndex)
        {
            totalAltObservations += obs.altObservations[altIndex];
        }

        sigTotal.ref += (obs.refObservations * obs.observationCount);
        sigTotal.alt += (totalAltObservations * obs.observationCount);
        sigTotal.locus += obs.observationCount;
    }
}



static
void
reportIndelErrorRateSet(
    const IndelErrorContext& context,
    const char* extendedContextTag,
    const SignalGroupTotal& sigTotal,
    const IndelErrorData& data,
    unsigned iter,
    const double loghood,
    const double indelErrorRate,
    const double theta,
    const double noisyLocusRate,
    std::ostream& os)
{
    static const std::string sep(", ");

    os << std::setprecision(10);
    os << context << "_" << extendedContextTag
       << sep << data.excludedRegionSkipped
       << sep << (sigTotal.locus + data.depthSkipped)
       << sep << sigTotal.locus
       << sep << sigTotal.ref
       << sep << sigTotal.alt
       << sep << iter
       << sep << loghood
       << sep << indelErrorRate
       << sep << theta
       << sep << noisyLocusRate
       << "\n";
}



static
void
reportExtendedContext(
    const bool isLockTheta,
    const IndelErrorContext& context,
    const std::vector<ExportedIndelObservations>& observations,
    const IndelErrorData& data,
    std::ostream& os)
{
    // Get summary counts for QC purposes. Note these are unrelated to minimization or model:
    SignalGroupTotal sigInsertTotal;
    getAltSigTotal(observations, INDEL_SIGNAL_TYPE::INSERT_1, INDEL_SIGNAL_TYPE::DELETE_1, sigInsertTotal);

    SignalGroupTotal sigDeleteTotal;
    getAltSigTotal(observations, INDEL_SIGNAL_TYPE::DELETE_1, INDEL_SIGNAL_TYPE::SIZE, sigDeleteTotal);


    // initialize conjugate direction minimizer settings and minimize lhood...
    //
    double minParams[MIN_PARAMS3::SIZE];

    unsigned iter;
    double x_all_loghood;
    {
        static const double line_tol(1e-10);
        static const double end_tol(1e-10);
        static const unsigned max_iter(40);

        // initialize parameter search
        minParams[MIN_PARAMS3::LN_INSERT_ERROR_RATE] = std::log(1e-3);
        minParams[MIN_PARAMS3::LN_DELETE_ERROR_RATE] = std::log(1e-3);
        minParams[MIN_PARAMS3::LN_NOISY_LOCUS_RATE] = std::log(0.4);
        minParams[MIN_PARAMS3::LN_THETA] = error_minfunc_model3::defaultLogTheta;

        static const unsigned SIZE2(MIN_PARAMS3::SIZE * MIN_PARAMS3::SIZE);
        double conjDir[SIZE2];

        std::fill(conjDir, conjDir + SIZE2, 0.);
        const unsigned dim(isLockTheta ? MIN_PARAMS3::SIZE - 1 : MIN_PARAMS3::SIZE);
        for (unsigned i(0); i < dim; ++i)
        {
            conjDir[i * (dim + 1)] = 0.0005;
        }

        double start_tol(end_tol);
        double final_dlh;
        error_minfunc_model3 errFunc(observations,isLockTheta);

        codemin::minimize_conj_direction(minParams, conjDir, errFunc, start_tol, end_tol, line_tol,
                                         x_all_loghood, iter, final_dlh, max_iter);
    }

    // report:
    {
        double normalizedParams[MIN_PARAMS3::SIZE];
        error_minfunc_model3::argToParameters(minParams, normalizedParams);

        const double theta(std::exp(normalizedParams[MIN_PARAMS3::LN_THETA]));

        const double insertErrorRate(std::exp(normalizedParams[MIN_PARAMS3::LN_INSERT_ERROR_RATE]));
        const double noisyLocusRate(std::exp(normalizedParams[MIN_PARAMS3::LN_NOISY_LOCUS_RATE]));
        reportIndelErrorRateSet(context, "I", sigInsertTotal, data, iter, -x_all_loghood, insertErrorRate, theta, noisyLocusRate, os);

        const double deleteErrorRate(std::exp(normalizedParams[MIN_PARAMS3::LN_DELETE_ERROR_RATE]));
        reportIndelErrorRateSet(context, "D", sigDeleteTotal, data, iter, -x_all_loghood, deleteErrorRate, theta, noisyLocusRate, os);
    }
}

}



void
indelModelVariantAndBinomialMixtureError(
    const SequenceErrorCounts& counts)
{
    const bool isLockTheta(false);

    std::ostream& ros(std::cout);

    ros << "context, excludedLoci, nonExcludedLoci, usedLoci, refReads, altReads, iter, lhood, errorRate, theta, noisyLocusRate\n";

    std::vector<ExportedIndelObservations> observations;
    for (const auto& contextInfo : counts.getIndelCounts())
    {
        const auto& context(contextInfo.first);
        const auto& data(contextInfo.second);

        data.exportObservations(observations);

        if (observations.empty()) continue;

        log_os << "INFO: computing rates for context: " << context << "\n";
        reportExtendedContext(isLockTheta, context, observations, data, ros);
    }
}
