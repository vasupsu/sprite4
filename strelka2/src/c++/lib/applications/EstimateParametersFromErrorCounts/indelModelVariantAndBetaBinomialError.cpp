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

#include "indelModelVariantAndBetaBinomialError.hh"

#include "blt_util/log.hh"
#include "blt_util/math_util.hh"
#include "blt_util/prob_util.hh"

//#define CODEMIN_DEBUG
#define CODEMIN_USE_BOOST
#include "minimize_conj_direction.h"

#include "boost/math/special_functions/beta.hpp"

#include <cmath>

#include <iomanip>
#include <iostream>


namespace
{

namespace MIN_PARAMS4
{
enum index_t
{
    LN_INDEL_ERROR_MEAN,
    LN_INDEL_ERROR_CONCENTRATION,
    LN_THETA,
    SIZE
};
}


//#define DEBUG_MODEL4



static
bool
isSaneVal(const double val)
{
    return (! (std::isnan(val) || std::isinf(val)));
}



static
void
checkSaneVal(const double val)
{
    assert(isSaneVal(val));
}



static
double
getObsLogLhood(
    const double logHomPrior,
    const double logHetPrior,
    const double logNoIndelPrior,
    const double indelErrorAlpha,
    const double indelErrorBeta,
    const double indelBetaDenom,
    const bool isInsert,
    /*    const double logNoIndelRefRate, */
    const ExportedIndelObservations& obs)
{
    static const double log0(-std::numeric_limits<double>::infinity());

    static const double homAltRate(0.99);
    static const double hetAltRate(0.5);

    static const double logHomAltRate(std::log(homAltRate));
    static const double logHomRefRate(std::log(1.-homAltRate));
    static const double logHetRate(std::log(hetAltRate));

    unsigned totalIndelObservations(0);

    if (isInsert)
    {
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::INSERT_1); altIndex<INDEL_SIGNAL_TYPE::DELETE_1; ++altIndex)
        {
            totalIndelObservations += obs.altObservations[altIndex];
        }
    }
    else
    {
        for (unsigned altIndex(INDEL_SIGNAL_TYPE::DELETE_1); altIndex<INDEL_SIGNAL_TYPE::SIZE; ++altIndex)
        {
            totalIndelObservations += obs.altObservations[altIndex];
        }
    }

    // get lhood of homref GT:
    double noindel(log0);
    {
        noindel = std::log(boost::math::beta((totalIndelObservations+indelErrorAlpha),
                                             (obs.refObservations+indelErrorBeta))
                           /indelBetaDenom);
    }

    // get lhood of het and hom GT:
    double het(log0);
    double hom(log0);
    {
        het = logHetRate*(obs.refObservations+totalIndelObservations);

        hom = (logHomAltRate*totalIndelObservations +
               logHomRefRate*obs.refObservations);
    }

    const double mix = log_sum( log_sum(logHomPrior+hom,logHetPrior+het), logNoIndelPrior+noindel);

    return mix;
}



static
double
contextLogLhood(
    const std::vector<ExportedIndelObservations>& observations,
    const double logIndelErrorMean,
    const double logIndelErrorConcentration,
    const bool isInsert,
    const double logTheta)
{
    const double indelErrorMean(std::exp(logIndelErrorMean));
    const double indelErrorConcentration(std::exp(logIndelErrorConcentration));

    checkSaneVal(indelErrorMean);
    checkSaneVal(indelErrorConcentration);

    const double indelErrorAlpha(indelErrorMean*indelErrorConcentration);
    const double indelErrorBeta(indelErrorConcentration*(1.-indelErrorMean));

#ifdef DEBUG_MODEL4
    log_os << "MODEL4: loghood input:"
           << " insert_ab: " << indelErrorAlpha << " " << indelErrorBeta
           << " insert_mc: " << indelErrorMean << " " << indelErrorConcentration
           << " isInsert: " << isInsert
           << " theta: " << std::exp(logTheta)
           << "\n";
#endif

    static const double log2(std::log(2));
    const double logHomPrior(logTheta-log2);
    const double logHetPrior(logTheta);
    const double theta(std::exp(logTheta));
    const double logNoIndelPrior(std::log(1-(theta*3./2.)));

    //const double insertErrorConcentration(insertErrorAlpha+insertErrorBeta);
    //const double insertErrorMean(insertErrorAlpha/insertErrorConcentration);

    const double indelBetaDenom(boost::math::beta(indelErrorAlpha, indelErrorBeta));

    // we haven't set this up for very good numerical stability, so the bounds on alpha/beta are fairly tight:
    assert((indelBetaDenom > 0.) && "Can't process proposed beta distribution parameters");

    double logLhood(0.);
    for (const auto& obs : observations)
    {
        const double mix(getObsLogLhood(logHomPrior, logHetPrior, logNoIndelPrior,
                                        indelErrorAlpha, indelErrorBeta, indelBetaDenom,
                                        isInsert, obs));

#ifdef DEBUG_MODEL4
        log_os << "MODEL4: loghood obs: mix/delta: " << mix << " " << (mix*obs.repeatCount) << "\n";
#endif

        logLhood += (mix*obs.observationCount);
    }

    checkSaneVal(logLhood);

#ifdef DEBUG_MODEL4
    log_os << "MODEL4: loghood output:" << logLhood << "\n";
#endif

    return logLhood;
}


static const double maxConcentration(2000);


struct error_minfunc_model4 : public codemin::minfunc_interface<double>
{
    explicit
    error_minfunc_model4(
        const std::vector<ExportedIndelObservations>& observations,
        const bool isInsert,
        const bool isLockTheta = false)
        : _obs(observations), _isInsert(isInsert), _isLockTheta(isLockTheta)
    {}

    unsigned dim() const override
    {
        return (_isLockTheta ? (MIN_PARAMS4::SIZE-1) : MIN_PARAMS4::SIZE);
    }

    double val(const double* in) override
    {
        minimizerParamsToModelParams(in,_params);
        return -contextLogLhood(_obs,
                                _params[MIN_PARAMS4::LN_INDEL_ERROR_MEAN],
                                _params[MIN_PARAMS4::LN_INDEL_ERROR_CONCENTRATION],
                                _isInsert,
                                (_isLockTheta ? defaultLogTheta : _params[MIN_PARAMS4::LN_THETA]));
    }

    /// normalize the minimization values back to usable parameters
    ///
    /// most values are not valid on [-inf,inf] -- the minimizer doesn't
    /// know this. here is where we fill in the gap:
    ///
    static
    void
    minimizerParamsToModelParams(
        const double* in,
        double* out)
    {
        auto rateSmoother = [](double a) -> double
        {
            static const double triggerVal(1e-3);
            static const double logTriggerVal(std::log(triggerVal));
            if (a>logTriggerVal)
            {
                a = std::log1p(a-logTriggerVal) + logTriggerVal;
            }
            return (a>maxLogRate ? maxLogRate-std::abs(a-maxLogRate) : a);
        };

#if 0
        auto locusRateSmoother = [](double a) -> double
        {
            static const double triggerVal(0.8);
            static const double logTriggerVal(std::log(triggerVal));
            if (a>logTriggerVal)
            {
                a = std::log(1+(a-logTriggerVal)) + logTriggerVal;
            }
            return (a>maxLogLocusRate ? maxLogLocusRate-std::abs(a-maxLogLocusRate) : a);
        };
#endif

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

            if (a>logTriggerVal)
            {
                a = std::log1p(a-logTriggerVal) + logTriggerVal;
            }
            return (a>maxLogTheta ? maxLogTheta-std::abs(a-maxLogTheta) : a);
        };

        out[MIN_PARAMS4::LN_INDEL_ERROR_MEAN] = rateSmoother(in[MIN_PARAMS4::LN_INDEL_ERROR_MEAN]);
        out[MIN_PARAMS4::LN_INDEL_ERROR_CONCENTRATION] = std::log(softMaxTransform(in[MIN_PARAMS4::LN_INDEL_ERROR_CONCENTRATION],0.0,maxConcentration));
        out[MIN_PARAMS4::LN_THETA] = thetaSmoother(in[MIN_PARAMS4::LN_THETA]);
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

private:
    const std::vector<ExportedIndelObservations>& _obs;
    bool _isInsert;
    bool _isLockTheta;
    double _params[MIN_PARAMS4::SIZE];
};

const double error_minfunc_model4::defaultLogTheta = std::log(1e-4);
const double error_minfunc_model4::maxLogTheta = std::log(0.4);
const double error_minfunc_model4::maxLogRate = std::log(0.35);



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
        for (unsigned altIndex(altBeginIndex); altIndex<altEndIndex; ++altIndex)
        {
            totalAltObservations += obs.altObservations[altIndex];
        }

        sigTotal.ref += (obs.refObservations*obs.observationCount);
        sigTotal.alt += (totalAltObservations*obs.observationCount);
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
    const double indelErrorMean,
    const double indelErrorConcentration,
    const double theta,
    std::ostream& os)
{
    static const std::string sep(", ");

//    const double indelErrorConcentration(indelErrorAlpha+indelErrorBeta);
//    const double indelErrorMean(indelErrorAlpha/indelErrorConcentration);

    const double indelErrorAlpha(indelErrorMean*indelErrorConcentration);
    const double indelErrorBeta(indelErrorConcentration*(1.-indelErrorMean));

    os << std::setprecision(10);
    os << context << "_" << extendedContextTag
       << sep << data.excludedRegionSkipped
       << sep << (sigTotal.locus + data.depthSkipped)
       << sep << sigTotal.locus
       << sep << sigTotal.ref
       << sep << sigTotal.alt
       << sep << iter
       << sep << loghood
       << sep << indelErrorAlpha
       << sep << indelErrorBeta
       << sep << indelErrorMean
       << sep << indelErrorConcentration
       << sep << theta
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
    for (unsigned indelTypeIndex(0); indelTypeIndex<2; ++indelTypeIndex)
    {
        const bool isInsert(indelTypeIndex==0);
        double minParams[MIN_PARAMS4::SIZE];

        unsigned iter;
        double x_all_loghood;
        {
            static const double line_tol(1e-10);
            static const double end_tol(1e-10);
            static const unsigned max_iter(40);

            // initialize parameter search
            minParams[MIN_PARAMS4::LN_INDEL_ERROR_MEAN] = std::log(1e-3);
            minParams[MIN_PARAMS4::LN_INDEL_ERROR_CONCENTRATION] = softMaxInverseTransform(100.0,0.0,maxConcentration);
            minParams[MIN_PARAMS4::LN_THETA] = error_minfunc_model4::defaultLogTheta;

            static const unsigned SIZE2(MIN_PARAMS4::SIZE*MIN_PARAMS4::SIZE);
            double conjDir[SIZE2];

            std::fill(conjDir,conjDir+SIZE2,0.);
            const unsigned dim(isLockTheta ? MIN_PARAMS4::SIZE-1 : MIN_PARAMS4::SIZE);
            for (unsigned dimIndex(0); dimIndex<dim; ++dimIndex)
            {
                conjDir[dimIndex*(dim+1)] = 0.0005;
            }

            double start_tol(end_tol);
            double final_dlh;
            error_minfunc_model4 errFunc(observations, isInsert, isLockTheta);

            codemin::minimize_conj_direction(minParams,conjDir,errFunc,start_tol,end_tol,line_tol,
                                             x_all_loghood,iter,final_dlh,max_iter);
        }

        // report:
        {
            double normalizedParams[MIN_PARAMS4::SIZE];
            error_minfunc_model4::minimizerParamsToModelParams(minParams,normalizedParams);

            const double theta(std::exp(normalizedParams[MIN_PARAMS4::LN_THETA]));
            const double indelErrorMean(std::exp(normalizedParams[MIN_PARAMS4::LN_INDEL_ERROR_MEAN]));
            const double indelErrorConcentration(std::exp(normalizedParams[MIN_PARAMS4::LN_INDEL_ERROR_CONCENTRATION]));
            const std::string tag(isInsert ? "I" : "D");
            reportIndelErrorRateSet(context, tag.c_str(), sigInsertTotal, data, iter, -x_all_loghood, indelErrorMean, indelErrorConcentration, theta, os);
        }
    }
}

}



void
indelModelVariantAndBetaBinomialError(
    const SequenceErrorCounts& counts)
{
    const bool isLockTheta(false);

    std::ostream& ros(std::cout);

    ros << "context, excludedLoci, nonExcludedLoci, usedLoci, refReads, altReads, iter, lhood, alpha, beta, mean, concentration, theta\n";

    std::vector<ExportedIndelObservations> observations;
    for (const auto& contextInfo : counts.getIndelCounts())
    {
        const auto& context(contextInfo.first);
        const auto& data(contextInfo.second);

        //if (context.repeatCount != 15) continue;

        data.exportObservations(observations);

        if (observations.empty()) continue;

        log_os << "INFO: computing rates for context: " << context << "\n";
        reportExtendedContext(isLockTheta, context, observations, data, ros);
    }
}
