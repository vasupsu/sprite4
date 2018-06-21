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
/// \author Mitch Bekritsky
///

#pragma once


/// \brief two-sided binomial exact probability
///
/// This is a two sided binomial exact pval wherein we find the
/// prob of n_success or more extreme number of successes and then
/// double it.
///
///
double
get_binomial_twosided_exact_pval(
    const double p,
    const unsigned n_success,
    const unsigned n_trials);

/// \brief two-sided binomial exact test
///
bool
is_reject_binomial_twosided_exact(
    const double alpha,
    const double p,
    const unsigned n_success,
    const unsigned n_trials);

/// \brief two-sided binomial test chi-sqr approximation
bool
is_reject_binomial_twosided_chi_sqr(
    const double alpha,
    const double p,
    const unsigned n_success,
    const unsigned n_trials);

/// \brief two-sided binomial test
///
/// Find the probability of n_success or more extreme success under B(n_trial,p)
///
/// This function chooses from the two testing methods above (exact/approx) based on trial size
///
bool
is_reject_binomial_twosided(
    const double alpha,
    const double p,
    const unsigned n_success,
    const unsigned n_trials);


/// \brief one-sided binomial exact probability
///
/// probability of n_success or more given B(n_trials,p)
///
/// matches R code: pbinom((n_success-1),n_trials,p,lower.tail=FALSE)
double
get_binomial_gte_n_success_exact_pval(
    const double p,
    const unsigned n_success,
    const unsigned n_trials);


/// \brief one-sided binomial exact test
///
/// tests whether n_success or greater can be rejected under
/// a null hypothesis of B(n_trials,p)
///
/// matches R code: binom.test(n_success, n_trials, p, "greater")$p.value <= alpha
bool
is_reject_binomial_gte_n_success_exact(
    const double alpha,
    const double p,
    const unsigned n_success,
    const unsigned n_trials);

/// returns the minimum number of successes to reject the null hypothesis
/// with a p-value of at most alpha for a given error rate and number of trials
///
/// matches R code 1 + qbinom(alpha, n_trials, p, lower.tail = FALSE)
double
min_count_binomial_gte_exact(
    const double alpha,
    const double p,
    const unsigned n_trials);
