r"""
Helper functions for modeling negative binomial on 
phip-seq data
"""


import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm


def nb_logsf(count, size, prob):
    # Assumes count, size, prob are greater than 0
    logsf = st.nbinom.logsf(count, size, prob)

    # value is not near limit of numerical precision, safe to return
    if logsf > -31:
        return logsf

    # Compute log(p) using logaddexp() for higher precision, terminating the
    # summation when contribution from the next term falls below a tolerance
    x = np.round(count + 1)
    accum = st.nbinom.logpmf(x, size, prob)
    while True:
        x = x + 1
        next_term = st.nbinom.logpmf(x, size, prob)
        new_accum = np.logaddexp(accum, next_term)
        if new_accum - accum < 0.0001:
            break
        accum = new_accum
    return new_accum


def nb_mlxp(counts, size, prob):
    mlxp = np.zeros(len(counts))  # default mlxp is 0
    if size > 0:
        for i in range(len(counts)):
            if counts[i] > 0:
                mlxp[i] = -nb_logsf(counts[i], size, prob) / np.log(10)

    return mlxp


def fit_neg_binom(x, p, outlier_reject_scale=10):
    """Fit a Type-I (p=1) or Type-II (p=2) negative binominal to each peptide across all samples
    Uses the "mu, alpha" parametrization in statsmodels
    If x is empty, all returned values are -1
    If the fit fails to converge, all returned values are -2
    Data lying [outlier_reject_scale]*[interquartile range] beyond the 75th percentile
    are dropped from the fit to improve fit convergence.
    """

    # Remove extreme outliers, which are data points that lie beyond
    # 10 times the interquartile range above the 75th percentile
    q75, q25 = np.percentile(x, [75, 25])
    iqr = q75 - q25
    if iqr > 1:
        outlier_cut = q75 + outlier_reject_scale * iqr
        x = np.array([xi for xi in x if xi < outlier_cut])

    mu = -1
    alpha = -1
    var = -1
    size = -1
    prob = -1

    exog = np.ones(len(x))
    if x.mean() > 0:
        start_param0 = np.log(x.mean())
        if p == 1:
            start_param1 = (x.var() - x.mean()) / x.mean()
        elif p == 2:
            start_param1 = (x.var() - x.mean()) / x.mean() / x.mean()

        res = sm.NegativeBinomialP(x, exog, p=p).fit(
            start_params=[start_param0, start_param1], disp=0
        )

        converged = not np.isnan(res.llr)
        if not converged:
            return (-2, -2, -2, -2, -2)

        mu = np.exp(res.params[0])
        alpha = res.params[1]

        # Compute "size, prob" parameters needed for SciPy nbinom
        # Conversion code based on: https://github.com/statsmodels/statsmodels/issues/106#issuecomment-43961704
        if p == 1:
            Q = 1
            var = mu + alpha * mu
        elif p == 2:
            Q = 0
            var = mu + alpha * mu * mu

        size = 1.0 / alpha * mu**Q
        prob = size / (size + mu)

    return (mu, alpha, var, size, prob)


def mlxp_neg_binom(counts, size, prob):
    """Compute -log10(pval) for negative binomial
    counts is DataFrame where each row is assumed to be drawn
    from a negative binomial with parameters size and prob.
    """

    mlxp = [
        nb_mlxp(counts.iloc[i].values, size[i], prob[i]) for i in range(counts.shape[0])
    ]
    mlxp = pd.DataFrame(data=mlxp, index=counts.index, columns=counts.columns)
    mlxp.replace(np.inf, -np.log10(np.finfo(np.float64).tiny), inplace=True)
    return mlxp
