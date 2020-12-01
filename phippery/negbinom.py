# Copyright 2017 Uri Laserson
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
import pandas as pd
import scipy.stats as st
import statsmodels.api as sm

def neg_binom_pvals(counts, size, prob):

    if size < 0:
        # no model for this peptide due to zero counts in mock-IPs,
        # assign p-value = 1
        pvals = np.ones(len(counts))
    else:
        pvals = 1. - st.nbinom.cdf(counts, size, prob)
    
    return pvals


def fit_neg_binom(x, p):
    """Fit a Type-I (p=1) or Type-II (p=2) negative binominal to each peptide across all samples
    Uses the "mu, alpha" parametrization in statsmodels
    """

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

        res = sm.NegativeBinomialP(x, exog, p=p).fit(start_params=[start_param0, start_param1], disp=0) 
        mu = np.exp(res.params[0])
        alpha = res.params[1]
        
        # Compute "size, prob" parameters needed for SciPy nbinom 
        # Conversion code based from: https://github.com/statsmodels/statsmodels/issues/106#issuecomment-43961704
        if p == 1:
            Q = 1
            var = mu + alpha * mu
        elif p == 2:
            Q = 0
            var = mu + alpha * mu * mu

        size = 1. / alpha * mu**Q
        prob = size / (size + mu)

    return (mu, alpha, var, size, prob)


def mlxp_neg_binom(counts, size, prob):
    """Compute -log10(pval) for negative binomial
    counts is DataFrame where each row is assumed to be drawn 
    from a negative binomial with parameters size and prob.
    """

    mlxp = [
        -np.log10(neg_binom_pvals(counts.iloc[i].values, size[i], prob[i]))
        for i in range(counts.shape[0])
    ]
    mlxp = pd.DataFrame(data=mlxp, index=counts.index, columns=counts.columns)
    mlxp.replace(np.inf, -np.log10(np.finfo(np.float64).tiny), inplace=True)
    return mlxp
