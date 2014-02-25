# Copyright (C) 2008 Simone Leo - CRS4. All Rights Reserved.
# 
# Permission to use, copy, modify, and distribute this software and its
# documentation for educational, research, and not-for-profit purposes, without
# fee and without a signed licensing agreement, is hereby granted, provided
# that the above copyright notice, this paragraph and the following two
# paragraphs appear in all copies, modifications, and distributions. Contact
# CRS4, Parco Scientifico e Tecnologico, Edificio 1, 09010 PULA (CA - Italy),
# +3907092501 for commercial licensing opportunities.
# 
# IN NO EVENT SHALL CRS4 BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL,
# INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF
# THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF CRS4 HAS BEEN ADVISED
# OF THE POSSIBILITY OF SUCH DAMAGE.
# 
# CRS4 SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
# HEREUNDER IS PROVIDED "AS IS". CRS4 HAS NO OBLIGATION TO PROVIDE MAINTENANCE,
# SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

"""
Functions for computing the False Discovery Rate (FDR) of a multiple test
procedure. The FDR is the expected proportion of falsely rejected hypotheses.

@see:
  Yoav Benjamini and Yosef Hochberg, I{Controlling the False Discovery Rate:
  A Practical and Powerful Approach to Multiple Testing}. Journal of the Royal
  Statistical Society, Series B (Methodological), Vol. 57, No. 1 (1995), pp.
  289-300.

  Yekutieli, D. and Benjamini, Y., I{Resampling-based false discovery rate
  controlling multiple test procedures for correlated test statistics}. J. of
  Statistical Planning and Inference, 82, pp. 171-196, 1999.

  John D. Storey (2002), I{A direct approach to false discovery rates}. Journal
  of the Royal Statistical Society: Series B (Statistical Methodology) 64 (3),
  479-498.
"""

__docformat__ = 'epytext en'


import operator
import random


def bh_rejected(pv, threshold):
    """
    Return the list of rejected p-values from C{pv} at FDR level C{threshold}.

    The downside of this approach is that the FDR level must be chosen in
    advance - the L{bh_qvalues} function yields a p-value-like output instead.

    @type pv: list
    @param pv: p-values from a multiple statistical test
    @type threshold: float
    @param threshold: the level at which FDR rate should be controlled

    @rtype: list
    @return: p-values of rejected null hypotheses
    """
    if threshold < 0 or threshold > 1:
        raise ValueError("the threshold must be between 0 and 1")
    if not pv:
        return []
    pv = sorted(pv)
    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")
    m = len(pv)
    for i in xrange(m-1, -1, -1):
        if pv[i] <= float(i+1)*threshold/m:
            return pv[:i+1]
    return []


def bh_qvalues(pv):
    """
    Return Benjamini-Hochberg FDR q-values corresponding to p-values C{pv}.

    This function implements an algorithm equivalent to L{bh_rejected} but
    yields a list of 'adjusted p-values', allowing for rejection decisions
    based on any given threshold.

    @type pv: list
    @param pv: p-values from a multiple statistical test

    @rtype: list
    @return: adjusted p-values to be compared directly with the desired FDR
      level
    """
    if not pv:
        return []
    m = len(pv)
    args, pv = zip(*sorted(enumerate(pv), None, operator.itemgetter(1)))
    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")
    qvalues = m * [0]
    mincoeff = pv[-1]
    qvalues[args[-1]] = mincoeff
    for j in xrange(m-2, -1, -1):
        coeff = m*pv[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[args[j]] = mincoeff
    return qvalues


def get_pi0(pv, lambdas):
    """
    Compute Storey's C{pi0} from p-values C{pv} and C{lambda}.

    this function is equivalent to::
    
        m = len(pv)
        return [sum(p >= l for p in pv)/((1.0-l) * m) for l in lambdas]
        
    but the above is C{O(m*n)}, while it needs only be C{O(m+n)
    (n = len(lambdas))}

    @type pv: list
    @param pv: B{SORTED} p-values vector
    @type lambdas: list
    @param lambdas: B{SORTED} lambda values vector

    @rtype: list
    @return: estimated proportion of null hypotheses C{pi0} for each lambda
    """
    m = len(pv)
    i = m - 1
    pi0 = []
    for l in reversed(lambdas):
        while i >= 0 and pv[i] >= l:
            i -= 1
        pi0.append((m-i-1)/((1.0-l)*m))
    pi0.reverse()
    return pi0


def storey_qvalues(pv, l=None):
    """
    Return Storey FDR q-values corresponding to p-values C{pv}.

    The main difference between B-H's and Storey's q-values is that the latter
    are weighted by the estimated proportion C{pi0} of true null hypotheses.

    @type pv: list
    @param pv: p-values from a multiple statistical test
    @type l: float
    @param l: lambda value for C{pi0} (fraction of null p-values) estimation

    @rtype: list
    @return: storey q-values corresponding to C{pv}
    """
    if not pv:
        return []
    m = len(pv)
    args, pv = zip(*sorted(enumerate(pv), None, operator.itemgetter(1)))
    if pv[0] < 0 or pv[-1] > 1:
        raise ValueError("p-values must be between 0 and 1")

    if l is None:
        # optimal lambda/pi0 estimation
        lambdas = [i/100.0 for i in xrange(0, 91, 5)]
        n = len(lambdas)
        pi0 = get_pi0(pv, lambdas)
        min_pi0 = min(pi0)
        mse = [0] * n        
        for i in xrange(1, 101):
            # compute bootstrap sample with replacement
            pv_boot = [pv[int(random.random()*m)] for j in xrange(m)]
            pi0_boot = get_pi0(sorted(pv_boot), lambdas)
            for j in xrange(n):
                mse[j] += (pi0_boot[j] - min_pi0) * (pi0_boot[j] - min_pi0)
        min_mse = min(mse)
        argmin_mse = [i for i, mse_i in enumerate(mse) if mse_i == min_mse]
        pi0 = min(pi0[i] for i in argmin_mse)
        pi0 = min(pi0, 1)
    else:
        try:
            l = float(l)
        except ValueError:
            raise TypeError("lambda must be a number")
        if l < 0 or l >= 1:
            raise ValueError("lambda must be within [0,1)")
        pi0 = get_pi0(pv, [l])
        pi0 = min(pi0[0], 1)

    qvalues = m * [0]
    mincoeff = pi0 * pv[-1]
    qvalues[args[-1]] = mincoeff
    for j in xrange(m-2, -1, -1):
        coeff = pi0*m*pv[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[args[j]] = mincoeff
    return qvalues
