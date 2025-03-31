LaTeX Examples in reStructuredText
================================

This page demonstrates how to use LaTeX math in reStructuredText documents for ReadTheDocs.

Inline Equations
---------------

You can use inline math with the :math: role. For example:

The negative binomial probability mass function is given by :math:`P(X=k) = \binom{k+r-1}{k} p^r (1-p)^k` for :math:`k = 0, 1, 2, \ldots`

The shape parameter :math:`r` controls the dispersion, with higher values leading to less overdispersion.

The mean of the distribution is :math:`\mu = \frac{r(1-p)}{p}` and the variance is :math:`\sigma^2 = \frac{r(1-p)}{p^2}`.

Display Equations
----------------

For display math, use the math directive:

.. math::

   P(X=k) = \binom{k+r-1}{k} p^r (1-p)^k

Or alternatively, for aligned equations:

.. math::
   :nowrap:

   \begin{align}
   P(X=k) &= \binom{k+r-1}{k} p^r (1-p)^k\\
   &= \frac{\Gamma(k+r)}{k!\Gamma(r)}p^r(1-p)^k
   \end{align}

The variance can be expressed as:

.. math::

   \text{Var}(X) = \mu + \frac{\mu^2}{r}

For the negative binomial distribution, the dispersion parameter is:

.. math::

   \phi = 1 + \frac{\mu}{r}

Equation Numbering
-----------------

You can also number equations:

.. math::
   :label: nb_pmf

   P(X=k) = \binom{k+r-1}{k} p^r (1-p)^k

You can reference equation :eq:`nb_pmf` in your text.

Tips for Complex Equations
------------------------

For complex equations or if you're experiencing rendering issues:

1. Break equations into smaller parts
2. Use simpler LaTeX commands when possible
3. Ensure proper escaping of special characters
4. Test your equations in a LaTeX editor before including them in documentation 