# robusTest

Implementation of corrected two sample tests for Pearson, Kendall and Spearman correlation tests, corrected Mann-Whitney (Wilcoxon) rank sum test, corrected Wilcoxon signed rank test. A robust variance test, a median test and a test for independence between two continuous variables based on the Kolmogorov-Smirnov's distance are also implemented. All the tests are asymptotically calibrated meaning that the probability of rejection under the null hypothesis is asymptotically equal to the level of the test. 

The package contains the functions:

- `cortest` that implements the corrected Pearson, Kendall's and Spearman's tests. As compared to the original tests in cor.test which all assume independence between the variables under the null hypothesis, the corrected tests assume that the correlation (of the different types) is equal to 0 under the null.
- `wilcoxtest` that implements the corrected version of the Wilcoxon test for two independent samples (which is equivalent to Mann-Whitney in that case) and for two paired samples. In the two independent sample case, as compared to the original test in wilcox.test which assumes that the variables have the same distribution under the null hypothesis, the corrected test assumes under the null that the probability that one variable exceeds the other is equal to 0.5.
- `vartest` that implements a variance test based on the Welch correction for the variables (x_i-mean(x))^2 and (y_i-mean(y))^2. As compared to the original test in var.test which only works under the Gaussian scenario, the corrected test works for any distribution of the two variables as long as the fourth order moments exist for both variables. The test can be applied to more than two groups. In the two sample case, the function also returns the confidence interval for the difference of the two variances.
- `mediantest` that tests if the median of the random variable (one sample case) or of the difference between two random variables (two sample case) is equal to 0.  The test is based on  asymptotic results on the rank statistics for the uniform distribution. Confidence intervals for the median (one sample case) or for the difference between two random variables (two sample case) are returned.
- `indeptest`that tests the independence between two continuous variables. The test is based on the maximum distance between the joint empirical cumulative distribution function and the product of the marginals. The distribution of this test has been numerically obtained, the test is exact for all n<=150 and approximated for n>150.
- `tiebreak` which randomly breaks ties in vectors, either inside the vector or between two vectors.

The dataset `Evans` can also be loaded from the **robusTest** package, this dataset was originally provided in the **lbreg** package.

