# STATA-ME-test 

Authors : Young Jun Lee and Daniel Wilhelm

This project provides the STATA command `dgmtest` which implements the test for significance by Delgado and Manteiga (2001) and can be used to test for the presence of measurement error as described in Wilhelm (2017).

Files contained in this package:

- The file `dgmtest.ado` contains the `dgmtest` command.
- The file `simul_dgmtest.ado` contains code to replicate the simulations in Delgado and Manteiga (2001).
- The file `example.do` contains the simple simulation example shown below.


## Installation
1. Download the package.
2. Change into the directory containing this package.
3. Use the command `dgmtest` as described below.

## Syntax
The command `dgmtest` tests the null hypothesis

```
H0:   E[Y | X,Z] = E[Y | X]
```

or

```
H0:   E[Y | X,XL,Z] = E[Y | X,XL]
```

against the alternative that the null does not hold, where

- Y is a scalar dependent variable (`depvar`)
- X is a vector of explanatory variables (`expvar`)
- Z is a vector of explanatory variables (`expvar`)
- XL is a vector of additively linear control variables (`expvar`)

Wilhelm (2017) shows that, under certain conditions, the null hypothesis H0 is equivalent to the null hypothesis of no measurement error in X. Most importantly, Z has to be an excluded variable for the outcome equation of Y. If we have additively linear control variables, we apply the test in Delgado and Manteiga (2001) after replacing Y with (Y - ahat * XL) where ahat is the Robinson (1988)'s estimator.

Syntax:

```
dgmtest depvar expvar [if] [in] [, q(integer) ql(integer) teststat(string) kernel(string) bootdist(string) bw(real) bootnum(integer) ngrid(integer) qgrid(real)]
```

where

- `depvar` is the outcome variable Y
- `expvar` is a list of variables containing all elements of X, Z, and XL (oder of `expvar`: X Z XL)

The options are as follows:

- `q` is the dimension of X (default = 1).
- `ql` is the dimension of XL (default = 0).
- `teststat` is the type of test statistic to be used: Cramer-van Mises (CvM, default) or Kolmogorov-Smirnov (KS).
- `kernel` is the kernel function: biweight, epanechnikov (default), epan2, epan4, normal, rectangle, triangular.
- `bw` is the bandwidth (default = n^(-1/3q), rule of thumb, where n is the sample size and q the dimension of X).
- `bootnum` is the number of bootstrap samples for the computation of the test's critical value (default = 500).
- `bootdist` is the distribution of the bootstrap multiplier variable: mammen (default), rademacher, uniform.
- `ngrid` is the number of equally spaced grid points used to compute the supremum of the KS statistic, if that statistic is chosen via the option `teststat`. The default is 0 which means that the sample serves as the grid.
- `qgrid` is a number between 0 and 1 to define the min and max values of the grid in the previous option. The min value is the `qgrid`-quantile and the max value is the (1-`qgrid`)-quantile. The default is 0 so that in that case the grid ranges from the min to the max value in the sample.

If options are left unspecified, the command runs on the default settings.


## Examples

Testing for measurement error in simulated data, using the default options:
```
set obs 200

// true regressor
generate Xstar = runiform()

// measurement error in X
generate etaX = runiform()

// mismeasured regressor
generate X = Xstar + 0.5*etaX

// additively linear control variable
generate XL = runiform()

// measurement error in Z
generate etaZ = runiform()

// second measurement of true regressor
generate Z = Xstar + 0.5*etaZ

// regression error
generate epsilon = runiform()

// outcome equation
generate Y1 = Xstar^2 + 0.2*Xstar + 0.5*epsilon
generate Y2 = 0.5*XL + Xstar^2 + 0.2*Xstar + 0.5*epsilon

// perform the test of the hypothesis of no measurement error in X
dgmtest Y1 X Z, kernel(epan2)
dgmtest Y2 X Z XL, ql(1) kernel(epan2)
```



# References
[Robinson, P. M. (1988), "Root-N-Consistent Semiparametric Regression", Econometrica, 56(4), p. 931-954](http://www.jstor.org/stable/1912705)

[Wilhelm, D. (2017), "Testing for the Presence of Measurement Error", working paper available soon](http://www.ucl.ac.uk/~uctpdwi)

[Delgado, M. and Manteiga, W. (2001), "Significance Testing in Nonparametric Regression Based on the Bootstrap", Annals of Statistics, 29(5), p. 1469-1507](http://www.jstor.org/stable/2699997)
