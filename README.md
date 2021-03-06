# STATA-ME-test 

Authors : Young Jun Lee and Daniel Wilhelm

This project provides the STATA command `dgmtest` which implements the test for significance by Delgado and Manteiga (2001) and can be used to test for the presence of measurement error as described in Wilhelm (2018) and Lee and Wilhelm (2018).

Files contained in this package:

- The file `dgmtest.ado` contains the `dgmtest` command.
- The file `dgmtest.sthlp` contains the Stata helpfile for the `dgmtest` command.
- The files `example_DGM2001.ado` and `simul_DGM2001.do` contain the code to replicate the simulations in Delgado and Manteiga (2001).
- The files `example_Wilhelm2018.ado` and `simul_Wilhelm2018.do` contain the code to replicate the simulations in Wilhelm (2018).
- The file `example.do` contains the simple simulation example shown below.


## Installation
1. Download the package.
2. Change into the directory containing this package.
3. Use the command `dgmtest` as described below.

## Syntax
The command `dgmtest` tests the null hypothesis

```
H0:   E[Y | X, W, Z] = E[Y | X, W]
```

against the alternative that the null does not hold, where

- Y is a scalar dependent variable
- X and W are vectors of explanatory variables
- Z is a vector of explanatory variables

The vector of explanatory variables, W, may contain elements that enter the conditional expectation in a linear, additively separable fashion. For example, decompose W=(W1,W2) where W1 enters nonseparably and W2 enters in a linear, additively separable fashion,

```
E[Y | X, W, Z] = f(X,W1,Z) + pi*W2
```

where f is some function and pi a row-vector of the same dimension as W2. In the presence of variables W2, we apply the test in Delgado and Manteiga (2001) after replacing Y with (Y - pihat*W2), where pihat is Robinson (1988)'s estimator of pi.

Syntax:

```
dgmtest depvar expvar [if] [in] [, qz(integer) qw2(integer) teststat(string) kernel(string) bootdist(string) bw(real) bootnum(integer) ngrid(integer) qgrid(real)]
```

where

- `depvar` is the outcome variable Y
- `expvar` is a list of variables containing all elements of X, W, and Z. The order of variables in the list should be: X, W, Z)

The options are as follows:

- `qz` is the dimension of Z (default = 1).
- `qw2` is the dimension of W2 (default = 0).
- `teststat` is the type of test statistic to be used: Cramer-van Mises (CvM, default) or Kolmogorov-Smirnov (KS).
- `kernel` is the kernel function: biweight, epanechnikov (default), epan2, epan4, normal, rectangle, triangular.
- `bw` is the bandwidth (default = n^(-1/3q), rule of thumb, where n is the sample size and q the dimension of X1).
- `bootnum` is the number of bootstrap samples for the computation of the test's critical value (default = 500).
- `bootdist` is the distribution of the bootstrap multiplier variable: mammen (default), rademacher, uniform.
- `ngrid` is the number of equally spaced grid points used to compute the supremum of the KS statistic, if that statistic is chosen via the option `teststat`. The default is 0 which means that the sample serves as the grid.
- `qgrid` is a number between 0 and 1 to define the min and max values of the grid in the previous option. The min value is the `qgrid`-quantile and the max value is the (1-`qgrid`)-quantile. The default is 0 so that in that case the grid ranges from the min to the max value in the sample.

If options are left unspecified, the command runs on the default settings.


## Testing for the presence of measurement error

Wilhelm (2018) shows that, under some conditions, the null hypothesis H0 is equivalent to the hypothesis of no measurement error in X. In this context, the variable Z must be excluded from the outcome equation. For example, it could be a second measurement or an instrumental variable. See Wilhelm (2018), Lee and Wilhelm (2018), and the examples below for more details.



## Examples

### Generate explanatory variables


```
set obs 200

// true regressor
generate Xstar = runiform()

// measurement error in X
generate etaX = runiform()

// mismeasured regressor
generate X1 = Xstar + 0.5*etaX

// additively linear control variable
generate X2 = runiform()

// measurement error in Z
generate etaZ = runiform()

// second measurement of true regressor
generate Z = Xstar + 0.5*etaZ

// regression error
generate epsilon = runiform()
```


### Generate outcome variable

We generate an outcome in two different ways, in a regression with and without additively separable, linear controls:

```
// outcome equation without controls
generate Y1 = Xstar^2 + 0.2*Xstar + 0.5*epsilon

// outcome equation with controls
generate Y2 = Xstar^2 + 0.2*Xstar + 0.5*X2 + 0.5*epsilon
```


### Perform the test of no measurement error

Perform the test using default options:

```
// perform the test of the hypothesis of no measurement error in X1
dgmtest Y1 X1 Z
dgmtest Y2 X1 X2 Z, qw2(1)
```

Perform the test, choosing the triangular kernel function:

```
// perform the test of the hypothesis of no measurement error in X1
dgmtest Y1 X1 Z, kernel(triangular)
dgmtest Y2 X1 X2 Z, qw2(1) kernel(triangular)
```



# References
[Delgado, M. and Manteiga, W. (2001), "Significance Testing in Nonparametric Regression Based on the Bootstrap", Annals of Statistics, 29(5), p. 1469-1507](http://www.jstor.org/stable/2699997)

[Robinson, P. M. (1988), "Root-N-Consistent Semiparametric Regression", Econometrica, 56(4), p. 931-954](http://www.jstor.org/stable/1912705)

[Wilhelm, D. (2018), "Testing for the Presence of Measurement Error", CeMMAP Working Paper CWP45/18](http://www.ucl.ac.uk/~uctpdwi/papers/cwp451818.pdf)

[Lee, Y. and Wilhelm, D. (2018), "Testing for the Presence of Measurement Error in Stata", working paper available soon](http://www.ucl.ac.uk/~uctpdwi)
