{smcl}
{* *! version 1.1.0 22aug2018}{...}
{cmd:help dgmtest} 
{hline}

{title:Title}

{p2colset 5 19 21 2}{...}
{p2col :{bf:dgmtest} {hline 2}}Nonparametric Testing for Significance of Regressors and for the Presence of Measurement Error{p_end}
{p2colreset}{...}


{title:Syntax}

{phang}
Nonparametric Testing for Significance of Regressors and for the Presence of Measurement Error

{p 8 16 2}{cmd:dgmtest} {depvar} [{expvars}] {ifin}
[, {it:{help fvw13##dgmoptions:DGM_options}}]

{marker dgmoptions}{...}
{synoptset 20 tabbed}{...}
{synopthdr :DGM_options}
{synoptline}
{syntab:Number of Variables}
{synopt :{opt qz(integer)}}dimension of the instrumental variables which are not
significant under the null hypothsis; default is {cmd:qz(1)}.{p_end}
{synopt :{opt qw2(integer)}}dimension of additively linear control variables;
default is {cmd:qw2(0)}.{p_end}

{syntab:Type of Test}
{synopt :{opt teststat(string)}}specify test statistic to be used among
Cramer-von Mises' (CvM) and Kolmogorov-Smirnov (KS) statistics; default is
{cmd:teststat(CvM)}.{p_end}
{synopt :{opt kernel(string)}}specify kernel function; default is
{cmd:kernel(epanechnikov)}.{p_end}
{synopt :{opt bootdist(string)}}specify distribution of the bootstrap multiplier
variables; default is {cmd:bootdist(mammen)}.{p_end}

{syntab:Options for Kernel Method and Bootstrap}
{synopt :{opt bw(real)}}half-width of kernel; default is n^(-1/3q), where
'n' is the number of observations and 'q' is the dimension of explanatory
variables which are nonparametrically significant under the null hypothesis.{p_end}
{synopt :{opt bootnum(integer)}}number of bootstrap samples for the computation
of the test's critical value; default is {cmd:bootnum(500)}.{p_end}

{syntab:Suboptions for Kolmogorov-Smirnov Statistic}
{synopt :{opt ngrid(integer)}}number of equally spaced grid points used to
compute the supremum of the Kolmogorov-Smirnov statistic, if that statistic is
chosen via the option teststat.{p_end}
{synopt :{opt qgrid(real)}}quantile probability between 0 and 1 to set the
min and max values of the grid points.{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}
{cmd:dgmtest} performs the nonparametric test by Delgado and Gonzalez Manteiga (2001) for the significance of regressors. It can also be used as a test for the presence of measurement error as described in Wilhelm (2018) and Lee and Wilhelm (2018).

{pstd}
The null hypothesis is H_0: E[Y|X,W,Z] = E[Y|X,W], i.e. that the variable Z is not significant for the conditional mean of the outcome Y ({it:depvar}) given X, W, and Z. {it:expvars} contains the list of variables X, W, Z. W can contain two types of explanatory variables, W=(W1,W2) where W1 is included in the conditional mean in a nonparametric fashion and W2 is included as additively separable. Therefore, {it:expvars} should contain at the very least two variables (X and Z) and possibly more (W), where the additional controls W are listed between X and Z.

{pstd}
Wilhelm (2018) shows that the null hypothesis can also be interpreted as the null of no measurement error in the explanatory variable X. In this case, Z is an instrument for X or a second measurement. More details can be found in Wilhelm (2018) and Lee and Wilhelm (2018).

{pstd}
{cmd:dgmtest} displays a test statistic and its bootstrap p-value with the
bootstrap critical values.


{title:Options for DGMTEST}

{dlgtab:Number of Variables}

{phang}
{opt qz(integer)} specifies the dimension of the instrumental variables Z which
are not significant in the null hypothesis. These are the last {cmd:{opt qz}}
variables listed in {it:expvars}.

{phang}
{opt qw2(integer)} specifies the dimension of additively linear control
variables W2. The default is 0, which means that there is no additively linear
control variable. If it is positive, given a root-n-consistent estimator for
those variables as proposed by Robinson (1988), we evaluate the test
statistic. Let M be the number of variables in {it:expvars}. Then, first set of
(M-{cmd:{opt qz}}-{cmd:{opt qw2}}) variables are the explanatory variables which
are nonparametrically significant under the null hypothesis. Second set of
{cmd:{opt qw2}} variables are additively linear control variables.

{dlgtab:Type of Test}

{phang}
{opt teststat(string)} specifies the type of test statistic. {cmd:CvM} and
{cmd:KS} represent the Cramer-von Mises' and Kolmogorov-Smirnov statistics,
respectively. The default is {cmd:CvM}.

{phang}
{opt kernel(string)} specifies the kernel function for use in calculating the
test statistic. The default kernel is the Epanechnikov kernel
({cmd:epanechnikov}).

{synoptset 20 tabbed}{...}
{synoptline}
{synopthdr :kernel}
{synoptline}
{synopt :{opt epanechnikov}} Epanechnikov kernel; the default{p_end}
{synopt :{opt epan2}} Epanechnikov kernel of order 2{p_end}
{synopt :{opt epan4}} Epanechnikov kernel of order 4{p_end}
{synopt :{opt biweight}} biweight kernel function{p_end}
{synopt :{opt normal}} Gaussian kernel function{p_end}
{synopt :{opt rectangle}} rectangle kernel function{p_end}
{synopt :{opt triangular}} triangle kernel function{p_end}
{synoptline}

{phang}
{opt bootdist(string)} specifies the distribution of the bootstrap multiplier
variables. Following Delgado and Gonzalez Manteiga (2001), it should have a zero mean and unit
variance. The default is {cmd:mammen}.

{synoptset 20 tabbed}{...}
{synoptline}
{synopthdr :bootdist}
{synoptline}
{synopt :{opt mammen}} Two point distribution attaching masses
(sqrt(5)+1)/2sqrt(5) and (sqrt(5)-1)/2sqrt(5) to the points -(sqrt(5)-1)/2 and
(sqrt(5)+1)/2, respectively{p_end}
{synopt :{opt rademacher}} Two point distribution attaching masses 1/2 to {-1,1}
{p_end}
{synopt :{opt uniform}} continuous uniform distribution on (-sqrt(3),sqrt(3))
{p_end}
{synoptline}

{dlgtab:Options for Kernel Method and Bootstrap}

{phang}
{opt bw(real)} specifies the half-width of kernel. The default is
n^(-1/3q), a rule of thumb in Delgado and Gonzalez Manteiga (2001), where 'n'
is the number of observations and 'q' is the dimension of explanatory variables
which are nonparametrically significant under the null hypothesis.

{phang}
{opt bootnum(integer)} specifies the number of bootstrap samples for the
computation of the test's critical values. We generate {cmd:{opt bootnum}}
bootstrap residual samples and compute their corresponding bootstrap statistics.
Then, we approximate the bootstrap critical value by using them.

{dlgtab:Suboptions for Kolmogorov-Smirnov Statistic}

{phang}
{opt ngrid(integer)} specifies the number of equally spaced grid points used to
compute the supremum of the Kolmogorov-Smirnov statistic, if that statistic is
chosen via the option {cmd:KS}. The default is 0, which means that the sample
serves as the grid. It is required for calculating the exact Kolmogorov-Smirnov
statistic, but it is a burden when we perform a simulation with a large sample.

{phang}
{opt qgrid(real)} specifies the quantile probability between 0 and 1 to set the
minimum or maximum values of the grid points in the previous option. If
{cmd:{opt qgrid}} is smaller than 0.5, the minimum value is the
{cmd:{opt qgrid}}-quantile and the max value is the
(1-{cmd:{opt qgrid}})-quantile. The default is 0, so that in that case the grid
ranges from the min to the max value in the sample.


{title:Examples}

{pstd}Setup{p_end}
{phang2}{cmd:. webuse 401k}{p_end}

{pstd}Cramer-von Mises' statistic with a kernel function {cmd:{opt epan2}}{p_end}
{phang2}{cmd:. dgmtest prate mrate ltotemp, kernel(epan2)}{p_end}
{phang2}{cmd:. dgmtest prate mrate totemp, kernel(epan2)}{p_end}

{pstd}Cramer-von Mises' statistic with a kernel function {cmd:{opt epan2}}{p_end}
{phang2}{cmd:. dgmtest prate mrate ltotemp age, kernel(epan2)}{p_end}

{pstd}Cramer-von Mises' statistic with additional linear control variable and a
kernel function {cmd:{opt epan2}}{p_end}
{phang2}{cmd:. dgmtest prate ltotemp mrate age, qw2(1) kernel(epan2)}{p_end}

{pstd}Kolmogorov-Smirnov statistic with a kernel function {cmd:{opt epan2}}{p_end}
{phang2}{cmd:. dgmtest prate mrate ltotemp, teststat(KS) kernel(epan2)}{p_end}
{phang2} In this case, there is an error message because the number of
observations is large. We need to use grid points with an option
{cmd:{opt ngrid}}.{p_end}

{pstd}Kolmogorov-Smirnov statistic with a kernel function {cmd:{opt epan2}} and
100 grid points{p_end}
{phang2}{cmd:. dgmtest prate mrate ltotemp, teststat(KS) kernel(epan2) ngrid(100)}{p_end}


{title:Saved results}

{pstd}
{cmd:dgmtest} saves the following in {cmd:e()}:

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(dimXW1)}}dimension of explanatory variables which are nonparametrically significant under the null hypothesis{p_end}
{synopt:{cmd:e(dimW2)}}dimension of additively linear control variables{p_end}
{synopt:{cmd:e(dimZ)}}dimension of explanatory variables which are not significant under the null hypothesis{p_end}
{synopt:{cmd:e(stat)}}test statiatic value{p_end}
{synopt:{cmd:e(btcv1)}}1% bootstrap critical value{p_end}
{synopt:{cmd:e(btcv5)}}5% bootstrap critical value{p_end}
{synopt:{cmd:e(btcv10)}}10% bootstrap critical value{p_end}
{synopt:{cmd:e(btpv)}}bootstrap p-value{p_end}
{synopt:{cmd:e(bw)}}bandwidth for kernel method{p_end}
{synopt:{cmd:e(bootnum)}}number of bootstrap samples{p_end}
{synopt:{cmd:e(ngrid)}}number of grid points{p_end}
{synopt:{cmd:e(qgrid)}}quantile probability for min or max value of grid points{p_end}


{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:dgmtest}{p_end}
{synopt:{cmd:e(title)}}{cmd:Nonparametric Significance Testing}{p_end}
{synopt:{cmd:e(teststat)}}type of test statistic{p_end}
{synopt:{cmd:e(kernel)}}type of kernel function{p_end}
{synopt:{cmd:e(bootdist)}}distribution of bootstrap multiplier variable{p_end}


{title:Reference}

{phang}
Delgado, M. and Gonzalez-Manteiga, W., 2001, Significance Testing in Nonparametric Regression Based on the Bootstrap. Annals of Statistics, 29(5), 1469-1507

{phang}
Lee, Y. J. and Wilhelm, D., 2018, Testing for the Presence of Measurement Error in Stata

{phang}
Robinson, P. M., 1988, Root-N-Consistent Semiparametric Regression, Econometrica, 56(4), 931-954

{phang}
Wilhelm, D., 2018, Testing for the Presence of Measurement Error, CeMMAP Working Paper CWP45/18


{title:Remarks}

{p 4 4}This is a first and preliminary version. Please feel free to share your comments, reports of bugs and
propositions for extensions.

{p 4 4}If you use this command in your work, please cite Delgado and Gonzalez Manteiga (2001) when testing for significance of regressors and Wilhelm (2018) when testing for the presence of measurement error.


{title:Disclaimer}

{p 4 4 2}THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT WARRANTY OF ANY KIND, EITHER EXPRESSED 
OR IMPLIED. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. 
SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.

{p 4 4 2}IN NO EVENT WILL THE COPYRIGHT HOLDERS OR THEIR EMPLOYERS, OR ANY OTHER PARTY WHO
MAY MODIFY AND/OR REDISTRIBUTE THIS SOFTWARE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY 
GENERAL, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM.


{title:Authors}

{p 4 6}Young Jun Lee and Daniel Wilhelm{p_end}
{p 4 6}University College London{p_end}
{p 4 6}young.lee.13@ucl.ac.uk and d.wilhelm@ucl.ac.uk{p_end}
