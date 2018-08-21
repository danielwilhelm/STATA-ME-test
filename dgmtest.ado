/*
	Nonparametric Testing for Significance and Presence of Measurement Error

21/08/2018

This program is developed for a testing methodology, in Delgado and Gonzalez-Manteiga (AoS, 2001),
for selecting explanatory variables in nonparametric regression.

	H0: E(Y|X,W,Z) = m(X,W1) + a*W2 a.s.,
	where	m(.) = E(Y|.),
			Y is a scalar, and
			X is R^qx valued, W1 is R^qw1 valued, W2 is R^qw2 valued, and Z is R^qz valued.

Syntax:
	dgmtest depvar expvar [if] [in]
		[, qz(integer) qw2(integer) teststat(string) kernel(string) bootdist(string) bw(real) bootnum(integer) ngrid(integer) qgrid(real)]
	where
	qz      : dimension of Z (default = 1)
	qw2     : dimension of additively linear control variables W2 (default = 0)
	teststat: Test statistic, options: CvM (default), KS
	kernel  : kernel function
	          options: biweight, epanechnikov (default), epan2 (DGM, 2001), epan4, normal, rectangle, triangular
	bw      : bandwidth of kernel function (default = 0 => bw = n^(-1/3q))
	bootdist: distribution with bounded support, zero mean, and unit variance
	          options: mammen (default), rademacher, and uniform
	bootnum : number of bootstrap samples (default = 500)
	ngrid   : number of equidistant grid points for each variable in KS statistic (default = 0 => we use exact observation values of W)
	qgrid   : quantile probability for maximum or minimum grid points in KS statistic(default = 0)

Outcome:
	e(N)      : number of observations
	e(dimXW1) : diminsion of (X,W1)
	e(dimZ)   : diminsion of Z
	e(dimW2)  : dimension of additively linear and correctly measured control variables W2
	e(stat)   : scalar value of the Cramer-von Mises statistic
	e(bootnum): number of bootstrap samples
	e(bw)     : bandwidth
	e(btpv)   : P[stat < statnstar]
	e(btcv1)  : bootstrap critical value at 1%
	e(btcv5)  : bootstrap critical value at 5%
	e(btcv10) : bootstrap critical value at 10%
	e(ngrid)  : number of grid points
	e(qgrid)  : quantile probability for min or max values of grid points

If unspecified, the command runs on a default setting.
In mata, W=(X,Z) with X=(X,W1).
*/

*************** MAIN DGMTEST CODE ********************************************
program define dgmtest, eclass
		version 12
		
		syntax varlist(min=2) [if] [in] [, qz(integer 1) qw2(integer 0) teststat(string) kernel(string) bw(real 0) bootdist(string) bootnum(integer 500) ngrid(integer 0) qgrid(real 0)]
		ereturn clear
		ereturn local cmd = "dgmtest"
		ereturn local title  = "Nonparametric Significance Test"
		marksample touse
				
		gettoken Y W : varlist
		local Y `Y'
		local W `W'
		
		preserve
		if `"`if'`in'"' != "" qui keep `if' `in'
		
		
		display "----------------------------------------------------- "
		display " Delgado and Manteiga test
		display	"----------------------------------------------------- "
		
		display " "
		if (`qw2' == 0) {
		display "H0: E[Y | X,W1,Z] = E[Y | X,W1]"
		}
		else {
		display "H0: E[Y | X,W,Z] = E[Y | X,W1] + a*W2 "
		}
		display " "
		
		display "----- parameter settings -----"
		display " "

		
	// Checking inputs
		if (`qz' < 1) {
		display "qz should be a positive integer"
		error 111
		}
		
		if (`qw2' < 0) {
		display "q2 should be a nonnegative integer"
		error 111
		}
		
		if ("`teststat'" == "") {
		display "Test statistic: CvM (default)"
		local teststat = "CvM"
		}
		else if (inlist("`teststat'","CvM","KS")) {
		display "Test statistic: `teststat'"
		}
		else {
		display "Choose a test statistic between Cramer-von Mises (CvM; default) and Kolmogorov-Smirnov (KS) statistics"
		error 197
		}
		
		if ("`kernel'" == "") {
		display "Kernel: epanechnikov (default)"
		local kernel = "epanechnikov"
		}
		else if (inlist("`kernel'","biweight","epanechnikov","epan2","epan4","normal","rectangle","triangular")) {
		display "Kernel: `kernel'"
		}
		else {
		display "Choose a kernel among biweight, epanechnikov (default), epan2, epan4, normal, rectangle, and triangular"
		error 197
		}
		
		if (`bw' == 0) {
		display "bw = n^(1/3q) (default)"
		}
		else if (`bw' < 0) {
		display "bandwidth should be greater than zero"
		error 111
		}
		
		if ("`bootdist'" == "") {
		local bootdist = "mammen"
		display "bootstrap multiplier distribution: mammen (default)"
		}
		else if (inlist("`bootdist'","mammen","uniform","rademacher")) {
		display "bootstrap multiplier distribution: `bootdist'"
		}
		else {
		display "Choose a distribution among mammen (default), rademacher and uniform"
		error 197
		}
		
		if (`bootnum' <= 0) {
		display "Choose a positive integer for our bootstrap sample"
		error 111
		}
		
		if (`ngrid' == 0 & "`teststat'" == "KS") {
		display "We use the observation values for KS statistic"
		}
		else if (`ngrid' < 0) {
		display "Choose a positive integer for the number of equi-distant grid points"
		error 111
		}

		if (`qgrid' < 0 | `qgrid' > 1) {
		display "The quantile for maximum or minimum grid points should be in [0,1]"
		error 111
		}

	// Running main program
		ereturn local teststat  = "`teststat'"
		ereturn local kernel 	= "`kernel'"
		ereturn local bootdist 	= "`bootdist'"
		
		mata: test("`Y'","`W'",`qz',`qw2',"`teststat'","`kernel'",`bw',"`bootdist'",`bootnum',`ngrid',`qgrid')
		
		display " "
		display as txt " number of observations: " as res e(N)
		display as txt " bandwidth: " as res e(bw)
		display as txt " dimension of (X,W1): " as res e(dimXW1)
		display as txt " dimension of W2: " as res e(dimW2)
		display as txt " dimension of Z: " as res e(dimZ)
		display as txt " number of bootstrap samples: " as res e(bootnum)
		
		display " "
		display "----- test results -----"
		display " "
				
		display as txt " `teststat' = " as res e(stat)
		display as txt " bootstrap critical value at 1%: " as res e(btcv1)
		display as txt " bootstrap critical value at 5%: " as res e(btcv5)
		display as txt " bootstrap critical value at 10%: " as res e(btcv10)
		display as txt " p(`teststat' < `teststat'*) = " as res e(btpv)
end

*************** Define a mata function computing test statistic ****************
mata:

// Significance test
void test(string scalar yname, string scalar wname, real scalar qz,
	  real scalar qw2, string scalar teststat, string scalar kernel,
	  real scalar bw, string scalar bootdist, real scalar bootnum,
	  real scalar ngrid, real scalar qgrid)
{	real matrix W, X
	real colvector Y, statst
	real scalar q, stat, pstat

	Y = st_data(., yname)
	W = st_data(., wname)
	q = cols(W)-qz-qw2
	if (hasmissing(Y)>0) {
	exit(_error(3351, "Y has missing values"))
	}
	else if (hasmissing(W)>0) {
	exit(_error(3351, "(X,W,Z) has missing values"))
	}
	else if (q < 1) {
	exit(_error(102, "dim(X,W1) should be a positive integer"))
	}
	X = W[|1,1 \ rows(Y),q|]
	
	if (bw == 0) {
	bw = rows(Y)^(-1/(3*q))
	}
	
	if (qw2 > 0) {
	W2  = W[|1,q+1 \ rows(Y),q+qw2|]
	eW2 = W2 - npreg(W2, X, kernel, bw)
	eY  = Y - npreg(Y, X, kernel, bw)
	Y   = Y - invsym(eW2'*eW2)*(eW2'*eY)*W2
	W   = X,W[|1,q+qw2+1 \ rows(Y),cols(W)|]
	}
	
	stat	= mstat(Y, W, X, teststat, kernel, bw, ngrid, qgrid)
	statst  = mstat(btrs(Y, X, kernel, bw, bootdist, bootnum), W, X, teststat, kernel, bw, ngrid, qgrid)
	
	pstat	= sum(statst :> stat)/bootnum
	st_numscalar("e(N)", rows(Y))
	st_numscalar("e(dimXW1)", q)
	st_numscalar("e(dimZ)", cols(W)-q)
	st_numscalar("e(dimW2)", qw2)
	st_numscalar("e(stat)", stat)
	st_numscalar("e(btpv)", pstat)
	st_numscalar("e(btcv1)", qtile(statst,0.99))
	st_numscalar("e(btcv5)", qtile(statst,0.95))
	st_numscalar("e(btcv10)", qtile(statst,0.9))
	st_numscalar("e(bw)", bw)
	st_numscalar("e(bootnum)", bootnum)
	st_numscalar("e(ngrid)", ngrid)
	st_numscalar("e(qgrid)", qgrid)
}

// Test statistic for conditional mean and its bootstrap versions
real colvector mstat(real matrix Y, real matrix W, real matrix X,
                     string scalar teststat, string scalar kernel,
                     real scalar bw, real scalar ngrid, real scalar qgrid)
{	real matrix K, dX, Ti, W1, dW, Tn
	real colvector iota, statV
	real scalar n, i
	
	n = rows(Y)
	K = J(n,n,1)
	iota = J(n,1,1)
	for (i = 1; i <= cols(X); i++) {
		dX = X[|1,i \ n,i|]*iota' - iota*X[|1,i \ n,i|]'
		K  = K :* Kij(dX:/bw, kernel)
		}
	Ti = rowsum(K):*Y - K*Y
	if (teststat == "CvM") {
	W1 = J(n,n,1)
	for (i = 1; i <= cols(W); i++) {
		dW = W[|1,i \ n,i|]*iota' - iota*W[|1,i \ n,i|]'
		W1 = W1 :* (dW:<=0)
		}
	Tn    = (Ti'*W1) :/ (n^2 * bw^cols(X))
	statV = rowsum(Tn:^2)
	}
	else if (teststat == "KS") {
	if (ngrid == 0) {
	W1 = iota
	for (i = 1; i <= cols(W); i++) {
		dW = W[|1,i \ n,i|]*iota' - iota*W[|1,i \ n,i|]'
		dW = uniqrows((dW:<=0)')'
		W1 = (W1#J(1,cols(dW),1)) :* (J(1,cols(W1),1)#dW)
		W1 = uniqrows(W1')'
		}
	}
	else {
	W1 = iota
	for (i = 1; i <= cols(W); i++) {
		dW = W[|1,i \ n,i|]*J(ngrid,1,1)' - iota*rangen(qtile(W[|1,i \ n,i|],qgrid),qtile(W[|1,i \ n,i|],1-qgrid),ngrid)'
		dW = uniqrows((dW:<=0)')'
		W1 = (W1#J(1,cols(dW),1)) :* (J(1,cols(W1),1)#dW)
		W1 = uniqrows(W1')'
		}
	}
	Tn    = (Ti'*W1) :/ (n^2 * bw^cols(X))
	statV = rowmax(abs((n^(1/2)):*Tn))
	}
	return(statV)
}

// Bootstrap resample
real matrix btrs(real matrix Y, real matrix X, string scalar kernel,
                 real scalar bw, string scalar bootdist, real scalar bootnum)
{	real matrix mX, e, V, Ystar
	real colvector iota
	
	mX= npreg(Y, X, kernel, bw)
	e = Y - mX
	V = uniform(rows(Y),cols(Y)*bootnum)
	if (bootdist == "mammen") {
	// P[V=(1-sqrt(5))/2] = (sqrt(5)+1)/(2*sqrt(5)) and P[V=(1+sqrt(5))/2] = (sqrt(5)-1)/(2*sqrt(5))
	V = (V :> ((sqrt(5)+1)/(2*sqrt(5)))) :* sqrt(5) :+ ((1-sqrt(5))/2)
	}
	else if (bootdist == "rademacher") {
	// Rademacher: P(V=-1) = P(V=1) = 1/2
	V = (V :> 1/2) :* 2 :- 1
	}
	else if (bootdist == "uniform") {
	// Uniform(-sqrt(3),sqrt(3))
	V = (V :- 1/2) :* 2 :* sqrt(3)
	}
	iota  = J(bootnum,1,1)
	Ystar = (iota'#mX) :+ ((iota'#e):*V)
	return(Ystar)
}

// Nonparametric regression
real matrix npreg(real matrix Y, real matrix X, string scalar kernel,
                  real scalar bw)
{	real matrix K, dX, mX
	real colvector iota, fX
	real scalar n, i
	
	n = rows(Y)
	K = J(n,n,1)
	iota = J(n,1,1)
	for (i = 1; i <= cols(X); i++) {
		dX = X[|1,i \ n,i|]*iota' - iota*X[|1,i \ n,i|]'
		K  = K :* Kij(dX:/bw, kernel)
		}
	fX = rowsum(K)
	mX = (K*Y) :/ fX
	return(mX)
}

// Kernels
real matrix Kij(real matrix Z, string scalar kernel)
{	real matrix K

	if (kernel == "biweight") {
	// BIWEIGHT kernel
	K = (1 :- (Z:^2)) :^ 2 :* (15/16)
	}
	else if (kernel == "epanechnikov") {
	// EPANECHNIKOV's asymptotically optimal kernel
	K = (abs(Z) :< sqrt(5)) :* 0.75 :* (1 :- (0.2:*(Z:^2))) :/ sqrt(5)
	}
	else if (kernel == "epan2") {
	// EPANECHNIKOV's kernel with support [-1,1]
	K = (abs(Z) :< 	1) :* 0.75 :* (1 :- (Z:^2))
	}
	else if (kernel == "epan4") {
	// 4th order EPANECHNIKOV's kernel
	K = (abs(Z) :< 	1) :* 0.75 :* (1 :- (Z:^2)) :* (15 :- (35:*(Z:^2))) :/ 8
	}
	else if (kernel == "normal") {
	// GAUSSIAN density kernel
	K = exp(-0.5 :* (Z:^2)) :/ sqrt(2*pi())
	}
	else if (kernel == "rectangle") {
	// RECTANGLE kernel
	K = (abs(Z):<1) :* (1/2)
	}
	else if (kernel == "triangular") {
	// TRIANGLE kernel
	K = (abs(Z):<1) :* (1 :- abs(Z))
	}
	return(K)
}

// Quantile
real scalar qtile(real colvector Z, real scalar qq)
{	real colvector Zsort, qind
	real scalar n, Zq

	n = rows(Z)
	Zsort = sort(Z,1)
	qind = (rowmax((floor(n*qq + 0.5),1)) \ rowmin((ceil(n*qq + 0.5),n)))
	Zq = Zsort[qind[1,1],1] + (n*qq + 0.5 - qind[1,1]) * (Zsort[qind[2,1],1] - Zsort[qind[1,1],1])
	return(Zq)
}
end
