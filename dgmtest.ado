/*
	Significance Testing in Nonparametric Regression Based on the Bootstrap

05/09/2017

This program is developed for a testing methodology, in Delgado and Gonzalez-Manteiga (AoS, 2001),
for selecting explanatory variables in nonparametric regression.

	H0: E(Y|W) = m(X) a.s.,
	where	m(.) = E(Y|X=.),
			Y is a scalar, and
			W = (X,Z), X is R^q valued and Z is R^p valued.

Syntax:
	dgmtest depvar expvar [if] [in]
		[, q(integer) teststat(string) kernel(string) bootdist(string) cbw(real) bootnum(integer)]
	where
	q       : dimension of X (default = 1)
	teststat: Test statistic, options: CvM (default), KS
	kernel  : kernel function
	          options: biweight, epanechnikov (default), epan2 (DGM, 2001), epan4, normal, rectangle, triangular
	cbw     : constant 'C' to be chosen for the bandwidth in DGM (AoS, 2001) (default = 1, bw = cbw*(n^(-1/3q)))
	bootdist: distribution with bounded support, zero mean, and unit variance
	          options: mammen (default), rademacher, and uniform
	bootnum : number of bootstrap samples (default = 500)

Outcome:
	r(stat) : scalar value of the Cramer-von Mises statistic
	r(pstat): P[stat < statnstar]

If unspecified, the command runs on a default setting.
*/

*************** MAIN DGMTEST CODE ********************************************
program define dgmtest
		version 12
		
		syntax varlist(min=2) [if] [in] [, q(integer 1) teststat(string) kernel(string) cbw(real 1) bootdist(string) bootnum(integer 500)]
		marksample touse
		
		gettoken Y W : varlist
		local Y `Y'
		local W `W'
	
	// Checking inputs
		// q should be a positive integer, but less than p+q, dimension of W
		//if (`q' < 1) {
		//display "q should be a positive integer"
		//error 111
		//}
		//else if (`q' > r(W)) {
		//display "q should be less than or equal to the dimension of W"
		//error 111
		//}
		
		if ("`teststat'" == "") {
		display "Test statistic is CvM statistic (default)"
		local teststat = "CvM"
		}
		else if (inlist("`teststat'","CvM","KS")) {
		display "Test statistic: `teststat'"
		}
		else {
		display "Choose a test statistic between Cramer-von Mises (CvM; default) and Kolmogorov-Smirnov (KS) statistics"
		error 111
		}
		
		if ("`kernel'" == "") {
		display "Kernel: epanechnikov, which is our default"
		local kernel = "epanechnikov"
		}
		else if (inlist("`kernel'","biweight","epanechnikov","epan2","epan4","normal","rectangle","triangular")) {
		display "Kernel: `kernel'"
		}
		else {
		display "Choose a kernel among biweight, epanechnikov (default), epan2, epan4, normal, rectangle, and triangular"
		error 111
		}
		
		// Constant 'cbw' for bandwidth should be a positive real number
		//if (`cbw' <= 0) {
		//display "bandwidth should greater than zero"
		//error 111
		//}
		
		if ("`bootdist'" == "") {
		local bootdist = "mammen"
		display "We generated r.v. in Mammen (AoS, 1993) for V, which is our default"
		}
		else if (inlist("`bootdist'","mammen","uniform","rademacher")) {
		display "We generated `bootdist' r.v. for V"
		}
		else {
		display "Choose a distribution among Mammen (default), Rademacher and Uniform for V"
		error 111
		}
		
		// Our program is based on at least one bootstrap sample
		//if (`bootnum' <= 0) {
		//display "Choose a positive integer for our bootstrap sample"
		//error 111
		//}

	// Running main program
		mata: test("`Y'","`W'",`q',"`teststat'","`kernel'",`cbw',"`bootdist'",`bootnum')
		display as txt " `teststat' = " as res r(stat)
		display as txt " p(`teststat' < `teststat'*) = " as res r(pstat)
end

*************** Define a mata function computing test statistic ****************
mata:

// Significance test
void test(string scalar yname, string scalar wname, real scalar q,
          string scalar teststat, string scalar kernel, real scalar cbw,
          string scalar bootdist, real scalar bootnum)
{	real matrix W, X
	real colvector Y, statst
	real scalar stat, pstat

	Y = st_data(., yname, 0)
	W = st_data(., wname, 0)
	X = W[|1,1 \ rows(Y),q|]
	bw = cbw*(rows(Y)^(-1/(3*q)))
	
	stat	= mstat(Y, W, X, teststat, kernel, bw)
	statst  = mstat(btrs(Y, X, kernel, bw, bootdist, bootnum), W, X, teststat, kernel, bw)
	
	pstat	= sum(statst :> stat)/bootnum
	st_numscalar("r(stat)", stat)
	st_numscalar("r(pstat)", pstat)
}

// Test statistic for conditional mean and its bootstrap versions
real colvector mstat(real matrix Y, real matrix W, real matrix X,
                     string scalar teststat, string scalar kernel,
                     real scalar bw)
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
	W1 = iota
	for (i = 1; i <= cols(W); i++) {
		dW = W[|1,i \ n,i|]*iota' - iota*W[|1,i \ n,i|]'
		dW = uniqrows((dW:<=0)')'
		W1 = (W1#J(1,cols(dW),1)) :* (J(1,cols(W1),1)#dW)
		W1 = uniqrows(W1')'
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
end
