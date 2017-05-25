/*
	Significance Testing in Nonparametric Regression Based on the Bootstrap

18/05/2017

This program is developed for a testing methodology, in Delgado and Gonzalez-Manteiga (AoS, 2001),
for selecting explanatory variables in nonparametric regression.

	H0: E(Y|W) = m(X) a.s.,
	where	m(.) = E(Y|X=.),
			Y is a scalar, and
			W = (X,Z), X is R^q valued and Z is R^p valued.

Syntax:
	dgmtest depvar expvar [if] [in]
		[, q(integer) kernel(string) bootdist(string) bw(real) bootnum(integer)]
	where
	q	  : dimension of X (default = 1)
	kernel: kernel function
			options: biweight, epanechnikov (default), epan2 (DGM, 2001), normal, rectangle, triangular
	bw	  : bandwidth of kernel function (default = 0.21544 = 100^(-1/3))
	bootdist: distribution with bounded support, zero mean, and unit variance
			options: mammen (default), rademacher, and uniform
	bootnum : number of bootstrap samples (default = 500)

Outcome:
	r(Cn) : scalar value of the Cramer-von Mises statistic
	r(pCn): P[Cn < Cnstar]

If unspecified, the command runs on a default setting.
*/

*************** MAIN DGMTEST CODE ********************************************
program define dgmtest
		version 12
		
		syntax varlist(min=2) [if] [in] [, q(integer 1) kernel(string) bw(real 0.21544) bootdist(string) bootnum(integer 500)]
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
		
		if ("`kernel'" == "") {
		display "Kernel: epanechnikov, which is our default"
		local kernel = "epanechnikov"
		}
		else if (inlist("`kernel'","biweight","epanechnikov","epan2","normal","rectangle","triangular")) {
		display "Kernel: `kernel'"
		}
		else {
		display "Choose a kernel among biweight, epanechnikov (default), epan2, normal, rectangle, and triangular"
		error 111
		}
		
		// bandwidth is a positive real number
		//if (`bw' <= 0) {
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
		mata: Cnstar("`Y'","`W'",`q',"`kernel'",`bw',"`bootdist'",`bootnum')
		display as txt " Cn = " as res r(Cn)
		display as txt " p(Cn < Cnstar) = " as res r(pCn)
end

*************** Define a mata function computing test statistic ****************
mata:

// Cn and Cnstar
void Cnstar(string scalar yname, string scalar wname, real scalar q,
			string scalar kernel, real scalar bw, string scalar bootdist,
			real scalar bootnum)
{	real matrix W, X
	real colvector Y, Cnst
	real scalar Cn0, pCn

	Y = st_data(., yname, 0)
	W = st_data(., wname, 0)
	X = W[|1,1 \ rows(Y),q|]
	
	Cn0	 = Cn(Y, W, X, kernel, bw)
	Cnst = Cn(btrs(Y, X, kernel, bw, bootdist, bootnum), W, X, kernel, bw)
	pCn  = sum(Cnst :> Cn0)/bootnum
	
	st_numscalar("r(Cn)", Cn0)
	st_numscalar("r(pCn)", pCn)
}

// Cn
real colvector Cn(real matrix Y, real matrix W, real matrix X,
				  string scalar kernel, real scalar bw)
{	real matrix K, dX, dY, Tij, Ti, W1, dW
	real colvector vecY, iotan, iotab, CnV
	real rowvector Tn
	real scalar n, i
	
	n = rows(Y)
	K = J(n,n,1)
	iotan = J(n,1,1)
	iotab = J(cols(Y),1,1)
	for (i = 1; i <= cols(X); i++) {
		dX = X[|1,i \ n,i|]#iotan' - X[|1,i \ n,i|]'#iotan
		K  = K :* Kij(dX:/bw, kernel)
		}
	vecY = vec(Y)
	dY   = vecY :- (Y'#iotan)
	Tij  = (iotab#K) :* dY
	Ti 	 = colshape(rowsum(Tij)',n)'
	W1   = J(n,n,1)
	for (i = 1; i <= cols(W); i++) {
		dW = W[|1,i \ n,i|]#iotan' - W[|1,i \ n,i|]'#iotan
		W1 = W1 :* (dW:<=0)
		}
	Tn	= colsum((Ti#iotan'):*(iotab'#W1)) :/ (n^2 * bw^cols(X))
	CnV = rowsum(colshape(Tn:^2,n))
	return(CnV)
}

// Bootstrap resample
real matrix btrs(real colvector Y, real matrix X, string scalar kernel,
				 real scalar bw, string scalar bootdist, real scalar bootnum)
{	real matrix V
	real colvector mX, e
	
	mX= npreg(Y, X, kernel, bw)
	e = Y - mX
	V = uniform(rows(Y),bootnum)
	if (bootdist == "mammen") {
	// P[V=(1-sqrt(5))/2] = (sqrt(5)+1)/(2*sqrt(5)) and P[V=(1+sqrt(5))/2] = (sqrt(5)-1)/(2*sqrt(5))
	V = (V :> ((sqrt(5)+1)/(2*sqrt(5)))) :* sqrt(5) :+ ((1-sqrt(5))/2)
	}
	else if (bootdist == "uniform") {
	// Uniform(-sqrt(3),sqrt(3))
	V = (V :- 1/2) :* 2 :* sqrt(3)
	}
	else if (bootdist == "rademacher") {
	// Rademacher: P(V=-1) = P(V=1) = 1/2
	V = (V :> 1/2) :* 2 :- 1
	}
	Ystar = mX :+ (e:*V)
	return(Ystar)
}

// Nonparametric regression
real colvector npreg(real colvector Y, real matrix X, string scalar kernel,
					 real scalar bw)
{	real matrix K, dX
	real colvector iota, fX, mX
	real scalar n, i
	
	n = rows(Y)
	K = J(n,n,1)
	iota = J(n,1,1)
	for (i = 1; i <= cols(X); i++) {
		dX = X[|1,i \ n,i|]#iota' - X[|1,i \ n,i|]'#iota
		K  = K :* Kij(dX:/bw, kernel)
		}
	fX = rowsum(K)
	mX = rowsum(K :* (Y'#iota)) :/ fX
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
