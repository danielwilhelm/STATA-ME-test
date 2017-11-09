/*
	Monte Carlo study in Wilhelm (2017, working)
	
25/10/2017

Samples are generated according to the model
	
	Y = Xstar^2 + 0.5*Xstar + sigME*uY
	where
	uY ~ N(0,1)
	Xstar ~ i.i.d. U(0,1) and independent of uY

	Modify (1) the number of observation,
	       (2) sigME and lamb to change stanard deviation of ME in X, or
	       (3) add specific option values in dgmtest
	           default: q(1) ql(0) teststat(CvM) kernel(epanechnikov) bw(0) bootdist(mammen) bootnum(500) ngrid(0) qgrid(0)
*/

*************** EXAMPLE_WILHELM_2017 CODE *******************************************
program define example_Wilhelm_2017, eclass

capture program drop dgmtest
clear

// Number of observations = sample size (Wilhelm (2017): n=200 or 500)
set obs 500

// Changing standard deviation of ME in X with lambda and sigma_ME
scalar sigME 	= 1
scalar lamb 	= 0.75

// Error u are generated from N(0,1)
generate uY 	= rnormal()
generate uME    = rnormal()
generate uZ 	= rnormal()

// Xstar are generated from U[0,1]
generate Xstar 	= runiform()

// Bernoulli random variable D with parameter 1-lambda
if (lamb>0 & lamb<1) {
generate D = rbinomial(1,1-lamb)
}
else {
scalar D = 1-lamb
}

// Model 1
generate Y1	= Xstar^2 + 0.5*Xstar + 0.5*uY
generate X1 = Xstar + D*sigME*uME
generate Z1 = Xstar + 0.3*uZ
// Model 2
generate Y2	= Xstar^2 + 0.5*Xstar + 0.5*uY
generate X2 = Xstar + D*sigME*uME*exp(-abs(Xstar-0.5))
generate Z2 = Xstar + 0.3*uZ
// Model 3
generate Y3	= Xstar^2 + 0.5*Xstar + 0.5*uY
generate X3 = Xstar + D*sigME*uME*exp(-abs(Xstar-0.5))
generate Z3 = Xstar + 0.3*uZ*exp(-abs(Xstar-0.5))
// Model 4
generate Y4 = Xstar^2 + 0.5*Xstar + 0.2*uY
generate X4 = Xstar + D*sigME*uME
generate Z4 = -(Xstar-1)^2 + 0.2*uZ

// Testing for the presence of measurement errors
// 0.1*[200^(-1/3)]=0.0171; 0.2*[200^(-1/3)]=0.0342; 0.5*[200^(-1/3)]=0.0855;
// 0.1*[500^(-1/3)]=0.0126; 0.2*[500^(-1/3)]=0.0252; 0.5*[500^(-1/3)]=0.0630;
//dgmtest Y1 X1 Z1, kernel(epan2) bootnum(100) ngrid(10) qgrid(0.05)
//dgmtest Y2 X2 Z2, kernel(epan2) bootnum(100) ngrid(10) qgrid(0.05)
//dgmtest Y3 X3 Z3, kernel(epan2) bootnum(100) ngrid(10) qgrid(0.05)
dgmtest Y4 X4 Z4, kernel(epan2) bootnum(100) ngrid(10) qgrid(0.05)
end
