// This do file provides Monte Carlo simulation results with dgmtest.ado file
// Modify the number of seed and reps, and run


capture program drop simul_dgmtest
clear

// Set the seed
set seed 12

// Simulation
// Cn	: reps by 1 vector of Cn values;
// pCn	: reps by 1 vector of P[Cn < Cnstar]
quietly timer clear 1
quietly timer on 1
simulate Cn = r(Cn) pCn = r(pCn), reps(100) nodots: simul_dgmtest
quietly timer off 1
quietly timer list 1
display "It took " r(t1) " seconds for 1000 Monte Carlo samples using 1000 bootstrap samples"

// Proportion of rejections in #reps Monte Carlo samples with different alpha
quietly generate pCn90 = (pCn<0.1)
quietly generate pCn95 = (pCn<0.05)
quietly generate pCn99 = (pCn<0.01)

// We just need to report the mean values in the table
quietly mean pCn90 pCn95 pCn99
matrix list e(b)
