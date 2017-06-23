// This do file provides Monte Carlo simulation results with dgmtest.ado file
// Modify the number of seed and reps, and run


capture program drop simul_dgmtest
clear

// Set the seed
set seed 12

// Simulation
// stat  : reps by 1 vector of test statistic values;
// pstat : reps by 1 vector of P[stat < stat*]
quietly timer clear 1
quietly timer on 1
simulate stat = r(stat) pstat = r(pstat), reps(2000) nodots: simul_dgmtest
quietly timer off 1
quietly timer list 1
display "It took " r(t1) " seconds for 2000 Monte Carlo samples using 2000 bootstrap samples"

// Proportion of rejections in #reps Monte Carlo samples with different alpha
quietly generate pstat90 = (pstat<0.1)
quietly generate pstat95 = (pstat<0.05)
quietly generate pstat99 = (pstat<0.01)

// We just need to report the mean values in the table
quietly mean pstat90 pstat95 pstat99
matrix list e(b)
