set obs 200

// true regressor
generate Xstar = runiform()

// measurement error in X
generate etaX = runiform()

// mismeasured regressor
generate X = Xstar + 0.5*etaX

// additively linear control variable
generate W2 = runiform()

// measurement error in Z
generate etaZ = runiform()

// second measurement of true regressor
generate Z = Xstar + 0.5*etaZ

// regression error
generate epsilon = runiform()

// outcome equation
generate Y1 = Xstar^2 + 0.2*Xstar + 0.5*epsilon
generate Y2 = 0.5*W2 + Xstar^2 + 0.2*Xstar + 0.5*epsilon

// perform the test of the hypothesis of no measurement error in X
dgmtest Y1 X Z, kernel(epan2)
dgmtest Y2 X W2 Z, qw2(1) kernel(epan2)
