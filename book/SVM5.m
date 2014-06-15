% This runs time tests on the SVM solver with the full rank and low rank
% representations.  It also compares the cost of the sparse SVM evaluations
% to ignorantly performing a full summation, which would be the associated
% cost of using an RBF network.  There is also an option to study the cost
% of different parameterizations for solving the quadratic program.

% To allow for the low-rank expansion parameter to be set
global GAUSSQR_PARAMETERS

% Choose which test you want to run
% 1 - Fixed parameterization, increasing size, low-rank vs. standard
% 2 - Range of SV amounts and associated summation cost
% 3 - Cost of solving optimization problem with various ep and bc values
test_opt = 1;