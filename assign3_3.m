% Jakob Horvath, u1092049
% Compares the linsolve and pvand functions in terms of timing and
% accuracy. Also, obtains the condition number of computed Vandermonde
% matrices.

maxRank = 300; % produces reasonable results
X1 = zeros(maxRank, maxRank/10); % holds the resulting x vectors of linsolve
X2 = zeros(maxRank, maxRank/10); % holds the resulting x vectors of pvand
T = zeros(maxRank/10, 2); % holds timing data of linsolve and pvand
C = zeros(maxRank/10, 1); % holds condition numbers resulting from vander

for N=10:10:maxRank
    x = linspace(0, 1, N);
    x = x.';
    b = func3_3(x); % e^x
    A = vander(x);
    C(N/10) = cond(A, 1); % store condition number of the current matrix
    
    tic
    X1(1:N, N/10) = linsolve(A, b);
    T(N/10, 1) = toc; % store timing data for calculating Ax=b using linsolve
    
    tic
    X2(1:N, N/10) = pvand(N, x, b).';
    T(N/10, 2) = toc; % store timing data for calculating Ax=b using pvand
end