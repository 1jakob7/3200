% Jakob Horvath, u1092049
% Computes the lower and upper triangular matrices for A and B from LU
% decomposition.

A = [4 1 -2; 4 4 -3; 8 4 2];
B = [2 1 -2; 4 4 -3; 8 4 4];

[L1,U1] = lu(A);
[L2,U2] = lu(B);

disp('Lower Triangular');
disp(L1);
disp(L2);
disp('Upper Triangular');
disp(U1);
disp(U2);