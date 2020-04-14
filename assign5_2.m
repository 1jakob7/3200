% Jakob Horvath, u1092049
% Finds the largest and/or smallest eigenvalues/eigenvectors of various
% matrices. In addition, applies this system to a biology example with
% birth and death rates. Makes use of the PowerMethod.m and PowerMethodInv.m functions.

tol = 1.e-7;

% a)
A = [1 2 3; 4 5 6; 7 8 9];
v1 = [1 2 1]';
vn = PowerMethod(A, v1, tol);
vlambda = vn'*A*vn / (vn'*vn);
eA = eig(A);
fprintf('\n Estimated largest eigenvalue of A: %f', vlambda);
fprintf('\n Actual largest eigenvalue of A: %f \n\n', eA(1));

% b) and c)
B = [2 3 2; 1 0 -2; -1 -3 -1];
g1 = [2 3 2]';
gn = PowerMethod(B, g1, tol); % evaluates to [0 1 0] and stops updating
glambda = gn'*B*gn / (gn'*gn);
eB= eig(B);
fprintf('\n Estimated largest eigenvalue of B: %f', glambda);
fprintf('\n Actual largest eigenvalue of B: %f \n\n', eB(1));
gnew = [1 -1 1]'; 
gnnew = PowerMethod(B, gnew, tol); % evaluates to itself
gnewlambda = gnnew'*B*gnnew / (gnnew'*gnnew);
fprintf('\n Estimated largest eigenvalue of B: %f', gnewlambda);
fprintf('\n Actual largest eigenvalue of B: %f \n\n', eB(1));

% d) 
x = [1, 2, 1]';
rn = PowerMethodInv(A, x, tol);
rn = linsolve(A, rn);
rlambda = rn'*A*rn / (rn'*rn);
fprintf('\n Estimated smallest eigenvalue of A: %d', rlambda);
fprintf('\n Actual smallest eigenvalue of A: %d \n\n', eA(3));

% e)
b = [0.3, 0.3, 0.3, 0.1];
d = [0.1, 0.2, 0.5, 0.9];
M = [b(1) b(2) b(3) b(4); 1-d(1) 0 0 0; 0 1-d(2) 0 0; 0 0 1-d(3) 1-d(4)];
q1 = [1, 1, 1, 1]';
qn = PowerMethod(M, q1, tol);
qlambda = qn'*M*qn / (qn'*qn);
eM = eig(M);
fprintf('\n Estimated largest eigenvalue of M: %f\n\n', qlambda); % < 1 so eventually the pop. will die out
q0 = [100, 200, 150, 75]';
L = (M^1000)*q0; % close to 0 -> the pop. has died out
Mnew = [b(1) b(2) b(3) b(4); 1-d(1) 0 0 0; 0 1-d(2) 0 0; 0 0 1-d(3) 1-0.01];
qnew = PowerMethod(Mnew, q1, tol);
qnewlambda = qnew'*Mnew*qnew / (qnew'*qnew);
fprintf('\n Estimated largest eigenvalue of Mnew: %f\n\n', qnewlambda); % > 1 so pop. growth is unbounded

