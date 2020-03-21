% Jakob Horvath, u1092049
% Estimates the condition number of various sparse/diagonal matrices. For
% n=161, solves the system of equations as well for some values of 'a'. 

warning('off', 'MATLAB:cond:SparseNotSupported'); % ignore sparse matrix warning when calculating condition number

n = [21, 41, 81, 161];
a = [1.0, 1.e-1, 1.e-3, 1.e-5, 1.e-7, 1.e-9, 1.e-11 1.e-13, 1.e-15];
C = zeros(9, 4);

for w=1:length(n)
    e = ones(n(w), 1);
    M = spdiags([e -2*e e], -1:1, n(w), n(w));
    for j=1:length(a)
        for i=round(n(w)/2):n(w)-1
           M(i, i-1) = a(j);
           M(i, i) = -2*a(j);
           M(i, i+1) = a(j);
        end
        % these matrix entries still need to be updated after previous loop
        M(round(n(w)/2), round(n(w)/2)) = -(1+a(j)); % absolute middle entry
        M(round(n(w)/2), round(n(w)/2+1)) = a(j); % ...
        M(round(n(w)/2), round(n(w)/2-1)) = 1; % ...
        M(n(w), n(w)-1) = a(j); % last column of the matrix
        M(n(w), n(w)) = -2*a(j); % last column of the matrix
        %
        C(j, w) = cond(M, 1); % a p-norm of 2 or Inf will produce the same condition number
    end
end

X = zeros(161, 3);
X2 = zeros(161, 3);
a = [1.0, 1.e-5, 1.e-15];
b = zeros(161, 1);
b(1) = -8;
for v=1:length(a)
    b(161) = -a(v)*4;
    [L,U,P] = lu(M);
    y = L\(P*b);
    x = U\y; % same as x = A \ b
    x2 = M*x; % store the matrix-vector multiplication using the estimated vector, 'x'
    for i=1:length(x)
       X(i, v) = x(i);
       X2(i, v) = x2(i);
    end
end