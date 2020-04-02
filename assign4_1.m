% Jakob Horvath, u1092049
% Finds an estimate for the vector 'b' in a linear equation of the
% form 'Ax=b' using the Jacobi iterative method.
% *** Does not give a reasonable estimate ***

% Construct matrix A %
N = 8;
e = ones(N, 1);
A = spdiags([e -4*e 6*e -4*e e], -2:2, N, N);
% Adjust main diagonal entries
A(1, 1) = 9;
A(7, 7) = 5;
A(8, 8) = 1;
% Adjust offset diagonal entries
A(7, 8) = 2;
A(8, 7) = -2;

% Construct vectors b and x %
b = 1/(N^4)*ones(8, 1);
x = zeros(8, 1);
x_new = zeros(8, 1); % needed for computing 'c' (convergence rate)

normVal = Inf;
tol = 1.e-5;
itr = 0;

% Compute estimate for 'b' using Jacobi method
while normVal > tol
    for i = 1:N
       sigma = 0;
       for j = 1:i-1
          sigma = sigma + A(i, j)*x(j); 
       end
       for j = i+1:N
          sigma = sigma + A(i, j)*x(j); 
       end
       x_new(i) = (1/A(i, i))*(b(i)-sigma);
    end
    itr = itr + 1;
    if (itr > 1)
        c = norm(x_new - x)/norm(x - x_old);
        normVal = (c/(1-c))*norm(x - x_old);
    end
    x_old = x;
    x = x_new;
end

fprintf('Solution of the system is : \n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f in %d iterations \n', x_new, itr);
b_est = A*x_new;